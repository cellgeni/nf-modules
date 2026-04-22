#!/usr/bin/env python3
"""
pe_to_ome_tif.py

Convert PerkinElmer Columbus/Harmony plate acquisitions to OME-TIFF.

Reads an Index.idx.xml (or Index.xml) file, locates the per-plane TIFF
images in the same directory, and writes one OME-TIFF per well-field.

Flat-field (vignetting) correction
------------------------------------
The script automatically looks for illumination-correction reference images
referenced in the index XML (<Maps> section, Harmony 4.x+ format, or the
older <IlluminationCorrectionImages> format).  When found, each raw plane is
corrected before writing:

  If flatfield (FF) + darkfield (DF) are both present:
      corrected = (raw - DF) / (FF - DF)  ×  mean(FF - DF)

  If only flatfield (FF) is present:
      corrected = raw / FF  ×  mean(FF)

The result is clipped to the original dtype range to avoid wrap-around.

The script raises a clear WARNING if no correction data is found and exits
with an error if correction files are referenced in the XML but cannot be
located on disk (use --no-correction to skip the check entirely).

Z-projection
------------
A max-intensity Z-projection is applied by default.
Use --no-project to keep all Z planes.

Usage examples
--------------
  # Max Z-project + flat-field correction (defaults):
  python pe_to_ome_tif.py /data/experiment/Index.idx.xml

  # Skip flat-field correction:
  python pe_to_ome_tif.py Index.idx.xml --no-correction

  # Keep all Z planes, write to a specific folder:
  python pe_to_ome_tif.py Index.idx.xml --no-project -o /output

  # Process specific wells only:
  python pe_to_ome_tif.py Index.idx.xml --wells A01,A02,B01

  # Mean projection with zstd compression:
  python pe_to_ome_tif.py Index.idx.xml -m mean -c zstd

OME companion file
------------------
After conversion a plate.companion.ome is written to the output directory.
Opening it in Fiji/Bio-Formats, napari, QuPath, or OMERO gives a single
HCS plate view that links all the individual well-field OME-TIFFs without
copying any pixel data.  Use --no-companion to suppress this step.

Requirements
------------
  pip install tifffile numpy
"""

import argparse
import ast
import re
import sys
import uuid
import xml.etree.ElementTree as ET
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np
import tifffile


OME_NS  = "http://www.openmicroscopy.org/Schemas/OME/2016-06"
METHODS = ("max", "mean", "min", "sum", "median")
M_TO_UM = 1e6   # metres → micrometres


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class PlateInfo:
    name: str
    plate_id: str
    measurement_date: str
    plate_type: str
    n_rows: int
    n_cols: int
    size_x: int
    size_y: int
    pixel_size_x_m: float
    pixel_size_y_m: float


@dataclass
class ImageEntry:
    """One plane/channel/timepoint entry from the PE index."""
    row: int           # 1-based plate row
    col: int           # 1-based plate column
    field: int         # 1-based field-of-view index
    plane: int         # 1-based Z-plane index (PlaneID)
    timepoint: int     # 1-based timepoint (TimepointID)
    channel: int       # 1-based channel index (ChannelID)
    channel_name: str
    filename: str      # filename relative to the index directory
    pos_x_m: float     # stage X in metres
    pos_y_m: float     # stage Y in metres
    z_offset_m: float  # Z offset in metres
    excitation_nm: float
    emission_nm: float
    exposure_ms: float
    time_offset_s: float
    abs_time: str
    ff_already_applied: bool  # True if PE already applied FF correction


@dataclass
class CorrectionMaps:
    """
    Per-channel illumination-correction arrays.

    flatfield  – 2-D float32 array representing the illumination profile.
                 Pixel values reflect relative brightness across the field.
    darkfield  – 2-D float32 array representing the camera offset / background.
                 May be None if not provided.
    channel_name – human-readable label for diagnostics.
    """
    channel_id: int
    channel_name: str
    flatfield: np.ndarray          # shape (Y, X), float32
    darkfield: Optional[np.ndarray]  # shape (Y, X), float32, or None


@dataclass
class TiffInfo:
    """Metadata for one output OME-TIFF — used to build the companion file."""
    path: Path
    filename: str        # basename only; used as FileName= in UUID elements
    well_row: int        # 1-based plate row;  0 = no well detected
    well_col: int        # 1-based plate col;  0 = no well detected
    field_idx: int       # 1-based field-of-view; 0 = not detected
    image_idx: int       # sequential 0-based index  → Image:N in OME-XML
    name: str            # display name
    size_x: int
    size_y: int
    size_z: int
    size_c: int
    size_t: int
    dim_order: str       # e.g. 'XYZCT'
    dtype_str: str       # OME type string, e.g. 'uint16'
    px_x_um: float
    px_y_um: float
    px_z_um: float
    stage_x_um: float
    stage_y_um: float
    stage_z_um: float
    channels: list[dict] # [{name, excitation_nm, emission_nm, …}, …]
    file_uuid: str       # per-file UUID used in companion TiffData/UUID


# ---------------------------------------------------------------------------
# PE XML parsing helpers
# ---------------------------------------------------------------------------

def _strip_namespaces(root: ET.Element) -> None:
    """
    Remove XML namespace prefixes from all tags in-place.

    ET.find(".//Tag") fails silently when the document has a default namespace
    because the actual stored tag is '{http://...}Tag'.  Stripping namespaces
    once after parsing lets all downstream code use plain tag names.
    """
    for elem in root.iter():
        if "}" in elem.tag:
            elem.tag = elem.tag.split("}", 1)[1]


def _float(elem, *tags, default=0.0) -> float:
    """Return the float text of the first matching child tag found."""
    for tag in tags:
        child = elem.find(tag)
        if child is not None and child.text:
            try:
                return float(child.text.strip())
            except ValueError:
                pass
    return default


def _float_or_none(elem, *tags) -> Optional[float]:
    """Return float for first matching tag, or None when tag is absent/invalid."""
    for tag in tags:
        child = elem.find(tag)
        if child is not None and child.text:
            try:
                return float(child.text.strip())
            except ValueError:
                pass
    return None


def _int(elem, *tags, default=0) -> int:
    """Return the int text of the first matching child tag found."""
    for tag in tags:
        child = elem.find(tag)
        if child is not None and child.text:
            try:
                return int(child.text.strip())
            except ValueError:
                pass
    return default


def _str(elem, *tags, default="") -> str:
    """Return the stripped text of the first matching child tag found."""
    for tag in tags:
        child = elem.find(tag)
        if child is not None and child.text:
            return child.text.strip()
    return default


def _xml_summary(root: ET.Element, depth: int = 2) -> str:
    """Return a short indented summary of the XML tree for diagnostics."""
    lines = []

    def _walk(elem, lvl):
        if lvl > depth:
            return
        children = list(elem)
        tag = elem.tag
        text = (elem.text or "").strip()[:40]
        suffix = f' = "{text}"' if text and not children else ""
        lines.append("  " * lvl + f"<{tag}>{suffix}")
        for child in children[:8]:   # cap to 8 children per level
            _walk(child, lvl + 1)
        if len(children) > 8:
            lines.append("  " * (lvl + 1) + f"… ({len(children)} children total)")

    _walk(root, 0)
    return "\n".join(lines)


def _parse_channel_map(root: ET.Element) -> dict[int, dict]:
    """
    Extract channel metadata from Harmony V6 <Maps>/<Map>/<Entry ChannelID> structure.

    Harmony V6 stores channel properties (name, pixel size, wavelengths, exposure)
    in <Entry ChannelID="N"> elements.  We identify the relevant entries by
    requiring a <ChannelName> child (the FlatfieldProfile and SkewcropParameters
    entries lack this tag).

    Returns dict: channel_id (int) → {
        name, px_x_m, px_y_m, size_x, size_y,
        excitation_nm, emission_nm, exposure_s
    }
    """
    channel_info: dict[int, dict] = {}
    for entry in root.findall(".//Entry"):
        ch_id_str = entry.get("ChannelID")
        if ch_id_str is None:
            continue
        ch_name = _str(entry, "ChannelName")
        if not ch_name:
            continue  # skip FlatfieldProfile / SkewcropParameters entries
        try:
            ch_id = int(ch_id_str)
        except ValueError:
            continue
        channel_info[ch_id] = dict(
            name=ch_name,
            px_x_m=_float(entry, "ImageResolutionX"),
            px_y_m=_float(entry, "ImageResolutionY"),
            size_x=_int(entry, "ImageSizeX"),
            size_y=_int(entry, "ImageSizeY"),
            excitation_nm=_float(entry, "MainExcitationWavelength"),
            emission_nm=_float(entry, "MainEmissionWavelength"),
            exposure_s=_float(entry, "ExposureTime"),
        )
    return channel_info


def parse_index(index_path: Path) -> tuple[PlateInfo, list[ImageEntry]]:
    """
    Parse a PerkinElmer Index.idx.xml / Index.xml.

    Handles multiple format variants:
      - Harmony 4.x+:  root=<EvaluationInputData>, plate under <Plates><Plate>
      - Older Columbus: root=<Plate> (plate IS the root element)
      - Namespaced XML: namespace prefixes are stripped before parsing

    Returns (PlateInfo, [ImageEntry, ...]).
    Raises ValueError with a diagnostic XML summary on failure.
    """
    tree = ET.parse(index_path)
    root = tree.getroot()
    _strip_namespaces(root)

    # Locate the plate container — may be root itself or a descendant
    if root.tag == "Plate":
        plate_elem = root
    else:
        plate_elem = root.find(".//Plate")

    if plate_elem is None:
        raise ValueError(
            f"No <Plate> element found in {index_path.name}.\n"
            f"XML structure (first 2 levels):\n{_xml_summary(root)}\n\n"
            f"If this is a valid PE index file, please open an issue and\n"
            f"include the first 30 lines of your index file."
        )

    # Plate-level metadata — try alternate field names used in older formats
    plate = PlateInfo(
        name=_str(plate_elem, "Name"),
        plate_id=_str(plate_elem, "PlateID", "ID"),
        measurement_date=_str(plate_elem, "MeasurementDate", "Date"),
        plate_type=_str(plate_elem, "PlateTypeName", "PlateType", "Type"),
        n_rows=_int(plate_elem, "PlateRows", "Rows", default=8),
        n_cols=_int(plate_elem, "PlateColumns", "Columns", "Cols", default=12),
        size_x=_int(plate_elem, "ImageSizeX", "SizeX", "Width"),
        size_y=_int(plate_elem, "ImageSizeY", "SizeY", "Height"),
        pixel_size_x_m=_float(plate_elem, "ImageResolutionX", "PixelSizeX",
                               "ResolutionX"),
        pixel_size_y_m=_float(plate_elem, "ImageResolutionY", "PixelSizeY",
                               "ResolutionY"),
    )

    # Image entries — Harmony uses <WellImage>, older formats use <Image>
    raw_images = root.findall(".//WellImage") or root.findall(".//Image")

    if not raw_images:
        raise ValueError(
            f"No image entries found in {index_path.name}.\n"
            f"Tried <WellImage> and <Image> elements.\n"
            f"XML structure:\n{_xml_summary(root)}"
        )

    entries: list[ImageEntry] = []
    for wi in raw_images:
        # Filename field: URL (Harmony) → Filename / File (older)
        url = _str(wi, "URL", "Filename", "File")
        if not url:
            continue

        # Channel ID: try several known field names
        ch_id = (
            _int(wi, "ChannelID")
            or _int(wi, "ChannelNr")
            or _int(wi, "Channel", default=1)
        )
        ff_applied = _str(wi, "FlatfieldCorrectionApplied", default="No").lower() == "yes"

        z_offset_m = _float_or_none(wi, "ZOffset", "ZPosition", "Z")
        if z_offset_m is None:
            # Filename fallback for datasets where Z is encoded as "...pNN..."
            # e.g. r01c03f66p14-ch1sk1fk1fl1.tiff -> 14
            m = re.search(r"[Pp](\d+)", Path(url).name)
            if m:
                z_offset_m = float(m.group(1))
            else:
                z_offset_m = 0.0
                print(
                    f"  WARNING: no Z metadata (ZOffset/ZPosition/Z) and no 'pNN' "
                    f"token in filename '{Path(url).name}'; using z_offset=0.0"
                )

        # Try to extract row/col/field from the filename as a fallback.
        # PE filenames encode all three: r02c03f01p01-ch1sk1fk1fl1.tiff
        m_rcf = re.search(r"[Rr](\d+)[Cc](\d+)[Ff](\d+)", Path(url).name)

        row_val   = _int(wi, "Row", default=0) or (int(m_rcf.group(1)) if m_rcf else 1)
        col_val   = _int(wi, "Col", "Column", default=0) or (int(m_rcf.group(2)) if m_rcf else 1)
        field_val = (
            _int(wi, "FieldID") or _int(wi, "FieldNr") or _int(wi, "Field", default=0)
            or (int(m_rcf.group(3)) if m_rcf else 1)
        )

        entries.append(ImageEntry(
            row=row_val,
            col=col_val,
            field=field_val,
            plane=_int(wi, "PlaneID") or _int(wi, "PlaneNr") or _int(wi, "Plane", default=1),
            timepoint=(
                _int(wi, "TimepointID")
                or _int(wi, "TimepointNr")
                or _int(wi, "Timepoint", default=1)
            ),
            channel=ch_id,
            channel_name=_str(wi, "ChannelName", default=f"Ch{ch_id}"),
            filename=url,
            pos_x_m=_float(wi, "PositionX", "StagePositionX"),
            pos_y_m=_float(wi, "PositionY", "StagePositionY"),
            z_offset_m=z_offset_m,
            excitation_nm=_float(wi, "MainExcitationWavelength", "ExcitationWavelength"),
            emission_nm=_float(wi, "MainEmissionWavelength", "EmissionWavelength"),
            exposure_ms=_float(wi, "ExposureTime"),
            time_offset_s=_float(wi, "MeasurementTimeOffset", "TimeOffset"),
            abs_time=_str(wi, "AbsTime", "DateTime"),
            ff_already_applied=ff_applied,
        ))

    if not entries:
        raise ValueError(
            f"No valid image entries (with a URL/Filename) found in "
            f"{index_path.name}."
        )

    # --- Harmony V6: enrich entries with channel metadata from <Maps> -------
    ch_map = _parse_channel_map(root)
    if ch_map:
        # Fill in plate spatial dimensions from the first channel entry if absent
        first_ch = next(iter(ch_map.values()))
        if not plate.size_x:
            plate.size_x = first_ch.get("size_x", 0)
        if not plate.size_y:
            plate.size_y = first_ch.get("size_y", 0)
        if not plate.pixel_size_x_m:
            plate.pixel_size_x_m = first_ch.get("px_x_m", 0.0)
        if not plate.pixel_size_y_m:
            plate.pixel_size_y_m = first_ch.get("px_y_m", 0.0)
        # Propagate channel name / wavelength / exposure to every ImageEntry
        for e in entries:
            info = ch_map.get(e.channel)
            if info is None:
                continue
            if not e.channel_name or e.channel_name.startswith("Ch"):
                e.channel_name = info["name"]
            if not e.excitation_nm:
                e.excitation_nm = info["excitation_nm"]
            if not e.emission_nm:
                e.emission_nm = info["emission_nm"]
            if not e.exposure_ms:
                e.exposure_ms = info["exposure_s"] * 1000.0  # s → ms

    return plate, entries


# ---------------------------------------------------------------------------
# Flat-field correction parsing
# ---------------------------------------------------------------------------

def parse_correction_info(
    index_path: Path,
) -> dict[int, tuple[Optional[str], Optional[str]]]:
    """
    Search an Index.idx.xml for illumination-correction file references.

    Handles two PE format variants:

    Harmony 4.x+  (<Maps><Map> with <Type>Flatfield/Darkfield</Type>):
      <Maps>
        <Map>
          <Type>Flatfield</Type>
          <ChannelID>1</ChannelID>
          <URL>r01c01-ch1ffcorr.tiff</URL>
        </Map>
        <Map>
          <Type>Darkfield</Type>
          <ChannelID>1</ChannelID>
          <URL>r01c01-ch1dfcorr.tiff</URL>
        </Map>
      </Maps>

    Older Columbus (<IlluminationCorrectionImage> with dedicated sub-elements):
      <IlluminationCorrectionImages>
        <IlluminationCorrectionImage>
          <ChannelID>1</ChannelID>
          <FlatfieldURL>ffcorr-ch1.tiff</FlatfieldURL>
          <DarkfieldURL>dfcorr-ch1.tiff</DarkfieldURL>
        </IlluminationCorrectionImage>
      </IlluminationCorrectionImages>

    Returns
    -------
    dict mapping channel_id (int, 1-based) to
        (flatfield_filename: str | None, darkfield_filename: str | None)
    Empty dict if nothing is found.
    """
    tree = ET.parse(index_path)
    root = tree.getroot()
    _strip_namespaces(root)

    result: dict[int, list[Optional[str]]] = {}  # ch_id → [ff, df]

    # --- Harmony 4.x: <Maps><Map> with <Type> --------------------------------
    for map_elem in root.findall(".//Map"):
        map_type = _str(map_elem, "Type", default="").lower()
        if map_type not in ("flatfield", "darkfield"):
            continue
        ch_id = _int(map_elem, "ChannelID", default=0)
        url   = _str(map_elem, "URL")
        if not ch_id or not url:
            continue
        if ch_id not in result:
            result[ch_id] = [None, None]  # [ff, df]
        if map_type == "flatfield":
            result[ch_id][0] = url
        else:
            result[ch_id][1] = url

    # --- Older Columbus: <IlluminationCorrectionImage> -----------------------
    for ici in root.findall(".//IlluminationCorrectionImage"):
        ch_id = _int(ici, "ChannelID", default=0)
        ff_url = _str(ici, "FlatfieldURL") or _str(ici, "URL")
        df_url = _str(ici, "DarkfieldURL")
        if not ch_id or not ff_url:
            continue
        if ch_id not in result:
            result[ch_id] = [None, None]
        result[ch_id][0] = result[ch_id][0] or ff_url
        result[ch_id][1] = result[ch_id][1] or (df_url or None)

    return {ch: (ff, df) for ch, (ff, df) in result.items()}


def load_correction_maps(
    corr_info: dict[int, tuple[Optional[str], Optional[str]]],
    image_dir: Path,
) -> dict[int, CorrectionMaps]:
    """
    Load the correction TIFF files referenced by parse_correction_info().

    Parameters
    ----------
    corr_info  : output of parse_correction_info()
    image_dir  : directory containing the index file (filenames are relative)

    Returns
    -------
    dict mapping channel_id → CorrectionMaps

    Raises
    ------
    FileNotFoundError if a referenced correction file cannot be located.
    """
    maps: dict[int, CorrectionMaps] = {}

    for ch_id, (ff_url, df_url) in corr_info.items():
        if ff_url is None:
            continue

        ff_path = _resolve_path(image_dir, ff_url)
        ff_arr  = _load_corr_image(ff_path)

        df_arr: Optional[np.ndarray] = None
        if df_url:
            df_path = _resolve_path(image_dir, df_url)
            df_arr  = _load_corr_image(df_path)

        maps[ch_id] = CorrectionMaps(
            channel_id=ch_id,
            channel_name=f"Ch{ch_id}",
            flatfield=ff_arr,
            darkfield=df_arr,
        )

    return maps


def _load_corr_image(path: Path) -> np.ndarray:
    """Load a correction TIFF as float32."""
    with tifffile.TiffFile(str(path)) as tif:
        arr = tif.pages[0].asarray().astype(np.float32)
    return arr


def _parse_flatfield_profile(text: str):
    """
    Parse the pseudo-JSON FlatfieldProfile string used in Harmony V6.

    The text uses unquoted keys so standard json.loads() fails.  We extract
    just the four fields we need with targeted regexes.

    Returns (coeffs, dims, origin, scale) or None on parse failure.
      coeffs – list-of-lists of floats (triangular polynomial coefficients)
      dims   – [H, W] image dimensions
      origin – [x0, y0] coordinate origin
      scale  – [sx, sy] normalisation scale factors
    """
    m = re.search(r"Coefficients:\s*(\[\s*\[.+?\]\s*\])", text, re.DOTALL)
    if not m:
        return None
    try:
        coeffs = ast.literal_eval(m.group(1))
    except (SyntaxError, ValueError):
        return None

    m = re.search(r"Dims:\s*\[(\d+),\s*(\d+)\]", text)
    if not m:
        return None
    dims = [int(m.group(1)), int(m.group(2))]

    m = re.search(r"Origin:\s*\[([0-9.]+),\s*([0-9.]+)\]", text)
    if not m:
        return None
    origin = [float(m.group(1)), float(m.group(2))]

    m = re.search(r"Scale:\s*\[([0-9.Ee+\-]+),\s*([0-9.Ee+\-]+)\]", text)
    if not m:
        return None
    scale = [float(m.group(1)), float(m.group(2))]

    return coeffs, dims, origin, scale


def _eval_poly2d(
    coeffs: list,
    dims: list[int],
    origin: list[float],
    scale: list[float],
) -> np.ndarray:
    """
    Evaluate a 2-D triangular polynomial on a pixel grid.

    coeffs[k][j] is the coefficient of  u^j * v^(k−j)
    where:
        u = (x_pixel − origin[0]) * scale[0]
        v = (y_pixel − origin[1]) * scale[1]

    Returns float32 array of shape (dims[0], dims[1]).
    """
    H, W = dims
    u = (np.arange(W, dtype=np.float32) - origin[0]) * scale[0]  # (W,)
    v = (np.arange(H, dtype=np.float32) - origin[1]) * scale[1]  # (H,)
    U, V = np.meshgrid(u, v)  # (H, W)

    result = np.zeros((H, W), dtype=np.float32)
    for k, row in enumerate(coeffs):
        for j, c in enumerate(row):
            result += float(c) * (U ** j) * (V ** (k - j))
    return result


def parse_flatfield_from_maps(
    index_path: Path,
) -> dict[int, CorrectionMaps]:
    """
    Parse Harmony V6 polynomial flat-field profiles from the index XML.

    Harmony V6 encodes illumination correction as a 2-D polynomial inside
    <Maps>/<Map>/<Entry ChannelID="N">/<FlatfieldProfile>.  This function
    evaluates the polynomial on a pixel grid and returns the same
    CorrectionMaps structure used by the file-based path.

    Returns dict: channel_id → CorrectionMaps  (darkfield is always None).
    Returns an empty dict if no FlatfieldProfile entries are found.
    """
    tree = ET.parse(index_path)
    root = tree.getroot()
    _strip_namespaces(root)

    result: dict[int, CorrectionMaps] = {}

    for entry in root.findall(".//Entry"):
        ch_id_str = entry.get("ChannelID")
        ff_elem = entry.find("FlatfieldProfile")
        if ch_id_str is None or ff_elem is None or not ff_elem.text:
            continue
        try:
            ch_id = int(ch_id_str)
        except ValueError:
            continue

        parsed = _parse_flatfield_profile(ff_elem.text.strip())
        if parsed is None:
            print(f"  WARNING: Could not parse FlatfieldProfile for channel {ch_id}.")
            continue

        coeffs, dims, origin, scale = parsed

        # Best-effort channel name from the pseudo-JSON text
        m = re.search(r"ChannelName:\s*([^,}]+?)(?:,|})", ff_elem.text)
        ch_name = m.group(1).strip() if m else f"Ch{ch_id}"

        result[ch_id] = CorrectionMaps(
            channel_id=ch_id,
            channel_name=ch_name,
            flatfield=_eval_poly2d(coeffs, dims, origin, scale),
            darkfield=None,
        )

    return result


# ---------------------------------------------------------------------------
# Flat-field correction application
# ---------------------------------------------------------------------------

def apply_flat_field(
    raw: np.ndarray,
    corr: CorrectionMaps,
) -> np.ndarray:
    """
    Apply flat-field (and optional dark-field) correction to a 2-D plane.

    With darkfield:
        corrected = (raw − DF) / (FF − DF)  ×  mean(FF − DF)

    Without darkfield:
        corrected = raw / FF  ×  mean(FF)

    The result is clipped to [0, dtype_max] and cast back to the input dtype
    so downstream code and file writers see no type change.
    """
    raw_f = raw.astype(np.float32)
    ff    = corr.flatfield

    if corr.darkfield is not None:
        df    = corr.darkfield
        denom = ff - df
        # Guard against zero/negative denominators (saturated or bad pixels)
        safe  = np.where(denom > 0, denom, 1.0)
        scale = float(np.mean(denom[denom > 0])) if np.any(denom > 0) else 1.0
        corrected = (raw_f - df) / safe * scale
    else:
        ff_mean = float(np.mean(ff[ff > 0])) if np.any(ff > 0) else 1.0
        safe    = np.where(ff > 0, ff, 1.0)
        corrected = raw_f / safe * ff_mean

    # Clip to the original integer range; preserve dtype
    if np.issubdtype(raw.dtype, np.integer):
        info = np.iinfo(raw.dtype)
        corrected = np.clip(corrected, info.min, info.max)
    else:
        corrected = np.clip(corrected, 0, np.finfo(raw.dtype).max)

    return corrected.astype(raw.dtype)


def apply_flat_field_stack(
    stack: np.ndarray,
    corr: CorrectionMaps,
) -> np.ndarray:
    """
    Apply flat-field (and optional dark-field) correction to a Z stack.

    stack shape: (Z, Y, X)
    """
    stack_f = stack.astype(np.float32)
    ff = corr.flatfield.astype(np.float32)

    if corr.darkfield is not None:
        df = corr.darkfield.astype(np.float32)
        denom = ff - df
        safe = np.where(denom > 0, denom, 1.0)
        scale = float(np.mean(denom[denom > 0])) if np.any(denom > 0) else 1.0
        corrected = (stack_f - df[None, :, :]) / safe[None, :, :] * scale
    else:
        ff_mean = float(np.mean(ff[ff > 0])) if np.any(ff > 0) else 1.0
        safe = np.where(ff > 0, ff, 1.0)
        corrected = stack_f / safe[None, :, :] * ff_mean

    if np.issubdtype(stack.dtype, np.integer):
        info = np.iinfo(stack.dtype)
        corrected = np.clip(corrected, info.min, info.max)
    else:
        corrected = np.clip(corrected, 0, np.finfo(stack.dtype).max)

    return corrected.astype(stack.dtype)


# ---------------------------------------------------------------------------
# Correction map export
# ---------------------------------------------------------------------------

def save_correction_maps(
    correction: dict[int, CorrectionMaps],
    output_dir: Path,
    prefix: str = "",
) -> None:
    """Save per-channel flatfield/darkfield maps as TIFF files."""
    if not correction:
        return

    stem = f"{prefix}_" if prefix else ""
    for ch_id in sorted(correction):
        corr = correction[ch_id]
        ff_path = output_dir / f"{stem}ch{ch_id:02d}_ffp.tiff"
        tifffile.imwrite(str(ff_path), corr.flatfield.astype(np.float32))

        if corr.darkfield is not None:
            df_path = output_dir / f"{stem}ch{ch_id:02d}_dfp.tiff"
            tifffile.imwrite(str(df_path), corr.darkfield.astype(np.float32))


# ---------------------------------------------------------------------------
# Grouping & well utilities
# ---------------------------------------------------------------------------

def well_label(row: int, col: int) -> str:
    """1-based (row, col) → label like 'A01'."""
    return f"{chr(ord('A') + row - 1)}{col:02d}"


def label_to_rowcol(label: str) -> tuple[int, int]:
    """'A01' → (1, 1).  Raises ValueError for bad labels."""
    label = label.strip().upper()
    m = re.fullmatch(r"([A-Z]+)(\d+)", label)
    if not m:
        raise ValueError(f"Cannot parse well label: {label!r}")
    row = sum((ord(c) - ord("A") + 1) * (26 ** i)
              for i, c in enumerate(reversed(m.group(1))))
    return row, int(m.group(2))


def parse_well_field(filename: str) -> tuple[int, int, int]:
    """
    Extract (row, col, field) from an output filename.

    Recognised patterns:
      A01_F001_maxproj.ome.tif   →  (1, 1, 1)
      B12_F003.ome.tif           →  (2, 12, 3)
      r01c02f003_max.ome.tif     →  (1, 2, 3)

    Returns (0, 0, 0) if no well pattern is found.
    """
    stem = filename
    for ext in (".ome.tif", ".ome.tiff", ".tif", ".tiff"):
        if stem.lower().endswith(ext):
            stem = stem[: -len(ext)]
            break

    # Pattern 1: {WELL}_F{field}  e.g. A01_F001, B12_F003
    # Use re.search so an optional prefix (e.g. "sampleA_A01_F001") still matches.
    m = re.search(r"([A-Za-z]{1,3}\d{1,3})_[Ff](\d+)", stem)
    if m:
        try:
            row, col = label_to_rowcol(m.group(1))
            return row, col, int(m.group(2))
        except ValueError:
            pass

    # Pattern 2: r{row}c{col}f{field}  e.g. r01c02f03
    m = re.search(r"[Rr](\d+)[Cc](\d+)[Ff](\d+)", stem)
    if m:
        return int(m.group(1)), int(m.group(2)), int(m.group(3))

    return 0, 0, 0


def _row_letter_to_int(letter: str) -> int:
    """'A' → 1, 'B' → 2, …, 'Z' → 26, 'AA' → 27, …"""
    result = 0
    for c in letter.strip().upper():
        if not c.isalpha():
            raise ValueError(f"Invalid row letter: {letter!r}")
        result = result * 26 + (ord(c) - ord("A") + 1)
    return result


def parse_row_range(spec: str) -> list[int]:
    """
    Parse a row-letter range spec into a list of 1-based row indices.

    Examples
    --------
    'B'    → [2]
    'B-D'  → [2, 3, 4]
    'AA'   → [27]
    """
    spec = spec.strip().upper()
    if "-" in spec:
        start_s, end_s = [s.strip() for s in spec.split("-", 1)]
        start, end = _row_letter_to_int(start_s), _row_letter_to_int(end_s)
        if start > end:
            raise ValueError(
                f"Row range start '{start_s}' is after end '{end_s}'."
            )
        return list(range(start, end + 1))
    return [_row_letter_to_int(spec)]


def parse_col_range(spec: str) -> list[int]:
    """
    Parse a column-number range spec into a list of 1-based column indices.

    Examples
    --------
    '4'    → [4]
    '4-6'  → [4, 5, 6]
    """
    spec = spec.strip()
    if "-" in spec:
        start_s, end_s = [s.strip() for s in spec.split("-", 1)]
        start, end = int(start_s), int(end_s)
        if start > end:
            raise ValueError(
                f"Column range start {start} is after end {end}."
            )
        return list(range(start, end + 1))
    return [int(spec)]


def group_by_well_field(
    entries: list[ImageEntry],
) -> dict[tuple[int, int, int], list[ImageEntry]]:
    """
    Returns {(row, col, field): [entries sorted by (timepoint, channel, z_offset)]}.
    """
    groups: dict[tuple, list] = defaultdict(list)
    for e in entries:
        groups[(e.row, e.col, e.field)].append(e)
    for key in groups:
        groups[key].sort(key=lambda e: (e.timepoint, e.channel, e.z_offset_m))
    return dict(groups)


# ---------------------------------------------------------------------------
# OME-XML construction
# ---------------------------------------------------------------------------

def _ome_dtype(dtype: np.dtype) -> str:
    return {
        "uint8": "uint8", "uint16": "uint16", "uint32": "uint32",
        "int8": "int8",   "int16": "int16",   "int32": "int32",
        "float32": "float", "float64": "double",
    }.get(dtype.name, "uint16")


def _sub(parent, local: str, **attribs) -> ET.Element:
    """Append a child element in the OME namespace."""
    el = ET.SubElement(parent, f"{{{OME_NS}}}{local}")
    for k, v in attribs.items():
        el.set(k, str(v))
    return el


def build_ome_xml(
    image_name: str,
    size_x: int,
    size_y: int,
    size_z: int,
    size_c: int,
    size_t: int,
    dtype: np.dtype,
    px_x_um: float,
    px_y_um: float,
    px_z_um: float,
    channels: list[dict],
    planes: list[dict],
    plate: PlateInfo,
    row: int,
    col: int,
    correction_applied: bool,
    correction_has_darkfield: bool,
) -> str:
    """
    Build an OME-XML string for one well-field.

    channels: list (ordered, 0-based) of dicts:
        name, excitation_nm, emission_nm, exposure_ms
    planes: list in IFD order (T outer → C → Z inner) of dicts:
        t, c, z, pos_x_um, pos_y_um, pos_z_um, delta_t_s, exposure_ms
    """
    ET.register_namespace("", OME_NS)
    ET.register_namespace("xsi", "http://www.w3.org/2001/XMLSchema-instance")

    ome = ET.Element(f"{{{OME_NS}}}OME")
    ome.set("xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance")
    ome.set("xsi:schemaLocation", f"{OME_NS} {OME_NS}/ome.xsd")

    # --- Instrument --------------------------------------------------------
    instr = _sub(ome, "Instrument", ID="Instrument:0")
    _sub(instr, "Microscope", Manufacturer="PerkinElmer", Model="Operetta/Opera")

    # --- Image + Pixels ----------------------------------------------------
    img = _sub(ome, "Image", ID="Image:0", Name=image_name)
    _sub(img, "InstrumentRef", ID="Instrument:0")

    pix = _sub(img, "Pixels",
               ID="Pixels:0",
               DimensionOrder="XYZCT",
               Type=_ome_dtype(dtype),
               SizeX=size_x, SizeY=size_y, SizeZ=size_z,
               SizeC=size_c, SizeT=size_t,
               PhysicalSizeX=f"{px_x_um:.6f}", PhysicalSizeXUnit="µm",
               PhysicalSizeY=f"{px_y_um:.6f}", PhysicalSizeYUnit="µm")
    if size_z > 1:
        pix.set("PhysicalSizeZ", f"{px_z_um:.6f}")
        pix.set("PhysicalSizeZUnit", "µm")

    # Channels
    for ci, ch in enumerate(channels):
        chan = _sub(pix, "Channel",
                   ID=f"Channel:0:{ci}",
                   Name=ch["name"],
                   SamplesPerPixel="1")
        if ch.get("excitation_nm"):
            chan.set("ExcitationWavelength", f"{ch['excitation_nm']:.1f}")
            chan.set("ExcitationWavelengthUnit", "nm")
        if ch.get("emission_nm"):
            chan.set("EmissionWavelength", f"{ch['emission_nm']:.1f}")
            chan.set("EmissionWavelengthUnit", "nm")
        _sub(chan, "LightPath")

    # TiffData — one entry per IFD
    for ifd, pl in enumerate(planes):
        _sub(pix, "TiffData",
             IFD=ifd,
             FirstT=pl["t"], FirstC=pl["c"], FirstZ=pl["z"],
             PlaneCount="1")

    # Plane — positions, exposure, timestamps
    for pl in planes:
        plane_el = _sub(pix, "Plane",
                        TheZ=pl["z"], TheT=pl["t"], TheC=pl["c"])
        if pl.get("exposure_ms") is not None:
            plane_el.set("ExposureTime", f"{pl['exposure_ms']:.3f}")
            plane_el.set("ExposureTimeUnit", "ms")
        if pl.get("delta_t_s") is not None:
            plane_el.set("DeltaT", f"{pl['delta_t_s']:.6f}")
            plane_el.set("DeltaTUnit", "s")
        # Always emit stage coordinates to avoid downstream "undefined" fallbacks.
        pos_x_um = float(pl.get("pos_x_um", 0.0))
        pos_y_um = float(pl.get("pos_y_um", 0.0))
        pos_z_um = float(pl.get("pos_z_um", 0.0))
        plane_el.set("PositionX", f"{pos_x_um:.3f}")
        plane_el.set("PositionXUnit", "µm")
        plane_el.set("PositionY", f"{pos_y_um:.3f}")
        plane_el.set("PositionYUnit", "µm")
        plane_el.set("PositionZ", f"{pos_z_um:.6f}")
        plane_el.set("PositionZUnit", "µm")

    # --- Plate (HCS metadata) ----------------------------------------------
    plate_el = _sub(ome, "Plate",
                    ID="Plate:0",
                    Name=plate.name or plate.plate_id,
                    Rows=plate.n_rows,
                    Columns=plate.n_cols)
    if plate.plate_type:
        plate_el.set("ExternalIdentifier", plate.plate_type)

    well_el = _sub(plate_el, "Well",
                   ID="Well:0",
                   Row=row - 1,    # OME Well rows/cols are 0-based
                   Column=col - 1)
    ws = _sub(well_el, "WellSample", ID="WellSample:0", Index="0")
    _sub(ws, "ImageRef", ID="Image:0")

    # --- StructuredAnnotations: record processing steps --------------------
    sa = _sub(ome, "StructuredAnnotations")
    ma = _sub(sa, "MapAnnotation",
              ID="Annotation:Processing:0",
              Namespace="pe_to_ome_tif.processing")
    mv = _sub(ma, "Value")
    for k, v in [
        ("FlatfieldCorrectionApplied", "Yes" if correction_applied else "No"),
        ("DarkfieldAvailable",         "Yes" if correction_has_darkfield else "No"),
        ("ConvertedBy",                "pe_to_ome_tif.py"),
    ]:
        m = _sub(mv, "M", K=k)
        m.text = v

    return ET.tostring(ome, encoding="unicode")


# ---------------------------------------------------------------------------
# OME companion file
# ---------------------------------------------------------------------------

def _read_tiff_info(path: Path, image_idx: int) -> TiffInfo:
    """
    Read OME-XML metadata from a written OME-TIFF (first IFD only — no pixels).
    Returns a TiffInfo ready for use in build_companion_xml().
    """
    with tifffile.TiffFile(str(path)) as tif:
        if not tif.is_ome:
            raise ValueError(f"{path.name} is not an OME-TIFF")
        ome_xml = tif.ome_metadata

    root = ET.fromstring(ome_xml)
    _strip_namespaces(root)
    pixels = root.find(".//Pixels")
    if pixels is None:
        raise ValueError(f"No <Pixels> in {path.name}")
    first_plane = pixels.find("Plane")
    stage_x_um = float(first_plane.get("PositionX", 0.0)) if first_plane is not None else 0.0
    stage_y_um = float(first_plane.get("PositionY", 0.0)) if first_plane is not None else 0.0
    stage_z_um = float(first_plane.get("PositionZ", 0.0)) if first_plane is not None else 0.0

    channels: list[dict] = []
    for ch_el in pixels.findall("Channel"):
        channels.append(dict(
            name=ch_el.get("Name", f"Ch{len(channels) + 1}"),
            excitation_nm=ch_el.get("ExcitationWavelength"),
            emission_nm=ch_el.get("EmissionWavelength"),
            excitation_unit=ch_el.get("ExcitationWavelengthUnit", "nm"),
            emission_unit=ch_el.get("EmissionWavelengthUnit", "nm"),
        ))
    size_c = int(pixels.get("SizeC", 1))
    if not channels:
        channels = [{"name": f"Ch{i + 1}"} for i in range(size_c)]

    row, col, fov = parse_well_field(path.name)
    name = re.sub(r"\.ome$", "", path.stem, flags=re.IGNORECASE)

    return TiffInfo(
        path=path,
        filename=path.name,
        well_row=row, well_col=col, field_idx=fov,
        image_idx=image_idx,
        name=name,
        size_x=int(pixels.get("SizeX", 1)),
        size_y=int(pixels.get("SizeY", 1)),
        size_z=int(pixels.get("SizeZ", 1)),
        size_c=size_c,
        size_t=int(pixels.get("SizeT", 1)),
        dim_order=pixels.get("DimensionOrder", "XYZCT"),
        dtype_str=pixels.get("Type", "uint16"),
        px_x_um=float(pixels.get("PhysicalSizeX", 1.0)),
        px_y_um=float(pixels.get("PhysicalSizeY", 1.0)),
        px_z_um=float(pixels.get("PhysicalSizeZ", 1.0)),
        stage_x_um=stage_x_um,
        stage_y_um=stage_y_um,
        stage_z_um=stage_z_um,
        channels=channels,
        file_uuid=str(uuid.uuid4()),
    )


def _iter_planes(dim_order: str, size_t: int, size_c: int, size_z: int):
    """Yield (t, c, z) in IFD order for the given DimensionOrder."""
    outer_to_inner = [d for d in reversed(dim_order) if d not in ("X", "Y")]
    sizes = {"T": size_t, "C": size_c, "Z": size_z}

    def _rec(dims):
        if not dims:
            yield {}
            return
        d = dims[0]
        for i in range(sizes.get(d, 1)):
            for rest in _rec(dims[1:]):
                yield {d: i, **rest}

    for combo in _rec(outer_to_inner):
        yield combo.get("T", 0), combo.get("C", 0), combo.get("Z", 0)


def build_companion_xml(
    tiff_infos: list[TiffInfo],
    plate: PlateInfo,
) -> str:
    """
    Return the OME-XML string for a companion file that links all TIFFs.

    Emits a <Plate>/<Well>/<WellSample> structure when well labels were
    detected in filenames; falls back to a plain multi-image XML otherwise.
    """
    ET.register_namespace("", OME_NS)
    ET.register_namespace("xsi", "http://www.w3.org/2001/XMLSchema-instance")

    ome = ET.Element(f"{{{OME_NS}}}OME")
    ome.set("xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance")
    ome.set("xsi:schemaLocation", f"{OME_NS} {OME_NS}/ome.xsd")

    # ---- One Image element per TIFF ----------------------------------------
    for ti in tiff_infos:
        img = _sub(ome, "Image", ID=f"Image:{ti.image_idx}", Name=ti.name)
        pix = _sub(img, "Pixels",
                   ID=f"Pixels:{ti.image_idx}",
                   DimensionOrder=ti.dim_order,
                   Type=ti.dtype_str,
                   SizeX=ti.size_x, SizeY=ti.size_y, SizeZ=ti.size_z,
                   SizeC=ti.size_c, SizeT=ti.size_t,
                   PhysicalSizeX=f"{ti.px_x_um:.6f}", PhysicalSizeXUnit="µm",
                   PhysicalSizeY=f"{ti.px_y_um:.6f}", PhysicalSizeYUnit="µm")
        if ti.size_z > 1 and ti.px_z_um > 0:
            pix.set("PhysicalSizeZ", f"{ti.px_z_um:.6f}")
            pix.set("PhysicalSizeZUnit", "µm")

        for ci, ch in enumerate(ti.channels):
            chan = _sub(pix, "Channel",
                        ID=f"Channel:{ti.image_idx}:{ci}",
                        Name=ch["name"], SamplesPerPixel="1")
            if ch.get("excitation_nm"):
                chan.set("ExcitationWavelength", str(ch["excitation_nm"]))
                chan.set("ExcitationWavelengthUnit", ch.get("excitation_unit", "nm"))
            if ch.get("emission_nm"):
                chan.set("EmissionWavelength", str(ch["emission_nm"]))
                chan.set("EmissionWavelengthUnit", ch.get("emission_unit", "nm"))
            _sub(chan, "LightPath")

        # TiffData: one entry per plane, each referencing this TIFF by filename
        for ifd, (t, c, z) in enumerate(
            _iter_planes(ti.dim_order, ti.size_t, ti.size_c, ti.size_z)
        ):
            td = _sub(pix, "TiffData",
                      IFD=ifd, FirstT=t, FirstC=c, FirstZ=z, PlaneCount="1")
            uuid_el = _sub(td, "UUID")
            uuid_el.set("FileName", ti.filename)
            uuid_el.text = f"urn:uuid:{ti.file_uuid}"
            # Preserve stage coordinates in companion metadata so readers like
            # Ashlar can use them as alignment priors.
            plane = _sub(pix, "Plane", TheZ=z, TheC=c, TheT=t)
            plane.set("PositionX", f"{ti.stage_x_um:.3f}")
            plane.set("PositionXUnit", "µm")
            plane.set("PositionY", f"{ti.stage_y_um:.3f}")
            plane.set("PositionYUnit", "µm")
            plane.set("PositionZ", f"{ti.stage_z_um:.6f}")
            plane.set("PositionZUnit", "µm")

    # ---- Plate / HCS structure ---------------------------------------------
    hcs_infos = [ti for ti in tiff_infos if ti.well_row > 0]
    if hcs_infos:
        n_rows = plate.n_rows or max(ti.well_row for ti in hcs_infos)
        n_cols = plate.n_cols or max(ti.well_col for ti in hcs_infos)
        plate_el = _sub(ome, "Plate",
                        ID="Plate:0",
                        Name=plate.name or plate.plate_id or "Plate",
                        Rows=n_rows, Columns=n_cols)
        if plate.plate_type:
            plate_el.set("ExternalIdentifier", plate.plate_type)

        wells: dict[tuple[int, int], list[TiffInfo]] = defaultdict(list)
        for ti in hcs_infos:
            wells[(ti.well_row, ti.well_col)].append(ti)

        for well_idx, ((row, col), fields) in enumerate(sorted(wells.items())):
            well_el = _sub(plate_el, "Well",
                           ID=f"Well:{well_idx}",
                           Row=row - 1, Column=col - 1)  # OME is 0-based
            for s_idx, ti in enumerate(
                sorted(fields, key=lambda t: t.field_idx or t.image_idx)
            ):
                ws = _sub(well_el, "WellSample",
                          ID=f"WellSample:{ti.image_idx}", Index=s_idx)
                _sub(ws, "ImageRef", ID=f"Image:{ti.image_idx}")

    return ET.tostring(ome, encoding="unicode")


def make_companion(
    output_dir: Path,
    companion_path: Path,
    plate: PlateInfo,
) -> None:
    """
    Scan *output_dir* for OME-TIFFs and write a companion file to
    *companion_path*.  Plate metadata comes from the already-parsed *plate*
    object so the index XML is not re-read.
    """
    tif_files = sorted(
        p for p in output_dir.iterdir()
        if p.name.lower().endswith(".ome.tif")
        or p.name.lower().endswith(".ome.tiff")
    )
    if not tif_files:
        print("  WARNING: no OME-TIFFs found in output dir — companion not written.")
        return

    tiff_infos: list[TiffInfo] = []
    for idx, path in enumerate(tif_files):
        try:
            tiff_infos.append(_read_tiff_info(path, idx))
        except Exception as exc:
            print(f"  WARNING: skipping {path.name} from companion: {exc}")

    if not tiff_infos:
        print("  WARNING: could not read metadata from any output file.")
        return

    xml_str = build_companion_xml(tiff_infos, plate)
    companion_path.write_text(
        '<?xml version="1.0" encoding="UTF-8"?>\n' + xml_str,
        encoding="utf-8",
    )

    hcs_count = sum(1 for ti in tiff_infos if ti.well_row > 0)
    n_wells    = len({(ti.well_row, ti.well_col)
                      for ti in tiff_infos if ti.well_row > 0})
    mode = "HCS plate" if hcs_count else "single"
    print(
        f"\nCompanion   : {companion_path.name}  ({mode})\n"
        + (f"  {n_wells} wells, {hcs_count} well-fields\n" if hcs_count else "")
        + f"  {len(tiff_infos)} Image series"
    )


# ---------------------------------------------------------------------------
# Per-well-field writer
# ---------------------------------------------------------------------------

def _project(frames: np.ndarray, method: str) -> np.ndarray:
    """frames shape: (Z, Y, X) → (Y, X)"""
    fn = {"max": np.max, "min": np.min, "sum": np.sum,
          "mean": np.mean, "median": np.median}[method]
    result = fn(frames, axis=0)
    if method in ("mean", "median"):
        result = result.astype(frames.dtype)
    return result


def write_well_field(
    entries: list[ImageEntry],
    plate: PlateInfo,
    image_dir: Path,
    output_path: Path,
    method: Optional[str],
    compression: str,
    bigtiff: bool,
    correction: dict[int, CorrectionMaps],  # channel_id → CorrectionMaps
) -> None:
    """
    Write one OME-TIFF for a single well-field.

    entries    : all ImageEntry objects for this (row, col, field), sorted by
                 (timepoint, channel, z_offset_m).
    method     : projection method string, or None to keep all Z planes.
    correction : per-channel CorrectionMaps (empty dict = no correction).
    """
    do_project = method is not None

    # --- Derive dimension extents ------------------------------------------
    timepoints = sorted({e.timepoint   for e in entries})
    channels   = sorted({e.channel     for e in entries})
    z_offsets  = sorted({e.z_offset_m  for e in entries})

    size_t = len(timepoints)
    size_c = len(channels)
    size_z = 1 if do_project else len(z_offsets)

    # Lookup: (timepoint_1based, channel_1based, z_offset_m) → entry
    entry_map: dict[tuple, ImageEntry] = {
        (e.timepoint, e.channel, e.z_offset_m): e for e in entries
    }

    # --- Peek at first file for dtype + spatial size -----------------------
    first_file = _resolve_path(image_dir, entries[0].filename)
    with tifffile.TiffFile(str(first_file)) as tif:
        page = tif.pages[0]
        dtype = page.dtype
        img_size_y, img_size_x = page.shape[:2]

    size_x = plate.size_x if plate.size_x else img_size_x
    size_y = plate.size_y if plate.size_y else img_size_y

    px_x_um = plate.pixel_size_x_m * M_TO_UM if plate.pixel_size_x_m else 1.0
    px_y_um = plate.pixel_size_y_m * M_TO_UM if plate.pixel_size_y_m else 1.0

    if len(z_offsets) > 1:
        diffs   = [z_offsets[i+1] - z_offsets[i] for i in range(len(z_offsets)-1)]
        px_z_um = float(np.median(diffs)) * M_TO_UM
    else:
        px_z_um = 1.0

    # --- Channel metadata --------------------------------------------------
    ch_meta: list[dict] = []
    for ch in channels:
        sample = next(e for e in entries if e.channel == ch)
        ch_meta.append(dict(
            name=sample.channel_name,
            excitation_nm=sample.excitation_nm,
            emission_nm=sample.emission_nm,
            exposure_ms=sample.exposure_ms,
        ))

    # --- Build plane list in IFD order: T outer → C → Z inner (XYZCT) -----
    plane_meta: list[dict] = []
    if do_project:
        for ti, t in enumerate(timepoints):
            for ci, ch in enumerate(channels):
                z0  = z_offsets[0]
                ref = entry_map.get((t, ch, z0)) or next(
                    e for e in entries if e.timepoint == t and e.channel == ch
                )
                plane_meta.append(dict(
                    t=ti, c=ci, z=0,
                    exposure_ms=ref.exposure_ms,
                    delta_t_s=ref.time_offset_s,
                    pos_x_um=ref.pos_x_m * M_TO_UM,
                    pos_y_um=ref.pos_y_m * M_TO_UM,
                    pos_z_um=ref.z_offset_m * M_TO_UM,
                ))
    else:
        for ti, t in enumerate(timepoints):
            for ci, ch in enumerate(channels):
                for zi, z in enumerate(z_offsets):
                    ref = entry_map.get((t, ch, z))
                    if ref is None:
                        # Preserve stage metadata for missing Z planes so
                        # downstream alignment does not see undefined positions.
                        ref = next(
                            (e for e in entries if e.timepoint == t and e.channel == ch),
                            entries[0],
                        )
                        plane_meta.append(dict(
                            t=ti, c=ci, z=zi,
                            exposure_ms=ref.exposure_ms,
                            delta_t_s=ref.time_offset_s,
                            pos_x_um=ref.pos_x_m * M_TO_UM,
                            pos_y_um=ref.pos_y_m * M_TO_UM,
                            pos_z_um=ref.z_offset_m * M_TO_UM,
                        ))
                    else:
                        plane_meta.append(dict(
                            t=ti, c=ci, z=zi,
                            exposure_ms=ref.exposure_ms,
                            delta_t_s=ref.time_offset_s,
                            pos_x_um=ref.pos_x_m * M_TO_UM,
                            pos_y_um=ref.pos_y_m * M_TO_UM,
                            pos_z_um=ref.z_offset_m * M_TO_UM,
                        ))

    # --- Build OME-XML -----------------------------------------------------
    row, col, fov = entries[0].row, entries[0].col, entries[0].field
    label       = well_label(row, col)
    image_name  = f"{label}_F{fov:03d}"
    corr_applied = bool(correction)
    has_dark     = any(cm.darkfield is not None for cm in correction.values())

    ome_xml = build_ome_xml(
        image_name=image_name,
        size_x=size_x, size_y=size_y, size_z=size_z,
        size_c=size_c, size_t=size_t,
        dtype=dtype,
        px_x_um=px_x_um, px_y_um=px_y_um, px_z_um=px_z_um,
        channels=ch_meta,
        planes=plane_meta,
        plate=plate,
        row=row, col=col,
        correction_applied=corr_applied,
        correction_has_darkfield=has_dark,
    )
    ome_xml_bytes = ome_xml.encode("utf-8")

    out_bytes = size_t * size_c * size_z * size_y * size_x * dtype.itemsize
    use_bigtiff = bigtiff or out_bytes > 3 * 2**30

    write_opts = dict(
        photometric="minisblack",
        compression=compression,
        metadata=None,
        contiguous=False,
    )

    # --- Write frames: T outer → C → Z inner ------------------------------
    with tifffile.TiffWriter(str(output_path), bigtiff=use_bigtiff) as writer:
        first_page = True

        for ti, t in enumerate(timepoints):
            for ci, ch in enumerate(channels):
                corr_map = correction.get(ch)  # may be None

                if do_project:
                    z_stack = []
                    for z in z_offsets:
                        e = entry_map.get((t, ch, z))
                        if e is None:
                            print(f"\n  WARNING: missing plane T={t} C={ch} "
                                  f"Z={z:.2e} — filling with zeros")
                            z_stack.append(np.zeros((size_y, size_x), dtype=dtype))
                        else:
                            img_path = _resolve_path(image_dir, e.filename)
                            plane = _load_plane(img_path, dtype, size_y, size_x)
                            z_stack.append(plane)
                    stack = np.stack(z_stack, axis=0)
                    if corr_map is not None:
                        stack = apply_flat_field_stack(stack, corr_map)
                    frame = _project(stack, method)
                else:
                    # written inside the Z loop; frame set per Z
                    pass

                if do_project:
                    if first_page:
                        writer.write(frame, description=ome_xml_bytes, **write_opts)
                        first_page = False
                    else:
                        writer.write(frame, **write_opts)
                else:
                    for zi, z in enumerate(z_offsets):
                        e = entry_map.get((t, ch, z))
                        if e is None:
                            print(f"\n  WARNING: missing T={t} C={ch} "
                                  f"Z={z:.2e} — zeros")
                            frame = np.zeros((size_y, size_x), dtype=dtype)
                        else:
                            img_path = _resolve_path(image_dir, e.filename)
                            frame = _load_plane(img_path, dtype, size_y, size_x)
                            if corr_map is not None:
                                frame = apply_flat_field(frame, corr_map)

                        if first_page:
                            writer.write(frame, description=ome_xml_bytes, **write_opts)
                            first_page = False
                        else:
                            writer.write(frame, **write_opts)

                print(f"  {label} F{fov:03d}  T={ti+1}/{size_t}  C={ci+1}/{size_c}",
                      end="\r", flush=True)


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def _resolve_path(image_dir: Path, filename: str) -> Path:
    """
    PE index filenames are usually bare names in the same folder,
    but some exports nest them in a subfolder.  Try several candidates.
    """
    candidates = [
        image_dir / filename,
        image_dir / Path(filename).name,
        image_dir.parent / filename,
    ]
    for p in candidates:
        if p.exists():
            return p
    raise FileNotFoundError(
        f"Image file not found: {filename!r}\n"
        f"  Searched in: {image_dir}"
    )


def _load_plane(path: Path, dtype: np.dtype, size_y: int, size_x: int) -> np.ndarray:
    """Load a single 2-D plane from a TIFF file."""
    with tifffile.TiffFile(str(path)) as tif:
        arr = tif.pages[0].asarray()
    if arr.dtype != dtype:
        arr = arr.astype(dtype)
    if arr.shape != (size_y, size_x):
        out = np.zeros((size_y, size_x), dtype=dtype)
        sy = min(arr.shape[0], size_y)
        sx = min(arr.shape[1], size_x)
        out[:sy, :sx] = arr[:sy, :sx]
        return out
    return arr


# ---------------------------------------------------------------------------
# Top-level converter
# ---------------------------------------------------------------------------

def convert(
    index_path: Path,
    output_dir: Path,
    method: Optional[str],
    compression: str,
    bigtiff: bool,
    wells_filter: Optional[set[tuple[int, int]]],
    apply_correction: bool,
    write_companion: bool,
    prefix: str = "",
    row_range: Optional[list[int]] = None,
    col_range: Optional[list[int]] = None,
) -> None:
    print(f"Parsing {index_path.name} …")
    plate, entries = parse_index(index_path)
    image_dir = index_path.parent

    print(
        f"  Plate : {plate.name or plate.plate_id}\n"
        f"  Date  : {plate.measurement_date}\n"
        f"  Layout: {plate.n_rows} rows × {plate.n_cols} cols\n"
        f"  Images: {len(entries)} entries"
    )

    # --- Flat-field correction handling ------------------------------------
    correction: dict[int, CorrectionMaps] = {}

    if apply_correction:
        # Check if PE already baked the correction into the exported images
        already_applied = {e.ff_already_applied for e in entries}
        if True in already_applied:
            n_applied = sum(1 for e in entries if e.ff_already_applied)
            print(
                f"\n  WARNING: {n_applied}/{len(entries)} WellImage entries report "
                f"FlatfieldCorrectionApplied=Yes.\n"
                f"  Applying correction again will double-correct those images.\n"
                f"  Use --no-correction to skip, or verify your export settings."
            )

        print("\nSearching for illumination-correction maps …")
        corr_info = parse_correction_info(index_path)

        if not corr_info:
            # No file-based references found — try Harmony V6 polynomial profiles
            print("  No file-based correction references found.  "
                  "Checking for V6 polynomial profiles …")
            correction = parse_flatfield_from_maps(index_path)
            if not correction:
                print(
                    "  WARNING: No flat-field correction data found in the index XML.\n"
                    "  Vignetting will NOT be corrected.  If correction images exist,\n"
                    "  check that the index XML contains <Maps> or\n"
                    "  <IlluminationCorrectionImages> sections.\n"
                    "  Use --no-correction to suppress this warning."
                )
            else:
                print(
                    f"  Found polynomial flat-field profiles for "
                    f"{len(correction)} channel(s) (Harmony V6 format).  "
                    f"Darkfield: no (polynomial profiles do not include darkfield)."
                )
        else:
            # Attempt to load the correction images; abort on missing files
            try:
                correction = load_correction_maps(corr_info, image_dir)
            except FileNotFoundError as exc:
                sys.exit(
                    f"\nERROR: Correction file referenced in the index XML "
                    f"could not be found:\n  {exc}\n"
                    f"Use --no-correction to skip flat-field correction."
                )

            # Verify all channels have correction data
            ch_ids_in_data  = {e.channel for e in entries}
            ch_ids_with_corr = set(correction.keys())
            missing = ch_ids_in_data - ch_ids_with_corr
            if missing:
                ch_names = {e.channel: e.channel_name for e in entries}
                names = ", ".join(f"{ch_names.get(c, c)}" for c in sorted(missing))
                print(
                    f"  WARNING: No flat-field correction found for "
                    f"channel(s): {names}\n"
                    f"  Those channels will be written without correction."
                )

            has_dark = any(cm.darkfield is not None for cm in correction.values())
            print(
                f"  Found correction maps for "
                f"{len(correction)}/{len(ch_ids_in_data)} channel(s).  "
                f"Darkfield: {'yes' if has_dark else 'no'}."
            )
    else:
        print("\n  Flat-field correction disabled (--no-correction).")

    # --- Group and filter wells -------------------------------------------
    groups = group_by_well_field(entries)
    if wells_filter:
        groups = {k: v for k, v in groups.items() if (k[0], k[1]) in wells_filter}
        if not groups:
            sys.exit("No matching wells found after applying --wells filter.")

    if row_range is not None or col_range is not None:
        allowed_rows = set(row_range) if row_range is not None else None
        allowed_cols = set(col_range) if col_range is not None else None
        groups = {
            k: v for k, v in groups.items()
            if (allowed_rows is None or k[0] in allowed_rows)
            and (allowed_cols is None or k[1] in allowed_cols)
        }
        range_desc = "".join([
            f" rows {chr(ord('A') + min(allowed_rows) - 1)}"
            f"-{chr(ord('A') + max(allowed_rows) - 1)}" if allowed_rows else "",
            f" cols {min(allowed_cols)}-{max(allowed_cols)}" if allowed_cols else "",
        ]).strip()
        if not groups:
            sys.exit(
                f"No wells found after applying range filter ({range_desc})."
            )

    output_dir.mkdir(parents=True, exist_ok=True)
    if correction:
        save_correction_maps(correction, output_dir, prefix=prefix)
        print(f"Saved correction maps (ffp/dfp) to: {output_dir}")

    proj_label = f"{method}-projection" if method else "all-Z-planes"
    corr_label = "with FF correction" if correction else "no FF correction"
    print(
        f"\nOutput dir  : {output_dir}\n"
        f"Mode        : {proj_label}, {corr_label}\n"
        f"Compression : {compression}\n"
        f"Well-fields : {len(groups)}\n"
    )

    for idx, ((row, col, fov), field_entries) in enumerate(
        sorted(groups.items()), start=1
    ):
        label    = well_label(row, col)
        base_name = f"{label}_F{fov:03d}"
        if method:
            base_name = f"{label}_F{fov:03d}_{method}proj"
        out_name = (f"{prefix}_{base_name}.ome.tif" if prefix
                    else f"{base_name}.ome.tif")
        output_path = output_dir / out_name

        n_t = len({e.timepoint   for e in field_entries})
        n_c = len({e.channel     for e in field_entries})
        n_z = len({e.z_offset_m  for e in field_entries})
        print(
            f"[{idx:4d}/{len(groups)}] {label} F{fov:03d}  "
            f"T={n_t} C={n_c} Z={n_z} → {out_name}"
        )

        try:
            write_well_field(
                entries=field_entries,
                plate=plate,
                image_dir=image_dir,
                output_path=output_path,
                method=method,
                compression=compression,
                bigtiff=bigtiff,
                correction=correction,
            )
        except FileNotFoundError as exc:
            print(f"\n  ERROR: {exc} — skipping.")
            continue

        print()  # newline after \r progress

    if write_companion:
        companion_name = f"{prefix}.companion.ome" if prefix else "plate.companion.ome"
        make_companion(
            output_dir=output_dir,
            companion_path=output_dir / companion_name,
            plate=plate,
        )

    print("Done.")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    ap = argparse.ArgumentParser(
        description=(
            "Convert PerkinElmer Columbus/Harmony plate data to OME-TIFF. "
            "Outputs one OME-TIFF per well-field with full OME metadata. "
            "Applies flat-field correction and optional Z-projection."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument(
        "index",
        help="Path to the PerkinElmer Index.idx.xml (or Index.xml) file.",
    )
    ap.add_argument(
        "-o", "--output-dir",
        default=None,
        help=(
            "Directory for output OME-TIFFs.  "
            "Defaults to a subfolder 'ome_tif' next to the index file."
        ),
    )
    ap.add_argument(
        "--no-project",
        action="store_true",
        help="Keep all Z planes instead of projecting.",
    )
    ap.add_argument(
        "-m", "--method",
        choices=METHODS,
        default="max",
        help="Z-projection method (ignored when --no-project is set).",
    )
    ap.add_argument(
        "--no-correction",
        action="store_true",
        help=(
            "Skip flat-field correction even if correction images are "
            "referenced in the index XML."
        ),
    )
    ap.add_argument(
        "-c", "--compression",
        default="zlib",
        help="TIFF compression: zlib, lzma, zstd, or 'none'.",
    )
    ap.add_argument(
        "--bigtiff",
        action="store_true",
        help="Force BigTIFF output (auto-enabled when a file exceeds 3 GB).",
    )
    ap.add_argument(
        "--wells",
        default=None,
        help=(
            "Comma-separated list of well labels to process, e.g. A01,A02,B03. "
            "Default: all wells."
        ),
    )
    ap.add_argument(
        "--row-range",
        default=None,
        metavar="START[-END]",
        help=(
            "Row letter range to process, inclusive.  "
            "Single row: 'B'.  Range: 'B-D' processes rows B, C, D.  "
            "Can be combined with --col-range.  Default: all rows."
        ),
    )
    ap.add_argument(
        "--col-range",
        default=None,
        metavar="START[-END]",
        help=(
            "Column number range to process, inclusive.  "
            "Single column: '4'.  Range: '4-6' processes columns 4, 5, 6.  "
            "Can be combined with --row-range.  Default: all columns."
        ),
    )
    ap.add_argument(
        "--no-companion",
        action="store_true",
        help=(
            "Do not write a plate.companion.ome file after conversion."
        ),
    )
    ap.add_argument(
        "--prefix",
        default="",
        help=(
            "Optional prefix prepended to every output filename, separated by '_'. "
            "E.g. --prefix mysample produces mysample_A01_F001_maxproj.ome.tif."
        ),
    )
    return ap.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    index_path = Path(args.index)
    if not index_path.exists():
        sys.exit(f"Error: {index_path} does not exist.")

    output_dir = (
        Path(args.output_dir)
        if args.output_dir
        else index_path.parent / "ome_tif"
    )

    method: Optional[str] = None if args.no_project else args.method
    compression = None if args.compression.lower() == "none" else args.compression

    wells_filter: Optional[set[tuple[int, int]]] = None
    if args.wells:
        wells_filter = set()
        for label in args.wells.split(","):
            try:
                wells_filter.add(label_to_rowcol(label.strip()))
            except ValueError as exc:
                sys.exit(f"Error: {exc}")

    row_range: Optional[list[int]] = None
    if args.row_range:
        try:
            row_range = parse_row_range(args.row_range)
        except ValueError as exc:
            sys.exit(f"Error in --row-range: {exc}")

    col_range: Optional[list[int]] = None
    if args.col_range:
        try:
            col_range = parse_col_range(args.col_range)
        except ValueError as exc:
            sys.exit(f"Error in --col-range: {exc}")

    convert(
        index_path=index_path,
        output_dir=output_dir,
        method=method,
        compression=compression,
        bigtiff=args.bigtiff,
        wells_filter=wells_filter,
        apply_correction=not args.no_correction,
        write_companion=not args.no_companion,
        prefix=args.prefix,
        row_range=row_range,
        col_range=col_range,
    )


if __name__ == "__main__":
    main()
