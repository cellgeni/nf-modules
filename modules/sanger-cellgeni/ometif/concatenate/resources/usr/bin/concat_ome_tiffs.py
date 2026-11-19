#!/usr/bin/env python3
import re

import click
import numpy as np
import tifffile

AXIS_CHOICES = ["C", "Z", "T"]


def _read_series(path):
    with tifffile.TiffFile(path) as tif:
        series = tif.series[0]
        array = series.asarray()
        axes = series.axes
        ome_metadata = tif.ome_metadata
    return array, axes, ome_metadata


def _physical_size_metadata(ome_xml):
    if not ome_xml:
        return {}
    metadata = {}
    for dim in ("X", "Y", "Z"):
        size = re.search(rf'PhysicalSize{dim}="([^"]+)"', ome_xml)
        unit = re.search(rf'PhysicalSize{dim}Unit="([^"]+)"', ome_xml)
        if size:
            metadata[f"PhysicalSize{dim}"] = float(size.group(1))
        if unit:
            metadata[f"PhysicalSize{dim}Unit"] = unit.group(1)
    return metadata


def _channel_names(ome_xml, n_channels):
    """Original per-channel Name attributes, or "C0", "C1", ... if absent/incomplete."""
    names = re.findall(r'<Channel\b[^>]*?\bName="([^"]*)"', ome_xml) if ome_xml else []
    if len(names) == n_channels and all(names):
        return names
    return [f"C{i}" for i in range(n_channels)]


@click.command()
@click.version_option(version="0.1.0")
@click.argument("input_files", nargs=-1, required=True, type=click.Path(exists=True))
@click.option("-o", "--output", required=True, type=click.Path(), help="Output OME-TIFF path")
@click.option(
    "--axis",
    type=click.Choice(AXIS_CHOICES),
    default="C",
    show_default=True,
    help="Axis to concatenate the input images along, growing a hyperstack",
)
def concat_ome_tiffs(input_files, output, axis):
    """Concatenate multiple OME-TIFFs into a single OME-TIFF hyperstack."""
    if len(input_files) < 2:
        raise click.UsageError("At least two input files are required")

    click.echo(f"Loading {len(input_files)} file(s)...")
    loaded = [_read_series(f) for f in input_files]
    arrays, axes_list, ome_list = zip(*loaded)

    ref_axes = axes_list[0]
    if any(a != ref_axes for a in axes_list):
        raise click.UsageError(f"All input images must share the same axes order, got: {set(axes_list)}")

    if axis not in ref_axes:
        # Every input is missing this axis: insert a new length-1 axis right
        # before Y (or at the front if there is no Y axis) so it can grow.
        insert_at = ref_axes.index("Y") if "Y" in ref_axes else 0
        ref_axes = ref_axes[:insert_at] + axis + ref_axes[insert_at:]
        arrays = [np.expand_dims(a, insert_at) for a in arrays]

    concat_axis = ref_axes.index(axis)

    click.echo(f"Concatenating along axis '{axis}' (array axis {concat_axis})...")
    stacked = np.concatenate(arrays, axis=concat_axis)

    metadata = {"axes": ref_axes}
    metadata.update(_physical_size_metadata(ome_list[0]))

    if axis == "C":
        # Stamp each channel with its 1-based position in the input series
        # (its "cycle"), so the series order stays traceable in the output
        # channel names even after everything has been merged into one file.
        channel_names = [
            f"cycle_{cycle}_{name}"
            for cycle, (arr, ome_xml) in enumerate(zip(arrays, ome_list), start=1)
            for name in _channel_names(ome_xml, arr.shape[concat_axis])
        ]
        metadata["Channel"] = {"Name": channel_names}

    click.echo(f"Writing output to {output}")
    tifffile.imwrite(output, stacked, ome=True, metadata=metadata, bigtiff=True)

    click.echo("Done!")


if __name__ == "__main__":
    concat_ome_tiffs()
