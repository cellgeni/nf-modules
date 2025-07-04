#!/usr/bin/env python3

import os
import glob
import pathlib
import logging
import argparse
import pandas as pd

logging.basicConfig(level="INFO", format="[%(asctime)s][%(levelname)s] %(message)s")

ACAPELLA_DEFAULTS = {
    # Supported values: max, sharpest, mean, none
    "acapella_zprojection": "max",
    # maximally allopwed gap in pixels
    "acapella_field_gap": "4000",
    # Well selection: ALL | A1,A2..E10
    "acapella_wells": "ALL",
    # Field selection: ALL | 1,2:4,7..9
    "acapella_fields": "ALL",
    # Plane selection: ALL | 1,5:2:16
    "acapella_planes": "ALL",
    # Timepoint selection:  ALL | 1,5:2:16
    "acapella_timepoints": "ALL",
    # Channel selection:  ALL | 1,2,4-6
    "acapella_channels": "ALL",
    # Note. ':' denotes exclusive range, '..' or  '-' inclusive range
}


def main(args):
    """Read manifest, clean it up, and convert it to JSON array"""
    # consider everything a string except 'Date' field
    # types = defaultdict(str, Date=datetime)
    # read manifest
    logging.info(f"Reading {args.manifest}")
    df = pd.read_excel(args.manifest, dtype=str)
    logging.info(f"File has {len(df)} rows")
    # remove empty rows
    logging.info("Removing empty rows...")
    df = df.dropna(how="all", axis="index")
    logging.info(f"{len(df)} rows left after removing empty")
    # cleanup cells and fill-in missing columns
    logging.info(f"Pre-processing...")
    df = preprocess(df)
    # write output json
    logging.info(f"Storing DataFrame as {args.out}")
    df.to_json(args.out, orient="records")


def preprocess(df):
    """Cleanup manifest records before processing"""
    # iterat using index so we can manipulate the dataframe
    logging.info(f"Update DataFrame...")
    for idx in df.index:
        logging.info(f"Preprocessing row {idx}")
        row = df.loc[idx]
        # clean any lead/tail spaces in cells
        for c in df.columns:
            if isinstance(row[c], str):
                row[c] = row[c].strip()

        # select id to use when searching (SlideID or AutomatedPlateID)
        if pd.isna(row.Automated_PlateID) or row.Automated_PlateID == "":
            _id = row.SlideID.strip()
            logging.info(f"Using SlideID {_id}")
        else:
            _id = row.Automated_PlateID.strip()
            logging.info(f"Using Automated_PlateID {_id}")
        df.loc[idx, "_id"] = _id

        # fix backslashes in paths
        Export_location = row.Export_location.replace("\\", "/")
        df.loc[idx, "Export_location"] = Export_location
        Team_dir = row.Team_dir.replace("\\", "/")
        df.loc[idx, "Team_dir"] = Team_dir

        # build search path for this measurement
        search_path = os.path.join(
            Team_dir,
            Export_location.replace("\\", "/"),
            f"{_id}__*Measurement {row.Measurement}",
        )

        # loosely search for both Index.xml (Phenix) and Index.idx.xml (Operetta)
        index_file = glob.glob(os.path.join(search_path, "Images", "Index*.xml"))
        if len(index_file) == 0:
            raise SystemExit(f"Index file not found in {search_path}.")
        df.loc[idx, "index_file"] = index_file[0]
        logging.info(f"Found index file {index_file}")

        # fill missing values on default columns ##############################
        if pd.isnull(row.SectionN) or row.SectionN == "":
            logging.info(f"'SectionN' column not provided, set value to 1")
            SectionN = 1
            df.loc[idx, "SectionN"] = SectionN
        else:
            SectionN = row.SectionN

        # decide to use MIP (max intesity projection) or full z-planes ########
        logging.info(f"Stitching_Z is '{row.Stitching_Z}'")
        if row.Stitching_Z not in ["max", "none"]:
            # this needs changing to support cherry pick planes
            Stitching_Z = ACAPELLA_DEFAULTS["acapella_zprojection"]
        else:
            Stitching_Z = row.Stitching_Z
        df.loc[idx, "acapella_zprojection"] = Stitching_Z

        # for naming the subfolder inside the stitched sample folder
        # this will allow to have mip, all_planes, 2-3_panes as subfolders
        if Stitching_Z == "max":
            zprojection_suffix = "mip"
        elif Stitching_Z == "none":
            zprojection_suffix = "all_planes"
        else:
            zprojection_suffix = (
                Stitching_Z.lower().replace(",", "_").replace("..", "-") + "_planes"
            )
        df.loc[idx, "zprojection_suffix"] = zprojection_suffix

        # geneare unique identifier for this
        df.loc[idx, "uuid"] = f"{_id}_{row.Measurement}_{Stitching_Z}"

        # fill/create acapella columns ########################################
        row = df.loc[idx]
        for acappela_field in ACAPELLA_DEFAULTS.keys():
            if (
                (acappela_field not in row)
                or pd.isnull(row[acappela_field])
                or row[acappela_field] == ""
            ):
                df.loc[idx, acappela_field] = ACAPELLA_DEFAULTS[acappela_field]
            logging.info(f"Using {acappela_field} = {df.loc[idx, acappela_field]}")

        # build output_dir path if mising #####################################
        if (
            ("output_dir" not in row)
            or pd.isnull(row["output_dir"])
            or row["output_dir"] == ""
        ):
            output_dir = os.path.join(
                row.Team_dir,
                "0HarmonyStitched",
                row.Project,
                pathlib.Path(index_file[0]).parents[1].stem,
                zprojection_suffix,
            )
            df.loc[idx, "output_dir"] = output_dir
        else:
            output_dir = row.output_dir
        logging.info(f"Using output_dir: {output_dir}")

        # fill/create omero columns ###########################################
        omero_columns = [
            "omero_group",
            "omero_username",
            "omero_project",
            "omero_dataset",
        ]
        for oc in omero_columns:
            if oc not in row:
                row[oc] = ""
                df.loc[idx, oc] = ""
        if pd.isnull(row["omero_group"]) or row["omero_group"] == "":
            df.loc[idx, "omero_group"] = row["Project"]
            logging.info(
                f"omero_group not provided using column 'Project'  with value '{row['Project']}'"
            )
        if pd.isnull(row["omero_project"]) or row["omero_project"] == "":
            df.loc[idx, "omero_project"] = row["Tissue_1"]
            logging.info(
                f"omero_project not provided using column 'Tissue_1' with value '{row['Tissue_1']}'"
            )
        if pd.isnull(row["omero_username"]) or row["omero_username"] == "":
            df.loc[idx, "omero_username"] = row["Researcher"]
            logging.info(
                f"omero_username not provided using column 'Researcher' with value '{row['Researcher']}'"
            )
        if pd.isnull(row["omero_dataset"]) or row["omero_dataset"] == "":
            # prepare yourself, this is pretty ugly
            try:
                sample_name = row[f"Sample_{SectionN}"]
            except:
                logging.warning(
                    f"SectionN value '{SectionN}' doesn match a column name with pattern 'Sample_N'"
                )
                sample_name = row[f"Sample_1"]
            sample_name = (
                ""
                if pd.isna(sample_name)
                else "-".join(str(sample_name).split("-")[:2])
            )
            sample_name = sample_name.replace("/", ".")
            targets = sorted(
                [
                    t
                    for t in row.keys()
                    if t.lower().startswith("target") and row[t] and not pd.isna(row[t])
                ]
            )
            dataset_list = [sample_name] + [
                row[t] for t in targets if not pd.isna(row[t]) and row[t]
            ]
            dataset_list = map(str, dataset_list)
            df.loc[idx, "omero_dataset"] = "_".join(
                [s for s in dataset_list if len(s) != 0]
            )
            logging.info(
                f"omero_dataset not provided using '{df.loc[idx, 'omero_dataset']}'"
            )

    return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=str, required=True)
    parser.add_argument("--out", type=str, default="manifest.json")
    args = parser.parse_args()
    main(args)
