#!/usr/bin/env python3

import os
import json
import argparse
import pandas as pd
import json
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(levelname)s: %(message)s'
)

def init_parser() -> argparse.ArgumentParser:
    """
    Initialise argument parser for the script
    """
    parser = argparse.ArgumentParser(
        description="Aggregates metadata for duplicated attributes and saves in .csv and .json format"
    )
    parser.add_argument(
        "--input",
        metavar="<file>",
        type=str,
        help="Specify a path to the .csv file with iRODS metadata",
    )
    parser.add_argument(
        "--dup-sep",
        metavar="<str>",
        type=str,
        default=",",
        help="Separator for duplicated metadata attributes (default: ',')",
    )
    parser.add_argument(
        "--index_name",
        metavar="<str>",
        type=str,
        default="id",
        help="Name of the index column for the output .csv file (default: 'id')",
    )
    parser.add_argument(
        "--id",
        metavar="<str>",
        type=str,
        default=None,
        help="Identifier to use for the index column in the output .csv file (default: None)",
    )
    return parser


def main():
    """
    Main function of the script
    """
    
    # parse script arguments
    parser = init_parser()
    args = parser.parse_args()

    # Check if file exists and not empty
    if not os.path.isfile(args.input) or os.path.getsize(args.input) == 0:
        logging.warning(f"Input file '{args.input}' is empty or does not exist - creating empty output files")
        
        # Create a DataFrame with only index column
        metadata = pd.DataFrame(index=[args.id] if args.id else None)
        metadata.index.name = args.index_name
    else:
        # read input metadata file
        irods_metadata = pd.read_csv(args.input, header=None, names=["attribute", "value", "unit"])

        # aggregate duplicated metadata attributes
        metadata = pd.pivot_table(irods_metadata, values="value", columns="attribute", aggfunc=lambda x: args.dup_sep.join(x))
        metadata.index = [args.id] if args.id else metadata.index
        metadata.index.name = args.index_name

    # save aggregated metadata
    metadata.to_csv("metadata.csv", index=True if args.id else False)
    metadata.to_json("metadata.json", orient="index" if args.id else "records", indent=4)

if __name__ == "__main__":
    main()