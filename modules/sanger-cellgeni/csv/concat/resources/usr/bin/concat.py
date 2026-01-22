#!/usr/bin/env python3

import argparse
import pandas as pd
import json

def init_parser() -> argparse.ArgumentParser:
    """
    Initialise argument parser for the script
    """
    parser = argparse.ArgumentParser(
        description="Concatenates .csv files and save the result in .csv and .json format"
    )
    parser.add_argument(
        "--input",
        metavar="<file>",
        nargs="+",
        type=str,
        help="Specify a path to the .csv files to concatenate",
    )
    parser.add_argument(
        "--axis",
        metavar="<str>",
        type=str,
        default="index",
        help="Axis to concatenate along ('columns' or 'index'; default: 'index')",
    )
    parser.add_argument(
        "--join",
        metavar="<str>",
        type=str,
        default="outer",
        help="How to handle indexes on other axis (or axes). Options are 'inner' and 'outer' (default: 'outer')",
    )
    parser.add_argument(
        "--prefix",
        metavar="<str>",
        type=str,
        default="output",
        help="Prefix for the output files (default: 'output')",
    )
    return parser


def main():
    """
    Main function of the script
    """
    # parse script arguments
    parser = init_parser()
    args = parser.parse_args()

    # read input files
    csv_files = [pd.read_csv(f) for f in args.input]

    # concatenate .csv files
    result = pd.concat(csv_files, axis=args.axis, join=args.join)

    # sort values and columns
    result = result.reindex(sorted(result.columns), axis=1)
    result = result.sort_values(by=result.columns.tolist())

    # save result
    result.to_csv(f"{args.prefix}.csv", index=False)
    result.to_json(f"{args.prefix}.json", orient="records", lines=True, indent=4)

if __name__ == "__main__":
    main()