#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) 2024 Wellcome Sanger Institute

import fire
from pandas import DataFrame
from shapely import from_wkt
from pathlib import Path


def load_wkts_as_table(wkt_files: list):
    names, ys, xs = [], [], []
    for wkt_file in wkt_files:
        current_stem = Path(wkt_file).stem
        with open(wkt_file, "r") as f:
            for g in from_wkt(f.read()).geoms:
                names.append(current_stem)
                ys.append(g.y)
                xs.append(g.x)
    return {"feature_name": names, "Y": ys, "X": xs}


def main(*wkt_files, output_name: str):
    tab = DataFrame(load_wkts_as_table(wkt_files))
    tab["instance_id"] = tab.index + 1
    tab.set_index("instance_id", inplace=True)
    tab["x_int"] = tab["X"].astype(int)
    tab["y_int"] = tab["Y"].astype(int)
    tab.to_csv(output_name)


if __name__ == "__main__":
    options = {"run": main, "version": "0.0.2"}
    fire.Fire(options)
