#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) 2025 Wellcome Sanger Institute

import pandas as pd
from starfish.core.codebook.codebook import Codebook
from starfish.types import Axes, Features
import numpy as np
import re


def to_starfish_codebook(
    codebook: pd.DataFrame, target_col: str, code_col: str, is_merfish: bool = False
) -> Codebook:
    mappings = []
    for _, row in codebook.iterrows():
        mapping = {}
        mapping[Features.TARGET] = row[target_col]
        codeward = []
        for r, c in enumerate(str(row[code_col])):
            if is_merfish:
                codeward.append(
                    {Axes.ROUND.value: r, Axes.CH.value: 0, Features.CODE_VALUE: c}
                )
            else:
                codeward.append(
                    {
                        Axes.ROUND.value: r,
                        Axes.CH.value: int(c) - 1,
                        Features.CODE_VALUE: 1,
                    }
                )
        mapping[Features.CODEWORD] = codeward
        mappings.append(mapping)
    if is_merfish:
        return Codebook.from_code_array(
            mappings, n_round=r + 1, n_channel=1
        )  # hardcoded to 4 channels since Cartana only uses 4 channels
    else:
        return Codebook.from_code_array(
            mappings, n_round=r + 1, n_channel=4
        )  # hardcoded to 4 channels since Cartana only uses 4 channels


def filter_columns_by_regex(df: pd.DataFrame, pattern: str) -> pd.DataFrame:
    regex = re.compile(pattern)
    return df.loc[:, df.columns.to_series().apply(lambda x: bool(regex.match(x)))]


def load_codebook(codebook_path: str, codebook_code_col: str) -> pd.DataFrame:
    """
    Read codebook from a file.

    Args:
        codebook_path (str): A file path to the codebook.

    Returns:
        pd.DataFrame: A pandas DataFrame containing the codebook.
    """
    if codebook_path.endswith(".xlsx"):
        codebook = pd.read_excel(codebook_path, dtype={codebook_code_col: str})
    elif codebook_path.endswith(".csv"):
        codebook = pd.read_csv(codebook_path, dtype={codebook_code_col: str})
    elif codebook_path.endswith(".tsv"):
        codebook = pd.read_csv(codebook_path, sep="\t", dtype={codebook_code_col: str})
    else:
        raise ValueError("Codebook file must be in .xlsx, .csv or .tsv format")
    return codebook.dropna()
    # codebook.to_json(f"{Path(codebook_p).stem}.json")


def filter_columns_by_regex(df: pd.DataFrame, pattern: str) -> pd.DataFrame:
    regex = re.compile(pattern)
    return df.loc[:, df.columns.to_series().apply(lambda x: bool(regex.match(x)))]


def qc_codebook(
    codebook, code_col: str, coding_col_prefix: str = "cycle\d_channel\d_*"
):
    df = filter_columns_by_regex(codebook, coding_col_prefix).astype(int)
    _, n_code = df.shape
    if "channel" not in coding_col_prefix:
        for i, code in enumerate(codebook[code_col]):
            assert code == "".join(
                df.iloc[i, :].values.astype(str)
            ), "ISS-like Codebook is not consistent with the coding columns"
    else:
        for i, code in enumerate(codebook[code_col]):
            current_code = np.zeros(n_code, dtype=int)
            for j, char in enumerate(code):
                current_code[j * n_code // len(code) + int(char)] = 1
            assert "".join(df.iloc[i, :].values.astype(str)) == "".join(
                current_code.astype(str)
            ), "Merifish-like Codebook is not consistent with the coding columns"
