#!/usr/bin/env python3

from __future__ import annotations
from pathlib import Path
import json
import importlib.util

import pytest

TRAIN_PATH = Path("../bin/nichecompass_train_sample_integration.py")


# Skip the suite if the script isn't where we expect it.
pytestmark = pytest.mark.skipif(
    not TRAIN_PATH.exists(),
    reason=f"Training script not found at {TRAIN_PATH}",
)


def _load_train_module():
    """Dynamically import the training script from usr/bin by file path."""
    spec = importlib.util.spec_from_file_location(
        "nichecompass_train_sample_integration", TRAIN_PATH.as_posix()
    )
    assert spec and spec.loader, "Failed to create import spec for training script"
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[attr-defined]
    return module


@pytest.fixture(scope="session")
def train():
    """Provide the imported training module to tests."""
    return _load_train_module()


def _parse_with_config_aware(train, argv: list[str]):
    """Mimic main(): pre-parse --config, load JSON, then full parse + merge (CLI wins)."""
    parser, pre = train.build_parser()
    pre_args, remaining = pre.parse_known_args(argv)
    cfg = {}
    if getattr(pre_args, "config", None):
        cfg = train.load_config_json(pre_args.config)
    args = parser.parse_args(remaining)
    return train.merge_config_and_args(args, cfg)


def test_minimal_config_only(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, train):
    # Put batches in config only, no CLI args.
    a = (tmp_path / "a.h5ad").as_posix()
    b = (tmp_path / "b.h5ad").as_posix()
    cfg_path = tmp_path / "config.json"
    cfg_path.write_text(json.dumps({"batches": [a, b]}), encoding="utf-8")

    # Ensure cwd-based default for outdir is predictable
    monkeypatch.chdir(tmp_path)

    rp = _parse_with_config_aware(train, ["--config", cfg_path.as_posix()])

    assert rp.batches == [Path(a), Path(b)]
    assert rp.outdir == tmp_path                     # default: Path.cwd()
    assert rp.prefix == "nichecompass"               # dataclass default
    assert rp.species == "human"                     # dataclass default
    assert rp.nichecompass_data_tag == "0.3.0"       # dataclass default
    assert rp.cat_covariates_keys == [rp.sample_key] # default to [sample_key]


def test_cli_overrides_config(tmp_path: Path, train):
    a = (tmp_path / "a.h5ad").as_posix()
    b = (tmp_path / "b.h5ad").as_posix()
    cfg_path = tmp_path / "config.json"
    cfg = {
        "batches": [a, b],
        "prefix": "cfg_prefix",
        "n_neighbors": 4,
        "use_cuda_if_available": True,
    }
    cfg_path.write_text(json.dumps(cfg), encoding="utf-8")

    rp = _parse_with_config_aware(
        train,
        [
            "--config", cfg_path.as_posix(),
            "--prefix", "cli_prefix",
            "--n_neighbors", "8",
            "--use_cuda_if_available", "false"
        ]
    )

    assert rp.prefix == "cli_prefix"
    assert rp.n_neighbors == 8
    assert rp.use_cuda_if_available is False
    assert rp.batches == [Path(a), Path(b)]


def test_missing_batches_error(tmp_path: Path, train):
    cfg_path = tmp_path / "config.json"
    cfg_path.write_text(json.dumps({"prefix": "x"}), encoding="utf-8")

    parser, pre = train.build_parser()
    pre_args, remaining = pre.parse_known_args(["--config", cfg_path.as_posix()])
    cfg = train.load_config_json(pre_args.config)
    args = parser.parse_args(remaining)

    with pytest.raises(ValueError, match="You must provide --batches"):
        train.merge_config_and_args(args, cfg)


def test_batches_space_separated_cli_only(tmp_path: Path, train):
    a = tmp_path / "a.h5ad"
    b = tmp_path / "b.h5ad"

    rp = _parse_with_config_aware(train, ["--batches", a.as_posix(), b.as_posix()])

    assert rp.batches == [a, b]
    assert isinstance(rp.batches[0], Path)


def test_outdir_defaults_to_cwd(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, train):
    a = tmp_path / "a.h5ad"
    monkeypatch.chdir(tmp_path)

    rp = _parse_with_config_aware(train, ["--batches", a.as_posix()])
    assert rp.outdir == tmp_path


def test_bool_list_conversion_and_keys_default(tmp_path: Path, train):
    a = tmp_path / "a.h5ad"

    rp = _parse_with_config_aware(
        train,
        [
            "--batches", a.as_posix(),
            "--sample_key", "donor",
            "--cat_covariates_no_edges", "true", "false",
            "--cat_covariates_embeds_nums", "2"
        ]
    )

    assert rp.cat_covariates_keys == ["donor"]   # default to [sample_key]
    assert rp.cat_covariates_no_edges == [True, False]  # string -> bool list
    assert rp.cat_covariates_embeds_nums == [2]         # string int -> int list


def test_use_cuda_if_available_string_to_bool(tmp_path: Path, train):
    a = tmp_path / "a.h5ad"

    rp = _parse_with_config_aware(train, ["--batches", a.as_posix(), "--use_cuda_if_available", "true"])
    assert rp.use_cuda_if_available is True

    rp2 = _parse_with_config_aware(train, ["--batches", a.as_posix(), "--use_cuda_if_available", "false"])
    assert rp2.use_cuda_if_available is False


def test_unknown_config_key_rejected(tmp_path: Path, train):
    a = (tmp_path / "a.h5ad").as_posix()
    cfg_path = tmp_path / "config.json"
    cfg_path.write_text(json.dumps({"batches": [a], "unknown_key": 123}), encoding="utf-8")

    parser, pre = train.build_parser()
    pre_args, remaining = pre.parse_known_args(["--config", cfg_path.as_posix()])
    cfg = train.load_config_json(pre_args.config)
    args = parser.parse_args(remaining)

    with pytest.raises(ValueError, match="Unknown keys in --config"):
        train.merge_config_and_args(args, cfg)


def test_derived_config_key_rejected(tmp_path: Path, train):
    a = (tmp_path / "a.h5ad").as_posix()
    cfg_path = tmp_path / "config.json"
    cfg_path.write_text(json.dumps({"batches": [a], "timestamp": "20250101_000000"}), encoding="utf-8")

    parser, pre = train.build_parser()
    pre_args, remaining = pre.parse_known_args(["--config", cfg_path.as_posix()])
    cfg = train.load_config_json(pre_args.config)
    args = parser.parse_args(remaining)

    with pytest.raises(ValueError, match="Unknown keys in --config"):
        train.merge_config_and_args(args, cfg)
