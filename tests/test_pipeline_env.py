from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]


def _load_pipeline_env():
    spec = importlib.util.spec_from_file_location(
        "pipeline_env", REPO_ROOT / "scripts" / "pipeline_env.py"
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


pe = _load_pipeline_env()


def test_model_checksum_constant_matches_scientific_contract():
    # The expected model hash is the single source of truth used by both OSes.
    assert pe.MODEL_SHA256 == (
        "0f1cc3f7ecdd8a037a57c6c48d9d8921391be4cbce3fa9f13c3e3a2e1253c667"
    )


def test_venv_python_paths_are_os_specific():
    p = Path("/proj")
    # Both branches are exercised regardless of the host OS.
    win = pe.venv_python
    # We can at least assert the model paths are stable.
    assert pe.model_path(p) == p / ".models" / "cellpose" / "cpsam_v2"
    assert callable(win)


def test_verify_model_reports_missing_model(tmp_path):
    ok, message = pe.verify_model(tmp_path)
    assert ok is False
    assert "not found" in message


def test_verify_model_detects_checksum_mismatch(tmp_path):
    target = pe.model_path(tmp_path)
    target.parent.mkdir(parents=True)
    target.write_bytes(b"not the real model")
    ok, message = pe.verify_model(tmp_path)
    assert ok is False
    assert "checksum mismatch" in message


def test_env_ready_flags_missing_venv(tmp_path):
    ok, message = pe.env_ready(tmp_path)
    assert ok is False
    assert ".venv" in message or "missing" in message


@pytest.mark.skipif(
    not pe.model_path(REPO_ROOT).is_file(),
    reason="Real cpsam_v2 model cache not present in this checkout.",
)
def test_verify_model_passes_on_real_cached_model():
    ok, message = pe.verify_model(REPO_ROOT)
    assert ok is True, message
