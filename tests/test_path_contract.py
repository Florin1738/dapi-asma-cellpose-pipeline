from pathlib import Path

import pytest

from dapi_norm.path_contract import to_posix_path


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        ("Plate_1\\Run_A\\XY01", "Plate_1/Run_A/XY01"),
        (Path("Plate_1") / "Run_A" / "XY01", "Plate_1/Run_A/XY01"),
    ],
)
def test_to_posix_path_preserves_portable_output_contract(value: str | Path, expected: str) -> None:
    assert to_posix_path(value) == expected
