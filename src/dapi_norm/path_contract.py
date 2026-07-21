"""Stable path text for output identifiers and browser-facing links.

Filesystem operations should keep using :class:`pathlib.Path`. Values that
become IDs, manifest fields, or HTML URLs use forward slashes on every OS so
that output contracts remain portable between Windows and Unix-like systems.
"""

from pathlib import Path


def to_posix_path(value: str | Path) -> str:
    """Return ``value`` with forward slashes, without changing filesystem I/O."""

    return str(value).replace("\\", "/")
