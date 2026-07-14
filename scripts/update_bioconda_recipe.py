#!/usr/bin/env python3
"""Update the official Bioconda RepliDec recipe for a PyPI release."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


VERSION_PATTERN = re.compile(r'({%\s*set\s+version\s*=\s*")[^"]+("\s*%})')
SHA256_PATTERN = re.compile(r"(^\s*sha256:\s*)[0-9a-fA-F]+\s*$", re.MULTILINE)
BUILD_PATTERN = re.compile(r"(^\s*number:\s*)[0-9]+\s*$", re.MULTILINE)
SEMVER_PATTERN = re.compile(r"^[0-9]+\.[0-9]+\.[0-9]+(?:[a-zA-Z0-9.-]+)?$")
SHA_PATTERN = re.compile(r"^[0-9a-f]{64}$")
HELP_COMMAND_PATTERN = re.compile(
    r"^(?P<indent>[ \t]*)-[ \t]+Replidec[ \t]+--help[ \t]*$",
    re.MULTILINE,
)
VERSION_COMMAND_PATTERN = re.compile(
    r"^[ \t]*-[ \t]+Replidec[ \t]+--version[ \t]*$",
    re.MULTILINE,
)


def _replace_once(pattern: re.Pattern[str], replacement: str, text: str, label: str) -> str:
    updated, count = pattern.subn(replacement, text, count=1)
    if count != 1:
        raise ValueError(f"Expected exactly one {label} field, found {count}")
    return updated


def update_recipe_text(text: str, version: str, sha256: str) -> str:
    if SEMVER_PATTERN.fullmatch(version) is None:
        raise ValueError(f"Invalid release version: {version}")
    sha256 = sha256.lower()
    if SHA_PATTERN.fullmatch(sha256) is None:
        raise ValueError("SHA256 must contain exactly 64 hexadecimal characters")

    updated = _replace_once(
        VERSION_PATTERN,
        rf"\g<1>{version}\g<2>",
        text,
        "version",
    )
    updated = _replace_once(
        SHA256_PATTERN,
        rf"\g<1>{sha256}",
        updated,
        "sha256",
    )
    updated = _replace_once(
        BUILD_PATTERN,
        r"\g<1>0",
        updated,
        "build number",
    )

    # RepliDec now uses vX.Y.Z tags; older recipes linked to v.X.Y.Z.
    updated = updated.replace("/blob/v.{{ version }}/", "/blob/v{{ version }}/")

    if VERSION_COMMAND_PATTERN.search(updated) is None:
        help_command = HELP_COMMAND_PATTERN.search(updated)
        if help_command is None:
            raise ValueError("Unable to locate the Bioconda command test section")
        version_command = f"{help_command.group('indent')}- Replidec --version"
        updated = (
            updated[: help_command.end()]
            + f"\n{version_command}"
            + updated[help_command.end() :]
        )

    return updated


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--recipe", type=Path, required=True)
    parser.add_argument("--version", required=True)
    parser.add_argument("--sha256", required=True)
    args = parser.parse_args()

    original = args.recipe.read_text(encoding="utf-8")
    updated = update_recipe_text(original, args.version, args.sha256)
    args.recipe.write_text(updated, encoding="utf-8")


if __name__ == "__main__":
    main()
