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
UPPER_ENTRY_POINT_PATTERN = re.compile(
    r"^(?P<indent>[ \t]*)-[ \t]+Replidec[ \t]*=[ \t]*Replidec\.Replidec_cmdline:main[ \t]*$",
    re.MULTILINE,
)
LOWER_ENTRY_POINT_PATTERN = re.compile(
    r"^[ \t]*-[ \t]+replidec[ \t]*=[ \t]*Replidec\.Replidec_cmdline:main[ \t]*$",
    re.MULTILINE,
)
HELP_COMMAND_PATTERN = re.compile(
    r"^(?P<indent>[ \t]*)-[ \t]+Replidec[ \t]+--help[ \t]*$",
    re.MULTILINE,
)
VERSION_COMMAND_PATTERN = re.compile(
    r"^(?P<indent>[ \t]*)-[ \t]+Replidec[ \t]+--version[ \t]*$",
    re.MULTILINE,
)
LOWER_HELP_COMMAND_PATTERN = re.compile(
    r"^(?P<indent>[ \t]*)-[ \t]+replidec[ \t]+--help[ \t]*$",
    re.MULTILINE,
)
LOWER_VERSION_COMMAND_PATTERN = re.compile(
    r"^[ \t]*-[ \t]+replidec[ \t]+--version[ \t]*$",
    re.MULTILINE,
)


def _replace_once(pattern: re.Pattern[str], replacement: str, text: str, label: str) -> str:
    updated, count = pattern.subn(replacement, text, count=1)
    if count != 1:
        raise ValueError(f"Expected exactly one {label} field, found {count}")
    return updated


def _insert_line_after(
    text: str,
    anchor_pattern: re.Pattern[str],
    existing_pattern: re.Pattern[str],
    line: str,
    section: str,
) -> str:
    if existing_pattern.search(text) is not None:
        return text

    anchor = anchor_pattern.search(text)
    if anchor is None:
        raise ValueError(f"Unable to locate the Bioconda {section} section")

    inserted = f"{anchor.group('indent')}- {line}"
    return text[: anchor.end()] + f"\n{inserted}" + text[anchor.end() :]


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

    updated = _insert_line_after(
        updated,
        UPPER_ENTRY_POINT_PATTERN,
        LOWER_ENTRY_POINT_PATTERN,
        "replidec = Replidec.Replidec_cmdline:main",
        "entry_points",
    )
    updated = _insert_line_after(
        updated,
        HELP_COMMAND_PATTERN,
        VERSION_COMMAND_PATTERN,
        "Replidec --version",
        "command test",
    )
    updated = _insert_line_after(
        updated,
        VERSION_COMMAND_PATTERN,
        LOWER_HELP_COMMAND_PATTERN,
        "replidec --help",
        "command test",
    )
    updated = _insert_line_after(
        updated,
        LOWER_HELP_COMMAND_PATTERN,
        LOWER_VERSION_COMMAND_PATTERN,
        "replidec --version",
        "command test",
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
