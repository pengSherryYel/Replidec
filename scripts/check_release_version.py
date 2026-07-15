#!/usr/bin/env python3
"""Validate a release tag against RepliDec's single source version."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


VERSION_PATTERN = re.compile(r'^__version__\s*=\s*["\']([^"\']+)["\']', re.MULTILINE)
TAG_PATTERN = re.compile(r"^v(?P<version>[0-9]+\.[0-9]+\.[0-9]+(?:[a-zA-Z0-9.-]+)?)$")


def read_project_version(project_root: Path) -> str:
    version_file = project_root / "Replidec" / "__init__.py"
    match = VERSION_PATTERN.search(version_file.read_text(encoding="utf-8"))
    if match is None:
        raise ValueError(f"Unable to read __version__ from {version_file}")
    return match.group(1)


def validate_release_tag(tag: str, project_version: str) -> str:
    match = TAG_PATTERN.fullmatch(tag)
    if match is None:
        raise ValueError(f"Release tag must use vX.Y.Z format, received: {tag}")
    tag_version = match.group("version")
    if tag_version != project_version:
        raise ValueError(
            f"Release tag {tag} does not match project version {project_version}"
        )
    return tag_version


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tag", required=True)
    parser.add_argument("--project-root", type=Path, default=Path.cwd())
    parser.add_argument("--github-output", type=Path)
    args = parser.parse_args()

    project_version = read_project_version(args.project_root)
    version = validate_release_tag(args.tag, project_version)

    if args.github_output:
        with args.github_output.open("a", encoding="utf-8") as output:
            output.write(f"version={version}\n")
            output.write(f"tag={args.tag}\n")

    print(version)


if __name__ == "__main__":
    main()
