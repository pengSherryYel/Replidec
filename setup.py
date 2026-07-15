#!/usr/bin/env python

import re
from pathlib import Path

from setuptools import setup, find_packages

ROOT = Path(__file__).resolve().parent

readme = (ROOT / "README.md").read_text(encoding="utf-8")
version_match = re.search(
    r'^__version__\s*=\s*["\']([^"\']+)["\']',
    (ROOT / "Replidec" / "__init__.py").read_text(encoding="utf-8"),
    re.MULTILINE,
)
if version_match is None:
    raise RuntimeError("Unable to determine RepliDec version")

version = version_match.group(1)

requirements = ["biopython>=1.77",
                "future>=0.18.2",]


setup(
    name='Replidec',
    version=version,
    author="Xue Peng",
    author_email='peng_sherry@outlook.com',
    keywords='Replidec',
    description='Replication Cycle Decipher for Phages',
    long_description=readme,
    long_description_content_type="text/markdown",
    url='https://github.com/pengSherryYel/Replidec',
    #packages=find_packages(['Replidec', 'Replidec.*']),
    packages=find_packages(),
    install_requires=requirements,
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Intended Audience :: Science/Research",],
    python_requires='>=3.8',
    license="MIT",
    zip_safe=True,
    entry_points={
        'console_scripts': [
            'Replidec=Replidec.Replidec_cmdline:main',
            'replidec=Replidec.Replidec_cmdline:main',
        ]
    },
)
