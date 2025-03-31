#!/usr/bin/env python

from setuptools import setup, find_packages

with open('README.rst', 'r') as readme_file:
    readme = readme_file.read()

requirements = ["biopython>=1.77",
                "future>=0.18.2",]


setup(
    name='Replidec',
    version='0.3.3',
    author="Xue Peng",
    author_email='xue.peng@helmholtz-muenchen.de',
    keywords='Replidec',
    description='Replication Cycle Detector for Phages',
    long_description=readme,
    long_description_content_type="text/markdown",
    url='https://github.com/deng-lab/Replidec',
    #packages=find_packages(['Replidec', 'Replidec.*']),
    packages=find_packages(),
    install_requires=requirements,
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Intended Audience :: Science/Research",],
    python_requires='>=3.8',
    license="MIT license",
    zip_safe=True,
    entry_points={ 'console_scripts': ['Replidec=Replidec.Replidec_cmdline:main']},
)
