from setuptools import setup,find_packages


with open("README.md", "r") as fh:
    long_description = fh.read()

INSTALL_REQUIRES=[
        "biopython>=1.77",
        "future>=0.18.2",
        ]

setup(name='Replidec',
      #version='0.2.1',
      version='0.2.0',
      description='Replication Cycle Detector for Phages',
      author='Xue Peng',
      author_email='xue.peng@helmholtz-muenchen.de',
      url='https://github.com/pengSherryYel/Replidec_v0.2.1',
      packages=find_packages(where='.',
                            include=['Replidec.*']),
      long_description=long_description,
      long_description_content_type="text/markdown",
      install_requires=INSTALL_REQUIRES,
      include_package_data=True,
      classifiers=[
          "Programming Language :: Python :: 3",
          "License :: OSI Approved :: MIT License",
          "Operating System :: OS Independent",
          "Topic :: Scientific/Engineering :: Bio-Informatics",
          "Intended Audience :: Science/Research",],
      python_requires='>=3.8',
      entry_points={ 'console_scripts': ['Replidec=Replidec.Replidec_cmdline:main']},
)
