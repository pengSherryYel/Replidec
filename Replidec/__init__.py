"""

Replidec package - bacteriophage lifestyle prediction (temperate|virulent|chronic) based on the naive bayes classifier.

Written by Xue Peng (Email: xue.peng@helmholtz-muenchen.de)

Python modules
----------------
The package consists of the following Python modules:
* bayes_classifier_batch
* bayes_classifier_contig
* bayes_classifier_genomes

"""
__version__ = "0.2.1"
__all__ = ["Replidec", "Replidec_cmdline","Replidec_multi","utility"]
from Replidec import *
