#!/usr/bin/env python3

"""
DGRPy is a python package for DGRP microbiome data analysis
The project is under early development and this is my initial release.
NOTE: This is an unstable early release so as soon as I publish stable version I will also add documentation.
Many functions are not at final form at all. For example, JadvalOTU can only load OTU data only in format as OTUs_sample.csv

@author: Farid MUSA(mmtechslv)
"""

__author__ = "Farid MUSA"
__copyright__ = "Copyright (C) 2019, DGRPy Project"
__credits__ = ["Farid MUSA"]
__license__ = "GPLv3"
__version__ = "1.1"
__maintainer__ = "Farid MUSA"
__email__ = "farid.musa.h@gmail.com"
__status__ = "Development"



from .jadval_otu import OtuTable as jOtuTable
from .jadval_tdb import TaxonomyGG as jGreengenes
from .jadval_tdb import TaxonomySILVA as jSILVA
from ._jadval_error import get_error_log as jGetErrorLog
from ._jadval_error import get_error as jGetError
from ._jadval_error import get_last_error as jGetLastError

