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


from ._jadval_shared import JadvalPrime
import jadval._jadval_shared as jShared

import jadval._jadval_constants as jConst

import jadval._jadval_error as jError
from ._jadval_error import raise_error as jRaiseError
from ._jadval_error import raise_warning as jRaiseWarning
from ._jadval_error import error_handler as jErrorHandler
