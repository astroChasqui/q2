"""q2 was created to facilitate standard 1D/LTE spectroscopic analyses
of stars. It helps determining fundamental parameters (Teff, logg,
[Fe/H], etc.) of solar-type stars given an observational data set. q2
requires a slightly modified version of the spectrum synthesis code
MOOG (http://www.as.utexas.edu/~chris/moog.html) and files containing
the model atmosphere grids.

This package is part of the q2 project: http://inti.as.utexas.edu/q2.
"""

from __future__ import print_function
from config import *
from star import *
import moog
import specpars
import errors
import abundances
import grids
import yyage
import logging


logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.ERROR)
#to change level interactively: q2.logger.setLevel('DEBUG')

__author__ = 'Ivan Ramirez (UT Austin)'
__email__ = 'ivan@astro.as.utexas.edu'

logger.info('q2 package successfully imported')
