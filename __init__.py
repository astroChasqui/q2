"""q2 was created to facilitate standard 1D/LTE spectroscopic analyses
of stars. It helps determining fundamental parameters (Teff, logg,
[Fe/H], etc.) of solar-type stars given an observational data set. q2
requires the 2014 version of the spectrum synthesis code MOOG
(http://www.as.utexas.edu/~chris/moog.html) and files containing the
model atmosphere and isochrone grids.
"""

from __future__ import print_function
import matplotlib
matplotlib.use('Agg')
from config import *
from star import *
import moog
import specpars
import errors
import abundances
import grids
import yypars
import logging

logger = logging.getLogger(__name__)

__author__ = 'Ivan Ramirez (UT Austin)'
__email__ = 'ivan@astro.as.utexas.edu'


logging.basicConfig(level=logging.ERROR)
logger.setLevel('WARNING')
moog_check = moog_is_available()
data_check = data_are_available()
logger.setLevel('ERROR')
