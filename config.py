import os
import logging

logger = logging.getLogger(__name__)


path = os.path.dirname(os.path.realpath(__file__))
path = os.path.join(path, 'Data')

COLORTEFF_PATH  = os.path.join(path, 'ColorTeff')
MODATM_PATH     = os.path.join(path, 'ModelAtmospheres')
ISOCHRONES_PATH = os.path.join(path, 'Isochrones')
OTHER_PATH      = os.path.join(path, 'Other')

def moog_is_available():
    """You should be able to run MOOGSILENT from the command line in order
    to use the MOOG features included in q2. This function checks if
    MOOG is available on your system. If False, you wont be able to
    connect q2 to MOOG and many things will fail.
    """
    if os.system('which MOOGSILENT >/dev/null'):
        logger.warning("MOOGSILENT is not available")
        return False
    else:
        logger.info("MOOGSILENT is available")
        return True

def data_are_available():
    """q2 needs data files with model atmosphere and isochrone grids.
    These files can be downloaded from:
    http://www.astrochasqui.com/projects/astro/share/q2Data.tar.gz
    They need to be extracted inside the q2 directory.
    'tar xvfz q2Data.tar.gz' will create the Data folder.
    """
    if os.path.exists(path):
        logger.info("Data folder exists")
        return True
    else:
        logger.warning("Data folder does not exist")
        return False
