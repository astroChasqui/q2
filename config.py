#import platform
import os


system_path = os.path.dirname(os.path.realpath(__file__))

COLORTEFF_PATH  = os.path.join(system_path, 'ColorTeff')
MODATM_PATH     = os.path.join(system_path, 'ModelAtmospheres')
ISOCHRONES_PATH = os.path.join(system_path, 'Isochrones')
OTHER_PATH      = os.path.join(system_path, 'Other')
