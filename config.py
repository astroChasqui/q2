import os


path = os.path.dirname(os.path.realpath(__file__))
path = os.path.join(path, 'Data')

COLORTEFF_PATH  = os.path.join(path, 'ColorTeff')
MODATM_PATH     = os.path.join(path, 'ModelAtmospheres')
ISOCHRONES_PATH = os.path.join(path, 'Isochrones')
OTHER_PATH      = os.path.join(path, 'Other')
