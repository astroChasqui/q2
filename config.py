import platform
import os


if platform.system() == 'Windows':
    system_path = 'C:\Users\Ivan\Documents\Dropbox\Code\python\q2\Data'
if platform.system() == 'Linux':
    system_path = '/home/ivan/Dropbox/Code/python/q2/Data'

COLORTEFF_PATH  = os.path.join(system_path, 'ColorTeff')
MODATM_PATH     = os.path.join(system_path, 'ModelAtmospheres')
ISOCHRONES_PATH = os.path.join(system_path, 'Isochrones')
OTHER_PATH      = os.path.join(system_path, 'Other')
