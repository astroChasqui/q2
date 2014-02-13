import platform


if platform.system() == 'Windows':
    MODATM_PATH = 'C:\Users\Ivan\Documents\Dropbox\Code\python\q2\Data\ModelAtmospheres\\'
    ISOCHRONES_PATH = 'C:\Users\Ivan\Documents\Dropbox\Code\python\q2\Data\Isochrones\\'
if platform.system() == 'Linux':
    MODATM_PATH = '/home/ivan/Dropbox/Code/python/q2/Data/ModelAtmospheres/'
    ISOCHRONES_PATH = '/home/ivan/Dropbox/Code/python/q2/Data/Isochrones/'
