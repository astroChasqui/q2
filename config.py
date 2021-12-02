import os
import logging
import matplotlib.pyplot as plt
import glob
import re


logger = logging.getLogger(__name__)

path = os.path.expanduser("~") + '/q2-tools'
path = os.path.join(path, 'Data')
moogpath = os.path.expanduser("~") + '/q2-tools/MOOG-for-q2'

COLORTEFF_PATH  = os.path.join(path, 'ColorTeff')
MODATM_PATH     = os.path.join(path, 'ModelAtmospheres')
ISOCHRONES_PATH = os.path.join(path, 'Isochrones')
OTHER_PATH      = os.path.join(path, 'Other')

plt.rc("font", family='serif', serif='Ubuntu', monospace='Ubuntu Mono', \
               size=14)
plt.rc("axes", labelsize=15, titlesize=12)
plt.rc("xtick", top=True, direction='in', labelsize=14)
plt.rc("xtick.major", size=8, width=1)
plt.rc("ytick", right=True, direction='in', labelsize=14)
plt.rc("ytick.major", size=8, width=1)
plt.rc("lines", markersize=10, markeredgewidth=2)
plt.rc("lines", linewidth=3)

def moog_is_available():
    """You should be able to run MOOGSILENT from the command line in order
    to use the MOOG features included in q2. This function checks if
    MOOG is available on your system. If False, you wont be able to
    connect q2 to MOOG and many things will fail.
    """

    #default MOOG-for-q2 path
    moogpath = os.path.expanduser("~") + '/q2-tools/MOOG-for-q2'

    if os.system('which '+moogpath+'/MOOGSILENT >/dev/null'):

        logger.warning("MOOGSILENT is not available. Don't worry, we'll take care of it. \n")

        download_moog()
        moogsilent_without_sm()
        installing_moog()
        MOOGq2_is_available()

    else:
        print("MOOGSILENT is available in " + moogpath)
        logger.info("MOOGSILENT is available in " + moogpath)
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
        logger.warning("Data folder don't exist. Don't worry, we'll take care of it.\n ")
        print('Data download requires internet and may take a few minutes...')
        download_data()


def download_data():
    """Download and extract Data Folder to a HOME directory"""
    dir = os.path.expanduser("~") + '/q2-tools/Data'

    # If directory already exists, we pass
    if os.path.exists(dir):
        pass
    else:
        os.system('wget https://www.dropbox.com/s/g8pabfpnzybaa46/Data.zip -P $HOME')
        os.system('unzip $HOME/Data.zip -d $HOME' )
        os.system('mv $HOME/Data $HOME/q2-tools')
        os.system('rm -rf $HOME/Data/')
        os.system('rm $HOME/Data.zip')


def download_moog():
    """Download and extract MOOG to a home directory"""

    dir = os.path.expanduser("~") + '/q2-tools/MOOG-for-q2'

    # If directory already exists, we pass
    if os.path.exists(dir):
        pass

    else:
        os.system('wget https://www.as.utexas.edu/~chris/MOOGNOV2019.tar.gz -P $HOME')
        os.system('tar -xzf $HOME/MOOGNOV2019.tar.gz -C $HOME' )

        os.makedirs(dir)
        os.system('mv $HOME/moognov2019/* $HOME/q2-tools/MOOG-for-q2')
        os.system('rm -rf $HOME/moognov2019/')
        os.system('rm $HOME/MOOGNOV2019.tar.gz')


def moogsilent_without_sm(moogpath = os.path.expanduser("~") + '/q2-tools/MOOG-for-q2'):
    """Editing files to remove sm dependencies"""

    files = glob.glob(moogpath + '/' + '*.f')

    #Exclude call plot lines

    # check : replace
    dict_replace = {
        "call pltabun": "\n         print *, 'No sm, no plots!'\nc          call pltabun",
        "call pltcog": "\n         print *, 'No sm, no plots!'\nc          call pltcog",
        "call pltflux": "\n         print *, 'No sm, no plots!'\nc          call pltflux",
        "call pltspec": "print *, 'No sm, no plots!'\nc call pltspec",
        "call makeplot": "\nc          call makeplot",
        "call drawcurs": "\nc          call drawcurs",   
        "call pointcurs": "\nc          call pointcurs"
    }

    for file in files:
        with open(file) as r:
            text = r.read()
            for check, replace in dict_replace.items():
                text = re.sub(check, replace, text)

        with open(file, "w") as w:
             w.write(text)

    # Edit Begin file
    with open(moogpath + '/' + 'Begin.f') as r:
            begin = r.read().replace("           write (*,*) 'maxline', maxline\n           pause", "").replace('num = 60', 'num = len(trim(moogpath))')

    with open(moogpath + '/' + "Begin.f", "w") as w:
             w.write(begin)

    dict_remove = {
    "Abunplot.o": " ",
    "Binplot.o": " ",
    "Cogplot.o": "",
    "Defcolor.o": "",
    "Drawcurs.o": "",
    "Fluxplot.o": "",
    "Makeplot.o": "",
    "Pointcurs.o": "",
    "Specplot.o": "",
    "Moog.o" : "Moogsilent.o",
    "X11LIB = .*" : "",
    "SMLIB = .*" : "",
    "FC = .*" : "FC = gfortran -std=legacy  -w",
    "CC = .*" : "CC = cc",
    "-L.*" : " "
    }

    makesilent= glob.glob(moogpath + '/' + '*silent')

    for file in makesilent:
        with open(file) as r:
            text = r.read()
            for check, remove in dict_remove.items():
                text = re.sub(check, remove, text)

        with open(file, "w") as w:
                w.write(text)


def installing_moog(moogpath = os.path.expanduser("~") + '/q2-tools/MOOG-for-q2', machine='pcl'):

    olddir = os.getcwd()
    dir = str(moogpath)
    os.chdir(dir)

    print('We need some information about your machine...\n')
    print(" 'rh64' for 64-bit linux system\n 'rh' for 32-bit linux system\n 'maclap' for mac laptop\n 'macdesk' for mac desktop\n")

    system = input(str("Please, enter the type of your system ('rh64', 'rh', 'maclap', 'macdesk') : "))

    while system not in ['rh64', 'rh', 'maclap', 'macdesk']:
        print('Sorry, that is not an option \n Try again!')
        system = input(str("Enter the type of your system ('rh64', 'rh', 'maclap', 'macdesk') : "))

    if system in ['maclap', 'macdesk']:
        machine = 'mac'
    else:
        machine = 'pcl'

    # Edit moogsilent file
    with open('Moogsilent.f') as r:
            file = r.read().replace("/Users/chris/CODES/moognov2019/", moogpath).replace('machine = "mac"', 'machine = ' + '"' + machine + '"')

    with open("Moogsilent.f", "w") as w:
             w.write(file)


    if system == "rh64":
        os.system('make -f Makefile.rh64silent')

    elif system == "rh":
        os.system('make -f Makefile.rhsilent')

    elif system == "maclap":
        os.system('make -f Makefile.maclapsilent')

    elif system == "macdesk":
        os.system('make -f Makefile.macdesksilent')

    os.chdir(olddir)


def MOOGq2_is_available(moogpath=os.path.expanduser("~") + '/q2-tools/MOOG-for-q2'):

    if os.system('which '+moogpath+'/MOOGSILENT >/dev/null'):
        print("MOOGSILENT is not available")
    else:
        print("MOOGSILENT without SM is available in " + str(moogpath))
