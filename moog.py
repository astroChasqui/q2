import numpy as np
import os
import logging
from config import *

logger = logging.getLogger(__name__)


class Driver:
    """Set the options for your MOOG driver."""
    def __init__(self, mode='abfind'):
        self.mode = mode
        self.standard_out = 'moog.std'
        self.summary_out = 'moog.sum'
        self.model_in = 'model.in'
        self.lines_in = 'lines.in'
        self.plot = 0
        self.hfs_species = None

    def create_file(self, file_name="batch.par"):
        """Creates the MOOG driver file."""
        self.file_name = file_name
        f = open(file_name, 'w')
        if self.mode == 'abfind':
            if self.hfs_species:
                f.write('blends\n')
            else:
                f.write('abfind\n')
        else:
            f.write(self.mode+'\n')
        f.write('standard_out '+self.standard_out+'\n')
        f.write('summary_out  '+self.summary_out+'\n')
        f.write('model_in     '+self.model_in+'\n')
        f.write('lines_in     '+self.lines_in+'\n')
        f.write('atmosphere   1\n')
        f.write('molecules    1\n')
        f.write('lines        1\n')
        f.write('flux/int     0\n')
        f.write('damping      1\n')
        f.write('freeform     1\n')
        f.write('plot         '+str(self.plot)+'\n')
        if self.hfs_species:
            f.write('blenlimits\n')
            f.write(' 2.0 0.01 '+self.hfs_species+'\n')
        if self.mode == 'cog':
            f.write('coglimits\n')
            f.write('  -6.5 -3.5 0.1 0 0\n')
        f.close()


def create_model_in(Star, file_name='model.in'):
    """Creates a model atmosphere file for MOOG from the model_atmosphere
    attribute of a Star object.
    """
    try:
        Star.vt
    except:
        logger.error('Moog model_in file requires a microturbulence (vt)')
        return None
    if hasattr(Star, 'model_atmosphere'):
        Star.moog_model_in_name = file_name
    else:
        logger.error('No model data to write to moog model_in file.')
        return None

    with open('head.tmp', 'w') as f:
        f.write('KURUCZ\n')
        f.write('TEFF='+str(Star.teff)+',LOGG='+str(Star.logg)+
                ',[FE/H]='+str(Star.feh)+','+Star.model_atmosphere_grid+'\n')
        nd = len(Star.model_atmosphere['T'])
        f.write('ND=       '+str(nd)+'\n')

    with open('body.tmp', 'w') as f:
        for idx in range(nd):
            f.write("{0:.8E} {1:.1F} {2:.3E} {3:.3E} {4:.3E}\n".format(\
                    Star.model_atmosphere['RHOX'][idx],\
                    Star.model_atmosphere['T'][idx],\
                    Star.model_atmosphere['P'][idx],\
                    Star.model_atmosphere['XNE'][idx],\
                    Star.model_atmosphere['ABROSS'][idx])
                   )

    with open('tail.tmp', 'w') as f:
        f.write('%5.2F\n' %Star.vt)
        if Star.model_atmosphere_grid != 'marcs':
            path = os.path.join(MODATM_PATH, 'kurucz')
            fabund = open(os.path.join(path, 'p00.'+Star.model_atmosphere_grid),
                          'r')
        else:
            path = os.path.join(MODATM_PATH, 'marcs')
            fabund = open(os.path.join(path, 'z+0.00'), 'r')

        line = fabund.readline()
        f.write(line[0:12]+' '+str(Star.feh)+'\n')
        line = fabund.readline()
        while line:
            species = line[0:2]
            if Star.model_atmosphere_grid == 'marcs':
                abund = float(line[3:9])+Star.feh
                #alpha-element enhancement
                if species==' 8' or species=='10' or species=='12' or \
                   species=='14' or species=='16' or species=='18' or \
                   species=='20' or species=='22':
                    afe = -0.4*Star.feh
                    if Star.feh >=  0: afe=0.0
                    if Star.feh <= -1: afe=0.4
                    abund = abund+afe
            else:
                abund = 12.+np.log10(np.power(10, float(line[3:9]))/0.92040)+ \
                        Star.feh
            abund = str('%5.2F' %abund)
            f.write(species+' '+abund+'\n')
            line = fabund.readline()
        fabund.close()
        f.write('NMOL      22\n')
        f.write('  101.0   106.0   107.0   108.0   112.0  126.0\n')
        f.write('  606.0   607.0   608.0\n')
        f.write('  707.0   708.0\n')
        f.write('  808.0   812.0   822.0\n')
        f.write('  10108.0 60808.0\n')
        f.write('  6.1     7.1     8.1   12.1  22.1  26.1\n')

    file_list = ['head.tmp', 'body.tmp', 'tail.tmp']
    with open(file_name, 'w') as outfile:
        for one_file in file_list:
            with open(one_file) as infile:
                outfile.write(infile.read())
    for one_file in file_list:
        os.unlink(one_file)
    logger.info('Moog infile model atmosphere created: '+file_name)


def create_lines_in(Star, species=0, file_name='lines.in'):
    """Creates a line list file for MOOG"""
    if species > 0:
        idx = np.where(np.logical_and(Star.linelist['species'] == species,\
                                       Star.linelist['ew'] >= 0))[0]
    else:
        #species = 0 means all species
        idx = np.where(np.logical_and(Star.linelist['species'] > species,\
                                       Star.linelist['ew'] >= 0))[0]

    nlines = len(idx)
    if nlines == 0:
        logger.warning('No lines found for '+Star.name)
        return False
    else:
        logger.info(str(nlines)+' lines found for '+Star.name)
    gf_values = Star.linelist['gf'][idx]
    gf10 = [10**gfx for gfx in Star.linelist['gf'][idx] if gfx >= 0]
    if len(gf10) == len(Star.linelist['gf'][idx]):
        logger.info('all gf values for this species are positive --> 10^gf')
        #gf_values = gf10
        Star.linelist['gf'][idx] = gf10
    #Star.linelist['gf'][idx] = gf_values

    with open(file_name, 'w') as f:
        f.write("MOOG linelist created by q2\n")
        for lidx in idx:
            f.write("{0:10.4f} {1:4.1f} {2:6.3f} {3:5.3f} 3 0 {4:5.1f}\n".format(\
                    Star.linelist['wavelength'][lidx],\
                    Star.linelist['species'][lidx],\
                    Star.linelist['ep'][lidx],\
                    Star.linelist['gf'][lidx],\
                    Star.linelist['ew'][lidx])
                   )

    Star.linelist['gf'][idx] = gf_values

    logger.info('Moog line list created: '+file_name)
    return True


def abfind(Star, species, species_id):
    """Runs MOOG with abfind driver for a given Star and species
    
    Star is a star object; must have all attributes in place
    species could be 26.0 for Fe I, for example
    species_id is a string that will become a new attribute for the Star object
    Example: abfind(s, 26.1, 'fe2')
    s.fe2 #shows result from abfind
    MD is the moog driver object
    """
    k = Star.linelist['species'] == species
    negs = [wx for wx in Star.linelist['wavelength'][k] if wx < 0]
    if len(negs) == 0:
        MD = Driver() #normal
    else:
        MD = Driver() #hfs
        MD.hfs_species = str(round(species))
    MD.create_file('batch.par')

    create_model_in(Star)
    found_lines = create_lines_in(Star, species=species)
    if not found_lines:
        logger.warning('Did not run abfind (no lines found)')
        return False
    os.system('MOOGSILENT > moog.log 2>&1')
    f = open(MD.summary_out, 'r')
    line=''
    while line[0:10] != 'wavelength':
        line = f.readline()
    if 'ID' in line:
        moogjul2014 = True
    else:
        moogjul2014 = False
    ww, ep, ew, rew, ab, difab = [], [], [], [], [], []
    while line:
        line = f.readline()
        if line[0:7] == 'average': break
        linesplit = line.split()
        if float(linesplit[6]) > 999.: #exclude dummies (hfs)
            continue
        ww.append(float(linesplit[0]))
        if moogjul2014: #MOOGJUL2014 adds a new column 'ID' to moog.sum
            ep.append(float(linesplit[2]))
            ew.append(float(linesplit[4]))
            rew.append(float(linesplit[5]))
            ab.append(float(linesplit[6]))
        else: #older versions of MOOG don't have 'ID' but 'EP' in 2nd col
            ep.append(float(linesplit[1]))
            ew.append(float(linesplit[3]))
            rew.append(float(linesplit[4]))
            ab.append(float(linesplit[5]))
        difab.append(None)
    f.close()
    os.unlink(MD.file_name)
    os.unlink(MD.model_in)
    os.unlink(MD.lines_in)
    os.unlink(MD.summary_out)
    os.unlink(MD.standard_out)
    os.unlink('moog.log')
    if os.path.isfile('fort.99'):
        os.unlink('fort.99')

    x = {'ww': np.array(ww), 'ep': np.array(ep), 'ew': np.array(ew),\
         'rew': np.array(rew), 'ab': np.array(ab), 'difab': np.array(difab)}
    setattr(Star, species_id, x)
    logger.info('Successfully ran abfind')
    return True


def cog(Star, species, cog_id):
    """Runs MOOG with cog driver for a given Star and species

    Star is a star object; must have all attributes need by MOOG set.
    species could be 26.0 for Fe I, for example. cog_id is a string that
    will become a new attribute for the Star object. For example:
    >>>cog(s, 26.1, 'cog_fe2')
    s.cog_fe2 #shows result from cog
    MD is the moog driver object
    """
    k = Star.linelist['species'] == species
    #negs = [wx for wx in Star.linelist['wavelength'][k] if wx < 0]
    MD = Driver(mode='cog')
    MD.create_file()
    create_model_in(Star)
    found_lines = create_lines_in(Star, species=species)
    if not found_lines:
        logger.warning('Did not run cog (no lines found)')
        return False
    os.system('MOOGSILENT > moog.log 2>&1')

    f = open(MD.summary_out, 'r')
    line = f.readline()
    cog_obj = {}
    while line:
        line = f.readline()
        if line.startswith('wavelength'):
            npt = int(line.split('=')[5]) #number of cog points
            #wavelength = round(float(line.split('=')[1].split()[0]), 1)
            wavelength = float(line.split('=')[1].split()[0])
            line = f.readline()
            x, y = [], []
            for i in range(int(np.ceil(npt/5.))):
                line = f.readline()
                for j in range(len(line.split())/2):
                    x.append(float(line.split()[2*j].replace(',', '')))
                    y.append(float(line.split()[2*j+1]))
            cog_obj[wavelength] = {'loggf': np.array(x), 'logrw': np.array(y)}
    f.close()

    os.unlink(MD.file_name)
    os.unlink(MD.model_in)
    os.unlink(MD.lines_in)
    os.unlink(MD.summary_out)
    os.unlink(MD.standard_out)
    os.unlink('moog.log')

    setattr(Star, cog_id, cog_obj)
