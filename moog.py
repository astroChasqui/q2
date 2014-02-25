from astropy.io import ascii
import numpy as np
import os
import logging
from config import *

logger = logging.getLogger(__name__)


class Driver:
    def __init__(self):
        self.standard_out = 'moog.std'
        self.summary_out = 'moog.sum'
        self.model_in = 'model.in'
        self.lines_in = 'lines.in'
        self.plot = 0
        self.hfs_species = None

    def create_file(self, file_name):
        self.file_name = file_name
        f = open(file_name, 'w')
        if self.hfs_species:
            f.write('blends\n')
        else:
            f.write('abfind\n')
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
        f.close()


def create_model_in(Star, file_name='model.in'):
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

    f = open('head.tmp', 'w')
    f.write('KURUCZ\n')
    f.write('TEFF='+str(Star.teff)+',LOGG='+str(Star.logg)+
            ',[FE/H]='+str(Star.feh)+','+Star.model_atmosphere_grid+'\n')
    nd = len(Star.model_atmosphere)
    f.write('ND=       '+str(nd)+'\n')
    f.close()

    ascii.write(Star.model_atmosphere, 'body.tmp', Writer=ascii.NoHeader,
                formats={'RHOX':'%.8E','T':'%.1F','P':'%.3E',
                'XNE':'%.3E','ABROSS':'%.3E'})

    f = open('tail.tmp', 'w')
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
    f.close()

    file_list = ['head.tmp', 'body.tmp', 'tail.tmp']
    with open(file_name, 'w') as outfile:
        for one_file in file_list:
            with open(one_file) as infile:
                outfile.write(infile.read())
    for one_file in file_list:
        os.unlink(one_file)
    logger.info('Moog infile model atmosphere created: '+file_name)


def create_lines_in(Star, species=0, file_name='lines.in'):
    #species = 0 means all species
    if species > 0:
        #idx = np.where(Star.linelist['species'] == species)
        idx = np.logical_and(Star.linelist['species'] == species,
                             Star.linelist['ew'] >= 0)
    else:
        idx = np.logical_and(Star.linelist['species'] > 0,
                             Star.linelist['ew'] >= 0)

    nlines = len(idx[idx==True])
    if nlines == 0:
        logger.warning('No lines found for '+Star.name)
        return False
    else:
        logger.info(str(nlines)+' lines found for '+Star.name)
    gf_values = Star.linelist['gf'][idx]
    gf10 = [10**gfx for gfx in Star.linelist['gf'][idx] if gfx > 0]
    if len(gf10) == len(Star.linelist['gf'][idx]):
        logger.info('all gf values for this species are positive --> 10^gf')
        gf_values = gf10
    ll_data = {'wavelength': Star.linelist['wavelength'][idx],
               'species': Star.linelist['species'][idx],
               'ep': Star.linelist['ep'][idx],
               'gf': gf_values,
               'damp': 3+np.zeros(nlines),
               'zero': np.zeros(nlines),
               'ew': Star.linelist['ew'][idx]}
    ascii.write(ll_data, file_name, Writer=ascii.CommentedHeader,
                names=['wavelength', 'species', 'ep', 'gf', 'damp', 'zero', 'ew'])
    logger.info('Moog line list created: '+file_name)
    return True


def abfind(Star, species, species_id,):
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
    os.system('MOOGSILENT > moog.log')
    f = open(MD.summary_out, 'r')
    line=''
    while line[0:20] != 'wavelength        EP':
        line = f.readline()
    ww, ep, ew, rew, ab, difab = [], [], [], [], [], []
    while line:
        line = f.readline()
        if line[0:3] == 'ave': break
        if float(line[50:60]) > 999.: #exclude dummies (hfs)
            continue
        ww.append(line[0:10])
        ep.append(float(line[10:20]))
        ew.append(float(line[30:40]))
        rew.append(np.log10(0.001*float(line[30:40])/float(line[0:10])))
        ab.append(float(line[50:60]))
        difab.append(None)
    f.close()
    os.unlink(MD.file_name)
    os.unlink(MD.model_in)
    os.unlink(MD.lines_in)
    os.unlink(MD.summary_out)
    os.unlink(MD.standard_out)
    os.unlink('moog.log')
    x = np.zeros((len(ww),), dtype=[('ww', 'f8'), ('ep', 'f8'), ('ew', 'f8'),
                                    ('rew', 'f8'), ('ab', 'f8'),
                                    ('difab', 'f8')])
    x['ww'] = ww
    x['ep'] = ep
    x['ew'] = ew
    x['rew'] = rew
    x['ab'] = ab
    x['difab'] = difab
    setattr(Star, species_id, x)
    logger.info('Successfully ran abfind')
    return True
