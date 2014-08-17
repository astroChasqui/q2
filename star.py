import numpy as np
import logging
import modatm
from config import *
from tools import read_csv

logger = logging.getLogger(__name__)


class Data:
    def __init__(self, fname_star_data, fname_lines=None):
        try:
            self.star_data = read_csv(fname_star_data)
            self.star_data_fname = fname_star_data
        except:
            logger.error('Star data file not found.')
            return None

        if fname_lines:
            try:
                self.lines = read_csv(fname_lines, is_lines_file=True)
                self.lines_fname = fname_lines
            except:
                logger.error('Lines file not found or could not be read.')
                return None
        else:
            logger.warning('No lines data. Wont be able to MOOG.')

        logger.info('Data object successfully created.')


class Star:
    def __init__(self, name='Unnamed star',
                       teff=None, logg=None, feh=None, vt=None):
        self.name = name
        self.teff = teff
        self.logg = logg
        self.feh = feh
        self.vt = vt
        logger.info('Star object successfully created.')

    def parameters(self):
        print(self.teff)
        print(self.logg)
        print(self.feh)
        print(self.vt)

    def get_data_from(self, Data):
        #idx must correspond to a unique id; hence the [0][0]
        try:
            idx = np.where(Data.star_data['id'] == self.name)[0][0]
            logger.info("Star '"+self.name+"' found in data object.")
        except:
            logger.error("Star '"+self.name+"' not found in data object.")
            return None
        try:
            self.teff = Data.star_data['teff_out'][idx]
            self.logg = Data.star_data['logg_out'][idx]
            self.feh = Data.star_data['feh_out'][idx]
            self.err_teff = Data.star_data['err_teff_out'][idx]
            self.err_logg = Data.star_data['err_logg_out'][idx]
            self.err_feh = Data.star_data['err_feh_out'][idx]
            try:
                self.vt = Data.star_data['vt_out'][idx]
                self.err_vt = Data.star_data['err_vt_out'][idx]
            except:
                logger.warning('No vt_out for this star.')
        except:
            self.teff = Data.star_data['teff_in'][idx]
            self.logg = Data.star_data['logg_in'][idx]
            self.feh = Data.star_data['feh_in'][idx]
            try:
                self.vt = Data.star_data['vt_in'][idx]
            except:
                logger.warning('No vt_in for this star.')
            try:
                self.err_teff = Data.star_data['err_teff_in'][idx]
                self.err_logg = Data.star_data['err_logg_in'][idx]
                self.err_feh = Data.star_data['err_feh_in'][idx]
            except:
                logger.info('No errors in _in parameters.')

        logger.info('Attributes teff, logg, feh added to star object.')
        if hasattr(self, 'err_teff'):
            logger.info('Attributes err_teff, err_logg, err_feh added to star object.')

        additional_parameters = ['v', 'err_v', 'plx', 'err_plx', 'converged']
        msg = []
        for ap in additional_parameters:
            if ap in Data.star_data.keys():
                if Data.star_data[ap][idx]:
                    setattr(self, ap, Data.star_data[ap][idx])
                    msg.append(ap)
        if msg:
            logger.info('Additional attribute(s) '+','.join(msg)+\
                        ' added to star object.')

        # gets line data excluding cells with no ew:
        if hasattr(Data, 'lines'):
            idx = np.where(Data.lines[self.name] != None and\
                           Data.lines[self.name] > -100000)
            self.linelist = {'wavelength': Data.lines['wavelength'][idx],
                             'species': Data.lines['species'][idx],
                             'ep': Data.lines['ep'][idx],
                             'gf': Data.lines['gf'][idx],
                             'ew': Data.lines[self.name][idx]}
            logger.info('Attribute linelist added to star object.')
        else:
            logger.warning('There is no line data to attach to Star object.')

    def get_model_atmosphere(self, grid):
        try:
            self.teff
            self.logg
            self.feh
        except:
            logger.error('To create model atmosphere, star must have all '+
                         'three fundamental parameters: Teff, logg, and '+
                         '[Fe/H].')
            return None
        x = modatm.interpolate(self.teff, self.logg,
                               self.feh, grid)
        if x != None:
            self.model_atmosphere = x
            self.model_atmosphere_grid = grid
