from astropy.io import ascii
import numpy as np
import logging
import modatm
from config import *

logger = logging.getLogger(__name__)


class Data:
    def __init__(self, fname_star_data, fname_lines=None):
        # need checks: star data must have id, teff_in, etc
        # also lines must have wave, ew, gf, etc.
        try:
            self.star_data = ascii.read(fname_star_data,
                                        fill_values=[('', '-9999')])
            self.star_data_fname = fname_star_data
        except:
            logger.error('Star data file not found.')
            return None
        if fname_lines:
            try:
                self.lines = ascii.read(fname_lines,
                                            fill_values=[('', '-9999')])
                self.lines_fname = fname_lines
            except:
                logger.error('Lines file not found.')
                return None
        else:
            logger.warning('No lines data. Wont be able to MOOG.')

        logger.info('Data object successfully created.')


class Star:
    def __init__(self, name='Unnamed star'):
        self.name = name
        logger.info('Star object successfully created.')

    def get_data_from(self, Data):
        idx = np.where(Data.star_data['id'] == self.name)
        try:
            idx[0][0]
            logger.info("Star '"+self.name+"' found in data object.")
        except:
            logger.error("Star '"+self.name+"' not found in data object.")
            return None
        try:
            self.teff = Data.star_data['teff_out'][idx[0][0]]
            self.logg = Data.star_data['logg_out'][idx[0][0]]
            self.feh = Data.star_data['feh_out'][idx[0][0]]
            self.err_teff = Data.star_data['err_teff_out'][idx[0][0]]
            self.err_logg = Data.star_data['err_logg_out'][idx[0][0]]
            self.err_feh = Data.star_data['err_feh_out'][idx[0][0]]
            try:
                self.vt = Data.star_data['vt_out'][idx[0][0]]
                self.vt = Data.star_data['err_vt_out'][idx[0][0]]
            except:
                logger.warning('No vt_out for this star.')
        except:
            self.teff = Data.star_data['teff_in'][idx[0][0]]
            self.logg = Data.star_data['logg_in'][idx[0][0]]
            self.feh = Data.star_data['feh_in'][idx[0][0]]
            try:
                self.vt = Data.star_data['vt_in'][idx[0][0]]
            except:
                logger.warning('No vt_in for this star.')
            try:
                self.err_teff = Data.star_data['err_teff_in'][idx[0][0]]
                self.err_logg = Data.star_data['err_logg_in'][idx[0][0]]
                self.err_feh = Data.star_data['err_feh_in'][idx[0][0]]
            except:
                logger.info('No errors in _in parameters.')

        logger.info('Attributes teff, logg, feh added to star object.')
        if hasattr(self, 'err_teff'):
            logger.info('Attributes err_teff, err_logg, err_feh added to star object.')

        additional_parameters = ['v', 'err_v', 'plx', 'err_plx', 'converged']
        msg = []
        for ap in additional_parameters:
            if ap in Data.star_data.keys():
                if Data.star_data[ap][idx[0][0]]:
                    setattr(self, ap, Data.star_data[ap][idx[0][0]])
                    msg.append(ap)
        if msg:
            logger.info('Additional attribute(s) '+','.join(msg)+\
                        ' added to star object.')

        # gets ews; excludes cells with no ew:
        if hasattr(Data, 'lines'):
            idx = np.where(((Data.lines['wavelength'] > -10000) &
                            (Data.lines[self.name] != '')))
            self.linelist = {'wavelength': Data.lines[idx]['wavelength'],
                             'species': Data.lines[idx]['species'],
                             'ep': Data.lines[idx]['ep'],
                             'gf': Data.lines[idx]['gf'],
                             'ew': Data.lines[idx][self.name]}
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
