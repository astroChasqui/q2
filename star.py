import numpy as np
import logging
import modatm
from config import *
from tools import read_csv

logger = logging.getLogger(__name__)


class Data:
    """q2 Data objects contain the information from the input star
    and lines CSV files as attributes 'star_data' and 'lines'.
    """
    def __init__(self, fname_star_data, fname_lines=None):
        """Takes as input a star data file (CSV, required) and a line-list
        data file (CSV, optional) to create a q2 data object.
        """
        try:
            self.star_data = read_csv(fname_star_data, file_type='stars')
            self.star_data_fname = fname_star_data
            if not self.star_data:
                logger.error('Star data file not read. Data.star_data '+\
                             'attribute set to None.')
        except:
            self.star_data = None
            self.star_data_fname = None
            logger.error('Star data file not found.')

        if fname_lines:
            try:
                self.lines = read_csv(fname_lines, file_type='lines')
                self.lines_fname = fname_lines
                if not self.lines:
                    logger.error('Lines data file not read. Data.lines '+\
                                 'attribute set to None.')
            except:
                self.lines = None
                self.lines_fname = None
                logger.error('Lines file not found.')
        else:
            self.lines = None
            self.lines_fname = None
            logger.warning('No lines data. Wont be able to MOOG.')

        if self.star_data:
            logger.info('Data object created with star_data attribute.')
        if self.lines:
            logger.info('lines_data attribute added to Data object.')

    def __repr__(self):
        if self.star_data:
            nstars = len(self.star_data['id'])
        else:
            nstars = 0
        if self.lines:
            nlines = len(np.where(self.lines['wavelength'] > 0)[0])
        else:
            nlines = 0
        return "Data object built from:\n"\
               "  stars file = {0} ({1} stars)\n"\
               "  lines file = {2} ({3} lines)".\
               format(self.star_data_fname, nstars, self.lines_fname, nlines)


class Star:
    """q2 Star objects contain information about a star (e.g., Teff, logg,
    etc). This information can be grabbed from a q2.Data object using the
    'get_data_from' method. If the Star object has parameters known, a
    model atmosphere can be computed using the 'get_model_atmosphere'
    method. This will attach a 'model_atmosphere' attribute to the Star
    object.
    """
    def __init__(self, name='Unnamed star',
                       teff=None, logg=None, feh=None, vt=None):
        self.name = name
        self.teff = teff
        self.logg = logg
        self.feh = feh
        self.vt = vt
        logger.info('Star object successfully created.')

    def __repr__(self):
        if hasattr(self, 'linelist'):
            nlines = len(self.linelist['wavelength'])
            species = ','.join([str(sp) for sp in\
                      set(self.linelist['species'])])
        else:
            nlines = species = None
        return "Star object named '{0}':\n"\
               "  Teff (K) = {1}, logg [cgs] = {2}, [Fe/H] = {3}, "\
                 "vt (km/s) = {4}\n"\
               "  Spectral lines = {5} (species: {6})".\
               format(self.name, self.teff, self.logg, self.feh, self.vt,
                      nlines, species)

    def get_data_from(self, Data):
        """If the Star object has a name that matches one of the id's in
        a Data object, the information from Data will be given to Star.
        """
        #idx must correspond to a unique id; hence the [0][0]
        try:
            idx = np.where(Data.star_data['id'] == self.name)[0][0]
            logger.info("Star '"+self.name+"' found in data object.")
        except:
            logger.error("Star '"+self.name+"' not found in data object.")
            return None

        parameters = ['teff', 'err_teff', 'logg', 'err_logg',
                      'feh', 'err_feh', 'vt', 'err_vt',
                      'v', 'err_v', 'plx', 'err_plx', 'converged']
        msg = []
        for par in parameters:
            if par in Data.star_data.keys():
                if Data.star_data[par][idx] != None:
                    setattr(self, par, Data.star_data[par][idx])
                    msg.append(par)
        if msg:
            logger.info('Attribute(s) '+','.join(msg)+\
                        ' added to star object.')

        # gets line data excluding cells with no ew:
        #if hasattr(Data, 'lines'):
        if Data.lines:
            idx = np.where(Data.lines[self.name] > 0)
            self.linelist = {'wavelength': Data.lines['wavelength'][idx],
                             'species': Data.lines['species'][idx],
                             'ep': Data.lines['ep'][idx],
                             'gf': Data.lines['gf'][idx],
                             'ew': Data.lines[self.name][idx]}
            logger.info('Attribute linelist added to star object.')
        else:
            logger.warning('There is no line data to attach to Star object.')

    def get_model_atmosphere(self, grid='odfnew'):
        """If teff, logg, and feh are set attributes for a Star object,
        a model atmosphere will be interpolated from one of the
        available grids: 'odfnew' (default), 'aodfnew', 'over', 'nover'
        (all Kurucz), or 'marcs'.
        """
        if self.teff == None or self.logg == None or self.feh == None:
            logger.error('To create model atmosphere, star must have all '+
                         'three fundamental parameters: Teff, logg, and '+
                         '[Fe/H].')
            return None
        x = modatm.interpolate(self.teff, self.logg,
                               self.feh, grid)
        if x != None:
            self.model_atmosphere = x
            self.model_atmosphere_grid = grid
