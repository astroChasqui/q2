import sqlite3
import numpy as np
import logging
import matplotlib.pyplot as plt
from scipy.integrate import simps
from scipy.interpolate import griddata
from config import *
import os


logger = logging.getLogger(__name__)

class SolvePars:
    def __init__(self, db='yy02.sql3', nsigma=5,
                 window_len_age=13):
        self.get_isochrone_points_db = db
        self.get_isochrone_points_nsigma = nsigma
        self.smooth_window_len_age = window_len_age
        self.smooth_window_len_mass = 0
        self.smooth_window_len_logl = 0
        self.smooth_window_len_mv = 0
        self.smooth_window_len_r = 0

class PlotPars:
    def __init__(self, figure_format='png', directory=""):
        self.age_xlim = [0, 14]
        self.mass_xlim = None
        self.logl_xlim = None
        self.mv_xlim = None
        self.r_xlim = None
        self.directory = directory
        self.figure_format = figure_format

def pdf(pdf_x, ips, prob, par, smooth_window_len):
    '''Calculates a probability distribution function (PDF) for parameter par
    given the x-values for the PDF, the isochrone points ips, and their 
    probability. Return PDF and smoothed PDF (using smooth_window_len) if
    possible (otherwise returns two non-smoothed PDFs), as well as a stats
    dictionary with mean, std, most probable value, etc.
    '''
    dx = 0.5*(pdf_x[1] - pdf_x[0])
    pdf_y = []
    for x in pdf_x:
        pdf_y.append(sum(prob[np.logical_and(ips[par] >= x-dx,
                                             ips[par] <  x+dx)]))
    pdf_y = np.array(pdf_y)
    pdf_y = pdf_y/simps(pdf_y, pdf_x)

    try:
        pdf_y_smooth = smooth(pdf_y, smooth_window_len)
        pdf_y_smooth = pdf_y_smooth/simps(pdf_y_smooth, pdf_x)
    except:
        pdf_y_smooth = pdf_y
        logger.warning('Unable to smooth '+par+' PDF.')

    stats = get_stats(pdf_x, pdf_y_smooth)

    return pdf_y, pdf_y_smooth, stats

def get_stats(pdf_x, pdf_y_smooth):
    stats = {}
    stats['most_probable'] = \
      np.mean(np.array(pdf_x)[pdf_y_smooth == max(pdf_y_smooth)])
    stats['mean'] = simps(pdf_y_smooth*pdf_x, pdf_x)
    stats['std'] = np.sqrt(simps(pdf_y_smooth*(pdf_x-stats['mean'])**2,\
                                 pdf_x))

    k = pdf_x <= stats['most_probable']
    pdf_y_left = 0.5*pdf_y_smooth[k]/simps(pdf_y_smooth[k], pdf_x[k])
    pdf_x_left = pdf_x[k]
    areas_left = []
    for x in pdf_x_left:
        areas_left.append(simps(pdf_y_left[pdf_x_left <= x],
                                pdf_x_left[pdf_x_left <= x]))
    areas_left = np.array(areas_left)

    k = pdf_x >= stats['most_probable']
    pdf_y_right = 0.5*pdf_y_smooth[k]/simps(pdf_y_smooth[k], pdf_x[k])
    pdf_x_right = pdf_x[k]
    areas_right = []
    for x in pdf_x_right:
        areas_right.append(simps(pdf_y_right[pdf_x_right <= x],
                                 pdf_x_right[pdf_x_right <= x]))
    areas_right = np.array(areas_right)

    stats['lower_limit_1sigma'] = \
      np.mean(griddata(areas_left, pdf_x_left, 0.158))
    stats['lower_limit_2sigma'] = \
      np.mean(griddata(areas_left, pdf_x_left, 0.022))
    stats['upper_limit_1sigma'] = \
      np.mean(griddata(areas_right, pdf_x_right, 0.341))
    stats['upper_limit_2sigma'] = \
      np.mean(griddata(areas_right, pdf_x_right, 0.477))

    return stats

def solve_one(Star, SolvePars, PlotPars):
    '''Calculates most likely parameters of Star using YY isochrone points
    '''
    ips = get_isochrone_points(Star,
                               SolvePars.get_isochrone_points_db,
                               SolvePars.get_isochrone_points_nsigma)
    logger.info('Using {0} Y2 isochrone points'.format(len(ips['age'])))
    Star.n_yy_points = len(ips['age'])
    ips['t'] = 10**ips['logt']
    ips['r'] = 10**(0.5*(np.log10(ips['mass'])-ips['logg']+4.437))
    prob = np.exp(-1*((ips['t']-Star.teff)/(1.414214*Star.err_teff))**2)*\
           np.exp(-1*((ips['logg']-Star.logg)/(1.414214*Star.err_logg))**2)*\
           np.exp(-1*((ips['feh']-Star.feh)/(1.414214*Star.err_feh))**2)

    #age
    ages = 0.1+np.arange(150)*0.1
    pdf_age_x = ages[np.logical_and(ages >= min(ips['age'])-0.2,
                                   ages <= max(ips['age'])+0.2)]
    pdf_age_y, pdf_age_y_smooth, Star.yyage = \
      pdf(pdf_age_x, ips, prob, 'age', SolvePars.smooth_window_len_age)
    Star.yypdf_age = {'x': pdf_age_x, 'y': pdf_age_y, 'ys': pdf_age_y_smooth}

    #mass
    masses = 0.4+np.arange(211)*0.01
    pdf_mass_x = masses[np.logical_and(masses >= min(ips['mass'])-0.02,
                                       masses <= max(ips['mass'])+0.02)]
    pdf_mass_y, pdf_mass_y_smooth, Star.yymass = \
      pdf(pdf_mass_x, ips, prob, 'mass', SolvePars.smooth_window_len_mass)

    #luminosity
    logls = -1.0+np.arange(401)*0.01
    pdf_logl_x = logls[np.logical_and(logls >= min(ips['logl'])-0.02,
                                      logls <= max(ips['logl'])+0.02)]
    pdf_logl_y, pdf_logl_y_smooth, Star.yylogl = \
      pdf(pdf_logl_x, ips, prob, 'logl', SolvePars.smooth_window_len_logl)

    #absolute magnitude
    mvs = -3.0+np.arange(1601)*0.01
    pdf_mv_x = mvs[np.logical_and(mvs >= min(ips['mv'])-0.02,
                                  mvs <= max(ips['mv'])+0.02)]
    pdf_mv_y, pdf_mv_y_smooth, Star.yymv = \
      pdf(pdf_mv_x, ips, prob, 'mv', SolvePars.smooth_window_len_mv)

    #radius
    rs = 0.4+np.arange(211)*0.01
    pdf_r_x = rs[np.logical_and(rs >= min(ips['r'])-0.02,
                                rs <= max(ips['r'])+0.02)]
    pdf_r_y, pdf_r_y_smooth, Star.yyr = \
      pdf(pdf_r_x, ips, prob, 'r', SolvePars.smooth_window_len_r)

    if not os.path.exists(PlotPars.directory) and PlotPars.directory != "":
        os.mkdir(PlotPars.directory)

    plt.figure(figsize=(7, 4))
    plt.rc("axes", labelsize=15, titlesize=12)
    plt.rc("xtick", labelsize=14)
    plt.rc("ytick", labelsize=14)
    plt.rc("lines", markersize=10, markeredgewidth=2)
    plt.rc("lines", linewidth=2)
    plt.rc("xtick.major", size=6, width=1)
    plt.rc("ytick.major", size=6, width=1)
    plt.xlim([0,15])
    plt.xlabel('Age (Gyr)')
    plt.ylabel('Probability density')
    k2 = np.logical_and(pdf_age_x >= Star.yyage['lower_limit_2sigma'],
                        pdf_age_x <= Star.yyage['upper_limit_2sigma'])
    k1 = np.logical_and(pdf_age_x >= Star.yyage['lower_limit_1sigma'],
                        pdf_age_x <= Star.yyage['upper_limit_1sigma'])
    plt.fill_between(pdf_age_x[k2], 0 , pdf_age_y_smooth[k2],
                     color='0.8', hatch="/")
    plt.fill_between(pdf_age_x[k1], 0 , pdf_age_y_smooth[k1],
                     color='0.6', hatch="X")
    plt.plot([Star.yyage['most_probable'], Star.yyage['most_probable']],
             [0, max(pdf_age_y_smooth)], 'g--')
    plt.plot(pdf_age_x, pdf_age_y_smooth, 'g')
    plt.text(14.2, 0.86*plt.ylim()[1], Star.name,
             horizontalalignment='right', size=16)
    fig_name = os.path.join(PlotPars.directory,
                            Star.name.replace(' ', '_')+'_yyage')
    plt.savefig(fig_name+'.'+PlotPars.figure_format, bbox_inches='tight')
    plt.close()

    plt.figure(figsize=(6, 12))
    plt.rc("axes", labelsize=15, titlesize=12)
    plt.rc("xtick", labelsize=14)
    plt.rc("ytick", labelsize=14)
    plt.rc("lines", markersize=10, markeredgewidth=2)
    plt.rc("lines", linewidth=2)
    plt.rc("xtick.major", size=6, width=1)
    plt.rc("ytick.major", size=6, width=1)
    plt.subplots_adjust(hspace=0.4)

    for panel in np.arange(5)+1:
        ax = plt.subplot(5,1,panel)
        ax.get_yaxis().set_visible(False)
        if panel == 1:
            pdf_x, pdf_y, pdf_y_smooth = \
              pdf_age_x, pdf_age_y, pdf_age_y_smooth
            par = Star.yyage
            ax.set_xlabel('Age (Gyr)')
            if PlotPars.age_xlim:
                ax.set_xlim(PlotPars.age_xlim)
        if panel == 2:
            pdf_x, pdf_y, pdf_y_smooth = \
              pdf_mass_x, pdf_mass_y, pdf_mass_y_smooth
            par = Star.yymass
            ax.set_xlabel('Mass ($M_\odot$)')
            if PlotPars.mass_xlim:
                ax.set_xlim(PlotPars.mass_xlim)
        if panel == 3:
            pdf_x, pdf_y, pdf_y_smooth = \
              pdf_logl_x, pdf_logl_y, pdf_logl_y_smooth
            par = Star.yylogl
            ax.set_xlabel('$\log\,(L/L_\odot)$')
            if PlotPars.logl_xlim:
                ax.set_xlim(PlotPars.logl_xlim)
        if panel == 4:
            pdf_x, pdf_y, pdf_y_smooth = \
              pdf_mv_x, pdf_mv_y, pdf_mv_y_smooth
            par = Star.yymv
            ax.set_xlabel('$M_V$')
            if PlotPars.mv_xlim:
                ax.set_xlim(PlotPars.mv_xlim)
        if panel == 5:
            pdf_x, pdf_y, pdf_y_smooth = \
              pdf_r_x, pdf_r_y, pdf_r_y_smooth
            par = Star.yyr
            ax.set_xlabel('Radius ($R_\odot$)')
            if PlotPars.r_xlim:
                ax.set_xlim(PlotPars.r_xlim)
        ax.plot(pdf_x, pdf_y, color='0.9')
        ax.plot([par['most_probable'], par['most_probable']],
                [0, max(pdf_y_smooth)], 'g--')
        ax.plot(pdf_x, pdf_y_smooth, 'g')

    fig_name = os.path.join(PlotPars.directory,
                            Star.name.replace(' ', '_')+'_yypars')
    plt.savefig(fig_name+'.'+PlotPars.figure_format, bbox_inches='tight')
    plt.close()

def get_isochrone_points(Star, db, nsigma):
    '''Looks in the db database for isochrone points within nsigma from
    the mean parameters of the Star and returns those values in a dict.
    '''
    conn = sqlite3.connect(os.path.join(ISOCHRONES_PATH, db))
    conn.row_factory = sqlite3.Row
    logtm = np.log10(Star.teff-nsigma*Star.err_teff)
    logtp = np.log10(Star.teff+nsigma*Star.err_teff)
    c = conn.cursor()
    x = c.execute('SELECT feh, age, mass, logt, logl, logg, mv ' +\
                  'FROM  fa, mtlgv ON fa.fa = mtlgv.fa WHERE '   +\
                  'logt >= ? AND logt <= ? AND '   +\
                  'feh  >= ? AND feh  <= ? AND '   +\
                  'logg >= ? AND logg <= ? ',
                  (logtm, logtp,
                   Star.feh-nsigma*Star.err_feh,
                   Star.feh+nsigma*Star.err_feh,
                   Star.logg-nsigma*Star.err_logg,
                   Star.logg+nsigma*Star.err_logg)
                 )
    feh, age = [], []
    mass, logt, logl, logg, mv = [], [], [], [], []
    for xx in x.fetchall():
        feh.append(xx['feh'])
        age.append(xx['age'])
        mass.append(xx['mass'])
        logt.append(xx['logt'])
        logl.append(xx['logl'])
        logg.append(xx['logg'])
        mv.append(xx['mv'])
    conn.close()
    return {'feh' : np.array(feh)+0.04,
            'age' : np.array(age),
            'mass': np.array(mass),
            'logt': np.array(logt),
            'logl': np.array(logl),
            'logg': np.array(logg),
            'mv'  : np.array(mv)
            }

def smooth(x, window_len=11, window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."

    if window_len<3:
        return x

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"

    s=np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    return y[(window_len/2):-(window_len/2)]
