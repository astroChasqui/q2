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
    def __init__(self, db='yy02.sql3', nsigma=5, window_len=13):
        self.get_isochrone_points_db = db
        self.get_isochrone_points_nsigma = nsigma
        self.smooth_window_len = window_len

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
      round(griddata(areas_left, pdf_x_left, 0.158),1)
    stats['lower_limit_2sigma'] = \
      round(griddata(areas_left, pdf_x_left, 0.022),1)
    stats['upper_limit_1sigma'] = \
      round(griddata(areas_right, pdf_x_right, 0.341),1)
    stats['upper_limit_2sigma'] = \
      round(griddata(areas_right, pdf_x_right, 0.477),1)

    return pdf_y, pdf_y_smooth, stats

def solve_one(Star, SolvePars):
    '''Calculates most likely parameters of Star using YY isochrone points
    '''
    ips = get_isochrone_points(Star,
                               SolvePars.get_isochrone_points_db,
                               SolvePars.get_isochrone_points_nsigma)
    logger.info('Using {0} Y2 isochrone points'.format(len(ips['age'])))
    Star.n_yy_points = len(ips['age'])
    ips_t = 10**ips['logt']
    prob = np.exp(-1*((ips_t-Star.teff)/(1.414214*Star.err_teff))**2)*\
           np.exp(-1*((ips['logg']-Star.logg)/(1.414214*Star.err_logg))**2)*\
           np.exp(-1*((ips['feh']-Star.feh)/(1.414214*Star.err_feh))**2)

    #age, mass, logl, mv

    ages = 0.1+np.arange(150)*0.1
    pdf_age_x = ages[np.logical_and(ages >= min(ips['age'])-0.2,
                                   ages <= max(ips['age'])+0.2)]
    pdf_age_y, pdf_age_y_smooth, Star.yyage = \
      pdf(pdf_age_x, ips, prob, 'age', SolvePars.smooth_window_len)

    masses = 0.4+np.arange(421)*0.01
    pdf_mass_x = masses[np.logical_and(masses >= min(ips['mass'])-0.02,
                                       masses <= max(ips['mass'])+0.02)]
    pdf_mass_y, pdf_mass_y_smooth, Star.yymass = \
      pdf(pdf_mass_x, ips, prob, 'mass', 1e10) #1e10=No smoothing

    plt.figure(figsize=(7, 4))
    plt.rc("axes", labelsize=15, titlesize=12)
    plt.rc("xtick", labelsize=14)
    plt.rc("ytick", labelsize=14)
    plt.rc("lines", markersize=10, markeredgewidth=2)
    plt.rc("lines", linewidth=2)
    plt.xlim([0,15])
    plt.xlabel('age (Gyr)')
    plt.ylabel('Probability density')
    #plt.plot(pdf_age_x, pdf_age_y, color='0.9')
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
    fig_name = Star.name.replace(' ', '_')+'_yyage'
    #plt.savefig(fig_name, bbox_inches='tight')
    plt.close()

    plt.figure(figsize=(7, 4))
    plt.rc("axes", labelsize=15, titlesize=12)
    plt.rc("xtick", labelsize=14)
    plt.rc("ytick", labelsize=14)
    plt.rc("lines", markersize=10, markeredgewidth=2)
    plt.rc("lines", linewidth=2)
    #plt.xlim([0,2])
    plt.xlabel('Mass ($M_\odot$)')
    plt.ylabel('Probability density')
    plt.plot(pdf_mass_x, pdf_mass_y, color='0.9')
    #k2 = np.logical_and(masses >= mass['lower_limit_2sigma'],
    #                    masses <= mass['upper_limit_2sigma'])
    #k1 = np.logical_and(masses >= mass['lower_limit_1sigma'],
    #                    masses <= mass['upper_limit_1sigma'])
    #plt.fill_between(masses[k2], 0 , probs[k2], color='0.8', hatch="/")
    #plt.fill_between(masses[k1], 0 , probs[k1], color='0.6', hatch="X")
    plt.plot([Star.yymass['most_probable'], Star.yymass['most_probable']],
             [0, max(pdf_mass_y_smooth)], 'g--')
    plt.plot(pdf_mass_x, pdf_mass_y_smooth, 'g')
    plt.text(plt.xlim()[1], 0.86*plt.ylim()[1], Star.name,
             horizontalalignment='right', size=16)
    fig_name = Star.name.replace(' ', '_')+'_yymass'
    plt.savefig(fig_name, bbox_inches='tight')
    plt.close()


    return

    ###

    plt.figure(figsize=(7, 5))
    plt.rc("axes", labelsize=15, titlesize=12)
    plt.rc("xtick", labelsize=14)
    plt.rc("ytick", labelsize=14)
    plt.rc("lines", markersize=10, markeredgewidth=2)
    plt.rc("lines", linewidth=3)
    plt.plot(ips_t, ips['logg'], 'y,')
    plt.xlim(Star.teff+6*Star.err_teff, Star.teff-6*Star.err_teff)
    plt.ylim(Star.logg+6*Star.err_logg, Star.logg-6*Star.err_logg)
    plt.xlabel('$T_\mathrm{eff}$ (K)')
    plt.ylabel('log g [cgs]')
    plt.errorbar(Star.teff, Star.logg, Star.err_logg, Star.err_teff, 'green')
    fig_name = Star.name.replace(' ', '_')+'_yyhr'
    plt.savefig(fig_name, bbox_inches='tight')
    plt.close()



    logls = -2.0+np.arange(501)*0.01 # logl between -2 and 3 in 0.01 dex steps
    logls = logls[np.logical_and(logls >= min(ips['logl'])-0.05,
                                 logls <= max(ips['logl'])+0.05)]

    plt.figure(figsize=(7, 4))
    plt.rc("axes", labelsize=15, titlesize=12)
    plt.rc("xtick", labelsize=14)
    plt.rc("ytick", labelsize=14)
    plt.rc("lines", markersize=10, markeredgewidth=2)
    plt.rc("lines", linewidth=2)
    #plt.xlim([0,2])
    plt.xlabel('$\log\,(L/L_\odot)$')
    plt.ylabel('Probability density')
    plt.plot([logl['most_probable'], logl['most_probable']],
             [0,max(probs)], 'g--')
    plt.plot(logls, probs, 'g')
    plt.text(0.1, 0.86*plt.ylim()[1], Star.name,
             horizontalalignment='right', size=16)
    fig_name = Star.name.replace(' ', '_')+'_yylogl'
    plt.savefig(fig_name, bbox_inches='tight')
    plt.close()

    #RADIUS
    #logR = 0.5* ( np.log10(mass['most_probable']) - Star.logg + 4.437 )
    #R = 10**logR
    #print (R)

    plt.figure(figsize=(7, 4))
    plt.rc("axes", labelsize=15, titlesize=12)
    plt.rc("xtick", labelsize=14)
    plt.rc("ytick", labelsize=14)
    plt.rc("lines", markersize=10, markeredgewidth=2)
    plt.rc("lines", linewidth=2)
    #plt.xlim([0,2])
    plt.xlabel('Mass ($M_\odot$)')
    plt.ylabel('Probability density')
    #k2 = np.logical_and(masses >= mass['lower_limit_2sigma'],
    #                    masses <= mass['upper_limit_2sigma'])
    #k1 = np.logical_and(masses >= mass['lower_limit_1sigma'],
    #                    masses <= mass['upper_limit_1sigma'])
    #plt.fill_between(masses[k2], 0 , probs[k2], color='0.8', hatch="/")
    #plt.fill_between(masses[k1], 0 , probs[k1], color='0.6', hatch="X")
    plt.plot([mass['most_probable'], mass['most_probable']],
             [0,max(probs)], 'g--')
    plt.plot(masses, probs, 'g')
    #plt.text(14.2, 0.86*plt.ylim()[1], Star.name,
    #         horizontalalignment='right', size=16)
    fig_name = Star.name.replace(' ', '_')+'_yymass'
    plt.savefig(fig_name, bbox_inches='tight')
    plt.close()

    #del(ips)
    #return None

    del(ips)


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
