import sqlite3
import numpy as np
import logging
import matplotlib.pyplot as plt
from scipy.integrate import simps
from scipy.interpolate import griddata
from config import *
import os


logger = logging.getLogger(__name__)

class YYagePars:
    def __init__(self, db='yy02.sql3', nsigma=5, window_len=11):
        self.get_isochrone_points_db = db
        self.get_isochrone_points_nsigma = nsigma
        self.smooth_window_len = window_len

def calculate_age(Star, YYagePars):
    '''Calculates most likely age of Star using YY isochrone points
    '''
    ips = get_isochrone_points(Star,
                               YYagePars.get_isochrone_points_db,
                               YYagePars.get_isochrone_points_nsigma)
    logger.info('Using {0} Y2 isochrone points'.format(len(ips['age'])))
    ips_t = 10**ips['logt']
    prob = np.exp(-1*((ips_t-Star.teff)/(1.414214*Star.err_teff))**2)*\
           np.exp(-1*((ips['logg']-Star.logg)/(1.414214*Star.err_logg))**2)*\
           np.exp(-1*((ips['feh']-Star.feh)/(1.414214*Star.err_feh))**2)

    ages = np.array(sorted(set(ips['age'])))
    probs = []
    for age in ages:
        probs.append(sum(prob[ips['age'] == age]))
    probs = np.array(probs)
    del(ips)
    try:
        probs_smooth = smooth(probs, YYagePars.smooth_window_len)
    except:
        logger.warning('Unable to smooth age PDF')
        return None
    age = {}
    age['age'] = round(np.array(ages)[probs_smooth == max(probs_smooth)],2)

    k = ages <= age['age']
    probs_sn_left = 0.5*probs_smooth[k]/simps(probs_smooth[k],ages[k])
    ages_left = ages[k]
    areas_left = []
    for agel in ages_left:
        areas_left.append(simps(probs_sn_left[ages_left <= agel],
                                ages_left[ages_left <= agel]))
    areas_left = np.array(areas_left)

    k = ages >= age['age']
    probs_sn_right = 0.5*probs_smooth[k]/simps(probs_smooth[k],ages[k])
    ages_right = ages[k]
    areas_right = []
    for ager in ages_right:
        areas_right.append(simps(probs_sn_right[ages_right <= ager],
                                 ages_right[ages_right <= ager]))
    areas_right = np.array(areas_right)

    age['lower_limit_1sigma'] = \
      round(griddata(areas_left, ages_left, 0.158),1)
    age['upper_limit_1sigma'] = \
      round(griddata(areas_right, ages_right, 0.341),1)
    age['lower_limit_2sigma'] = \
      round(griddata(areas_left, ages_left, 0.022),1)
    age['upper_limit_2sigma'] = \
      round(griddata(areas_right, ages_right, 0.477),1)

    Star.age = age
    print(age)

    plt.figure(figsize=(7, 4))
    plt.rc("axes", labelsize=15, titlesize=12)
    plt.rc("xtick", labelsize=14)
    plt.rc("ytick", labelsize=14)
    plt.rc("lines", markersize=10, markeredgewidth=2)
    plt.rc("lines", linewidth=3)
    plt.xlim([0,15])
    plt.xlabel('age (Gyr)')
    plt.ylabel('Relative probability')
    #plt.plot(ages, probs, color='0.9')
    k2 = np.logical_and(ages >= age['lower_limit_2sigma'],
                        ages <= age['upper_limit_2sigma'])
    k1 = np.logical_and(ages >= age['lower_limit_1sigma'],
                        ages <= age['upper_limit_1sigma'])
    plt.fill_between(ages[k2], 0 , probs_smooth[k2], color='0.9', hatch="/")
    plt.fill_between(ages[k1], 0 , probs_smooth[k1], color='0.7', hatch="X")
    plt.plot([age['age'], age['age']], [0,max(probs_smooth)], 'g--')
    plt.plot(ages, probs_smooth, 'g')
    plt.text(14.2, 0.86*plt.ylim()[1], Star.name,
             horizontalalignment='right', size=16)
    plt.savefig(Star.name+"_age.eps", bbox_inches='tight')
    plt.close()


    #print(min(ips['mass']))
    #print(max(ips['mass']))

    '''
    plt.figure(figsize=(7, 5))
    plt.rc("axes", labelsize=15, titlesize=12)
    plt.rc("xtick", labelsize=14)
    plt.rc("ytick", labelsize=14)
    plt.rc("lines", markersize=10, markeredgewidth=2)
    plt.rc("lines", linewidth=3)
    plt.plot(ips_t, ips['logg'], 'y,')
    plt.xlim(Star.teff+4*Star.err_teff, Star.teff-4*Star.err_teff)
    plt.ylim(Star.logg+4*Star.err_logg, Star.logg-4*Star.err_logg)
    plt.xlabel('$T_\mathrm{eff}$ (K)')
    plt.ylabel('log g [cgs]')
    plt.errorbar(Star.teff, Star.logg, Star.err_logg, Star.err_teff, 'green')
    plt.savefig(Star.name, bbox_inches='tight')
    plt.close()
    '''

def get_isochrone_points(Star, db, nsigma):
    '''Looks in a database db for isochrone points within nsigma from
    the mean parameters of the Star and returns those values in a dict.
    '''
    path = ISOCHRONES_PATH
    conn = sqlite3.connect(os.path.join(path, db))
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
