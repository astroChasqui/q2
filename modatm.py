import matplotlib.pyplot as plt
import logging
import tools
import os
from config import *
from scipy.integrate import simps
import numpy as np
from scipy.interpolate import griddata

logger = logging.getLogger(__name__)


def get_from_file(teff, logg, feh, grid):
    """Extracts a Kurucz or MARCS model atmosphere

    If model file does not exist, the program returns None, None.
    Output is model, tauRoss. The latter is necessary for proper
    model interpolation.
    """

    # location of model files
    if grid == 'marcs':
        path = os.path.join(MODATM_PATH, 'marcs')
    else:
        path = os.path.join(MODATM_PATH, 'kurucz')

    # determine filename from input parameters
    logg = float(logg)
    g = 'g'+str(logg).replace('.', '')
    feh = float(feh)
    if feh >= 0:
        f='p'+str(abs(feh)).replace('.', '')
    else:
        f='m'+str(abs(feh)).replace('.', '')
    file_name = 't'+str(teff)+g+f+'.'+grid

    # if file doesn't exist return None
    try:
        with open(os.path.join(path, file_name)): pass
    except IOError:
        logger.warning('Model file not found: '+file_name)
        return None, None

    logger.info('Model found: '+file_name)

    with open(os.path.join(path, file_name), 'r') as f:
        xf = f.readlines()
    keys = ['RHOX', 'T', 'P', 'XNE', 'ABROSS']
    x = {'RHOX': [], 'T': [], 'P': [], 'XNE': [], 'ABROSS': []}
    for xfi in xf:
        for key, xfij in zip(keys,xfi.split(",")):
            x[key].append(float(xfij.strip("\n")))
    for key in keys:
        x[key] = np.array(x[key])

    # calculate tauRoss scale
    tau0 = x['RHOX'][0]*x['ABROSS'][0]
    xdtau, ydtau, tau = [], [], []
    for rhox, abross in zip(x['RHOX'], x['ABROSS']):
        xdtau.append(rhox)
        ydtau.append(abross)
        tau.append(simps(ydtau, xdtau))
    tau[0] = tau0
    tau = np.array(tau)

    return x, tau


def interpolate(teff, logg, feh, grid):
    """Interpolates within Kurucz or MARCS model atmosphere grids

    Returns original (not interpolated) model if it exists (grid node).
    Returns None if for some reason the program cannot interpolate.
    """

    model, tau = get_from_file(teff, logg, feh, grid)
    if model != None:
        return model

    avail_teff = range(3500, 7501, 250)
    avail_logg = [i*0.5 for i in range(11)]
    avail_feh_over = [+1.0,+0.5,+0.2,+0.1,+0.0,-0.1,-0.2,-0.3,-0.5,
                      -1.0,-1.5,-2.0,-2.5,-3.0,-3.5,-4.0,-4.5,-5.0]
    avail_feh_nover = [+0.5,+0.0,-0.5,-1.0,-1.5,-2.0,-2.5]
    avail_feh_odfnew = [+0.5,+0.2,+0.0,-0.5,-1.0,-1.5,-2.0,-2.5]
    avail_feh_aodfnew = [+0.5,+0.0,-0.5,-1.0,-1.5,-2.0,-2.5,-4.0]
    avail_feh_marcs = [+0.5,+0.25,+0.0,-0.25,-0.5,-0.75,
                       -1.0,-1.5,-2.0,-2.5,-3.0,-4.0]
    if grid == 'over':
        avail_feh = avail_feh_over
    if grid == 'nover':
        avail_feh = avail_feh_nover
    if grid == 'odfnew':
        avail_feh = avail_feh_odfnew
    if grid == 'aodfnew':
        avail_feh = avail_feh_aodfnew
    if grid == 'marcs':
        avail_feh = avail_feh_marcs
        avail_teff = [3500, 3600, 3700, 3800, 3900]
        [avail_teff.append(t) for t in range(4000, 7001, 250)]

    try:
        avail_feh
    except NameError:
        logger.error('The type of model requested is not available.')
        return None

    logger.info('Interpolating in the '+grid+' model grid:')

    #dt = [abs(t-teff) for t in avail_teff]
    #dt_idx = sorted(range(len(dt)), key=lambda k: dt[k])
    #tx = sorted([avail_teff[dt_idx[0]], avail_teff[dt_idx[1]]])

    dt = [t-teff for t in avail_teff]
    for idx, dtx in enumerate(dt[:-1]):
        if dt[idx]*dt[idx+1] <= 0:
            break
    tx = avail_teff[idx], avail_teff[idx+1]

    logger.info('Searching for '+str(teff)+' in range: '+str(tx))
    if teff < tx[0] or teff > tx[1]:
        logger.error('Cannot interpolate in Teff.')
        return None

    dg = [abs(g-logg) for g in avail_logg]
    dg_idx = sorted(range(len(dg)), key=lambda k: dg[k])
    gx = sorted([avail_logg[dg_idx[0]], avail_logg[dg_idx[1]]])
    logger.info('Searching for '+str(logg)+' in range: '+str(gx))
    if logg < gx[0] or logg > gx[1]:
        logger.error('Cannot interpolate in log g.')
        return None

    df = [abs(f-feh) for f in avail_feh]
    df_idx = sorted(range(len(df)), key=lambda k: df[k])
    if (df_idx[0] == len(df_idx)-1 and feh < avail_feh[df_idx[0]]) or \
       (df_idx[0] == 0 and feh > avail_feh[df_idx[0]]):
        logger.error('Too low or high feh')
        return None
    else:
        if (feh > 0 and feh <= avail_feh[df_idx[0]]) or \
           (feh < 0 and feh < avail_feh[df_idx[0]]):
            fx = [avail_feh[df_idx[0]+1], avail_feh[df_idx[0]]]
        else:
            fx = [avail_feh[df_idx[0]], avail_feh[df_idx[0]-1]]
    logger.info('Searching for '+str(feh)+' in range: '+str(fx))
    if feh < fx[0] or feh > fx[1]:
        logger.error('Cannot interpolate in [Fe/H].')
        return None

    logger.info('Reading available models:')
    T0G0F0, tau000 = get_from_file(tx[0], gx[0], fx[0], grid)
    T0G0F1, tau001 = get_from_file(tx[0], gx[0], fx[1], grid)
    T0G1F0, tau010 = get_from_file(tx[0], gx[1], fx[0], grid)
    T0G1F1, tau011 = get_from_file(tx[0], gx[1], fx[1], grid)
    T1G0F0, tau100 = get_from_file(tx[1], gx[0], fx[0], grid)
    T1G0F1, tau101 = get_from_file(tx[1], gx[0], fx[1], grid)
    T1G1F0, tau110 = get_from_file(tx[1], gx[1], fx[0], grid)
    T1G1F1, tau111 = get_from_file(tx[1], gx[1], fx[1], grid)

    tau_min = max(tau000[0], tau001[0], tau010[0], tau011[0],
                  tau100[0], tau101[0], tau110[0], tau111[0] )
    n = len(tau000) - 1
    tau_max = min(tau000[n], tau001[n], tau010[n], tau011[n],
                  tau100[n], tau101[n], tau110[n], tau111[n] )
    tau = tau000[(tau000 >= tau_min) & (tau000 <= tau_max)]
    tau_new = griddata(np.array(range(len(tau))), tau,
                       np.array(range(n+1))*np.array(float(len(tau)-1))/n)
    #for name in T0G0F0.dtype.names:
    for name in T0G0F0.keys():
        T0G0F0[name] = griddata(tau000, T0G0F0[name], tau_new)
        T0G0F1[name] = griddata(tau001, T0G0F1[name], tau_new)
        T0G1F0[name] = griddata(tau010, T0G1F0[name], tau_new)
        T0G1F1[name] = griddata(tau011, T0G1F1[name], tau_new)
        T1G0F0[name] = griddata(tau100, T1G0F0[name], tau_new)
        T1G0F1[name] = griddata(tau101, T1G0F1[name], tau_new)
        T1G1F0[name] = griddata(tau110, T1G1F0[name], tau_new)
        T1G1F1[name] = griddata(tau111, T1G1F1[name], tau_new)
        T0G0F0[name][0] = T0G0F0[name][1]*0.999
        T0G0F1[name][0] = T0G0F1[name][1]*0.999
        T0G1F0[name][0] = T0G1F0[name][1]*0.999
        T0G1F1[name][0] = T0G1F1[name][1]*0.999
        T1G0F0[name][0] = T1G0F0[name][1]*0.999
        T1G0F1[name][0] = T1G0F1[name][1]*0.999
        T1G1F0[name][0] = T1G1F0[name][1]*0.999
        T1G1F1[name][0] = T1G1F1[name][1]*0.999

    #feh-interpolate
    T0G0 = tools.linterp(T0G0F0, T0G0F1, fx[0], fx[1], feh)
    T0G1 = tools.linterp(T0G1F0, T0G1F1, fx[0], fx[1], feh)
    T1G0 = tools.linterp(T1G0F0, T1G0F1, fx[0], fx[1], feh)
    T1G1 = tools.linterp(T1G1F0, T1G1F1, fx[0], fx[1], feh)
    #logg-interpolate
    T0 = tools.linterp(T0G0, T0G1, gx[0], gx[1], logg)
    T1 = tools.linterp(T1G0, T1G1, gx[0], gx[1], logg)
    #teff-interpolate
    model = tools.linterp(T0, T1, tx[0], tx[1], teff)

    logger.info('Model interpolation was successful.')
    return model
