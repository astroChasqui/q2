import numpy as np
from numpy import matrix
from star import Star
import specpars
import logging

logger = logging.getLogger(__name__)


def error_one(Star_in, SolvePars, Ref=object):

    s = Star()
    s.__dict__ = Star_in.__dict__.copy()
    try:
        Ref.get_model_atmosphere(SolvePars.grid)
        logger.info('Relative abundances')
    except:
        logger.info('Absolute abundances')

    dteff = 20
    dvt = 0.02
    dlogg = 0.02

    s.teff = s.teff + dteff
    s.get_model_atmosphere(SolvePars.grid)
    specpars.iron_stats(s, Ref=Ref)
    fx1  = s.iron_stats['slope_ep']
    efx1 = s.iron_stats['err_slope_ep']
    gx1  = s.iron_stats['slope_rew']
    egx1 = s.iron_stats['err_slope_rew']
    hx1  = s.iron_stats['afe1'] - s.iron_stats['afe2']
    ehx1 = np.sqrt(s.iron_stats['err_afe1']**2+
                   s.iron_stats['err_afe2']**2)/\
           np.sqrt(s.iron_stats['nfe1']+s.iron_stats['nfe2'])
    s.teff = s.teff - 2*dteff
    s.get_model_atmosphere(SolvePars.grid)
    specpars.iron_stats(s, Ref=Ref)
    fx2  = s.iron_stats['slope_ep']
    efx2 = s.iron_stats['err_slope_ep']
    gx2  = s.iron_stats['slope_rew']
    egx2 = s.iron_stats['err_slope_rew']
    hx2  = s.iron_stats['afe1'] - s.iron_stats['afe2']
    ehx2 = np.sqrt(s.iron_stats['err_afe1']**2+
                   s.iron_stats['err_afe2']**2)/\
           np.sqrt(s.iron_stats['nfe1']+s.iron_stats['nfe2'])
    s.teff = s.teff + dteff
    dfx = 1e0*(fx2-fx1)/(2*dteff)
    efx = 1e0*(efx1+efx2)/2
    dgx = 1e0*(gx2-gx1)/(2*dteff)
    egx = 1e0*(egx1+egx2)/2
    dhx = 1e0*(hx2-hx1)/(2*dteff)
    ehx = 1e0*(ehx1+ehx2)/2
    dfdt, dgdt, dhdt = dfx, dgx, dhx
    eft, egt, eht = efx, egx, ehx

    s.vt = s.vt + dvt
    s.get_model_atmosphere(SolvePars.grid)
    specpars.iron_stats(s, Ref=Ref)
    fx1  = s.iron_stats['slope_ep']
    efx1 = s.iron_stats['err_slope_ep']
    gx1  = s.iron_stats['slope_rew']
    egx1 = s.iron_stats['err_slope_rew']
    hx1  = s.iron_stats['afe1'] - s.iron_stats['afe2']
    ehx1 = np.sqrt(s.iron_stats['err_afe1']**2+
                   s.iron_stats['err_afe2']**2)/\
           np.sqrt(s.iron_stats['nfe1']+s.iron_stats['nfe2'])
    s.vt = s.vt - 2*dvt
    s.get_model_atmosphere(SolvePars.grid)
    specpars.iron_stats(s, Ref=Ref)
    fx2  = s.iron_stats['slope_ep']
    efx2 = s.iron_stats['err_slope_ep']
    gx2  = s.iron_stats['slope_rew']
    egx2 = s.iron_stats['err_slope_rew']
    hx2  = s.iron_stats['afe1'] - s.iron_stats['afe2']
    ehx2 = np.sqrt(s.iron_stats['err_afe1']**2+
                   s.iron_stats['err_afe2']**2)/\
           np.sqrt(s.iron_stats['nfe1']+s.iron_stats['nfe2'])
    s.vt = s.vt + dvt
    dfx = 1e0*(fx2-fx1)/(2*dvt)
    efx = 1e0*(efx1+efx2)/2
    dgx = 1e0*(gx2-gx1)/(2*dvt)
    egx = 1e0*(egx1+egx2)/2
    dhx = 1e0*(hx2-hx1)/(2*dvt)
    ehx = 1e0*(ehx1+ehx2)/2
    dfdv, dgdv, dhdv = dfx, dgx, dhx
    efv, egv, ehv = efx, egx, ehx

    s.logg = s.logg + dlogg
    s.get_model_atmosphere(SolvePars.grid)
    specpars.iron_stats(s, Ref=Ref)
    fx1  = s.iron_stats['slope_ep']
    efx1 = s.iron_stats['err_slope_ep']
    gx1  = s.iron_stats['slope_rew']
    egx1 = s.iron_stats['err_slope_rew']
    hx1  = s.iron_stats['afe1'] - s.iron_stats['afe2']
    ehx1 = np.sqrt(s.iron_stats['err_afe1']**2+
                   s.iron_stats['err_afe2']**2)/\
           np.sqrt(s.iron_stats['nfe1']+s.iron_stats['nfe2'])
    s.logg = s.logg - 2*dlogg
    s.get_model_atmosphere(SolvePars.grid)
    specpars.iron_stats(s, Ref=Ref)
    fx2  = s.iron_stats['slope_ep']
    efx2 = s.iron_stats['err_slope_ep']
    gx2  = s.iron_stats['slope_rew']
    egx2 = s.iron_stats['err_slope_rew']
    hx2  = s.iron_stats['afe1'] - s.iron_stats['afe2']
    ehx2 = np.sqrt(s.iron_stats['err_afe1']**2+
                   s.iron_stats['err_afe2']**2)/\
           np.sqrt(s.iron_stats['nfe1']+s.iron_stats['nfe2'])
    s.logg = s.logg + dlogg
    dfx = 1e0*(fx2-fx1)/(2*dlogg)
    efx = 1e0*(efx1+efx2)/2
    dgx = 1e0*(gx2-gx1)/(2*dlogg)
    egx = 1e0*(egx1+egx2)/2
    dhx = 1e0*(hx2-hx1)/(2*dlogg)
    ehx = 1e0*(ehx1+ehx2)/2
    dfdg, dgdg, dhdg = dfx, dgx, dhx
    efg, egg, ehg = efx, egx, ehx
    
    d = matrix( [ [dfdt, dfdv, dfdg],
                  [dgdt, dgdv, dgdg],
                  [dhdt, dhdv, dhdg] ] )
    di = d.I

    s0 = np.mean([eft,efv,efg])
    s1 = np.mean([egt,egv,egg])
    s2 = np.mean([eht,ehv,ehg])

    eteff = np.sqrt ( (s0*di[0, 0])**2 +
                      (s1*di[0, 1])**2 +
                      (s2*di[0, 2])**2  )

    evt   = np.sqrt ( (s0*di[1, 0])**2 +
                      (s1*di[1, 1])**2 +
                      (s2*di[1, 2])**2  )

    elogg = np.sqrt ( (s0*di[2, 0])**2 +
                      (s1*di[2, 1])**2 +
                      (s2*di[2, 2])**2  )

    s.teff = s.teff + eteff
    s.get_model_atmosphere(SolvePars.grid)
    specpars.iron_stats(s, Ref=Ref)
    ap = s.iron_stats['afe']
    s.teff = s.teff - 2*eteff
    s.get_model_atmosphere(SolvePars.grid)
    specpars.iron_stats(s, Ref=Ref)
    am = s.iron_stats['afe']
    s.teff = s.teff + eteff
    eat = (ap-am)/2

    s.logg = s.logg + elogg
    s.get_model_atmosphere(SolvePars.grid)
    specpars.iron_stats(s, Ref=Ref)
    ap = s.iron_stats['afe']
    s.logg = s.logg - 2*elogg
    s.get_model_atmosphere(SolvePars.grid)
    specpars.iron_stats(s, Ref=Ref)
    am = s.iron_stats['afe']
    s.logg = s.logg + elogg
    eag = (ap-am)/2

    s.vt = s.vt + evt
    s.get_model_atmosphere(SolvePars.grid)
    specpars.iron_stats(s, Ref=Ref)
    ap = s.iron_stats['afe']
    s.vt = s.vt - 2*evt
    s.get_model_atmosphere(SolvePars.grid)
    specpars.iron_stats(s, Ref=Ref)
    am = s.iron_stats['afe']
    s.vt = s.vt + evt
    eav = (ap-am)/2

    ea = np.sqrt(eat**2+eag**2+eav**2+s2**2)

    Star_in.sp_err = {'teff': int(eteff), 'logg': elogg, 'afe': ea, 'vt': evt}
