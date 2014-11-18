import specpars
import numpy as np
from star import *
import sys

stdout = sys.stdout

def one(StarIn, RefIn):
    ngp = 20
    dteff = 10
    dlogg = 0.01
    tg = np.arange(StarIn.teff-ngp*dteff, StarIn.teff+ngp*dteff, dteff)
    gg = np.arange(StarIn.logg-ngp*dlogg, StarIn.logg+ngp*dlogg, dlogg)
    sp = specpars.SolvePars()
    sp.step_teff = 0
    sp.step_logg = 0
    sp.step_vt = 0.0
    sp.niter = 1
    s = Star(StarIn.name)
    s.feh = StarIn.feh
    s.vt = StarIn.vt
    s.vt = 1.30
    s.linelist = StarIn.linelist
    print "slep,eslep,slrew,eslrew,dfe,edfe,t,g,f,v"
    for t in tg:
        for g in gg:
            s.teff = t
            s.logg = g
            sys.stdout = open('grids.tmp', 'w')
            specpars.solve_one(s, sp, RefIn)
            sys.stdout = stdout
            afe1 = s.iron_stats['afe1']
            eafe1 = s.iron_stats['err_afe1']
            nfe1 = s.iron_stats['nfe1']
            afe2 = s.iron_stats['afe2']
            eafe2 = s.iron_stats['err_afe2']
            nfe2 = s.iron_stats['nfe2']
            dfe = afe1 - afe2
            edfe = np.sqrt(eafe1**2/nfe1+eafe2**2/nfe2)
            print "{0:.5f},{1:.5f},{2:.6f},{3:.6f},{4:.3f},{5:.3f},"\
                  "{6:4d},{7:.2f},{8:.3f},{9:.2f}".\
                   format(
                          s.iron_stats['slope_ep'],
                          s.iron_stats['err_slope_ep'],
                          s.iron_stats['slope_rew'],
                          s.iron_stats['err_slope_rew'],
                          dfe, edfe,
                          t, g, s.feh, s.vt,
                          )
