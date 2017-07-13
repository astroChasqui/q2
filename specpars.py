import numpy as np
import os
import logging
import matplotlib.pyplot as plt
from . import moog, errors
from .tools import linfit
from .star import Star
import datetime
from scipy import ma
from collections import OrderedDict
from bokeh.plotting import *
from bokeh.models import HoverTool

logger = logging.getLogger(__name__)


class SolvePars:
    def __init__(self, grid='odfnew'):
        self.step_teff = 32
        self.step_logg = 0.32
        self.step_vt = 0.32
        self.niter = 50
        self.grid = grid
        self.solar_afe = 7.45
        self.errors = False
        self.check_converged = True
        self.ignore = []

class PlotPars:
    def __init__(self):
        self.afe = None
        self.wavelength_range = None
        self.make_figure = True
        self.figure_format = 'png'
        self.title = None
        self.title_inside = None

def iron_stats(Star, Ref=object, plot=None, PlotPars=object, silent=True):
    if hasattr(Ref, 'name'):
        if Star.name == Ref.name:
            x = {'afe': 0, 'err_afe': 0,
                 'afe1': 0, 'err_afe1': 0, 'nfe1': 0,
                 'afe2': 0, 'err_afe2': 0, 'nfe2': 0,
                 'slope_ep': 0,
                 'err_slope_ep': 0,
                 'slope_rew': 0,
                 'err_slope_rew': 0,
                 'rank': 0,
                 'reference': Ref.name}
            Star.iron_stats = x
            return None
    logger.info('Begin iron_stats for '+Star.name)
    logger.info('Calculating abundances for '+Star.name)
    fe1_done = moog.abfind(Star, 26.0, 'fe1')
    fe2_done = moog.abfind(Star, 26.1, 'fe2')
    if not fe1_done and not fe2_done:
        logger.warning('No fe1/fe2 attribute(s) added to '+Star.name)
        return None
    if hasattr(Ref, 'name'):
        logger.info('Differential analysis. Reference star is '+Ref.name)
        if not (hasattr(Ref, 'fe1')):
            logger.info('Reference star does not have abundances as '+\
                        'attributes')
            logger.info('Calculating abundances for reference star')
            moog.abfind(Ref, 26.0, 'fe1')
            moog.abfind(Ref, 26.1, 'fe2')
        ww1, ww2 = Star.fe1['ww'], Star.fe2['ww']
        ww1r, ww2r = Ref.fe1['ww'], Ref.fe2['ww']
        w1, w2 = np.intersect1d(ww1, ww1r), np.intersect1d(ww2, ww2r)
        k1 = [i for i, w in zip(range(len(ww1)), ww1) if w in w1]
        k1r = [i for i, w in zip(range(len(ww1r)), ww1r) if w in w1]
        k2 = [i for i, w in zip(range(len(ww2)), ww2) if w in w2]
        k2r = [i for i, w in zip(range(len(ww2r)), ww2r) if w in w2]
        afe1 = Star.fe1['ab'][k1] - Ref.fe1['ab'][k1r]
        afe2 = Star.fe2['ab'][k2] - Ref.fe2['ab'][k2r]
        rew1 = np.log10(1e-3*Star.fe1['ew'][k1]/Star.fe1['ww'][k1])
        rew2 = np.log10(1e-3*Star.fe2['ew'][k2]/Star.fe2['ww'][k2])
        ep1, ep2 = Star.fe1['ep'][k1], Star.fe2['ep'][k2]
        w1 = Star.fe1['ww'][k1]
        w2 = Star.fe2['ww'][k2]
        #
        Star.fe1['ww'], Star.fe2['ww'] = w1, w2
        Star.fe1['ep'], Star.fe2['ep'] = ep1, ep2
        Star.fe1['ew'], Star.fe2['ew'] = Star.fe1['ew'][k1], Star.fe2['ew'][k2]
        Star.fe1['rew'], Star.fe2['rew'] =rew1, rew2
        Star.fe1['ab'], Star.fe2['ab'] = Star.fe1['ab'][k1], Star.fe2['ab'][k2]
        Star.fe1['difab'], Star.fe2['difab'] = afe1, afe2
        #
        if plot:
            #ylabel = '$\Delta$[Fe/H]'
            ylabel = '[Fe/H]'
    else:
        logger.info('Working with absolute abundances')
        w1, w2 = Star.fe1['ww'], Star.fe2['ww']
        afe1 = Star.fe1['ab']
        afe2 = Star.fe2['ab']
        rew1 = np.log10(1e-3*Star.fe1['ew']/w1)
        rew2 = np.log10(1e-3*Star.fe2['ew']/w2)
        ep1, ep2 = Star.fe1['ep'], Star.fe2['ep']
        if plot:
            ylabel = 'A(Fe)'

    mfe1, efe1 = np.mean(afe1), np.std(afe1, ddof=1)
    mfe2, efe2 = np.mean(afe2), np.std(afe2, ddof=1)
    mafe = np.mean(list(afe1)+list(afe2))
    eafe = np.std(list(afe1)+list(afe2))
    nfe1, nfe2 = len(afe1), len(afe2)

    zero_ep, slope_ep, err_slope_ep = linfit(ep1, afe1)
    zero_rew, slope_rew, err_slope_rew = linfit(rew1, afe1)
    x_epfit = np.array([min(ep1), max(ep1)])
    y_epfit = zero_ep + slope_ep*x_epfit
    x_rewfit = np.array([min(rew1), max(rew1)])
    y_rewfit = zero_rew + slope_rew*x_rewfit

    if plot:
        logger.info('Making figure')
        plt.figure(figsize=(7, 9))
        title = Star.name+' : '+str(int(Star.teff))+', '+str(Star.logg)+', ' \
                +str(round(Star.feh,3))+', '+str(Star.vt)
        if hasattr(Ref, 'name'):
            title += ' ['+Ref.name+']'
        if hasattr(PlotPars, 'title'):
            if PlotPars.title != None:
                title = PlotPars.title
        plt.suptitle(title, fontsize=16)
        plt.subplots_adjust(hspace=0.35, top=0.93, left=0.2)

        try:
            if PlotPars.afe[0] != -1000:
                ylim = [PlotPars.afe[0], PlotPars.afe[1]]
            else:
                ylim = [mafe-abs(PlotPars.afe[1]),
                        mafe+abs(PlotPars.afe[1])]
        except:
            ylim = [mafe-4*eafe, mafe+4*eafe]

        panel_a = plt.subplot(311)
        plt.xlabel('EP = $\chi$ (eV)')
        plt.ylabel(ylabel)
        plt.xlim(-0.2, 5.2)
        plt.ylim(ylim)
        if hasattr(PlotPars, 'title_inside'):
            if PlotPars.title_inside != None:
                plt.text(plt.xlim()[0]+0.50*(plt.xlim()[1]-plt.xlim()[0]),
                         plt.ylim()[0]+0.85*(plt.ylim()[1]-plt.ylim()[0]),
                         PlotPars.title_inside,
                         horizontalalignment='center',
                         size=16)
        panel_b = plt.subplot(312)
        plt.xlabel('REW = log (EW/$\lambda$)')
        plt.ylabel(ylabel)
        plt.xlim(1.02*min(list(rew1)+list(rew2)),
                 0.98*max(list(rew1)+list(rew2)))
        plt.ylim(ylim)

        panel_c = plt.subplot(313)
        plt.xlabel('Wavelength ($\mathrm{\AA}$)')
        plt.ylabel(ylabel)
        try:
            plt.xlim(PlotPars.wavelength_range[0], PlotPars.wavelength_range[1])
        except:
            plt.xlim(4100, 7900)
        plt.ylim(ylim)

        panel_a.plot(ep1, afe1, 'b+')
        panel_a.plot(ep2, afe2, 'go')
        panel_a.plot(x_epfit, y_epfit, 'b')

        panel_b.plot(rew1, afe1, 'b+')
        panel_b.plot(rew2, afe2, 'go')
        panel_b.plot(x_rewfit, y_rewfit, 'b')

        panel_c.plot(w1, afe1, 'b+')
        panel_c.plot(w2, afe2, 'go')
        panel_c.plot([4000, 8000], [mafe, mafe], 'black')

        if hasattr(PlotPars, 'directory'):
            if not os.path.exists(PlotPars.directory):
                os.mkdir(PlotPars.directory)
            plot = PlotPars.directory+'/'+plot
        if hasattr(PlotPars, 'figure_format'):
            plot = plot+'.'+PlotPars.figure_format
        plt.savefig(plot, bbox_inches='tight')
        #plt.close()

    if hasattr(Ref, 'name'):
        ref_star = Ref.name
    else:
        ref_star = None

    dfe = mfe1 - mfe2
    edfe = np.sqrt(efe1**2/nfe1+efe2**2/nfe2)

    x = {'afe': round(mafe, 3), 'err_afe': round(eafe, 3),
         'afe1': round(mfe1, 3), 'err_afe1': round(efe1, 3), 'nfe1': nfe1,
         'afe2': round(mfe2, 3), 'err_afe2': round(efe2, 3), 'nfe2': nfe2,
         'slope_ep': slope_ep,
         'err_slope_ep': err_slope_ep,
         'slope_rew': slope_rew,
         'err_slope_rew': err_slope_rew,
         'reference': ref_star}
    Star.iron_stats = x

    if not silent:
        print("FeI  : {0:6.3f} +/- {1:5.3f} (n={2:3.0f})".\
              format(mfe1, efe1, nfe1))
        print("FeII : {0:6.3f} +/- {1:5.3f} (n={2:3.0f})".\
              format(mfe2, efe2, nfe2))


def solve_one(Star, SolveParsInit, Ref=object, PlotPars=object):
    sp = SolvePars()
    sp.__dict__ = SolveParsInit.__dict__.copy()
    if not hasattr(Star, 'model_atmosphere_grid'):
        logger.info('Star has no model yet. Calculating.')
        Star.get_model_atmosphere(sp.grid)
    if not hasattr(Star, 'model_atmosphere'):
        print('Unable to find a starting model atmosphere for this star')
        return None
    if Star.model_atmosphere_grid != sp.grid:
        logger.info('Inconsistent model atmosphere grids '+
                     '(Star and SolvePars). '+
                     'Fixing problem now.')
        Star.get_model_atmosphere(sp.grid)

    if hasattr(Ref, 'name'):
        if not hasattr(Ref, 'model_atmosphere_grid'):
            logger.info('Ref star has no model yet. Calculating.')
            Ref.get_model_atmosphere(sp.grid)
        if Ref.model_atmosphere_grid != sp.grid:
            logger.info('Inconsistent model atmosphere grids '+
                         '(Ref star and SolvePars). '+
                         'Fixing problem now.')
            Ref.get_model_atmosphere(sp.grid)

    dtv, dgv, dvv, stop_iter = [], [], [], False
    if hasattr(Star, 'converged'):
        if not Star.converged:
            Star.converged = False
    else:
        Star.converged = False
    Star.stop_iter = sp.niter
    if sp.niter == 0:
        Star.converged = True

    print('it Teff logg [Fe/H]  vt           [Fe/H]')
    print('-- ---- ---- ------ ----      --------------')

    for i in range(sp.niter+1):
        if sp.step_teff <= 1 and sp.step_logg <= 0.01 \
           and sp.step_vt <= 0.01:
            if not stop_iter:
                Star.converged = False
                if SolveParsInit.niter > 0:
                    print('-- Begin final loop')
            stop_iter = True

        if i > 0:
            if Star.iron_stats['slope_ep'] > 0:
                Star.teff += sp.step_teff
            else:
                Star.teff -= sp.step_teff
            if Star.teff > 7000:
                Star.teff = 7000
            if Star.iron_stats['slope_rew'] > 0:
                Star.vt += sp.step_vt
            else:
                Star.vt -= sp.step_vt
            if Star.vt < 0:
                Star.vt = 0
            dfe = Star.iron_stats['afe1'] - Star.iron_stats['afe2']
            if dfe > 0:
                Star.logg += sp.step_logg
            else:
                Star.logg -= sp.step_logg
            if Star.logg > 5.0:
                Star.logg = 5.0

            if hasattr(Ref, 'name'):
                Star.feh = Ref.feh + Star.iron_stats['afe']
            else:
                Star.feh = Star.iron_stats['afe'] - sp.solar_afe
            if Star.feh > 1.0:
                Star.feh = 1.0
            if Star.feh > 0.5 and sp.grid != 'over':
                Star.feh = 0.5

            Star.get_model_atmosphere(sp.grid)

        if i+1 == sp.niter or sp.niter == 0:
            plot = Star.name
            if hasattr(Ref, 'name'):
                plot = Star.name+'-'+Ref.name
                if Star.name == Ref.name:
                    plot = None
                    Star.converged = ''
        else:
            plot = None

        is_done = iron_stats(Star, Ref=Ref, plot=plot, PlotPars=PlotPars)

        print("{0:2.0f} {1:4.0f} {2:4.2f} {3:6.3f} {4:4.2f}"\
              " ---> {5:6.3f}+/-{6:5.3f}".\
                format(i, Star.teff, Star.logg, Star.feh, Star.vt,
                          Star.iron_stats['afe'], Star.iron_stats['err_afe']))

        dtv.append(Star.teff)
        dgv.append(Star.logg)
        dvv.append(Star.vt)

        if i >= 4:
            if np.std(dtv[-5:]) <= 0.8*sp.step_teff and \
               np.std(dgv[-5:]) <= 0.8*sp.step_logg and \
               np.std(dvv[-5:]) <= 0.8*sp.step_vt:
                print('-- Converged at iteration '+str(i)+ \
                      ' of '+str(sp.niter))
                if stop_iter:
                    plot = Star.name
                    if hasattr(Ref, 'name'):
                        plot = Star.name+'-'+Ref.name
                    iron_stats(Star, Ref=Ref, plot=plot, PlotPars=PlotPars)
                    Star.converged = True
                    Star.stop_iter = i
                    break
                sp.step_teff = sp.step_teff/2
                sp.step_logg = sp.step_logg/2
                sp.step_vt = sp.step_vt/2
                if sp.step_teff < 1 and sp.step_teff > 0:
                    sp.step_teff = 1
                if sp.step_logg < 0.01 and sp.step_logg > 0:
                    sp.step_logg = 0.01
                if sp.step_vt < 0.01 and sp.step_vt > 0:
                    sp.step_vt = 0.01

    if not Star.converged:
        if hasattr(Ref, 'name'):
            if Star.name == Ref.name or SolveParsInit.niter == 0:
                print('--')
            else:
                print('-- Did not achieve final convergence.')
        else:
            print('-- Did not achieve final convergence.')

    print('------------------------------------------------------')

    if hasattr(Ref, 'name'):
        print('   D[Fe/H]    ||    D[Fe/H] Fe I   |   D[Fe/H] Fe II')
    else:
        print('    A(Fe)     ||      A(Fe I)      |     A(Fe II)   ')

    print("{0:6.3f} {1:6.3f} || {2:6.3f} {3:6.3f} {4:3d} "\
          "| {5:6.3f} {6:6.3f} {7:3d}".\
            format(Star.iron_stats['afe'], Star.iron_stats['err_afe'],
                   Star.iron_stats['afe1'], Star.iron_stats['err_afe1'],
                   Star.iron_stats['nfe1'],
                   Star.iron_stats['afe2'], Star.iron_stats['err_afe2'],
                   Star.iron_stats['nfe2']))
    print('------------------------------------------------------')

    Star.sp_err = {'teff': 0, 'logg': 0, 'afe': 0, 'vt': 0}
    if ((Star.converged and sp.errors == True) or \
        (sp.niter == 0 and sp.errors == True and Star.converged != '')):
        errors.error_one(Star, sp, Ref)
        Star.err_teff = int(Star.sp_err['teff'])
        Star.err_logg = Star.sp_err['logg']
        Star.err_feh = Star.sp_err['afe']
        Star.err_vt = Star.sp_err['vt']
        print("Solution with formal errors:")
        print("Teff    = {0:6d} +/- {1:5d}".\
              format(int(Star.teff), int(Star.sp_err['teff'])))
        print("log g   = {0:6.3f} +/- {1:5.3f}".\
              format(Star.logg, Star.sp_err['logg']))
        if hasattr(Ref, 'name'):
            print("D[Fe/H] = {0:6.3f} +/- {1:5.3f}".\
                  format(Star.iron_stats['afe'], Star.sp_err['afe']))
        else:
            print("A(Fe)   = {0:6.3f} +/- {1:5.3f}".\
                  format(Star.iron_stats['afe'], Star.sp_err['afe']))
        print("vt      = {0:6.2f} +/- {1:5.2f}".\
              format(Star.vt, Star.sp_err['vt']))
        print('------------------------------------------------------')


def solve_all(Data, SolveParsInit, output_file, reference_star=None,
              PlotPars=object):
    print('------------------------------------------------------')
    print('Initializing ...')
    start_time = datetime.datetime.now()
    print('- Date and time: '+start_time.strftime('%d-%b-%Y, %H:%M:%S'))
    print('- Model atmospheres: '+SolveParsInit.grid)
    print('- Star data: '+Data.star_data_fname)
    print('- Line list: '+Data.lines_fname)
    print('------------------------------------------------------')
    if reference_star:
        Ref = Star(reference_star)
        Ref.get_data_from(Data)
    else:
        Ref = None
    fout = open(output_file, 'w')
    if SolveParsInit.errors:
        fout.write('id,teff,logg,feh_model,vt,feh,err_feh_,'+
                   'feh1,err_feh1,nfe1,feh2,err_feh2,nfe2,'
                   'slope_ep,err_slope_ep,slope_rew,err_slope_rew,'
                   'stop_iter,converged,'
                   'err_teff,err_logg,err_feh,err_vt\n')
    else:
        fout.write('id,teff,logg,feh_model,vt,feh,err_feh,'+
                   'feh1,err_feh1,nfe1,feh2,err_feh2,nfe2,'
                   'slope_ep,err_slope_ep,slope_rew,err_slope_rew,'
                   'stop_iter,converged,'
                   'err_teff,err_logg,err_feh_,err_vt\n')
    for star_id in Data.star_data['id']:
        print('')
        print('*'*len(star_id))
        print(star_id)
        print('*'*len(star_id))
        s = Star(star_id)
        try:
            s.get_data_from(Data)
        except:
            logger.warning('No data found for '+s.name+\
                        '. Excluded from output file.')
            print('Data not found.')
            #fout.write("{0},,,,,,,,,,"\
            #           ",,,,,,,,,,,,\n".\
            #           format(s.name))
            continue
        if ma.count(Data.lines[star_id]) == 0:
            print('Line data not found.')
            continue
        sp = SolvePars()
        sp.__dict__ = SolveParsInit.__dict__.copy()
        if reference_star:
            if s.name == Ref.name:
                sp.niter = 0
                print('Reference star. No calculations needed.')
                #continue
        if hasattr(s, 'converged') and sp.check_converged:
            if s.converged == 'True':
                print('Already converged.')
                continue
                #sp.niter = 0
                #s.converged = True
        if s.name in sp.ignore:
            print('Asked to ignore.')
            continue

        solve_one(s, sp, Ref, PlotPars=PlotPars)

        if sp.niter == 0:
            s.converged = ''

        fout.write("{0},{1:4.0f},{2:5.3f},{3},{4:4.2f},{5},{6:5.3f},"\
                   "{7},{8:5.3f},{9},"\
                   "{10},{11:5.3f},{12},{13:.6f},{14:.6f},"\
                   "{15:.6f},{16:.6f},{17},{18},"\
                   "{19:3.0f},{20:5.3f},{21:5.3f},{22:4.2f}\n".\
                   format(s.name, s.teff, s.logg, str(round(s.feh,3)), s.vt,
                          str(round(s.iron_stats['afe'],3)),
                          s.iron_stats['err_afe'],
                          str(round(s.iron_stats['afe1'],3)),
                          s.iron_stats['err_afe1'],
                          s.iron_stats['nfe1'],
                          str(round(s.iron_stats['afe2'],3)),
                          s.iron_stats['err_afe2'],
                          s.iron_stats['nfe2'],
                          s.iron_stats['slope_ep'],
                          s.iron_stats['err_slope_ep'],
                          s.iron_stats['slope_rew'],
                          s.iron_stats['err_slope_rew'],
                          s.stop_iter,
                          s.converged,
                          s.sp_err['teff'], s.sp_err['logg'],
                          s.sp_err['afe'], s.sp_err['vt']
                          ))
    fout.close()

    print('')
    print('------------------------------------------------------')
    end_time = datetime.datetime.now()
    print('- Date and time: '+end_time.strftime('%d-%b-%Y, %H:%M:%S'))
    delta_t = (end_time - start_time).seconds
    hours, remainder = divmod(delta_t, 3600)
    minutes, seconds = divmod(remainder, 60)
    print('- Time elapsed: %sH %sM %sS' % (hours, minutes, seconds))
    print('Done!')
    print('------------------------------------------------------')
    print('')


def make_single_solution_table(solution_files, single_solution_file):
    """Takes q2.specpars.solve_all outputs and creates a single final one

    Files must be in the order in which they were computed!
    """
    #solution_files = ['starsDec_solution1.csv', 'starsDec_solution2.csv']
    #single_solution_file = 'starsDec_solution.csv'
    fout = open(single_solution_file, 'w')
    with open(solution_files[0], 'r') as f:
        lines = f.readlines()
    for line in lines:
        sid = line[0:line.index(',')]
        if 'True' in line or 'id,teff' in line:
            #nline = line[0:line.rfind('\n')]
            #linew = nline[0:nline.rfind('\n')]
            #fout.write(linew+'\n')
            #print line
            fout.write(line)
        else:
            for i in range(1, len(solution_files)):
                with open(solution_files[i], 'r') as f2:
                    lines2 = f2.readlines()
                for line2 in lines2:
                    sid2 = line2[0:line2.index(',')]
                    #nline2 = line2[0:line2.rfind(',')]
                    #line2w = nline2[0:nline2.rfind(',')]
                    if 'True' in line2 and sid == sid2:
                        #fout.write(line2w+'\n')
                        fout.write(line2)
    fout.close()


def fancy_ironstats_plot(Star):
    """Makes bokeh hover-ing plots

    Function written to look for outliers and investigate line-to-line scatter
    """
    if not hasattr(Star, 'iron_stats'):
        logger.error('Star object ('+Star.name+') has no ironstats attribute.')
        return None
    ww = np.concatenate((Star.fe1['ww'], Star.fe2['ww']))
    ep = np.concatenate((Star.fe1['ep'], Star.fe2['ep']))
    ew = np.concatenate((Star.fe1['ew'], Star.fe2['ew']))
    rew = np.concatenate((Star.fe1['rew'], Star.fe2['rew']))
    if Star.iron_stats['reference']:
        ab = np.concatenate((Star.fe1['difab'], Star.fe2['difab']))
        y_axis_label = '[Fe/H]'
    else:
        ab = np.concatenate((Star.fe1['ab'], Star.fe2['ab']))
        y_axis_label = 'A(Fe)'
    ws = [str(round(w, 1)) for w in ww]
    colors = np.concatenate((["blue"]*len(Star.fe1['ww']),
                             ["green"]*len(Star.fe2['ww'])))

    TOOLS="pan,wheel_zoom,box_zoom,reset,hover"
    output_notebook()

    title = Star.name
    if getattr(Star, 'iron_stats')['reference']:
        title += ' - '+getattr(Star, 'iron_stats')['reference']

    p1 = figure(title=title, plot_width=650, plot_height=300,
                x_axis_label='EP (eV)',
                y_axis_label=y_axis_label,
                tools=TOOLS, active_scroll = 'wheel_zoom')
    p1.xaxis.axis_label_text_font_style = "normal"
    p1.yaxis.axis_label_text_font_style = "normal"

    abst = [str(round(xab, 3)) for xab in ab]
    source = ColumnDataSource(
        data=dict(
            ws = ws,
            ep = ep,
            rew = rew,
            ab = ab,
            abst = abst,
            ew = ew,
            colors = colors,
        )
    )

    p1.scatter('ep', 'ab', size=10, color='colors',
            source=source, marker='circle')

    hover = p1.select(dict(type=HoverTool))
    hover.tooltips = OrderedDict([
        ("Wavelength, EP", "@ws A, @ep eV"),
        ("EW, REW", "@ew mA, @rew"),
        ("Abundance", "@abst"),
    ])

    show(p1)

    p2 = figure(title='', plot_width=650, plot_height=300,
                x_axis_label='REW',
                y_axis_label=y_axis_label,
                tools=TOOLS, active_scroll = 'wheel_zoom')
    p2.xaxis.axis_label_text_font_style = "normal"
    p2.yaxis.axis_label_text_font_style = "normal"

    p2.scatter('rew', 'ab', size=10, color='colors',
            source=source, marker='circle')

    hover = p2.select(dict(type=HoverTool))
    hover.tooltips = OrderedDict([
        ("Wavelength, EP", "@ws A, @ep eV"),
        ("EW, REW", "@ew mA, @rew"),
        ("Abundance", "@abst"),
    ])

    show(p2)
