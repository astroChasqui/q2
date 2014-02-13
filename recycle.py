"""Code that is no longer used by q2

Saved just for future reference.
"""

import inspect
from config import *


def q2print(message, err=0):
    """Prefixes a message with the calling module and function name

    err = 0 is an INFO message
    err = 1 is an ERROR message
    """

    if err == 0 and Messages.info == 0:
        return None
    if err == 1 and Messages.error == 0:
        return None

    maxc = 70
    if err == 0:
        err_strg = 'INFO'
    else:
        err_strg = 'ERROR'
    x = inspect.stack()[1]
    mname = inspect.getmodule(x[0]).__name__
    fname = inspect.stack()[1][3]
    msg = '['+mname+'.'+fname+': '+err_strg+'] '+message
    lmsg = len(msg)
    if lmsg <= maxc:
        print(msg)
    else:
        msga =  msg.split(' ')
        words = ''
        for word in msga:
            words += word+' '
            if len(words) >= maxc:
                print(words)
                words = '    '
        if words.replace(' ', '') != '':
            print(words)

def create_moog_model_in(Star, file_name='model.in'):
    try:
        Star.vt
    except:
        logger.error('Moog model_in file requires a microturbulence (vt)')
        return None
    if hasattr(Star, 'model_atmosphere'):
        Star.moog_model_in_name = file_name
    else:
        logger.error('No model data to write to moog model_in file.')
        return None

    f = open('head.tmp', 'w')
    f.write('KURUCZ\n')
    f.write('TEFF='+str(Star.teff)+',LOGG='+str(Star.logg)+
            ',[FE/H]='+str(Star.feh)+','+Star.model_atmosphere_grid+'\n')
    nd = len(Star.model_atmosphere)
    f.write('ND=       '+str(nd)+'\n')
    f.close()

    ascii.write(Star.model_atmosphere, 'body.tmp', Writer=ascii.NoHeader,
                formats={'RHOX':'%.8E','T':'%.1F','P':'%.3E',
                'XNE':'%.3E','ABROSS':'%.3E'})

    f = open('tail.tmp', 'w')
    f.write('%5.2F\n' %Star.vt)
    if Star.model_atmosphere_grid != 'marcs':
        path = MODATM_PATH+'kurucz/'
        fabund = open(path+'p00.'+Star.model_atmosphere_grid, 'r')
    else:
        path = MODATM_PATH+'marcs/'
        fabund = open(path+'z+0.00', 'r')

    line = fabund.readline()
    f.write(line[0:12]+' '+str(Star.feh)+'\n')
    line = fabund.readline()
    while line:
        species = line[0:2]
        if Star.model_atmosphere_grid == 'marcs':
            abund = float(line[3:9])+Star.feh
            #alpha-element enhancement
            if species==' 8' or species=='10' or species=='12' or \
               species=='14' or species=='16' or species=='18' or \
               species=='20' or species=='22':
                afe = -0.4*Star.feh
                if Star.feh >=  0: afe=0.0
                if Star.feh <= -1: afe=0.4
                abund = abund+afe
        else:
            abund = 12.+np.log10(np.power(10, float(line[3:9]))/0.92040)+ \
                    Star.feh
        abund = str('%5.2F' %abund)
        f.write(species+' '+abund+'\n')
        line = fabund.readline()
    fabund.close()
    f.write('NMOL      22\n')
    f.write('  101.0   106.0   107.0   108.0   112.0  126.0\n')
    f.write('  606.0   607.0   608.0\n')
    f.write('  707.0   708.0\n')
    f.write('  808.0   812.0   822.0\n')
    f.write('  10108.0 60808.0\n')
    f.write('  6.1     7.1     8.1   12.1  22.1  26.1\n')
    f.close()

    file_list = ['head.tmp', 'body.tmp', 'tail.tmp']
    with open(file_name, 'w') as outfile:
        for one_file in file_list:
            with open(one_file) as infile:
                outfile.write(infile.read())
    for one_file in file_list:
        os.unlink(one_file)
    logger.info('Moog infile model atmosphere created: '+file_name)

def create_moog_lines_in(Star, species=0, file_name='lines.in'):
    #species = 0 means all species
    if species > 0:
        idx = np.where(Star.linelist['species'] == species)
    else:
        idx = np.where(Star.linelist['species'] > 0)
    nlines = len(idx[0])
    ll_data = [Star.linelist['wavelength'][idx],
               Star.linelist['species'][idx],
               Star.linelist['ep'][idx],
               Star.linelist['gf'][idx],
               3.+np.zeros(nlines),
               np.zeros(nlines),
               Star.linelist['ew'][idx]]
    ascii.write(ll_data, file_name)





'''Tests ran to see if speed was improved (q2.yyage)
Loading isochrones to memory does make it run
a little bit faster, but not enough to go
through the trouble of filling up the memory
'''

from StringIO import StringIO

def db_in_memory(db='yy01.sql3'):
    # Read database to tempfile
    con = sqlite3.connect(db)
    tempfile = StringIO()
    for line in con.iterdump():
        tempfile.write('%s\n' % line)
    con.close()
    tempfile.seek(0)
    # Create a database in memory and import from tempfile
    con = sqlite3.connect(":memory:")
    con.cursor().executescript(tempfile.read())
    con.commit()
    return con.cursor()

def gipmem(Star, c, nsigma=3):
    logtm = np.log10(Star.teff-nsigma*Star.err_teff)
    logtp = np.log10(Star.teff+nsigma*Star.err_teff)
    x = c.execute('SELECT feh, age, mass, logt, logl, logg, mv ' +\
                  'FROM  fa, mtlgv ON fa.fa = mtlgv.fa WHERE '   +\
                  'logt >= ? AND logt <= ? AND '   +\
                  'feh  >= ? AND feh  <= ? AND '   +\
                  'logg >= ? AND logg <= ? ',
                  (logtm, logtp,
                   Star.feh-nsigma*Star.err_feh,
                   Star.feh+nsigma*Star.err_feh,
                   Star.logg-nsigma*Star.err_logg,
                   Star.logg+nsigma*Star.err_logg
                  )
                 )

    feh, age = [], []
    mass, logt, logl, logg, mv = [], [], [], [], []
    for xx in x.fetchall():
        feh.append(xx[0])
        age.append(xx[1])
        mass.append(xx[2])
        logt.append(xx[3])
        logl.append(xx[4])
        logg.append(xx[5])
        mv.append(xx[6])
    return {'feh' : np.array(feh),
            'age' : np.array(age),
            'mass': np.array(mass),
            'logt': np.array(logt),
            'logl': np.array(logl),
            'logg': np.array(logg),
            'mv'  : np.array(mv)
            }
