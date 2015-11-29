"""
Provides the appropriate hierarchy of the as-sorted tree below
to facilitate the analysis.

`-- bySample
    |-- Asteel
    |   |-- A_BB
    |   |   |-- 2012Jul25-26
    |   |   |   `-- lin
    |   |   `-- 2012July24
    |   |-- A_PS_RD
    |   |   `-- 2012July18
    |   `-- A_PS_TD
    |       |-- 2012July18
    |       `-- 2012July19
    |-- Bsteel
    |   |-- B_BB
    |   |   |-- 2012July11
    |   |   |-- 2012July20
    |   |   `-- 2012July23
    |   |-- B_PS_RD
    |   |   `-- 2012July13
    |   `-- B_PS_TD
    |       `-- 2012July17
    `-- etc
        |-- 2011Dec02
        |-- 2011Dec10
        `-- 2012July10
"""
def main(sample = 'Asteel'):
    """
    Takes a certain sample's data through the hierchy of the directory.
    Returns the data as a dictionary to facilitate any further analysis.
    """
    import os
    from __reader__ import __logreader__ as logreader
    bpath = os.path.join(os.getcwd(), 'bySample')

    fdict = dict(filenames=[], linfn = [] ,BB = [], PS_RD = [], PS_TD = [])

    cwd1 = os.path.join(bpath, sample)
    nspots1 = [o for o in os.listdir(cwd1) \
                   if os.path.isdir(os.path.join(cwd1,o))]

    for j in xrange(len(nspots1)): # BB, PS
        print 'nspots1:', nspots1
        cwd2 = os.path.join(cwd1, nspots1[j])
        nspots2 = [o for o in os.listdir(cwd2) \
                       if os.path.isdir(os.path.join(cwd2,o))]

        for k in xrange(len(nspots2)): # 2012July25-26 ...
            print 'nspots2:', nspots2[k]
            cwd3 = os.path.join(cwd2, nspots2[k])
            fn = os.path.join(cwd3, 'log') # log file.
            trfns, lnfns, images = logreader(fn)
            print 'nfns', lnfns
            print 'trfns', trfns

            for i in xrange(len(trfns)):
                fdict['filenames'].append(os.path.join(cwd3, trfns[i]))
                for l in xrange(len(lnfns[i])):
                    fdict['linfn'].append(os.path.join(cwd3, lnfns[i][l]))

            if nspots1[j].split('_')[-1]=='BB':
                fdict['BB'].append(fn)
            elif nspots1[j].split('_')[-2]+'_'\
                    +nspots1[j].split('_')[-1]=='PS_RD':
                fdict['PS_RD'].append(fn)
            elif nspots1[j].split('_')[-2]+'_'\
                    +nspots1[j].split('_')[-1]=='PS_TD':
                fdict['PS_TD'].append(fn)
            else:
                raise IOError, 'Unexpected case.'
    return fdict

def copy(sample='Asteel'):
    """
    """
    import numpy as np
    import shutil
    import os
    import __reader__
    pwd = os.getcwd()
    fdict = main(sample)

    # for i in xrange(len(fdict['filenames'])):
    #     fn = fdict['filenames'][i]
    #     filename = os.path.split(fn)[-1]

    #     try: shutil.copy(fn, os.path.join(pwd, sample, filename))
    #     except: print 'Failed copying %s'%fn
    #     else: pass #print os.path.join(pwd, sample, filename)
    print np.array(fdict['linfn']).shape

    missing = []
    for j in xrange(len(fdict['linfn'])):

        fn = fdict['linfn'][j]
        filename = os.path.split(fn)[-1]
        print filename

        try:shutil.copy(fn, os.path.join(pwd, sample, filename))
        except:
            print '%s is missing'%fn
            missing.append(fn)
    print 'missing file names:'
    print missing
