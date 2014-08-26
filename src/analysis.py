""" Analysis script for the X-ray results """
import __sort__
import __reader__

def psidsp(ids=None, phi=-90):
    """ psi sine square vs d-spacing """
    import matplotlib.pyplot as plt
    fnlin = __sort__.main(mode='lin', i=ids)
    fn = None
    for i in range(len(fnlin)):
        p = float(fnlin[i].split('Phi')[1]\
                      .split('.')[0])
        if p==phi:
            fn = fnlin[i]
            break
    if fn==None: raise IOError, \
            '%Phi is not found'%phi

    d1, d2 = __reader__.main(fn = fn, mode='lin')
    sin2psi1 = d1.get('sin2psi')
    print 'sin2psi1', sin2psi1
    dspc1 = d1.get('dspacing')
    sin2psi2 = d2.get('sin2psi')
    dspc2 = d2.get('dspacing')

    ax = plt.gca()
    ax.plot(sin2psi1, dspc1, marker='x',
            label='detector1  %i at phi=%f'%(ids,phi))
    ax.plot(sin2psi2, dspc2, marker='+',
            label='detector2  %i at phi=%f'%(ids,phi))
    ax.legend(loc='best', fancybox=True).get_frame().set_alpha(0.5)

def psidsp_phi(ids=4, phi=[-90, 0, 45, 135]):
    """ sin2psi vs dspacing plots for many phis for a certain ids """
    for i in range(len(phi)):
        psidsp(ids=ids, phi=phi[i])

def psidsp_ids(phi=-90, ids=[4,5,6]):
    """ sin2psi vs dspacing plots for many ids for a certain phi angle"""
    for i in range(len(ids)):
        psidsp(ids=ids[i], phi=phi)

def readstress(i1=0,i2=0,path='.'):
    """ stress and its error bar for the given
    (i1, i2) component """
    fntr, dummy_ids = __sort__.main(mode='tr',path=path)
    sig = []
    err = []
    for i in range(len(fntr)):
        cfn = fntr[i]
        s, er = __reader__.main(cfn, mode='tr')
        sig.append(s[i1,i2])
        err.append(er[i1,i2])
    return sig, err

def readstressfn(i1=0,i2=0, fn=None):
    s, er = __reader__.main(fn, mode='tr')
    return s[i1, i2], er[i1,i2]

def ys(ix=0, iy=1):
    """ plot s11-s22 yield points """
    import matplotlib.pyplot as plt
    s11, s11e = readstress(ix,ix)
    s22, s22e = readstress(iy,iy)
    plt.errorbar(s11, s22, s11e, s22e, 'x-',
                 ecolor='red',)

def sigeps(ix=0, iy=1, path='.'):
    """ stress-strain curve ... """
    import numpy as np
    import matplotlib.pyplot as plt
    s11, s11e = readstress(ix, ix, path=path)
    s22, s22e = readstress(iy, iy, path=path)

    print '# of data points: ', len(s11)

    # An virtual strain
    ind = np.arange(len(s11))

    ax = plt.gca()
    ax.errorbar(ind, s11, yerr=s11e,
                marker='o', ls='--',
                label='%i%i'%(ix+1,ix+1))
    ax.errorbar(ind, s22, yerr=s22e,
                marker='o', ls='--',
                label='%i%i'%(iy+1,iy+1))
    ax.legend(loc='best')

    ax.set_ylabel(r'$\bar{\Sigma}$', fontsize=20)
    ax.set_xlabel('Index', fontsize=20)
