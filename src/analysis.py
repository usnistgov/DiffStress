""" Analysis script for the X-ray results """
import __sort__
import __reader__

def psidsp(ids=None, phi=-90):
    """ psi sine square vs d-spacing """
    import matplotlib.pyplot as plt
    fnlin = __sort__.main(mode='lin', i=ids)
    fn = None
    for i in xrange(len(fnlin)):
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
    for i in xrange(len(phi)):
        psidsp(ids=ids, phi=phi[i])

def psidsp_ids(phi=-90, ids=[4,5,6]):
    """ sin2psi vs dspacing plots for many ids for a certain phi angle"""
    for i in xrange(len(ids)):
        psidsp(ids=ids[i], phi=phi)

def readstress(i1=0,i2=0,path='.'):
    """ stress and its error bar for the given
    (i1, i2) component """
    fntr, dummy_ids = __sort__.main(mode='tr',path=path)
    sig = []
    err = []
    for i in xrange(len(fntr)):
        cfn = fntr[i]
        s, er = __reader__.main(cfn, mode='tr')
        sig.append(s[i1,i2])
        err.append(er[i1,i2])
    return sig, err, fntr

def readstressfn(i1=0,i2=0, fn=None):
    s, er = __reader__.main(fn, mode='tr')
    return s[i1, i2], er[i1,i2]

def ys(ix=0, iy=1):
    """ plot s11-s22 yield points """
    import matplotlib.pyplot as plt
    s11, s11e, dum = readstress(ix,ix)
    s22, s22e, dum = readstress(iy,iy)
    plt.errorbar(s11, s22, s11e, s22e, 'x-',
                 ecolor='red',)

def sigeps(ix=0, iy=1, path='.',dicfn=None,fig=None,iax=0,label=None,
           ieps_flip = False, isig_flip=False):
    """ stress-strain curve ... """
    import numpy as np
    from MP.lib import mpl_lib, axes_label
    eqv = axes_label.__eqv__

    from __reader__ import __logreader__ as read_log
    from os import sep
    from MP.mat import mech
    fc=mech.FlowCurve()
    
    s11, s11e, fntr = readstress(ix, ix, path=path)
    s22, s22e, dum  = readstress(iy, iy, path=path)

    print '# of data points: ', len(s11)

    fn_trs=[]
    for i in xrange(len(fntr)):
        fn_trs.append(fntr[i].split(sep)[-1])

    # An virtual strain
    if dicfn==None: 
        xs = np.arange(len(s11))
    else:
        ind, ex, ey, exy, ex_e, ey_e, exy_e =np.loadtxt(dicfn).T
        fn, ln_fn, images = read_log(fn='%s%slog'%(path,sep))
        S1=[];S1e=[]; S2=[];S2e=[]

        if isig_flip:
            dum = s11[::]
            s11 = s22[::]
            s22 = dum[::]

            dum = s11e[::]
            s11e = s22e[::]
            s22e = dum[::]

        for i in xrange(len(images)):
            fn_triaxial = fn[i]
            idx = np.array(fn_trs).searchsorted(fn_triaxial)
            S1.append(s11[idx])
            S2.append(s22[idx])
            S1e.append(s11e[idx])
            S2e.append(s22e[idx])

        ## Stress states
        fc.get_stress(S1,i=0,j=0)
        fc.get_stress(S2,i=1,j=1)
        fc.set_zero_sigma_ij(2,2)
        fc.set_zero_shear_stress()

        e11=[];e22=[];e33=[]
        for i in xrange(len(images)):
            i0s=images[i]
            i0 = i0s[0]-1

            e11.append(ex[i0])
            e22.append(ey[i0])
            e33.append(-ex[i0]-ey[i0])

        ## x-y flip for the strain
        if ieps_flip:
            dum = e11[::]
            e11 = e22[::]
            e22 = dum[::]
        fc.get_strain(e11,0,0)
        fc.get_strain(e22,1,1)
        fc.get_strain(e33,2,2)
        fc.set_zero_shear_strain()

    fc.get_eqv()
    if fig==None:
        wf=mpl_lib.wide_fig
        fig=wf(nw=1,nh=1)
        ax=fig.axes[0]
    if fig!=None:
        ax=fig.axes[iax]

    ax.plot(fc.epsilon_vm, fc.sigma_vm, 'o', label=label)
    eqv(ax=ax,ft=12)
    return fig,fc,fn_trs

    
