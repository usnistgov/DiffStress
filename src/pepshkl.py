"""
pepshkl.out (int_eps_ph?.out) file reader and plotter
===============================

# def 1
 - reader(fn='pephkl.out', step=0, phi=0, skiprows=1, iopt=1)
 - returns:
   psi1, psi2, e1, e2, phkl1, phkl2 (iopt=1)
   psi1, psi2, e1, e2, phkl1, phkl2, vol1, vol2,
                 s11, s22, e11, e22 (iopt=2)

# def 2
 - ex01(fn='pepshkl.out',istep=0,iphi=0,ax01=None,ax01x=None)
 - plots ps1 vs eps1 ps2 vs eps2 curves in the given ax01 and ax01x

# def 3
 - ex02(fn='int_eps_ph1.out',istep=0,iphi=0,ax01=None,
     ax01x=None,ax02=None,
     label=None,ls=['-x','--o'])
 - plots 1) psi1 vs eps1; 2) psi2 vs eps2; 3) psi1 vs sig1;
         4) psi2 vs sig2 curves

# def 4
 - uniq(a)
 - Return only unique elements contained in a numpy array
    with its index contained also in another numpy array
"""
#print __doc__

import numpy as np
#from rearrange_sort import main as ss
from MP import ssort
ss=ssort.ind_swap
import matplotlib.pyplot as plt

def reader(fn='pepshkl.out',step=0,phi=0,skiprows=1,iopt=1):
    """
    Arguments
    =========
    fn   = 'pepshkl.out'
    step = 0
    phi  = 0
    skiprows = 1
    iopt     = 1(pepshkl.out) 2(int_eps_ph1.out) 3(
    """
    dat = np.loadtxt(fn, skiprows=skiprows, dtype='float').T
    # steps
    istep = dat[0]
    steps, isteps = uniq(istep)
    print 'isteps:', isteps
    nsteps = len(steps)
    # phis
    phis = map(int,dat[1])
    phis, iphis = uniq(phis)
    nphis = len(phis)
    print 'nphis:', nphis
    # betas
    betas = dat[2]
    betas, ibetas = uniq(betas)
    nbetas = len(betas)
    print 'betas:', betas
    # psis1
    psis1 = dat[3]
    psis1, ipsis1 = uniq(psis1)
    # psis2
    psis2 = dat[4]
    psis2, ipsis2 = uniq(psis2)
    # eps1
    eps1 = dat[5]
    # eps2
    eps2 = dat[6]
    # phkl1
    phkls1 = dat[7]
    # phkl2
    phkls2 = dat[8]

    if iopt==2 or iopt==3:
        ngr1 = dat[9]
        ngr2 = dat[10]
        vol1 = dat[11]
        vol2 = dat[12]

        s11, idum  = uniq(dat[13])
        s22, idum  = uniq(dat[14])
        s33, idum  = uniq(dat[15])
        s23, idum  = uniq(dat[16])
        s13, idum  = uniq(dat[17])
        s12, idum  = uniq(dat[18])
        e11, idum  = uniq(dat[19])
        e22, idum  = uniq(dat[20])
        e33, idum  = uniq(dat[21])
        e23, idum  = uniq(dat[22])
        e13, idum  = uniq(dat[23])
        e12, idum  = uniq(dat[24])
    #raise IOError

    # --
    istart = isteps[step]
    istart = istart + iphis[phi]
    iend = istart + nbetas

    psi1 = dat[3][istart:iend]
    psi2 = dat[4][istart:iend]
    e1 = dat[5][istart:iend]
    e2 = dat[6][istart:iend]
    phkl1 = dat[7][istart:iend]
    phkl2 = dat[8][istart:iend]
    vol1  = vol1[istart:iend]
    vol2  = vol2[istart:iend]

    if iopt==3:
        ngr1 = ngr1[istart:iend]
        ngr2 = ngr2[istart:iend]

    if len(psi1)==0 or len(psi2)==0 or len(e1)==0 or \
            len(e2)==0 or len(phkl1)==0 or len(phkl2)==0:
        raise IOError

    if iopt==1:
        return psi1,psi2,e1,e2,phkl1,phkl2
    elif iopt==2:
        return psi1,psi2,e1,e2,phkl1,phkl2,vol1,vol2,s11,s22,e11,e22
    elif iopt==3:
        return psi1,psi2,e1,e2,phkl1,phkl2,vol1,vol2,s11,s22,e11,e22,ngr1,ngr2


def reader2(fn='igstrain_fbulk_ph1.out',
            iopt=1,isort=True,verbose=False):
    """
    Arguments
    =========
    fn = 'igstrain_fbulk_ph1.out style'
    iopt 1 or 2
      iopt1: suitable for igstrain_fbulk_ph1.out, in which
            diffraction results are repeatedly obtained for
            individual elastic loading for stress factor probing
      iopt2: ??
    verbose=False
    """
    import MP.ssort as sort
    shsort = sort.shellSort
    ss     = sort.ind_swap

    if verbose:
        print '\n\n###########################################'
        print '# Reader for igstrain_fbulk_ph?.out files #'
        print '###########################################\n'

    dstr = open(fn,'r').read()
    if iopt==1: dstr = dstr.split('-- F')[1::]
    if iopt==2: dstr = dstr.split('--')[1::]

    # analyse the structure of the file.
    isf = []; blocks = []
    nblock = 0
    incomplete = False
    for i in xrange(len(dstr)):
        rst = dstr[i].split('\n')
        ablock = rst[1:-1]
        blocks.append(ablock)
        if i==0:
            npb = len(ablock)
            phis = []
            psis = []
            for j in xrange(npb):
                phi = float(ablock[j].split()[0])
                psi = float(ablock[j].split()[1])
                if not(phi in phis): phis.append(phi)
                if not(psi in psis): psis.append(psi)

            nphi = len(phis); npsi = len(psis)

        if len(ablock)!=npb:
            incomplete=True
            #break

        nblock = nblock + 1
        if iopt==1:
            i1, i2 = rst[0].split()[0],rst[0].split()[1]; i1 = i1[0]
            i1, i2 = map(int, [i1,i2])
            isf.append([i1,i2])

    #print 'total number of data block', nblock
    if iopt==1:
        usf = []
        for i in xrange(len(isf)):
            i1,i2 = isf[i]
            if not(isf[i] in usf):
                usf.append(isf[i])
            else: break
        nsf = len(usf)
        if verbose: print nsf, 'kinds of stress factor are probed:',

        for i in xrange(nsf):
            i1,i2 = usf[i]
            if verbose: print 'F%1i%1i'%(i1,i2),

    elif iopt==2: nsf = 1

    if iopt==1:nst=nblock / nsf
    elif iopt==2: nst=nblock

    if verbose: print '\nIG strains were measured at', nst, 'steps'

    if iopt==1: tdat = np.zeros((10,nst,nsf,nphi,npsi))
    if iopt==2: tdat = np.zeros((9,nst,nphi,npsi))

    for i in xrange(nst):
        for j in xrange(nsf):
            n = nsf * i + j
            ablock = blocks[n]
            for k in xrange(nphi):
                for l in xrange(npsi):
                    ipb = npsi * k + l

                    if iopt==1:
                        try:
                            dd = map(float, ablock[ipb].split())
                            nd = len(dd)
                            if nd==11:
                                phi,psi,sin2psi,ehkl,e,ehkle,\
                                    fhkl,fbulk,ige,sij,rsq\
                                    = dd
                            elif nd==10:
                                rsq = np.nan
                                phi,psi,sin2psi,ehkl,e,ehkle,\
                                    fhkl,fbulk,ige,sij\
                                    = dd


                        except IndexError:
                            phi,psi,sin2psi,ehkl,e,ehkle,\
                                fhkl,fbulk,ige,sij,rsq\
                                = np.nan,np.nan,np.nan,\
                                np.nan,np.nan,np.nan,\
                                np.nan,np.nan,\
                                np.nan,np.nan,np.nan

                        tdat[0,i,j,k,l] = ehkl       #e(hkl,phi,psi)
                        tdat[1,i,j,k,l] = e          #macro
                        tdat[2,i,j,k,l] = ehkle      #e-macro
                        tdat[3,i,j,k,l] = fhkl       #F(hkl,phi,psi)
                        tdat[4,i,j,k,l] = fbulk      #F(bulk,phi,psi)
                        tdat[5,i,j,k,l] = ige        #e(hkl,phi,psi) - f_hkl * Sij
                        tdat[6,i,j,k,l] = sij        #Sij
                        tdat[7,i,j,k,l] = phi        #Phi
                        tdat[8,i,j,k,l] = psi        #Psi
                        tdat[9,i,j,k,l] = rsq        #rsq (goodness of the 'linear' fit)

                    elif iopt==2:
                        #return ablock
                        try:
                            phi,psi,sin2psi,ehkl,e,ehkle,ige,\
                                f11,f22,s11,s22\
                                = map(float, ablock[ipb].split())
                        except IndexError:
                            phi,psi,sin2psi,ehkl,e,ehkle,ige,f11,f22,s11,s22\
                                = np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan

                        tdat[0,i,k,l] = ehkl       #e(hkl,phi,psi)
                        tdat[1,i,k,l] = e          #macro
                        tdat[2,i,k,l] = ehkle      #e-macro
                        tdat[3,i,k,l] = ige        #e(hkl) - F^hkl_ij * Sij (ij=1,1 and 2,2)
                        tdat[4,i,k,l] = f11        #F^hkl_11
                        tdat[5,i,k,l] = f22        #F^hkl_22
                        tdat[6,i,k,l] = s11        #S11
                        tdat[7,i,k,l] = s22        #S22
                        tdat[8,i,k,l] = psi        #S22

    if isort:
        if iopt==1: ndat = 10
        elif iopt==2: raise IOError,\
             'Unexpected ioption for isort==True'

        psi = tdat[8,0,0,0,:]
        _psi_,ind = shsort(psi)
        _tdat_ = tdat.swapaxes(0,-1)[::]
        tdat = ss(_tdat_[::],ind)[::]
        tdat = tdat.swapaxes(0,-1)[::]

    if verbose: print '\ntdat(property, ist, isf, iphi, ipsi)'
    if iopt==1: return tdat,usf
    if iopt==2: return tdat
    raise IOError, 'Unexpected iopt was given'
    ##

def reader3(fn='igstrain_unloads_avg.out',isort=False):
    """
    Reader for files in the 'igstrain_unloads_avg.out' template.

    Arguments
    =========
    fn    = 'igstrain_unloads_avg.out'
    isort = False
    """
    from MP.ssort import sh as sort
    # from sff_converter import condition
    # difl, nphi, phis, nbeta, neps, eps = condition(difile)
    ds = open(fn,'r').read()
    blocks = ds.split('--')[1::]
    nstp = len(blocks)
    print 'nstp:', nstp
    # file structure
    phis = []; psis = []
    lines = blocks[0].split('\n')[1:-1]
    for il in xrange(len(lines)):
        phi,psi,eps,f11,f22,f33,f23,f13,f12 = \
            map(float,lines[il].split())

        if not(phi in phis): phis.append(phi)
        if len(phis)==1: psis.append(psi)

    dat = np.zeros((nstp,len(phis),len(psis)))
    fij = np.zeros((nstp,len(phis),len(psis),6))
    for ist in xrange(nstp):
        ablock = blocks[ist]
        lines = ablock.split('\n')[1:-1]
        for iphi in xrange(len(phis)):
            x=[]; y=[]; z=[]
            for ipsi in xrange(len(psis)):
                l = len(psis) * iphi + ipsi
                ph,ps,ep,f1,f2,f3,f4,f5,f6 = \
                    map(float, lines[l].split())
                x.append(ps); y.append(ep); z.append([f1,f2,f3,f4,f5,f6])
                if not(isort):
                    dat[ist,iphi,ipsi] = ep
                    fij[ist,iphi,ipsi,:] = z[0][:]
            if isort:
                sx, sy, sz = sort(x,y,z)
                for ix in xrange(len(sx)):
                    dat[ist,iphi,ix] = sy[ix]
                    fij[ist,iphi,ix,:] = sz[ix][:]

    if isort: psis = sort(psis)
    psis=np.array(psis)
    return dat, psis, fij

def reader4(fn='igstrain_unload_ph1.out',
            ndetector=2,iopt=0):
    """
    Arguments
    =========
    fn        = 'igstrain_unloads_ph1.out'
    ndetector = 2
    iopt      = 0 (default return) (1: return vdat,
                                   ngrd in addition)
    """
    from MP.ssort import sh as sort
    npb = int(open(fn,'r').readlines()[1].split()[0])
    ds  = np.loadtxt(fn,skiprows=2).T
    # steps
    steps = np.unique(ds[0])
    nst   = len(steps)
    # phi
    phis = np.unique(ds[1])
    nphi = len(phis)
    # npsi
    nbeta = npb/nphi
    npsi  = nbeta * ndetector
    # lines
    nl  = len(ds[0])
    dat = ds.T

    #tdat = np.zeros((len(steps),
    # collect psi1 and psi2
    psis = []; psis1 = []; psis2 = []
    for il in xrange(nbeta):
        d = dat[il]
        istep, phi, beta, psi1, psi2, eps1,\
            eps2, sig1, sig2 = d[:9]
        psis.append(psi1); psis.append(psi2)

    psis = np.array(psis)

    # IG strain array
    tdat = np.zeros((nst,nphi,npsi))
    if iopt!=0:
        # volume/number of grains
        vdat = np.zeros((nst,nphi,npsi))
        ngrd = np.zeros((nst,nphi,npsi))

    il = 0
    for ist in xrange(nst):
        for iphi in xrange(nphi):
            for ib in xrange(nbeta):
                d = dat[il]
                istep, phi, beta, psi1, psi2, \
                    eps1, eps2, sig1, sig2 = d[:9]
                tdat[ist,iphi,ib*2]   = eps1
                tdat[ist,iphi,ib*2+1] = eps2
                ## volume/ngr
                if iopt!=0:
                    n1,n2,v1,v2 = d[9:13]
                    vdat[ist,iphi,ib*2]   = v1
                    vdat[ist,iphi,ib*2+1] = v2
                    ngrd[ist,iphi,ib*2]   = n1
                    ngrd[ist,iphi,ib*2+1] = n2
                il = il + 1
    if   iopt==0: return tdat,psis
    elif iopt==1: return tdat,psis,vdat,ngrd
    else:         raise IOError, 'Unexpected iopt given'

def ex01(fn='pepshkl.out',istep=0,iphi=0,ax01=None,ax01x=None):
    """
    Arguments
    =========
    fn    = 'pepshkl.out'
    istep = 0
    iphi  = 0
    ax01  = None
    ax01x = None
    """
    psi1,psi2,eps1,eps2,phkl1,phkl2 = reader(
        fn,istep,iphi,skiprows=1,iopt=1)
    if ax01==None or ax01x==None:
        fig01 = plt.figure(1)
        ax01  = fig01.add_subplot(111)
        ax01x = ax01.twinx()
    ax01.plot(psi1,eps1,'r-x')
    l0, = ax01.plot(psi2,eps2,'b-x')
    ax01x.plot(psi1,phkl1,'r--o')
    l1, = ax01x.plot(psi2,phkl2,'b--o')
    ax01.legend([l0,l1],[r'$\varepsilon^{hkl}$',r'$p^{hkl}_{\phi\psi}$'])\
        .get_frame().set_alpha(0.5)
    ax01.set_xlabel(r'$\psi$',dict(fontsize=28))

def ex02(fn='int_eps_ph1.out',istep=0,iphi=0,ax01=None,
         ax01x=None,ax02=None,
         label=None,ls=['-x','--o'],color='r'):
    """
    Arguments
    =========
    fn    = 'int_eps_ph1.out'
    istep = 0
    iphi  = 0
    ax01  = None
    ax01x = None
    ax02  = None
    label = None
    ls    = ['-x','--o']
    """
    import MP.ssort as sort
    sort = sort.shellSort
    psi1,psi2,eps1,eps2,sig1,sig2,\
        vol1,vol2,s11,s22,e11,e22,ngr1,ngr2=\
        reader(fn,istep,iphi,skiprows=2,iopt=3)

    sin2psi1 = np.sin(psi1*np.pi/180.)**2
    sin2psi2 = np.sin(psi2*np.pi/180.)**2

    if ax01==None or ax01x==None:
        fig01 = plt.figure(1)
        ax01  = fig01.add_subplot(211)
        ax02  = fig01.add_subplot(212)
        ax01x = ax01.twinx()
    if label==None:
        sin2psi1 = np.sin(psi1*np.pi/180.)**2
        sin2psi2 = np.sin(psi2*np.pi/180.)**2

        l0, = ax01.plot(sin2psi1,eps1,'r-x')
        ax01.plot(sin2psi2,eps2,'b-x')
        l1, = ax01x.plot(sin2psi1,sig1,'r--o')
        ax01x.plot(sin2psi2,sig2,'b--o')

#        l0, = ax01.plot(psi1,eps1,'r-x')
#        ax01.plot(psi2,eps2,'b-x')
#        l1, = ax01x.plot(psi1,sig1,'r--o')
#        ax01x.plot(psi2,sig2,'b--o')
    else:
        x = sin2psi1
        x,ind = sort(x)

        eps1 = ss(ind, eps1)

        ax01.plot(sin2psi2,eps2,'b'+ls[0],label=label)

        ngr1 = ss(ind, ngr1)
        l1, = ax01x.plot(
            x,ngr1,color+ls[0],label=label)
        #ax01x.plot(sin2psi2,sig2,'b'+ls[0],label=r'')
        vol1 = ss(ind, vol1)
        l2, = ax02.plot(
            x,vol1,color+ls[0],label=label)
        #ax02.plot(sin2psi2,vol2,'b'+ls[0],label=r'')

#        l0, = ax01.plot(psi1,eps1,'r'+ls[0],label=label)
#        ax01.plot(psi2,eps2,'b'+ls[0],label=r'')
#        l1, = ax01x.plot(psi1,sig1,'r'+ls[1],label=label)
#        ax01x.plot(psi2,sig2,'b'+ls[0],label=r'')
#        l2, = ax02.plot(psi1,vol1,'r'+ls[0],label=r'')
#        ax02.plot(psi2,vol2,'b'+ls[0],label=r'')

    if ax01==None or ax01x==None:
        ax01.legend([l0,l1],[r'$\varepsilon^2$',r'$\sigma^2$'])\
            .get_frame().set_alpha(0.5)
        ax01.set_xlabel(r'$\sin^2{\psi}$',dict(fontsize=28))
    # ax02.legend([l0,l1],[r'$\varepsilon^2$',r'$\sigma^2$'])\
    #     .get_frame().set_alpha(0.5)

def ex03(fn='int_eps_ph1.out',iphi=0,ifig=1):
    """
    Arguments
    =========
    fn   = 'int_eps_ph1.out'
    iphi = 0
    ifig = 1
    """
    # fig01 = plt.figure(ifig,figsize=(10,10))
    # ax01  = fig01.add_subplot(321)
    # ax01x = fig01.add_subplot(323)
    # ax02  = fig01.add_subplot(325)
    fig01=plt.figure(ifig, figsize=(6.5,5));  ax01 =fig01.add_subplot(111)
    fig02=plt.figure(ifig+1, figsize=(6.5,5));ax02 =fig02.add_subplot(111)
    fig03=plt.figure(ifig+2, figsize=(6.5,5));ax03 =fig03.add_subplot(111)

    steps = np.loadtxt(fn,skiprows=2).T[0]
    u, idum = uniq(steps)
    s11   = np.loadtxt(fn,skiprows=2).T[13]
    s11, idum = uniq(s11)
    s22   = np.loadtxt(fn,skiprows=2).T[14]
    s22, idum = uniq(s22)
    e11   = np.loadtxt(fn,skiprows=2).T[19]
    e11, idum = uniq(e11)
    e22   = np.loadtxt(fn,skiprows=2).T[20]
    e22, idum = uniq(e22)

    l00 = []
    l11 = []

    ls=[['-x','-x'],['-+','-+'],['-,','-,'],['-d','-d'],['-v','-v'],
        ['->','->'],['-2','-2'],['-4','-4'],['--o','--o'],['--.','-+'],
        ['-^','-,'],[':','-d'],['-<','-v'],['-1','->'],['-3','-2'],
        ['-4','-4']]
    colors = ['k','k','k','k','k','k','k',
              'b','k','m','r','g','b','k',
              'm','r','g','b','k','m']

    for i in xrange(len(u)):
#        ex02(fn,i,iphi,ax01,ax01x,ax02,
        dum = e11[i]+e22[i]#)#*2
        if i==0: dum = +0.
        if dum<0: dum = 0
        ex02(fn,i,iphi,ax01,ax02,ax03,
             label=r'$\bar{E}^{\mathrm{eff}}=$%4.2f'%dum,
             ls=ls[i],color=colors[i])
        ax03.legend(loc='best')

#    ax01.legend(bbox_to_anchor=(1.00,1),loc=2,
#                fancybox=True).get_frame().set_alpha(0.5)
    #ax01.legend(fancybox=True).get_frame().set_alpha(0.5)

    ax01.set_xlabel(r'$\sin^2{\psi}$',dict(fontsize=28))
    ax02.set_xlabel(r'$\sin^2{\psi}$',dict(fontsize=28))
    ax03.set_xlabel(r'$\sin^2{\psi}$',dict(fontsize=28))
    ax01.set_ylabel(r'$\varepsilon^{hkl} [\mu\varepsilon]$',dict(fontsize=28))
    ax02.set_ylabel('Number of grains',dict(fontsize=28))
    ax03.set_ylabel('Volume fraction',dict(fontsize=28))
    ax01.grid('on');ax02.grid('on');ax03.grid('on');
    plt.tight_layout()
    plt.show()

def uniq(a):
    """
    Arguemtns
    =========
    a

    Return only unique elements contained in a numpy array
    with its index contained also in another numpy array
    """
    u   = []
    ind = []
    for i in xrange(len(a)):
        if a[i] in u: pass
        else:
            u.append(a[i])
            ind.append(i)
    return np.array(u), np.array(ind)

def ex_igb(fn='igstrain_fbulk_ph1.out',ifig=2,iphi=0,isf=0,
           mxnst=None,flow=None):
    """
    Read ig strain file and analyze... (igstrain_fbulk_ph1.out)

    igstrain_fbulk_ph1.out contains diffraction for each and every
    reloading from 'unloaded' polycrystal.

    Arguments
    =========
    fn = 'igstrain_fbulk_ph1.out'
    ifig = 2
    iphi=0
    isf=0
    mxnst=None
    """
    import matplotlib.pyplot as plt
    from sff_converter import condition
    from MP.ssort import sh as sort
    from MP.mat import voigt
    from MP.lib.mpl_lib import wide_fig as wf

    difl, nphi, phis, nbeta, dum1, dum2 = condition(fn=None)

    if flow==None: raise IOError, 'Flow stress is missing'
    eps = flow.epsilon_vm[::]
    neps = flow.nstp

    nw = 4
    remainder = np.mod(neps,nw)
    if remainder==0:
        nh = neps/nw
    elif remainder>0:
        nh = neps/nw + 1

    fig  =wf(nw=2, nh=2,ifig=ifig,iarange=True);  axes  =fig.axes
    ax01, ax02, ax03, ax04 = axes
    fig01=wf(nw=nw,nh=nh,ifig=ifig+1,iarange=True,nfig=neps);axes01=fig01.axes
    fig02=wf(nw=nw,nh=nh,ifig=ifig+2,iarange=True,nfig=neps);axes02=fig02.axes
    fig03=wf(nw=nw,nh=nh,ifig=ifig+3,iarange=True,nfig=neps);axes03=fig03.axes
    fig04=wf(nw=nw,nh=nh,ifig=ifig+4,iarange=True,nfig=neps);axes04=fig04.axes

    tdat, usf = reader2(fn,iopt=1)
    print tdat.shape
    nsf = len(usf); nst = len(tdat[0]); npsi = len(tdat[0,0,0,0,:])
    if nst!=neps: raise IOError, 'Inconsistency between neps and nst'
    markers = ['o','x','+','^','d','*','o','x','+','^','d','*''o','x','+','^','d','*']
    colors  = ['r','g','b','k','m','y','r','g','b','k','m','y','r','g','b','k','m','y']
    from string import ascii_lowercase as alphabet
    axe_lab=[]
    for ab in alphabet: axe_lab.append('(%s)'%ab)
    yl = 0; yh = 0
    avg = np.zeros((2,nst,npsi)) #
    sft  = np.zeros((6,nst,npsi)) # f^{hkl}
    rsqt = np.zeros((6,nst,npsi)) # f^{hkl}

    if mxnst!=None: nst=mxnst
    for i in xrange(nst):
        ehkl  = tdat[0,i,isf,iphi,:] # e(hkl,phi,psi)
        e     = tdat[1,i,isf,iphi,:] # macro
        ehkle = tdat[2,i,isf,iphi,:] # e - macro
        fhkl  = tdat[3,i,isf,iphi,:] # Fhkl
        fbulk = tdat[4,i,isf,iphi,:] # Fbulk
        ige   = tdat[5,i,isf,iphi,:] # e - F_ij *Sij
        sij   = tdat[6,i,isf,iphi,:] # Sij
        psi   = tdat[8,i,isf,iphi,:]
        rsq   = tdat[9,i,isf,iphi,:]
        sin2psi = np.sin(psi*np.pi/180.)**2
        axes01[i].set_title(r'%s $\bar{E}^{\mathrm{eff}} = %4.2f$'%(
                axe_lab[i],eps[i]))
        axes02[i].set_title(r'%s $\bar{E}^{\mathrm{eff}} = %4.2f$'%(
                axe_lab[i],eps[i]))
        axes03[i].set_title(r'%s $\bar{E}^{\mathrm{eff}} = %4.2f$'%(
                axe_lab[i],eps[i]))
        axes04[i].set_title(r'%s $\bar{E}^{\mathrm{eff}} = %4.2f$'%(
                axe_lab[i],eps[i]))
        igys = []; igos = []; igts = []
        for j in xrange(nsf): ## Number of elastic loading for stress factor calc.
            ehkl_ = tdat[0,i,j,iphi,:] # e(hkl,phi,psi)
            e_    = tdat[1,i,j,iphi,:] # macro
            fhkl_ = tdat[3,i,j,iphi,:] # Fhkl
            fbulk_= tdat[4,i,j,iphi,:] # Fbulk
            sij_  = tdat[6,i,j,iphi,:] # Sij
            rsq_  = tdat[6,i,j,iphi,:] # R^2

            sft[j,i,:]  = fhkl_[::] #
            rsqt[j,i,:] = rsq_[::]

            i1,i2 = usf[j]
            igo = ehkl_ - e_ ## old way
            igt = []
            m = 0
            for k in xrange(len(fbulk_)):
                if fbulk_[k]==0:
                    m=m+1
                    igt.append(0.)
                else: igt.append(
                        ehkl_[k] - \
                        fhkl_[k] / fbulk_[k] * e_[k])

            if m!=0: print m, 'case of fbulk'+\
               '(phi,psi)=0 happened'

            igt = np.array(igt)
            igy = ehkl_ - fhkl_ * sij_ ## new way (YJ)
            if j==0:
                igy_avg = igy/nsf
                igo_avg = igo/nsf
                igt_avg = igt/nsf
            else:
                igy_avg = igy_avg + igy/nsf
                igo_avg = igo_avg + igo/nsf
                igt_avg = igt_avg + igt/nsf
            igys.append(igy);igos.append(igo);igts.append(igt)

            # The old method
            if j==4: m = markers[j]
            else: m= '--'
            y = igo
            x = np.sin(psi[:]*np.pi/180.)**2
            newx,newy = sort(x,y) # shell sort
            axes01[i].plot(
                newx,newy*1e6,m,color=colors[j],
                label=r'$\Sigma_{%1i%1i}$'%(i1,i2))
            if i==0:axes01[i].legend(loc='best',
                fancybox=True).get_frame().set_alpha(0.5)
            if j==0:axes01[i].grid('on')
            yl0, yh0 = axes01[i].get_ylim()
            if yl0<yl: yl = yl0
            if yh0>yh: yh = yh0

            # GHT
            y = igt
            x = np.sin(psi[:]*np.pi/180.)**2
            newx,newy = sort(x,y) # shell sort
            axes02[i].plot(
                newx,newy*1e6,m,color=colors[j],
                label=r'$\Sigma_{%1i%1i}$'%(i1,i2))
            if i==0:axes02[i].legend(loc='best',
                fancybox=True).get_frame().set_alpha(0.5)
            if j==0:axes02[i].grid('on')
            yl0, yh0 = axes02[i].get_ylim()
            if yl0<yl: yl = yl0
            if yh0>yh: yh = yh0

            # YJ
            y = igy
            x = np.sin(psi[:]*np.pi/180.)**2
            newx,newy = sort(x,y) # shell sort
            axes03[i].plot(
                newx,newy*1e6,m,color=colors[j],
                label=r'$\Sigma_{%1i%1i}$'%(i1,i2))
            if i==0:axes03[i].legend(loc='best',
                fancybox=True).get_frame().set_alpha(0.5)
            if j==0:axes03[i].grid('on')
            yl0, yh0 = axes03[i].get_ylim()
            if yl0<yl: yl = yl0
            if yh0>yh: yh = yh0

            pass # over nsf

        # errors for 6 cases of SF loading
        igos = np.array(igos).swapaxes(0,1)
        igys = np.array(igys).swapaxes(0,1)
        igts = np.array(igts).swapaxes(0,1)
        igoe=[]; igye=[]; igte=[]
        for m in xrange(len(igo_avg)):
            igoe.append(igos[m].std())
            igte.append(igts[m].std())
            igye.append(igys[m].std())
        igoe,igte,igye = np.array(igoe),np.array(igte),np.array(igye)
        ##

        # average old
        y = igo_avg
        ye= igoe
        x = np.sin(psi[:]*np.pi/180.)**2
        newx, newy, newye = sort(x,y,ye)
        axes04[i].errorbar(
            x=newx,y=newy*1e6,yerr=newye*1e6,
            ls='-',color='r',label='old')

        # average YU
        y  = igy_avg
        avg[0,i,:] = psi[:]; avg[1,i,:] = y[:]

        ye = igye
        x  = np.sin(psi[:]*np.pi/180.)**2

        newx, newy, newye = sort(x,y,ye)

        axes04[i].errorbar(
            x=newx,y=newy*1e6,yerr=newye*1e6,
            ls='-',ms=3,color='b',label='YJ')
        axes04[i].grid('on')
        axes04[0].legend(loc='best',fancybox=True)\
            .get_frame().set_alpha(0.5)

        i1,i2 = voigt.ijv[:,isf]
        sin2psi = np.sin(psi*np.pi/180.)**2

        ax01.plot(sin2psi,ige*1e6,markers[i]#,color='k',
           ,label=r'$\bar{E}^{\mathrm{eff}} = %4.2f$'%eps[i])
        ax02.plot(sin2psi,ehkle*1e6,markers[i]#,color='k'
           ,label=r'$\bar{E}^{\mathrm{eff}} = %4.2f$'%eps[i])
        l, = ax04.plot(sin2psi,fhkl*1e6,markers[i]#,color='k'
           ,label=r'$\bar{E}^{\mathrm{eff}} = %4.2f$'%eps[i])
        if i==0: ax04.plot(sin2psi,fbulk*1e6,'-',alpha=0.2,color='k')
        # ax03.plot(sin2psi,e*1e6,markers[i],
        #           label=r'$\bar{E}^{\mathrm{eff}} = %4.2f$'%eps[i])
        ax03.plot(sin2psi,(ehkl-fhkl/fbulk*e)*1e6,markers[i]
           ,label=r'$\bar{E}^{\mathrm{eff}} = %4.2f$'%eps[i])
        sin2psi = np.sin(psi*np.pi/180.)**2

    ## Deco the axes
    for i in xrange(nst):
        axes01[i].set_ylim(yl,yh)
        axes02[i].set_ylim(yl,yh)
        axes03[i].set_ylim(yl,yh)
        #axes04[i].set_ylim(yl,yh)
        axes01[i].set_xlim(0.,); axes02[i].set_xlim(0.,)
        axes03[i].set_xlim(0.,); axes04[i].set_xlim(0.,)
    for i in xrange(len(axes)):
        axes[i].set_xlabel(r'$\sin^2{\psi}$', dict(fontsize=20))
        axes01[i].set_xlabel(r'$\sin^2{\psi}$', dict(fontsize=20))
        axes02[i].set_xlabel(r'$\sin^2{\psi}$', dict(fontsize=20))
        axes03[i].set_xlabel(r'$\sin^2{\psi}$', dict(fontsize=20))
        axes04[i].set_xlabel(r'$\sin^2{\psi}$', dict(fontsize=20))
    for i in xrange(len(axes)): axes[i].grid('on')
    for i in xrange(len(axes)): axes[i].set_xlim(0.,0.5)
    ax01.set_ylabel(
        r'$\varepsilon^{hkl} - F^{hkl}_{ij}\bar{\Sigma}_{ij}$'+\
        r'  $[\mu\varepsilon]$',dict(fontsize=15))
    ax02.set_ylabel(
        r'$\varepsilon^{hkl} - E$  $[\mu\varepsilon]$',dict(fontsize=15))
    ax03.set_ylabel(r'$\varepsilon^{hkl} - $'+
        r'$F_{ij}^{hkl}/F^{\mathrm{bulk}}_{ij}$'+
        r'$ \bar{E} $  $[\mu \varepsilon]$',dict(fontsize=15))
    ax04.set_ylabel(r'$F_{ij}^{hkl}, F^{\mathrm{bulk}}_{ij}$  $[TPa^{-1}]$',
        dict(fontsize=15))
    ax01.set_title(r'(a) $\varepsilon(hkl,%3.0f^\circ,\psi)$'\
        r'$ - F_{%1i%1i}(hkl,%3.0f^\circ,\psi)$'\
        r'$\bar{\Sigma}_{%1i%1i}$'%(phis[iphi],i1,i2,phis[iphi],i1,i2),
        dict(fontsize=13))
    ax02.set_title(
        r'(b) $\varepsilon(hkl,%3.0f^\circ,\psi)$'%(phis[iphi])+\
            r'$ - \bar{E}(%3.0f^\circ,\psi)$'%\
            (phis[iphi]),dict(fontsize=13))
    ax03.set_title(
        r'(c) $\varepsilon(hkl,%3.0f^\circ,\psi)$'%(phis[iphi])+
        r'$-F_{%1i%1i}(hkl,%3.0f^\circ,\psi)/$'%(i1,i2,phis[iphi])+
        r'$F_{%1i%1i}^{\mathrm{bulk}}(%3.0f^\circ,\psi)$'%(
            i1,i2,phis[iphi])+
        r'$\bar{E}(%3.0f^\circ,\psi)$'%(phis[iphi]),dict(fontsize=13))
    ax04.set_title(r'(d) $F^{hkl}_{%1i%1i}$ and '%(i1,i2)+\
                       r'$F^{\mathrm{bulk}}_{%1i%1i}$'%
                   (i1,i2), dict(fontsize=13))
    ax04.legend(loc='lower right',fancybox=True).get_frame().set_alpha(0.5)
    # fig.tight_layout(); fig01.tight_layout(); fig02.tight_layout()
    # fig03.tight_layout(); fig04.tight_layout()
    fig.savefig('ig_bulk.pdf');fig01.savefig('ig_bulk_Old.pdf')
    fig02.savefig('ig_bulk_GHT.pdf'); fig03.savefig('ig_bulk_YJ.pdf')
    fig04.savefig('ig_bulk_avg.pdf')
    ## save the averaged IG strain

    return avg, sft # avg[2,nst,npsi], sft[6,nst,npsi]

def ex_igb_t(fn='igstrain_fbulk_ph1.out',
             fnout='igstrain_unloads_avg.out',
             mxnst=None,flow=None):
    """

    Arguments
    =========
    fn   ='igstrain_fbulk_ph1.out'
    fnout='igstrain_unloads_avg.out'
    mxnst=None
    """
    import matplotlib.pyplot as plt
    from sff_converter import condition
    import numpy as np
    difl, nphi, phis, nbeta, dum1, dum2 = condition(fn=None)
    if flow==None: raise IOError, 'Flow should be given'
    elif flow!=None:
        flow.get_eqv()
        eps = flow.epsilon_vm[::]

    avg = []; sfs = []; rsqs = [];
    for i in xrange(nphi):
        a, sf, rsq = ex_igb(
            fn=fn,ifig=i*7+1,iphi=i,isf=0,
            mxnst=mxnst,flow=flow)# a[2,nst,npsi], sf[6,nst,npsi]
        avg.append(a) # average ieps0
        sfs.append(sf)
        rsqs.append(rsq)

    avg = np.array(avg) #avg[nphi,2,nst,npsi]
    avg = avg.swapaxes(0,1)    #avg[2,nphi,nst, npsi]
    avg = avg.swapaxes(1,2)    #avg[2,nst, nphi,npsi]

    sfs = np.array(sfs) #sfs[nphi,6,nst,npsi]
    sfs = sfs.swapaxes(1,3) # [nphi,npsi,nst, 6]
    sfs = sfs.swapaxes(0,2) # [nst, npsi,nphi,6]
    sfs = sfs.swapaxes(1,2) # [nst, nphi,npsi,6]

    rsqs = np.array(rsqs)     #rsqs[nphi,6,nst,npsi]
    rsqs = rsqs.swapaxes(1,3) # [nphi,npsi,nst, 6]
    rsqs = rsqs.swapaxes(0,2) # [nst, npsi,nphi,6]
    rsqs = rsqs.swapaxes(1,2) # [nst, nphi,npsi,6]

    f = open(fnout,'w')
    f.writelines('%3s %8s %13s %14s %14s %14s %14s %14s %14s\n'
                 %('phi','psi','eps0','F11','F22','F33','F23','F13','F12'))

    if mxnst==None: nst=len(sfs)
    else: nst=mxnst

    for ist in xrange(nst):
        f.writelines('-- eps^{eff}: %5.3f\n'%eps[ist])
        for iphi in xrange(len(avg[0][ist])):
            phi = phis[iphi]
            psi = avg[0,ist,iphi,:]
            ep  = avg[1,ist,iphi,:]
            for i in xrange(len(psi)):
                ps  = psi[i]
                e   = ep[i]
                fij = sfs[ist,iphi,i,:]
                f.writelines('%3.0f %+7.2f %+12.7e'%(phi,ps,e))
                for j in xrange(len(fij)):
                    f.writelines(' %+12.7e'%fij[j])
                f.writelines('\n')
    f.close()
    return avg, sfs, rsqs  #avg[2,nst,nphi,npsi], sfs[

def ex_igb_bix(fn='igstrain_bix_ph1.out',ifig=1,iphi=0,
               mxnst=None,flow=None):
    """
    """
    import matplotlib.pyplot as plt
    from sff_converter import condition
    from MP.ssort import sh as sort
    difl, nphi, phis, nbeta, neps, eps = condition(fn=None)

    if flow!=None: eps = flow.epsilon_vm[::]
    else:          eps = eps * 2. ## Under condition balanced biaxial

    fig01 = plt.figure(ifig,figsize=(15.5,8.5))
    ax01=fig01.add_subplot(231)
    ax02=fig01.add_subplot(232)
    ax03=fig01.add_subplot(233)
    ax04=fig01.add_subplot(234)
    ax05=fig01.add_subplot(235)
    axes01 = [ax01,ax02,ax03,ax04,ax05]

    fig02 = plt.figure(ifig+1,figsize=(15.5,8.5))
    ax01=fig02.add_subplot(231)
    ax02=fig02.add_subplot(232)
    ax03=fig02.add_subplot(233)
    ax04=fig02.add_subplot(234)
    ax05=fig02.add_subplot(235)
    axes02 = [ax01,ax02,ax03,ax04,ax05]

    tdat = reader2(fn,iopt=2,isort=False)
    nst = len(tdat[0])
    npsi=len(tdat[0,0,0])

    rst = np.zeros((2,nst,nphi,npsi))
    sf2 = np.zeros((2,nst,nphi,npsi)) # f11 and f22

    markers=['o','x','+','^','d','*','o','x','+','^','d','*']
    colors =['r','g','b','k','m','y','r','g','b','k','m','y']
    axe_lab = [r'(a) $\varepsilon-F_{ij}\bar{\Sigma}_{ij}$',
               r'(b) $\bar{E}^{el}(\phi,\psi)$',
               r'(c) $\varepsilon^{hkl}$ and $F_{ij}\bar{\Sigma}_{ij}$',
               r'(d) $F_{ij}$',
               r'(e) $\bar{\Sigma}_{ij}$','(f)','(g)']

    if mxnst!=None: nst=mxnst

    for i in xrange(nst):
        ehkl  = tdat[0,i,iphi,:] #e(hkl)
        e     = tdat[1,i,iphi,:] #macro
        ehkle = tdat[2,i,iphi,:] #e-macro
        ige   = tdat[3,i,iphi,:] #e(hkl) - F^hkl_ij * Sij (ij=1,1 and 2,2)
        f11   = tdat[4,i,iphi,:] #F^hkl_11
        f22   = tdat[5,i,iphi,:] #F^hkl_22
        s11   = tdat[6,i,iphi,:] #S11
        s22   = tdat[7,i,iphi,:] #S22
        psi   = tdat[8,i,iphi,:] #psi
        sin2psi = np.sin(psi*np.pi/180.)**2

        x = np.sin(psi*np.pi/180.)**2
        y0 = ehkl-e; y1=e; y2=f11*s11+f22*s22; y3=f11; y4=f22;
        y5 = s11; y6= s22; y7=ehkl
        newx, Y0,Y1,Y2,Y3,Y4,Y5,Y6,Y7 = sort(
            x,y0,y1,y2,y3,y4,y5,y6,y7)
        ## ax01
        l, = axes01[0].plot(newx,Y0*1e6,markers[i])
        dum = ehkl - (f11*s11 + f22*s22)
        x = np.sin(psi*np.pi/180.)**2
        newx, newy = sort(x, dum)
        axes01[0].plot(newx,newy*1e6,'-',color=l.get_color())
        ## ax02
        axes01[1].plot(newx,Y1*1e6,markers[i],
           label=r'$\bar{E}^{\mathrm{eff}}$=%4.2f'%eps[i])
        ## ax03
        l, = axes01[2].plot(newx,Y7*1e6,markers[i])
        axes01[2].plot(newx,Y2*1e6,'-',color=l.get_color())
        ## ax04
        l, = axes01[3].plot(newx,Y3*1e6,markers[i])
        axes01[3].plot(newx,Y4*1e6,'-',color=l.get_color())
        ## ax05
        l, = axes01[4].plot(newx,Y5,markers[i])
        axes01[4].plot(newx,Y6,'-',color=l.get_color())
        if i==0:
            axes01[0].text(0.25,0.1,
                r'lines: $\varepsilon-(F_{11}\bar{\Sigma}_{11}$'+\
                r'$+F_{22}\bar{\Sigma}_{22})$',
                transform=axes01[0].transAxes)
            axes01[0].text(0.25,0.2,r'symbols: $\varepsilon-\bar{E}^{el} $',
                           transform = axes01[0].transAxes)
            axes01[2].text(0.6,0.3,r'symbols: $\varepsilon^{hkl}$',
                transform=axes01[2].transAxes)
            axes01[2].text(0.6,0.2,r'lines: $F_{ij}\bar{\Sigma_{ij}}$',
                 transform=axes01[2].transAxes)
            axes01[3].text(0.5,0.5,r'symbols: $F_{11}$',transform=
                           axes01[3].transAxes)
            axes01[3].text(0.5,0.4,r'lines: $F_{22}$',transform=
                           axes01[3].transAxes)
            axes01[4].text(0.5,0.5,r'symbols: $\bar{\Sigma}_{11}$',
                           transform=axes01[4].transAxes)
            axes01[4].text(0.5,0.4,r'lines: $\bar{\Sigma}_{22}$',
                           transform=axes01[4].transAxes)
            axes01[0].set_ylabel(r'IG strain $[\mu\varepsilon]$',
                                 dict(fontsize=15))
            axes01[1].set_ylabel(r'Macro strain $[\mu\varepsilon]$',
                                 dict(fontsize=15))
            axes01[2].set_ylabel(r'$\varepsilon^{hkl}(\phi,\psi)$'+\
                                     r' and $F_{ij}\bar{\Sigma}_{ij}$'+\
                                     r' $[\mu\varepsilon]$',
                                 dict(fontsize=15))
            axes01[3].set_ylabel(r'$F_{11}$ and $F_{22}$'+\
                                     r' $ [TPa^{-1}]$',
                                 dict(fontsize=15))
            axes01[4].set_ylabel(r'$\bar{\Sigma}_{11}$'+\
                                     r' and $\bar{\Sigma}_{22}$ [MPa]',
                                 dict(fontsize=15))

        for j in xrange(nphi):
            ehkl  = tdat[0,i,j,:] #e(hkl)
            e     = tdat[1,i,j,:] #macro
            ehkle = tdat[2,i,j,:] #e-macro
            ige   = tdat[3,i,j,:] #e(hkl) - F^hkl_ij * Sij (ij=1,1 and 2,2)
            f11   = tdat[4,i,j,:] #F^hkl_11
            f22   = tdat[5,i,j,:] #F^hkl_22
            s11   = tdat[6,i,j,:] #S11
            s22   = tdat[7,i,j,:] #S22
            psi   = tdat[8,i,j,:] #psi

            dum = ehkl - (f11*s11 + f22*s22)
            rst[0,i,j,:] = psi[:]
            rst[1,i,j,:] = dum[:]

            sf2[0,i,j,:]=f11[:]
            sf2[1,i,j,:]=f22[:]

            sin2psi = np.sin(psi*np.pi/180.)**2
            x = sin2psi
            y = ige*1e6

            newx, newy = sort(x,y)
            axes02[j].plot(newx,newy,'-')

            if i==0:
                axes02[j].set_title(r'$\phi=%3.0f^\circ$'%phis[j])
#    axes01[0].legend(loc='best',fancybox=True).get_frame().set_alpha(0.5)
    axes01[1].legend(loc='best',fancybox=True).get_frame().set_alpha(0.5)

    for i in xrange(len(axes01)):
        #axes01[i].set_title(axe_lab[i])
        axes01[i].set_xlabel(r'$\sin^2{\psi}$', dict(fontsize=20))
        axes01[i].grid('on')

    for i in xrange(len(axes02)):
        axes02[i].set_xlabel(r'$\sin^2{\psi}$', dict(fontsize=20))
        axes02[i].set_ylabel(r'IG strain $[\mu\varepsilon]$', dict(fontsize=20))
        axes02[i].grid('on')
    fig01.tight_layout()
    fig02.tight_layout()
    fig01.savefig('bix_analysis.pdf')
    fig02.savefig('bix_eps0_at_phis.pdf')
    return rst, sf2 #rst[2,nst,nphi,npsi], sf2[2,nst,nphi,npsi]


def pub_plot(ifig=10):
    """
    plot experimental and model-predicted IG strains
    """
    import sff_plot
    import matplotlib.pyplot as plt
    from MP.ssort import sh as sort

    fij,e0,phis,psis,exx = sff_plot.reader(fn='exp_dat/Bsteel/YJ_Bsteel_BB.sff')
    print phis
    # fij[nstr,nphi,npsi,6)
    # e0[nstr,nphi,npsi)
    # exx[nst]


    fig=plt.figure(ifig,figsize=(15.5,8.5/2.))
    ax=fig.add_subplot(131)
    ax1=fig.add_subplot(132)
    ax2=fig.add_subplot(133)

    ## Experimental IG strain
    for ist in xrange(len(exx)-1):
        if ist==0: pass
        else:
            ax.plot(np.sin(psis*np.pi/180.)**2,
                    e0[ist,0,:]*1e6,'-x',
                    label=r'$\bar{E}^{\mathrm{eff}}$=%3.1f'%(
                    exx[ist]*2))
            ax1.plot(np.sin(psis*np.pi/180.)**2,
                     e0[ist,1,:]*1e6,'-x',
                     label=r'$\bar{E}^{\mathrm{eff}}$=%3.1f'%(
                    exx[ist]*2))
            ax2.plot(np.sin(psis*np.pi/180.)**2,
                     e0[ist,2,:]*1e6,'-x',
                     label=r'$\bar{E}^{\mathrm{eff}}$=%3.1f'%(
                    exx[ist]*2))


    ax.legend(loc='best',fancybox=True,fontsize=12)
    ax.set_ylabel(r'IG strain $[\mu\varepsilon]$',
                  dict(fontsize=15))
    ax.set_xlabel(r'$\sin^2{\psi}$', dict(fontsize=18))
    ax.grid('on')
    fig.tight_layout()


    # ## Model IG strain
    # rst,sf2 = ex_igb_bix(mxnst=4,fn='igstrain_bix_ph1.out')
    # nst = len(rst[0][0])
    # for ist in xrange(nst):
    #     psi = rst[0,ist,0,:]
    #     x = np.sin(psi*np.pi/180.)**2
    #     y = rst[1,ist,0,:]
    #     x,y = sort(x,y)
    #     ax2.plot(x,y*1e6)

    # ax2.grid('on')
    # ax2.set_ylabel(r'IG strain $[\mu\varepsilon]$',
    #                dict(fontsize=15))
    # ax2.set_xlabel(r'$\sin^2{\psi}$', dict(fontsize=18))
    # fig.tight_layout()


def ex_igb_bix_t(fn='igstrain_bix_ph1.out',
                 fnout='igstrain_loads_avg.out',
                 flow=None):
    """
    Default input fn is 'igstrain_bix_ph1.out':
     this output file contains diffraction data measured
     prior to the 'interruption' by unloading. Note that
     this unloading is to elastically and uniaixally
     deform the polycrystal in order to obtain 'stress factor'
    """
    import matplotlib.pyplot as plt
    from sff_converter import condition
    import numpy as np
    difl, nphi, phis, nbeta, dum, eps = condition(fn=None)

    if flow!=None: eps = flow.epsilon_vm[::]
    else: eps = eps * 2

    #rst[2,nst,nphi,npsi], sf2[2,nst,nphi,npsi]
    rst, sf2 = ex_igb_bix(fn=fn,ifig=123,flow=flow)

    f   = open(fnout,'w')
    f.writelines('%3s %8s %13s %14s %14s %14s %14s %14s %14s\n'
                 %('phi','psi','eps0','F11','F22',
                   'F33','F23','F13','F12'))

    for ist in xrange(len(rst[0])):
        f.writelines('-- eps^{eff}: %5.3f\n'%eps[ist])
        for iphi in xrange(nphi):
            phi = phis[iphi]
            psi = rst[0,ist,iphi,:]
            ep = rst[1,ist,iphi,:]
            for i in xrange(len(psi)):
                ps = psi[i]
                e  = ep[i]
                f11 = sf2[0,ist,iphi,i]
                f22 = sf2[1,ist,iphi,i]
                f.writelines('%3.0f %+7.2f %+12.7e'%(phi,ps,e))
                f.writelines(' %+12.7e %+12.7e'%(f11,f22))
                for j in xrange(4): f.writelines(' %+12.7e'%0.)
                f.writelines('\n')
    f.close()
    plt.close('all')
    return rst

def ex_igb_cf(fnu='igstrain_fbulk_ph1.out',
              fnl='igstrain_bix_ph1.out'):
    """
    compare the ig strain at load and unloads.
    """
    import matplotlib.pyplot as plt
    from sff_converter import condition
    import numpy as np
    from MP.ssort import sh as sort
    difl, nphi, phis, nbeta, neps, eps = condition(fn=None)

    plt.ioff()
    avgu,sfs,rsqs=ex_igb_t(fn=fnu,fnout='igstrain_unloads_avg.out')   # at unlaods
    avgl=ex_igb_bix_t(fn=fnl,fnout='igstrain_loads_avg.out') # at loads

    print avgu.shape#[2,nst,nphi,nspi]
    print avgl.shape#[2,nst,iphi,npsi]

    plt.close('all'); plt.ion()

    fig = plt.figure(2,figsize=(15.5,8.5))
    fig02 = plt.figure(203); ax_dum=fig02.add_subplot(111)
    at = ['(a)','(b)','(c)','(d)','(e)','(f)']
    markers=['x','o','.','d','*']
    ls     =['--x','--o','--.','--d','--*']
    axes=[]
    for i in xrange(nphi):
        axes.append(fig.add_subplot(2,3,i+1))
    for ist in xrange(len(avgu[0])):
        for iphi in xrange(len(avgu[0][ist])):
            ax = axes[iphi]
            phi = phis[iphi]
            psi = avgu[0,ist,iphi,:]
            eu  = avgu[1,ist,iphi,:]
            el  = avgl[1,ist,iphi,:]
            sin2psi = np.sin(psi*np.pi/180.)**2

            x = sin2psi
            y0 = eu*1e6
            y1 = el*1e6
            newx,Y0,Y1 = sort(x,y0,y1)
            if iphi==0:
                l, = ax.plot(newx,Y0,'x',
                    label=r'$\bar{E}^{\mathrm{eff}}=%4.2f$'%(eps[ist]*2.0))
                ax_dum.plot(newx,Y0,markers[ist],color='k',
                            label=r'$\bar{E}^{\mathrm{eff}}=%4.2f$'%(eps[ist]*2.0))
                ax_dum.plot(newx,Y1,ls[ist],color='k')
            elif iphi!=0:
                l, = ax.plot(newx,Y0,'x')
            ax.plot(newx,Y1,'-',color=l.get_color())

    ax_dum.set_xlabel(r'$\sin^2{\psi}$', dict(fontsize=20))
    ax_dum.set_ylabel(r'IG strains $[\mu\varepsilon]$', dict(fontsize=20))
    ax_dum.legend(loc='best',fancybox=True).get_frame().set_alpha(0.5)
    ax_dum.text(0.1, 0.85, 'Lines: at loads',transform=ax_dum.transAxes)
    ax_dum.text(0.1, 0.92, 'Symbols: at unloads',transform=ax_dum.transAxes)
    ax_dum.grid('on')
    fig02.tight_layout()
    fig02.savefig('ige_cf_phi1.pdf')

    for i in xrange(len(axes)):
        ax = axes[i]
        ph = phis[i]
        #ax.set_title(r'%s $\phi=%3.0f^\circ$'%(at[i],ph))
        ax.set_xlabel(r'$\sin^2{\psi}$', dict(fontsize=20))
        ax.set_ylabel(r'IG strains $[\mu\varepsilon]$', dict(fontsize=20))
        if i==0:
            ax.legend(loc='best',fancybox=True).get_frame().set_alpha(0.5)
            ax.text(0.1, 0.85, 'Lines: at loads',transform=ax.transAxes)
            ax.text(0.1, 0.92, 'Symbols: at unloads',transform=ax.transAxes)
        ax.grid('on')
    fig.tight_layout()
    fig.savefig('ige_cf.pdf')


# def igb_reader(fns=['igstrain_unlaods_avg.out','igstrain_loads_avg.out'],
#                ifig=None):
#     if ifig!=None:import matplotlib.pyplot as plt
#     import numpy as np
# #    from ssort import sh as sort
#     fu = open(fns[0],'r').read()
