"""
Diffraction module for NIST X-ray configuration
"""
import os
import numpy as np
import matplotlib.pyplot as plt

isfile = os.path.isfile
def fc(fn):
    """ file check
    """
    if not(isfile(fn)): raise IOError, '%s file not existing'%fn

"""
# format
# 1 npb (phi beta combination number)   [0]


#     step #  ipb and 'detector 1'           [0]
#        wgt                                 [1]
#        phi_bragg                           [2]
#        step # and dtector 2                [3]
#        wgt                                 [4]
#        phi_bragg                           [5]

# - repeat [0]~[5] for ipb+1 times.
"""

def arcscanout(istep=1, ipb=1, iph=1):
    """
    Read arcscan_ph?.out and returns weights, phi^{Bragg},
    and average of phi^{Bragg} for each detector

    wgt1, phi1, wgt2, phi2, avg1, avg2 = out(istep=2, ipb=28)

    Arguments
    =========
    istep = 1
    ipb   = 1
    iph   = 1
    """
    fn  = 'arcscan_ph%i.out'%iph
    fc(fn) # filecheck
    dl  = open(fn, 'r').readlines()
    npb = int(dl[0].split()[0])
    dl  = dl[1:]

    # print 'len(dl)', len(dl)
    dat = dl[6*(istep-1)*npb+ 6*(ipb-1):6*(istep-1)*npb+6*ipb]
    print 'len(dat)', len(dat)
    wgt1 = np.array(map(float, dat[1].split()[1:]))
    phi1 = np.array(map(float, dat[2].split()[1:]))
    ang1 = float(dat[2].split()[0])
    # print 'ang1', ang1

    tp1 = phi1 + ang1
    avg1 = np.average(tp1)

    wgt2 = np.array(map(float, dat[4].split()[1:]))
    phi2 = np.array(map(float, dat[5].split()[1:]))
    ang2 = float(dat[5].split()[0])
    # print 'ang2', ang2
    tp2 = phi2 + ang2
    avg2 = np.average(tp2)

    # print 'avg1, avg2', avg1, avg2

    return wgt1, tp1, wgt2, tp2, avg1, avg2

def binning(phi, wgt, n=20):
    """
    Binning phi-wgt relationship by n number of blocks.

    Arguments
    =========
    phi
    wgt
    n   = 20

    # example:
    x, y = binning(phi1, wgt1, n=10)
    plot(x, y, '+-')
    """
    dwgt = np.zeros((n-1,))
    dphi = (max(phi) - min(phi))/n
    p0s  = np.linspace(min(phi), max(phi), n)
    p0s  = p0s[::-1]

    wbar = np.zeros((len(p0s),))
    for i in range(len(phi)):
        p = phi[i]
        for j in range(len(p0s)):
            if p>=p0s[j]:
                wbar[j] = wbar[j] + wgt[i]
                break
    return p0s, wbar

def plotpeak(steps=[0,1,2], n=20, ipb=1, iph=1):
    """
    Plot intensity vs Bragg angle peak profile

    Arguments
    =========
    steps = [0,1,2]
    n     = 20     : number of bins
    ipb   = 1      : phi and beta combination id
    iph   = 1      : phase


    # example:
    plotpeak(steps=[4, 5, 10], n=50, ipb=10, iph=1)
    """
    # -- binning figures
    plt.figure(1); ax1=plt.gca()
    plt.figure(2); ax2=plt.gca()

    # -- average figures
    plt.figure(3); ax3=plt.gca()
#    plt.figure(4); ax4=plt.gca()

    a1 = []
    a2 = []
    for i in range(len(steps)):
        wgt1, phi1, wgt2, phi2, avg1, avg2 \
            = arcscanout(istep=steps[i], ipb=ipb, iph=iph)
        a1.append(avg1)
        a2.append(avg2)

        # -- detector 1
        p0s, wbar = binning(n=n, phi=phi1, wgt=wgt1)
        p0s = p0s[::-1]
        wbar = wbar[::-1]
        ax1.plot(p0s, wbar, '-o', label='%i'%steps[i])

        # -- detector 2
        p0s, wbar = binning(n=n, phi=phi2, wgt=wgt2)
        p0s = p0s[::-1]
        wbar = wbar[::-1]
        ax2.plot(p0s, wbar, '-o', label='%i'%steps[i])


    ax3.plot(steps, a1, '-o', mfc='None', mec='red')
    ax3.plot(steps, a2, '-x', mec='blue')

    ax1.set_xlabel(r'$\theta^{B}$', dict(fontsize=25))
    ax2.set_xlabel(r'$\theta^{B}$', dict(fontsize=25))
    ax3.set_xlabel(r'steps', dict(fontsize=30))
    ax1.set_ylabel('Intensity')
    ax2.set_ylabel('Intensity')
    ax3.set_ylabel('Intensity')

#    ax1.legend(loc='best', fancybox=True).get_frame().set_alpha(0.5)
#    ax2.legend(loc='best', fancybox=True).get_frame().set_alpha(0.5)
#    ax3.legend(loc='best', fancybox=True).get_frame().set_alpha(0.5)

    plt.show()

    #    Below is for effective modulus and stress factor
"""
# eff_mod_ph?.out's structure:
# 28 npb  (header)         [-1]
    # --- for each ipb
# 'sigma bar'              [0]
# sig1i                    [1]
# sig2i                    [2]
# sig3i                    [3]
# istep   ipb   phi, beta  [4]     [0]*
# -- detector 1            [5]
# eps^{hkl}   ...          [6]
# 'DEC'                    [7]
#  c1j                     [8]
#  c2j                     [9]
#  c3j                     [10]
#  c4j                     [11]
#  c5j                     [12]
#  c6j                     [13]
# -- detector 2            [14]
# eps^{hkl}                [15]
#  c1j                     [16]
#           .
#           .
#           .
#                          [22]    [18]*

#   [0]*~[17]* repeats as many as npb

# one step has 'one block' of (4 + 17*npb) lines
"""

def effmodout(istep=1, ipb=1, iph=1):
    """
    istep, ipb and iph start from 1

    Arguments
    =========
    istep = 1
    ipb = 1
    iph = 1

    # example
    phi, beta, eps1, eps2, cij1, cij2, sbar\
       = effmodout(istep=1, ipb=1, iph=1)
    """
    istep = istep - 1
    ipb = ipb - 1

    filename = 'eff_mod_ph%i.out'%iph
    fc(filename)
    dl = open(filename, 'r').readlines()
    npb = int(dl[0].split()[0])
    #    print 'npb:', npb

    if ipb>=npb: raise IOError,\
            'ipb should be smaller than npb: %i'%npb
    dat = dl[1:]
    #print 'len(dat)',len(dat)
    dat = dat[istep*(4+19*npb):]

    # macroscopic stress reads
    sbar = np.zeros((3,3))
    for i in range(3):
        sbar[i,0],sbar[i,1],sbar[i,2] = map(float, dat[1+i].split())


    dat = dat[4+ipb*19:]
    #    print 'starting index:', 4+ipb*19
    #    print dat[0]
    dum = dat[0].split()
    step = int(dum[0])
    #if istep+1!=step:
    #    raise IOError, 'STEP is not equal'
    ipb = int(dum[1])
    phi = float(dum[2])
    beta = float(dum[3])

    eps1 = float(dat[2].split()[1])
    dec1 = np.zeros((6,6))
    for i in range(6):
        dec1[i,0], dec1[i,1], dec1[i,2], dec1[i,3], \
            dec1[i,4], dec1[i,5] = map(float, dat[4+i].split())

    eps2 = float(dat[11].split()[1])

    dec2 = np.zeros((6,6))
    for i in range(6):
        dec2[i,0], dec2[i,1], dec2[i,2], dec2[i,3],\
            dec2[i,4],dec2[i,5] =  map(float, dat[13+i].split())

    return phi, beta, eps1, eps2, dec1, dec2, sbar


#   Stress factor ##

def eps_sig_fij(i=1, j=1, ipb=1, iph=1):
    """
    Read 'eff_mod_ph?.out' and returns diffraction strain and sigma_{ij}

    Arguments
    =========
    i   = 1
    j   = 1
    ipb = 1
    iph = 1

    # example
    phi, beta, cij1, cij2 = eps_sig_fij(i=1, j=1, ipb=1, iph=1)
    """
    i = i - 1; j = j - 1
    eps1 = []; eps2 = []
    sigij = []
    phi = None
    beta = None
    iquit = False
    istep = 1
    while(not(iquit)):
        try: p, b, e1, e2, c1, c2, s\
                = effmodout(istep=istep, ipb=ipb, iph=iph)
        except:
            iquit=True
            #print 'phi:', phi
            #print 'beta:', beta
            return phi, beta, np.array(eps1), \
                np.array(eps2), np.array(sigij)
        else:
            if phi==None and beta==None: phi = p; beta = b
            else:
                if phi!=p or beta!=b: raise IOError, \
                        'phi and bet not matched.'
            eps1.append(e1)
            eps2.append(e2)
            sigij.append(s[i,j])
        istep = istep + 1

        if istep>100000: raise IOError, 'Abnormally many iter happend.'

def fijfit(eps,sigij):
    """
    Linearly fitting fij

    Arguments
    =========
    eps, sigij
    """
    rst = np.polyfit(sigij,eps,1)
    return rst

def fijout(i=1, j=1, ipb=1, istep=1, iph=1):
    """
    f_ij the stress factor for a certain phi and beta set.
    istep, i and j start from 1

    Arguments
    =========
    i     = 1 (it starts from 1)
    j     = 1 (it starts from 1)
    ipb   = 1
    istep = 1
    iph   = 1
    """
    i = i - 1; j = j - 1
    # 40 npb
    # 1 1 -90 23. 0 [0]
    # fij1          [1]
    # psi 11.2      [2]
    # f1i           [3]
    # f2i           [4]
    # f3i           [5]
    # fij2          [6]
    # psi  34.80    [7]
    # f1i           [8]
    # f2i           [9]
    # f3i           [10]   -- upto for one ipb case (11 lines)
    # 1 2 -90 20. 0 [11]
    #  ...

    if istep<=0: raise IOError, 'istep should start from 1'

    filename = 'sig_fac_ph%i.out'%iph
    fc(filename)
    dl       = open(filename, 'r').readlines()
    npb      = int(dl[0].split()[0])
    dat      = dl[1:]
    dat      = dat[11 * npb * (istep-1) + 11*(ipb-1):\
                       11 * npb * (istep-1) + 11*ipb]
    #print 'len(dat)', len(dat)
    #print 'dat[0]', dat[0]

    try: step = int(dat[0].split()[0])
    except IndexError: return False

    ip   = int(dat[0].split()[1])
    phi  = float(dat[0].split()[2])
    beta = float(dat[0].split()[3])

    #if step!=istep:
    #    print 'Given istep is %i but found step is %i'%(istep, step)
    #    raise IOError, 'Step is not matched'
    if ip!=ipb: raise IOError, 'ipb is not matched'

    psi1 = float(dat[2].split()[1])
    psi2 = float(dat[7].split()[1])
    fij1 = np.zeros((3,3)); fij2 = np.zeros((3,3))

    for k in range(3):
        fij1[k] = np.array(map(float, dat[3+k].split()))
        fij2[k] = np.array(map(float, dat[8+k].split()))

    return phi, beta, psi1, psi2, fij1[i,j], fij2[i,j]


def fijstep(i=1, j=1, istep=1, iph=1, difn=None):
    """
    Arguments
    =========
    i = 1
    j = 1
    istep = 1
    iph = 1
    difn = None
    """
    import numpy as np
    fc(difn)
    difl = open(difn,'r').readlines()
    nphi = int(difl[5].split()[0])
    phis = map(float,difl[6].split()[:nphi])
    nbeta = int(difl[7].split()[0])
    betas = map(float,difl[8].split()[:nbeta])

    npsi = nbeta * 2 # two detectors

    psis = []
    for k in range(nbeta):
        psis.append(betas[k]+11.8)
        psis.append(betas[k]-11.8)

    # print 'phis:',  phis
    # print 'betas:', betas
    # print 'psis:',  psis

    f = np.zeros((nphi, npsi))

    ipb = 1
    for k in range(nphi):
        for l in range(nbeta):
            phi  = phis[k]
            beta = betas[l]

            PHI, BETA, PSI1, PSI2, F1, F2 = \
                fijout(i, j, ipb, istep, iph)

            if abs(phi-PHI)>1.0: raise IOError
            if abs(beta-BETA)>1.0: raise IOError

            ipb = ipb + 1

            f[k,l*2]   = F1
            f[k,l*2+1] = F2

    return np.array(psis), np.array(f) # phi, psi

def fij(i=1, j=1, iph=1, difn=None):
    """
    Arguments
    =========
    i   = 1
    j   = 1
    iph = 1

    Fij (step, phi, psi)
    psi (step, phi)
    """
    fc(difn)
    from ssort import shellSort as ssort
    istep = 1
    Psi = []
    F   = []
    while True:
        print 'istep:', istep
        try: psi, f = fijstep(i, j, istep, iph, difn=difn)
        except: break
        else:
            ## -- sorting block starts --
            #f(phi, psi)
            ft = f.T
            x = psi
            sortedpsi, ind = ssort(x)
            temp = []
            for k in range(len(sortedpsi)):
                temp.append(ft[ind[k]-1])
            ft = np.array(temp)
            f = ft.T
            psi = sortedpsi
            ## -- sorting block ends --

            F.append(f)
            Psi.append(psi)

            istep = istep + 1

    # F(step, phi, psi)
    # Psi(step, phi)

    if len(F)==0 or len(Psi)==0:
        print "Something fishy happened: len(F)==0 or len(Psi)==0"
        print len(F), len(Psi)

    return np.array(F), np.array(Psi)

def internalepsplot(iph=1, phi=1, difn=None):
    """
    Arguments
    =========
    iph  = 1
    phi  = 1
    difn = None
    """
    fc(difn)
    difl = open(difn, 'r').readlines()
    nphi = int(difl[5].split()[0])
    phis = map(float,difl[6].split()[:nphi])
    neps = int(difl[11].split()[0])
    acueps = map(float, difl[12].split()[:neps])
    acueps = acueps[::-1]#; acueps.append(0)
    acueps = np.array(acueps[::-1])

    nstep = len(acueps)

    for istep in range(nstep):
        for p in phis:
            dat = intepsphiout(istep=istep, iph=iph, phi=p,
                               isort=True)
            psi, eps, sig, ms11, ms22, ms33, ms23, ms13, ms12= \
                dat[:9]

def eps_sin2psi(difile=None, istep=0, phi=0, iph=1, ax=None,iopt=0):
    """
    Arguments
    =========
    difile : diffraction input file name
    istep  : i-th step      (starts from 0)
    phi    : phi angle (0, 45, 90, 135, 180 ... )
    iph    : i-th phase     (starts from 1)
    ax     : None
    iopt   : 0->int_eps_ph%i.out;  1-> int_els_ph%i.out
    """
    from sff_converter import condition
    difl, nphi, phis, nbeta, neps, eps = condition(difile)
    dat = intepsphiout(istep=istep,phi=phi,iph=iph,isort=True,iopt=iopt)
    psi, eps, sig, ms11, ms22, ms33, ms23, ms13, ms12 = dat[0:9]
    #me11, me22, me33, me23, me13, me12 = dat[9:15]
    sin2psi = np.sin(psi*np.pi/180.)**2
    print 'psi', psi
    print 'sin2psi', sin2psi
    print 'eps', eps
    if ax==None: ax = plt.gca()

#    if iopt==1: label=r'$\sigma_{11}=$%7.3f'%ms11
#    elif iopt==1: label=r'$\sigma_{11}=$%7.3f'%ms11
    label=None

    ax.plot(sin2psi,eps*10**3,'o-',label=label)
    ax.set_xlabel(r'$\sin^2{\psi}$', dict(fontsize=28))
    ax.set_ylabel(r'$\varepsilon^{hkl}\times10^{3}$', dict(fontsize=28))
    if label!=None:
        ax.legend(loc='lower right',fancybox=True).get_frame().set_alpha(0.5)
    plt.tight_layout()
    plt.show()

def fijplot(i=1,j=1,iph=1,ifig=0, ired=False, difn=None):
    """
    Arguments
    =========
    i    = 1    i and j reprsents first and second
    j    = 1    components of the stress factor F_ij
    iph  = 1    Phase
    ifig = 0    figure numbering
    ired = False  ! reduction flag
                 (refer to source code for more details)
    difn = None
    """
    fc(difn)
    f, psi = fij(i,j,iph,difn=difn)
    print f.shape
    print psi.shape

    difl = open(difn, 'r').readlines()
    nphi = int(difl[5].split()[0])
    phis = map(float,difl[6].split()[:nphi])
    neps = int(difl[11].split()[0])
    eps = map(float, difl[12].split()[:neps])
    eps = eps[::-1]; eps.append(0)
    eps = np.array(eps[::-1])

    # F   (step, phi, psi)
    # Psi (step, psi)
    f   =   f.swapaxes(0,1)
    # F   (phi, step, psi)

    import matplotlib.pyplot as plt
    markers = ['x-','+-','d-','.-','*-','^-']

    fact=5.
    ncol=3
    nrow = 1 + int(nphi/(ncol+0.0001))

    figsize=(ncol*fact, nrow*fact)
    fig = plt.figure(ifig, figsize=figsize)

    axes = []
    for k in range(len(f)): #  phi
        fig.add_subplot(nrow, ncol, k + 1)
        ax = plt.gca()
        axes.append(ax)
        for l in range(len(f[i])): # step
            label=r'$\varepsilon$ %6.3f'%eps[l]
            # Reducing of psi range -- as requested by Thomas - 10/1/2012
            # - ~ +
            if ired:
                for m in range(len(psi[l])):
                    if psi[l][m]>-45:
                        m0 = m
                        break
                for m in range(len(psi[l])):
                    if psi[l][m]>45:
                        m1 = m
                        break
                if k==0:
                    ax.plot(psi[l][m0:m1],
                            f[k][l][m0:m1]*(10**6),
                            markers[l],label=label)
                else:
                    ax.plot(psi[l][m0:m1],
                            f[k][l][m0:m1]*(10**6),
                            markers[l])
            else:
                if k==0: ax.plot(psi[l], f[k][l]*(10**6),\
                                     markers[l],label=label)
                else: ax.plot(psi[l], f[k][l]*(10**6), \
                                  markers[l])
        __fijplotdeco__(ax=ax, i=i,j=j, phi=phis[k])
    plt.tight_layout()

def __fijplotdeco__(ax,i,j,phi):
    ax.set_xlabel(r'$\psi$', dict(fontsize=28))
    ax.set_ylabel(r'$F_{%i%i}$ [TPA$^{-1}$]'%(i,j), dict(fontsize=28))
    ax.set_title(r'$\phi$: %5.3f'%phi)
    try: ax.legend(loc='best').get_frame().set_alpha(0.5)
    except: pass

def __inteps_plotdeco__(ax, phi):
    ax.set_xlabel(r'$\psi$', dict(fontsize=28))
    ax.set_ylabel(r'$\varepsilon^{hkl}$ [$\mu$]', dict(fontsize=28))
    ax.set_title(r'$\phi$: %5.3f'%phi)
    ax.legend(loc='best').get_frame().set_alpha(0.5)

def __sigeps_plotdeco__(ax, phi):
    ax.set_xlabel(r'$\sigma$', dict(fontsize=28))
    ax.set_ylabel(r'$\varepsilon^{hkl}$', dict(fontsize=28))
    ax.set_title(r'$\phi$: %5.3f'%phi)

def cijout(i=1,j=1,ipb=1,iph=1):
    """
    c_ij, the elastic constant for a certain phi and beta set

    Arguments
    =========
    i   = 1
    j   = 1
    ipb = 1
    iph = 1

    # examples
    phi, beta, cij1, cij2 = cijout(i=1, j=1, ipb=1, iph=1)
    """
    iquit = False
    istep = 1
    cij1  = []
    cij2  = []
    sbar  = []
    phi = None; beta = None
    while(not(iquit)):
        try: p, b, f1, f2, c1, c2, s \
                = effmodout(istep=istep, ipb=ipb, iph=iph)
        except:
            iquit = True
            return phi, beta, np.array(cij1),\
                np.array(cij2), np.array(sbar)
        else:
            if phi==None and beta==None: phi=p; beta = b
            if phi!=p or beta!=b: raise IOError, \
                    'Unexpected phi and beta'
            cij1.append(c1[i,j]); cij2.append(c2[i,j])
            sbar.append(s)
        istep = istep + 1
    return phi, beta, np.array(cij1), np.array(cij2), np.array(sbar)

def eps_sig_fij_sort(i=1, j=1, ipb=[1,2,3,4], opt=1, iph=1, difn=None):
    """
    Sorting the result from fij for ipb based on psi
    regardless of the detector it comes from.

    Arguments
    =========
    i   = 1
    j   = 1
    ipb = [1,2,3,4]
    opt = 1 (1: sort based on sin2psi,  2: sort based on psi)
    iph = 1

    # examples
    psi, sin2psi, eps, sigma = \
       eps_sig_fij_sort(i=1,j=1,ipb=[1,2,3], opt=1, iph=1)
    """
    # import the shell sort f2py-wrapped binary module
    fc(difn)
    from sort import shellsort as ssort

    phi   = []; beta  = []
    psi   = []; sigma = []
    eps   = []

    for k in range(len(ipb)):
        p, b, e1, e2, s = eps_sig_fij(i=i,j=j,ipb=ipb[k], iph=iph, difn=difn)
        phi.append(p); beta.append(beta)
        psi1 = b - 11.8; psi2 = b + 11.8
        psi.append(psi1); psi.append(psi2)
        sigma.append(s); sigma.append(s)
        eps.append(e1); eps.append(e2)

    # sort
    sin2psi = sin(np.array(psi)*pi/180.)**2
    if opt==1:
        sorted, ind = ssort(sin2psi, len(sin2psi))
    elif opt==2:
        sorted, ind = ssort(psi, len(psi))
    else: raise IOError, 'Wrong option'

    temp1 = []
    temp2 = []
    temp3 = []
    for k in range(len(sorted)):
        temp1.append(sigma[ind[k]-1])
        temp2.append(eps[ind[k]-1])
        temp3.append(psi[ind[k]-1])
    sigma = temp1
    eps = temp2
    psi = temp3

    return psi, sin2psi, eps, sigma

#   internal strain
def intepsout(istep=1, ipb=1, iph=1):
    filename='int_eps_ph%i.out'%iph
    fc(filename)
    dl = open(filename,'r').readlines()
    npb = int(dl[1].split()[0])
    if ipb>npb: raise IOError, 'ipb is exceeding npb'
    dat = dl[2:]
    dat = dat[(istep-1)*npb + (ipb-1)]
    #print len(dat)
    dat = map(float, dat.split())
    step, phi, beta, psi1, psi2, eps1, eps2,\
        sig1, sig2, n1, n2, v1, v2,\
        m11, m22, m33, m23, m13, m12,\
        e11, e22, e33, e23, e13, e12 = dat
    #print 'step:', step
    #print 'phi:', phi
    #print 'beta:', beta
    return phi, beta, eps1, eps2, sig1, sig2

def intepsphiout(istep=0, phi=0, iphi=0, iph=1, isort=False, iopt=0):
    """
    Arguments
    =========
    istep = 0 (starts from 0)
    phi   = -90 or iphi as index for iopt 4 and 5.
    iph   = 1 (phase starts from 1)
    isort = False
    iopt  =
            0: ('int_eps_ph%i.out') - at unloads,
            1: ('int_els_ph%i.out') - at loads
            2: ('igstrain_load_ph%i.out') - at
            3: ('igstrain_unload_ph%i.out') - at
            4: ('igstrain_unloads_avg.out') at unloads
            5: ('igstrain_loads_avg.out')   at loads
    """
    from ssort import shellSort as ssort
    from ssort import sh
    from ssort import ind_swap

    eps, sig, psi = [], [], []
    if iopt==0:   filename='int_eps_ph%i.out'%iph
    elif iopt==1: filename='int_els_ph%i.out'%iph
    elif iopt==2: filename='igstrain_load_ph%i.out'%iph
    elif iopt==3: filename='igstrain_unload_ph%i.out'%iph
    elif iopt==4: filename='igstrain_unloads_avg.out'
    elif iopt==5: filename='igstrain_loads_avg.out'
    fc(filename)

    if iopt<4:
        dl = open(filename, 'r').readlines()
        npb = int(dl[1].split()[0])
        dat = dl[2:]
        dat = dat[(istep)*npb: (istep)*npb + npb]
        if len(dat)==0: raise IOError, 'empty dat is returned ..'
        for i in range(len(dat)):
            temp = map(float, dat[i].split())
            step, p, beta, psi1, psi2, eps1, eps2,\
                sig1, sig2, n1, n2, v1, v2 = temp[:13]

            ms11, ms22, ms33, ms23, ms13, ms12 = temp[13:19]
            me11, me22, me33, me23, me13, me12 = temp[19:25]

            if p==phi:
                eps.append(eps1); eps.append(eps2)
                sig.append(sig1); sig.append(sig2)
                psi.append(psi1) # detector 1
                psi.append(psi2) # detector 2
        if isort:
            temp1, temp2 = [], []
            sortedarray, ind = ssort(psi)#, len(psi))
            for k in range(len(sortedarray)):
                temp1.append(eps[ind[k]-1])
                temp2.append(sig[ind[k]-1])

            eps = temp1
            sigma = temp2
            psi = sortedarray

        return np.array(psi), np.array(eps), np.array(sig),\
            ms11, ms22, ms33, ms23, ms13, ms12

    if iopt==4 or iopt==5:
        from sff_converter import condition
        from pepshkl import reader3 as reader
        difl, nphi, phis, nbeta, neps, eps = condition()
        tdat, psi = reader(filename,isort) #[ist,iphi,ipsi]
        sin2psi = np.sin(psi*np.pi/180.)**2
        npsi = len(psi)
        eps = tdat[istep, iphi]
        return psi, eps

    raise IOError, 'Unexpected ioption give in intepsphiout'

def internal_eps_plot(
    iph=1, ifig=1, ired=False,
    ppsi=[1,2], difn=None):
    """
    Arguments
    =========
    iph  = 1
    ifig = 1
    ired = False
    ppsi = [1,2]
    difn = None
    """
    fc(difn)
    difl = open(difn, 'r').readlines()
    nphi = int(difl[5].split()[0])
    phis = map(float,difl[6].split()[:nphi])
    NBETA = int(difl[7].split()[0])
    NPSI = NBETA * 2
    BETAs = map(float, difl[8].split()[:NBETA])
    neps = int(difl[11].split()[0])
    acueps = map(float, difl[12].split()[:neps])
    acueps = acueps[::-1]#; acueps.append(0)
    acueps = np.array(acueps[::-1])

#    psi_eps_dat = np.zeros((nphi,NPSI)) # phi and psi
#    psi_sig_dat = np.zeros((nphi,NPSI))

    print '-- Beta\n',BETAs

    # alternating colors
    colors=['r','g','b','c','m','y','k','r','g','b','c','m','y','k',
            'r','g','b','c','m','y','k','r','g','b','c','m','y','k',
            'r','g','b','c','m','y','k','r','g','b','c','m','y','k',
            'r','g','b','c','m','y','k','r','g','b','c','m','y','k']

    nstep = len(acueps)
    axes = []
    axes01 = []
    import matplotlib.pyplot as plt
    markers = ['x-','+-','d-','.-','*-','^-','x-','+-','d-',
               '.-','*-','^-','x-','+-','d-','.-','*-','^-']

    fact = 5.
    ncol = 3
    nrow = 1 + int(nphi/(ncol+0.0001))

    figsize = (ncol*fact, nrow*fact)
    if ifig!=None:
        fig = plt.figure(ifig,figsize=figsize)
        fig01 = plt.figure(ifig+1,figsize=figsize)

    psi_eps_dat = []
    psi_sig_dat = []

    for i in range(nphi): # phi
        if ifig!=None:
            fig.add_subplot(nrow,ncol,i+1)
            ax = fig.axes[-1]
            fig01.add_subplot(nrow,ncol,i+1)
            ax01 = fig01.axes[-1]
            axes.append(ax)
            axes01.append(ax01)
            ax.set_title(r'$\phi$: %5.3f'%phis[i])
            ax01.set_title(r'$\phi$: %5.3f'%phis[i])

        temp1 = np.zeros((NPSI,neps+1)) # internal strain
        temp2 = np.zeros((6,NPSI, neps+1)) # macro stress

        for j in range(neps):
            psi, eps, sig, ms11, ms22, ms33, ms23, ms13, ms12, \
                me11, me22, me33, me23, me13, me12 \
                = intepsphiout(istep=j, phi=phis[i],
                               iph=iph, isort=True)

            if i==0 and j==0: npsi=len(psi)

            temp1[:,j]   = eps[:]
            temp2[0,:,j] = ms11
            temp2[1,:,j] = ms22
            temp2[2,:,j] = ms33
            temp2[3,:,j] = ms23
            temp2[4,:,j] = ms13
            temp2[5,:,j] = ms12

            if ifig!=None:
                if ired:
                    for m in range(len(psi)):
                        if psi[m]>-45:
                            m0 = m
                            break
                    for m in range(len(psi)):
                        if psi[m]>45:
                            m1 = m
                            break
                    ax.plot(
                        psi[m0:m1],np.array(eps)[m0:m1]*10**6,
                        markers[j],
                        label=r'$\varepsilon$ %6.3f'%acueps[j])
                else:
                    ax.plot(
                        psi,np.array(eps)*10**6,markers[j],
                        label=r'$\varepsilon$ %6.3f'%acueps[j])

                    iprob = 0
                    for k in range(len(psi)):
                        if k in ppsi:
                            iprob = iprob + 1
                            ax01.plot(ms11, np.array(eps[k])*10**6,'.',
                                      color=colors[iprob])

                            ax01.plot(ms22, np.array(eps[k])*10**6,'x',
                                      color=colors[iprob])

            if ifig!=None:
                __inteps_plotdeco__(ax, phi=phis[i])
                __sigeps_plotdeco__(ax01,phi=phis[i])

        psi_eps_dat.append(temp1) # phi, psi, neps
        psi_sig_dat.append(temp2) # phi, 6, psi, neps

    psi_eps_dat = np.array(psi_eps_dat)
    psi_sig_dat = np.array(psi_sig_dat).swapaxes(0,1)
    # psi_sig_dat # 6, phi, psi, neps

    if ifig!=None: plt.tight_layout()
    return psi_eps_dat, psi_sig_dat, psi, phis

def lin_sff(ax=None,ipsi=0,iph=0,isig=0,ifig=30,difn=None):
    """
    Arguments
    =========
    ax   = None (optional)
    ipsi = 0
    iph  = 0
    isig = 0 (0:sig11, 1:sig22, 2:sig33, 3:sig23, 4:sig13, 5:sig12)
    ifig = 30
    difn = None
    """
    plt.figure(ifig)
    ax=plt.gca()
    # eps [phi, psi, neps]
    # sig [6, phi, psi, neps]
    eps, sig, psi, phis = internal_eps_plot(iph=1, ifig=None, difn=difn)
    ax.plot(sig[isig,iph,ipsi], eps[iph,ipsi], '-+',
            label=r'$\phi: %5.1f$  $\psi: %5.1f$'%(
            phis[iph], psi[ipsi]
            ))
    ax.legend(loc='best',fancybox=True).get_frame().set_alpha(0.5)

def sff_calc(str_fn = 'STR_STR.OUT', hist='ex02_NIST_diff/BB',
             difn=None):
    """
    Arguments
    =========
    str_fn = 'STR_STR.OUT'
    hist   = 'ex02_NIST_diff/BB'
    difn   = None
    """
    fc(str_fn)
    fc(hist)
    fc(difn)
    difl = open(difn, 'r').readlines()
    nphi = int(difl[5].split()[0])
    phis = map(float,difl[6].split()[:nphi])
    neps = int(difl[11].split()[0])
    acueps = map(float, difl[12].split()[:neps])
    acueps = acueps[::-1]#; acueps.append(0)
    acueps = np.array(acueps[::-1])
