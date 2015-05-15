"""
Structure factor (*.sff) plotting scripts.
* sff files are processed file from EVPSC output
  for Thomas Gnaupel's PF software


# main: f1122
 - f1122(fn='debug.sff',ifig=1,iphi=0,i=1,j=1,ilab=True):
 - plots
     sin^2(psi) versus F
           psi  versus F
"""

#print __doc__

import numpy as np
from matplotlib.ticker import MaxNLocator
def readpsi(fn='debug.sff',lb='\r'):
    """
    Argument
    ========
    fn = 'debug.sff'
    lb = '\r'  line breaker?

    return based on psi
    """
    #datl = open(fn,'r').readlines()
    datl = open(fn,'r').read().split(lb)
    header = datl[:5]
    npsi = int(header[3].split()[-1])
    dat = np.loadtxt(fn,skiprows=6).T

    rst = dat[1][:npsi]
    print len(rst)
    return rst

def readeps(fn='debug.sff',lb='\r'):
    """
    Argument
    ========
    fn = 'debug.sff'

    return based on eps
    """
    # datl = open(fn,'r').readlines()
    datl = open(fn,'r').read().split(lb)
    header = datl[:5]
    neps = int(header[1].split()[-1])
    rst = map(float, header[4].split()[1:])
    return rst

def readphi(fn='debug.sff',lb='\r'):
    datl = open(fn,'r').read().split(lb)
    header = datl[:5]
    print header[2]
    nphi = int(header[2].split()[-1])
    print 'nphi:', nphi
    phis = np.loadtxt(fn,skiprows=6,dtype='str').T[0]
    phis.sort()
    phis = map(float,phis)
    phis = np.unique(np.array(phis,dtype='float'))
    return phis[:nphi:]

def reader(fn='debug.sff',lb='\r'):
    """
    Arguments
    =========
    fn = 'debug.sff'
    """
    dat=np.loadtxt(fn,skiprows=6).T
    phis=np.unique(dat[0])
    datl = open(fn,'r').readlines()
    if len(datl)==1:
        datl = datl[0].split(lb)
    # header
    header = datl[:5]
    nstr = int(header[1].split()[-1])
    nphi = int(header[2].split()[-1])
    npsi = int(header[3].split()[-1])
    exx  = map(float,header[4].split()[1:])

    dat = np.loadtxt(fn,skiprows=6).T
    fij = np.zeros((nstr,nphi,npsi,6)) # stress factor
    e0  = np.zeros((nstr,nphi,npsi))   # strain

    # (x,y) = ? (x = phi,psi,f11,f22,f33,f23,f12,eps0), (y=phi,psi)
    # each eps_xx takes a block that has (8,nphi*npsi) elements.

    # ix,iy = 0,0
    for istr in range(nstr):
        ix = istr * 9
        for iphi in range(nphi):
            iy0 = iphi * npsi
            if istr==0 and iphi==0:
                psis = []
            for ipsi in range(npsi):
                iy = iy0 + ipsi
                e0[istr,iphi,ipsi] = dat[ix+8,iy]
                for ij in range(6):
                    fij[istr,iphi,ipsi,ij] = dat[ix+2+ij,iy]

                if istr==0 and iphi==0:
                    psis.append(dat[1,iy])

    psis=np.array(psis)
    return fij,e0,phis,psis,exx

def f1122(fn='debug.sff',ifig=1,iphi=0,i=1,j=1,ilab=True,
          title='EVPSC',ls='-',ieps=None):

    """
    Arguments
    =========
    fn   = 'debug.sff'
    ifig = 1
    iphi = 1
    i    = 1 (starts from 1)
    j    = 1 (starts from 1)
    ilab = True
    ls   = '-'
    """
    import matplotlib.pyplot as plt


    try:
        readpsi(fn,lb='\n')
        readeps(fn,lb='\n')
        readphi(fn,lb='\n')
    except:
        lb='\r'
    else:
        lb='\n'

    fig01  = plt.figure(ifig,figsize=(4.,3))
    #fig02  = plt.figure(ifig+1,figsize=(6.,5.5))
    fig03  = plt.figure(ifig+2,figsize=(4.,3))
    fig01.clf();fig03.clf()
    ax01   = fig01.add_subplot(111)
    #ax02   = fig02.add_subplot(111)
    ax03   = fig03.add_subplot(111)
    fij,e0,dum1,dum2,dum3 = reader(fn) # nstr, nphi, npsi, 6
    axes   = [ax01,ax03]
    nstr, nphi, npsi, dum = fij.shape
    psi    = readpsi(fn,lb=lb)
    eps    = readeps(fn,lb=lb)#[::2]
    phi    = readphi(fn,lb=lb)

    color = ['r','g','b','k','m','r','g',
             'b','k','m','r','g','b','k',
             'm','r','g','b','k','m']
    markers = ['-o','-x','-+','-^','-d','-*']

    ijs = [[1,1],[2,2]]

    istr0 = -1
    for istr in range(len(eps)):
        if istr==0 and title=='EXP': pass
        elif ieps==None or istr in ieps:
            istr0 = istr0 + 1
            for i in range(2):
                f  = fij[istr][iphi]
                i0,j0 = ijs[i]

                iv = ij2voigt(i0,j0)
                y  = f.T[iv-1]

                if i==0: ilab=True
                if i==1: ilab=False

                if istr==0 and i==0:
                    ax01.text(0.2,0.1,r'Bold lines: $F_{11}$',
                              transform=ax01.transAxes)
                    ax01.text(0.2,0.2,r'Broken lines : $F_{22}$',
                              transform=ax01.transAxes)

                if ilab:
                    label=r'$\bar{E} = %4.2f$'%(eps[istr])
                    ax01.plot(np.sin(psi*np.pi/180.)**2,y,#markers[istr0],
                              label=label,color=color[istr0],alpha=0.7)
                    # ax02.plot(psi,y,markers[istr0],label=label,color='k',
                    #           alpha=0.7)
                    y = e0[istr0][iphi] * (1e6)
                    ax03.plot(np.sin(psi*np.pi/180.)**2,y,markers[istr0],
                              ms=8,label=label,color=color[istr0], alpha=0.7)
                else:
                    ax01.plot(np.sin(psi*np.pi/180.)**2,y,ls='--',#markers[istr0],
                              lw=1.5,color=color[istr0])#,alpha=0.9)
                    # ax02.plot(psi,y,'--',color='k',
                    #           alpha=0.7)
                    y = e0[istr][iphi] * (1e6)
                    ax03.plot(np.sin(psi*np.pi/180.)**2,y,markers[istr0],
                              ms=8,color=color[istr0], alpha=0.5)


    for ax in axes:
        ax.tick_params(axis='both', which='major',
                       labelsize=8)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=3))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=3))
        ax.minorticks_on()

    ft = 15
    ax01.set_ylabel(r'$\mathrm{F}_{11}$, $ \mathrm{F}_{22}$ $[TPa^{-1}]$',dict(fontsize=ft))
    ax01.set_xlabel(r'$\sin^2{\psi}$',dict(fontsize=ft))
    ax03.set_ylabel(r'$\varepsilon^\mathrm{IG}$ [$\mu\varepsilon$]',dict(fontsize=ft))
    ax03.set_xlabel(r'$\sin^2{\psi}$',dict(fontsize=ft))

    for i in range(len(axes)):
        ax = axes[i]
        ax.grid('on')
        ax.legend(loc='best',fancybox=True).\
            get_frame().set_alpha(0.5)
        if i==0: ax.set_title(r'%s at $\phi=%3.0f^\circ$'%(
                title,phi[iphi]),dict(fontsize=ft*1.2))
        if i==0: ax.set_ylim(-2.0,2.0)

    fig01.tight_layout()#pad=4)
    #fig02.tight_layout()#pad=4)
    fig03.tight_layout()#pad=4)
    fig01.savefig('F11_F22_%s.pdf'%(title))
    fig03.savefig('IGstrain_%s.pdf'%(title))
    # plt.close(fig02)
    # plt.close(fig03)
    plt.show()

def voigt2ij(i):
    """
    Convert Voigt into ij
    i starts from 1
    """
    vo = np.zeros((6,2))
    # 1 -> 1,1
    vo[0,0] = 1
    vo[0,1] = 1
    # 2 -> 2,2
    vo[1,0] = 2
    vo[1,1] = 2
    # 3 -> 3,3
    vo[2,0] = 3
    vo[2,1] = 3
    # 4 -> 2,3
    vo[3,0] = 2
    vo[3,1] = 3
    # 5 -> 1,3
    vo[4,0] = 1
    vo[4,1] = 3
    # 6 -> 1,2
    vo[5,0] = 1
    vo[5,1] = 2
    return vo[i-1]

def ij2voigt(i=1,j=1):
    """
    Convert ij into Voigt notation
    i and j start from 1
    """
    ivo = np.zeros((3,3))
    #  1,1 -> 1
    ivo[0,0] = 1
    #  2,2 -> 2
    ivo[1,1] = 2
    #  3,3 -> 3
    ivo[2,2] = 3
    #  2,3 -> 4
    ivo[1,2] = 4
    ivo[2,1] = 4
    #  1,3 -> 5
    ivo[0,2] = 5
    ivo[2,0] = 5
    #  1,2 -> 6
    ivo[0,1] = 6
    ivo[1,0] = 6
    return ivo[i-1,j-1]

def f1122_3d(fn='debug.sff', ifig=1, iphi=0, i=1, j=1, ilab=True,
             ls='x-'):
    """
    3d plotter.

    x: psi
    y: eps_xx
    z: fij(psi)
    """
    from mpl_toolkits.mplot3d import Axes3D as a3
    import matplotlib.pyplot as plt



    iv    = ij2voigt(i,j)
    # fig01 = plt.figure(ifig)
    # fig02 = plt.figure(ifig+1)
    fig03 = plt.figure(ifig+2)
    # ax01  = fig01.add_subplot(111)
    # ax02  = fig02.add_subplot(111)
    ax03  = fig03.add_subplot(111, projection='3d')
    fij,dum1,dum2,dum3,dum4   = reader(fn) # nstr, nphi, npsi, 6
    nstr, nphi, npsi, dum = fij.shape
    psi   = readpsi(fn)
    eps   = readeps(fn)#[::2]

    color = ['r','g','b','k','m',
             'r','g','b','k','m',
             'r','g','b','k','m']



    for istr in range(len(eps)):
        f = fij[istr][iphi]
        print len(f)
        y = f.T[iv-1]
        print 'iv:', iv
        print 'psi:',psi
        print 'f:', y

        ys = np.zeros(psi.shape)
        ys[:] = eps[istr]

#        if ilab:
#            label='%7.3f'%eps[istr]
            # ax01.plot(psi,y,ls,label=label)
            # ax01.plot(np.sin(psi*np.pi/180.)**2,y,ls,
            #           label=label,color=color[istr],alpha=0.7)
            # ax02.plot(psi,y,ls,label=label,color=color[istr],
            #           alpha=0.7)
        ax03.plot(xs=np.sin(psi*np.pi/180.)**2,ys=ys,zs=y,
                  color='k',marker='+')
        # else:
        #     ax01.plot(np.sin(psi*np.pi/180.)**2,y,ls,
        #               color=color[istr], alpha=0.7)
        #     ax02.plot(psi,y,ls,color=color[istr], alpha=0.7)

    # if ilab:
    #     ax01.legend(loc='best',fancybox=True).\
    #         get_frame().set_alpha(0.5)
    #     ax02.legend(loc='best',fancybox=True).\
    #         get_frame().set_alpha(0.5)


    # ax01.set_ylabel(r'F$_{%i%i}$'%(i,j),dict(fontsize=28))
    # ax01.set_xlabel(r'$\sin^2{\psi}$',dict(fontsize=28))
    # plt.figure(ifig)
    plt.tight_layout(pad=4)
    # ax02.set_ylabel(r'F$_{%i%i}$'%(i,j),dict(fontsize=28))
    # ax02.set_xlabel(r'$\psi$',dict(fontsize=28))
    ax03.set_xlabel(r'$\sin^2{\psi}$',dict(fontsize=28))
    ax03.set_ylabel(r'$\varepsilon_{xx}$',dict(fontsize=28))
    ax03.set_zlabel(r'F$_{%i%i}$'%(i,j),dict(fontsize=28))

    # plt.figure(ifig+1)
    # plt.tight_layout(pad=4)

def fij_epsxx(fn='debug.sff', iphi=0, ipsi=3, i=1, j=1, ifig=1):
    import matplotlib.pyplot as plt

    iv    = ij2voigt(i,j)

    fig03 = plt.figure(ifig+2)
    ax03  = fig03.add_subplot(111)
    fij,dum1,dum2,dum3,dum4   = reader(fn) # nstr, nphi, npsi, 6

    nstr, nphi, npsi, dum = fij.shape
    psi   = readpsi(fn)
    eps   = readeps(fn)#[::2]
    phis  = readphi(fn)

    color = ['r','g','b','k','m',
             'r','g','b','k','m',
             'r','g','b','k','m']

    F = []


    for istr in range(len(eps)):
        f = fij[istr][iphi]
        y = f.T[iv-1]
        F.append(y[ipsi])

    print 'eps:', eps
    print 'F:', F

    ax03.plot(eps, F,'-x',label=r'$psi=%7.3f$'%psi[ipsi])
    #ax03.legend(loc='best',fancybox=True).get_frame().set_alpha(0.5)

    # plt.tight_layout(pad=4)
    # ax03.set_xlabel(r'$\varepsilon_{xx}$', dict(fontsize=28))
    # ax03.set_ylabel(r'F$_{%i%i}^{\phi=%5.1f}$'%(i,j,phis[iphi]),
    #                 dict(fontsize=28))

    # ax03.plot(xs=np.sin(psi*np.pi/180.)**2,ys=ys,zs=y,
    #           color='k',marker='+')

    # plt.tight_layout(pad=4)
    # ax03.set_xlabel(r'$\sin^2{\psi}$',dict(fontsize=28))
    # ax03.set_ylabel(r'$\varepsilon_{xx}$',dict(fontsize=28))
    # ax03.set_zlabel(r'F$_{%i%i}$'%(i,j),dict(fontsize=28))

def fij_epsxxs(fn='YJ_Bsteel_BB.sff',iphi=0,i=1,j=1,psi=3):
    import numpy as np
    #psi = np.arange(10)
    psi = readpsi(fn)
    print 'psi:', psi
    for p in range(len(psi)):
        fij_epsxx(fn=fn,iphi=iphi,i=i,j=j,ipsi=p)
