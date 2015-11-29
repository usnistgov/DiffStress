#from rearrange_sort import main as ss
#from ssort import shellSort as sort
## from ssort import sh as sort
from MP.ssort import sh as sort
import matplotlib.pyplot as plt
import numpy as np
def main(fn=None,ref='Bsteel_BB_opt0.txt'):
    """   """
    markers = ['-o','--x','-+','-^','-d','-*']
    datl = open(ref,'r').readlines()
    datl = datl[2:]
    ind = [0,2,11,16,18]
    fn = []; strain = []

    fig = plt.figure(3)
#    fig = plt.figure(3,figsize=(12,5))
    ax1 = fig.add_subplot(111); ax1.grid('on')
#    ax1 = fig.add_subplot(121); ax1.grid('on')
#    ax2 = fig.add_subplot(122); ax2.grid('on')

    for i in xrange(len(ind)):
        f = datl[ind[i]].split()[0]
        epsx = float(datl[ind[i]].split()[1])
        epsy = float(datl[ind[i]].split()[2])
        eps = epsx+epsy
        strain.append(eps)
        f = f+'Data2Phi0.txt'
        fn.append(f)

    for i in xrange(len(fn)):
        sin2psi, ehkl,dspacing,psi = read(fn[i],iopt=1) # d: strain * 10^3
        x = np.array(sin2psi).copy()
        x,ehkl = sort(x,ehkl)
        #x,ind = sort(x,ehkl)
        #ehkl = ss(ind, ehkl)
        # dspacing = ss(ind, dspacing)
        ax1.plot(x, ehkl,markers[i],ms=8,color='k',
                 label=r'$\bar{E}^{\mathrm{eff}}=%5.2f$'%strain[i])
        # ax2.plot(x, dspacing,markers[i],ms=8,
        #          label=r'$\bar{E}^{\mathrm{eff}}=%5.2f$'%strain[i])

    ax1.legend(loc='best',fancybox=True).get_frame().set_alpha(0.5)
    #ax2.legend(loc='best',fancybox=True).get_frame().set_alpha(0.5)
    ax1.set_ylabel(r'$\varepsilon(hkl,\phi,\psi) \times 10^3 $',
                   dict(fontsize=28))
    #ax2.set_ylabel(r'$d(hkl,\phi,\psi)$',
    #dict(fontsize=28))
    ax1.set_xlabel(r'$\sin^2{\psi}$',dict(fontsize=28))
    #ax2.set_xlabel(r'$\sin^2{\psi}$',dict(fontsize=28))
    fig.tight_layout()
    fig.savefig('sin2psi_ex.pdf')

def read(fn,iopt=0,isort=False):
    datl = open(fn,'r').readlines()

    i=13
    psi1      = map(float, datl[i].split()[2:])
    sin2psi1  = map(float, datl[i+1].split()[2:])
    dspacing1 = map(float, datl[i+2].split()[1:])
    ehkl1     = map(float, datl[i+4].split()[1:])  # mu strain * 10^3

    i=25
    psi2      = map(float, datl[i].split()[2:])
    sin2psi2  = map(float, datl[i+1].split()[2:])
    dspacing2 = map(float, datl[i+2].split()[1:])
    ehkl2     = map(float, datl[i+4].split()[1:])  # mu strain * 10^3

    psi      = np.array(psi1 + psi2)

    # sin2psi  = sin2psi1+sin2psi2
    sin2psi = np.sign(psi) * np.sin(psi*np.pi/180.)**2
    ehkl     = ehkl1+ehkl2
    ehkl = np.array(ehkl) / 1e3
    dspacing = dspacing1 +dspacing2

    if iopt==1:
        psi=psi1
        sin2psi = sin2psi1
        ehkl = ehkl1
        dspacing=dspacing1

    if isort:
        sin2psi, dspacing, ehkl, psi = sort(sin2psi, dspacing, ehkl, psi)

    return sin2psi, ehkl, dspacing, psi

def ex_01(fns=[
        'powder_collection_25JUL12/25JUL12_0002Data1Phi-90.txt',
        'powder_collection_25JUL12/25JUL12_0002Data2Phi0.txt',
        'powder_collection_25JUL12/25JUL12_0002Data3Phi45.txt',
        'powder_collection_25JUL12/25JUL12_0002Data4Phi135.txt']):
    from MP.lib.mpl_lib import wide_fig as wf
    from MP.lib.mpl_lib import fancy_legend
    from MP.lib.axes_label import __deco__ as deco

    fig = wf(nh=1,nw=1,uw=4,left=0.2); ax=fig.axes[0]
    ehkl_dat = []
    for i in xrange(len(fns)):
        sign_sin2psi,ehkl,dspacing,psi = read(fns[i],iopt=0,isort=True)
        ax.plot(sign_sin2psi,ehkl,'-x',label=fns[i].split('Phi')[1].split('.')[0])
        ehkl_dat.append(ehkl)


    print 'mean:', np.array(ehkl_dat).mean()
    print 'std:', np.array(ehkl_dat).std()
    deco(iopt=0,ft=15,ax=ax,ipsi_opt=1)
    fancy_legend(ax=ax,size=10,ncol=2)
