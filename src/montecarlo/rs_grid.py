"""
An application for DiffStress package:
Conduct Monte Carlo experiments to quantify uncertainties
pertaining to multiaxial flow stress measurements.

ex01:
combination of N(psi) and N(F). In each cell,
a multiple number of precision error and bias are calculated.
"""
from MP.lib import whichcomp
whereami = whichcomp.guessWhereami()
def determineEnvironment(whereami):
    if whereami=='palmetto':
        submitCommand = 'qsub'
    else:
        submitCommand = None
    from MP.lib import checkX
    if checkX.main()!=0:
        availX = False
    else:
        availX = True
    return submitCommand, availX

submitCommand, availX = determineEnvironment(whereami=whereami)

## dependency
if not(availX):
    import matplotlib as mpl
    mpl.use('Agg') ## In case X-window is not available.

from numpy import sin, arcsin, pi, sqrt
from MP.lib import mpl_lib, axes_label
import tempfile,time,os
from glob import glob
from jac_par_ex1 import plot as reader
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
# from detect_ques import find_availability as fa

try:
    import cPickle as pickle
except:
    import pickle

## convert rules between sin2psi and psi
## psi to sin2psi
fpsi = lambda psi: sin(psi*pi/180.)**2
fx   = lambda   x: arcsin(sqrt(x)) * 180./ pi

def tick_bin(ax,nbin=3):
    from matplotlib.ticker import MaxNLocator
    ax.xaxis.set_major_locator(MaxNLocator(nbins=nbin))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=nbin))
    ax.minorticks_on()

def fpsis(*args):
    lst = []
    for arg in args: lst.append(fpsi(arg))
    return lst

def main_ex01(
    ## main arguments (2d plane)
    N_psi    = [3,5,9,13,30], ## x-axis
    # N_F      = [30,18,6],    ## y-axis, ## N_F      = [9,4,3],    ## y-axis
    N_F      = [98,49,33,24],

    ## secondary arguments
    ensemble = 5,
    sigmas   = [5e-5,2.5e-5,1e-5], ## CSE sigma
    bragg    = round(78.2*np.pi/180.,4),    ## {211} Cr
    # bragg    = round(81.1*np.pi/180.,4),    ## {310} Co

    ird      = 0.0457,             ## Random distribution Intensity with w=5 degree for {211, 310} (multiplicity)
    ## original bounds (up to the real experimental tilting angle maximum)
    # bounds   = [0, 0.329],
    bounds   = [0, 0.750], ## original bounds (up to 60 degree)

    nfrq     = 50,
    isub     = True,

    ## title written in PBS job file through write_bash_rs
    title='MC-rs_grid',
    ##
    ncpu     = False,
    noti=True):
    """
    Choice of axes: N(psi) and N(F).
    Dynamically renew the resources if iwait==True.

    Prepare/Run Monte Carlos virtual experiments for
    a 3D matrix of (N_psi, N_F, S^CSE).

    First of all, prepare bash scripts <rs_par_i_j.sh>
    and then submit jobs if <isub>==True.

    The feature that dynamically assess the availability of
    computing resource has been deprecated since the method
    was based on CTCMS and is not suitable for Palmetto in
    Clemson - 20160527

    Arguments
    =========
    ## main arguments
    N_psi = [7,11,15,29]
    N_F   = [9,6,3,1]

    ## secondary arguments

    ## number of ensembles to be created
    ensemble = 10

    ## Array includes levels of count. Stats. Error under examination
    sigmas   = [1e-5,1e-4,5e-4]

    ## tilting angle boundary in terms of +/- sin2psi value
    bounds   = [0, 0.5]

    ## Flag to determine if the job will be submitted
    isub     = True

    ## Flag used in CTCMS cluster (deprecated)
    noti     = False
    """
    ## sorting the CSE standard deviation according
    ## to the ascending order
    sigmas.sort()

    njobs = len(N_psi) * len(N_F)

    if type(ncpu).__name__=='NoneType':
        if whereami=='palmetto':
            ncpu=8 ## at least 8 cpus
        elif whereami=='mac':
            ncpu=4
        else:
            print 'whereami:',whereami
            raise IOError, 'Could not determine where am I'
    else:
        ## the number of ncpu is as-given
        pass

    bash_fns=[]
    for i in range(len(N_F)):
        for j in range(len(N_psi)):
            ### Begin sh scripting
            fn = 'rs_par_%1i_%1i.sh'%(i+1,j+1)
            bash_fns.append(fn)
            with open(fn,'w') as f:
                write_bash_rs(
                    f,i,j,noti,
                    'unknown',
                    title=title,
                    ## Below: **kwargs
                    sigmas=sigmas,
                    bounds=bounds,
                    bragg=bragg,
                    ird=ird,
                    ss=N_F[i],
                    nbins=N_psi[j],
                    nsample=ensemble,
                    ncpu=ncpu,
                    nfrq=nfrq)
            ## end of sh scripting

            ## submitting the job right away
            if isub:
                cmd = '%s %s'%(submitCommand, fn)
                os.popen(cmd)
    return bash_fns


## list of ploting functions
# plot_ex01
# plot_ex02

def hist(ifig=2):
    """
    Plot an examplary histogram
    to demonstrate the behavior of the scatter
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    from jac_par_ex1 import plot as reader

    fns = glob('rs_par_?_?.dat')
    fns.sort()
    dat, dicts = reader(fns[0],irst=True)
    sigmas = dicts['sigmas']

    hist=np.load('hist.npy'); bins=np.load('bins.npy')

    nsig,nstpx,nb = hist.shape
    ## analyze up to 5 initial plastic increments.
    if nsig != len(sigmas):
        raise IOError, \
            'incompatible sigmas found between hist and *.dat'

    nstp = 1
    fig = plt.figure(ifig,figsize=(6,3))
    gs = gridspec.GridSpec(
        nstp,nsig,wspace=0.0,hspace=0.0,
        left=0.2, top=0.80, bottom=0.2)

    for isig in range(nsig):
        for istp in range(nstp):
            h = hist[isig,istp,:]
            b = bins[isig,istp,:]
            width = 1.00 * (b[1]-b[0])
            center = (b[:-1]+b[1:])/2.
            ax = fig.add_subplot(gs[istp,isig])

            ax.bar(center,h,align='center',width=width,color='gray')
            xl = ax.get_xticks()[0]
            xh = ax.get_xticks()[-1]

            plt.locator_params(nbins=3)
            ax.set_xlim(-0.20,0.20)

            if True: #isig==0:
                if isig==0: s='$s^{CSE} =  $%i'%(
                    sigmas[isig]*1e6),
                elif isig!=0 and isig!=nsig-1: s='%i'%(
                    sigmas[isig]*1e6)
                elif isig==nsig-1:s='%i $\mu$ strain'%(
                    sigmas[isig]*1e6)
                ax.annotate(
                    s=r'%s'%s,
                    size=12,
                    horizontalalignment='center',
                    verticalalignment='center',
                    # rotation='vertical',
                    xy=(0.5,1.1),
                    xycoords='axes fraction')

            if isig==0 and istp==0:
                ax.set_xlabel(r'$\sigma^e$')
                ax.set_ylabel('Frequency')

            else:
                mpl_lib.rm_lab(ax,axis='x')
                mpl_lib.rm_lab(ax,axis='y')

    mpl_lib.tune_x_lim(fig.axes,axis='y')

    fig.savefig('hist.pdf')

def plot_ex01(iplot=True,isave=True):
    """
    1. Plots the data resulting from main_ex01
    2. Save dat from main_ex01

    Variables: N(DEC):   x axis (horizontal);
               N(psi): y axis (vertical);

    Plot in matric of (N_DEC, N_PSI)

    Arguments
    =========
    iplot = True
    isave = True

    Returns
    =======
    nfs, npsis, scse, bias, precision (if iplot==False)
    """
    wc = 'rs_par_?_?.dat' ## can be an argument
    # data from rs_par_y_x.dat
    fns = glob(wc)
    fns.sort()

    if len(fns)==0:
        print 'Not much to write home about'
        return

    if iplot:
        fig1=plt.figure(figsize=(6,3))
        fig1.add_axes((0.12,0.20,0.33,0.75))
        fig1.add_axes((0.62,0.20,0.33,0.75))
        ax1=fig1.axes[0]
        ax2=fig1.axes[1]
        fig=plt.figure(figsize=(12,12))

    day,m,date,t,y = time.asctime().split()
    _fn_pdf_ = y+m+date+'_'+t
    print _fn_pdf_,
    if isave:
        fn_tar = tempfile.mkstemp(
            dir=os.getcwd(),prefix='rs_grid_%s__'%(y+m+date),
            suffix='__.tar')[1]
        print fn_tar
        cmd = 'tar -cvf %s '%fn_tar
        for i in xrange(len(fns)):
            cmd = '%s %s'%(cmd,fns[i])
        os.system(cmd)
    else:
        print

    N_DEC=int(fns[-1].split('.dat')[0].split('_')[2])
    N_PSI=int(fns[-1].split('.dat')[0].split('_')[3])
    nx = N_DEC
    ny = N_PSI
    dum, dum1 = reader(fn=fns[0], irst=True)
    nz=len(dum1['sigmas'])
    if iplot:
        gs=gridspec.GridSpec(
            ncols=N_DEC,nrows=N_PSI,wspace=0.0,
            hspace=0.0,left=0.22,
            right=0.8,top=0.80)
        ls=['^','s','o','d','+']
        fc=['k','k','k','k','k']
    nfs=[];npsis=[];scse=dum1['sigmas']
    ## x, y, z (NF, Npsi, sCSE)

    precision = np.zeros((nx,ny,nz))
    bias      = np.zeros((nx,ny,nz))
    axes_xy   = []
    ixys=[]
    for i in xrange(len(fns)): ## In each cell of (nx,ny) matric
        ix=int(fns[i].split('.dat')[0].split('_')[2]) # N_DEC
        iy=int(fns[i].split('.dat')[0].split('_')[3]) # N_PSI
        ixys.append([ix,iy])
        if True:
            try:
                dat,dicts = reader(fn=fns[i], irst=True)
            except:
                pass
            else:
                ## flow_raw
                fn_fc = fns[i].split('.dat')[0]+'.pck'
                if not(os.path.isfile(fn_fc)):
                    raise IOError,\
                        'Could not find the pickled FC object'

                with open(fn_fc,'rb') as pck_fo:
                    fc_raw = pickle.load(pck_fo)

                sigmas = dicts['sigmas']
                ss     = dicts['ss']
                nbins  = dicts['nbins']
                irow = iy-1; ## N_psi
                icol = ix-1; ## N_DEC
                EPS = fc_raw.epsilon_vm   # dat[0]
                if icol==0: npsis.append(nbins)
                if irow==0: nfs.append(len(EPS[::ss]))
                for j in xrange(len(sigmas)): ## sDEC
                    mu = dat[1+j]
                    sigma1 = dat[1+j+len(sigmas)]
                    precision[icol,irow,j] \
                        = np.average(sigma1)*100.*2
                    bias[icol,irow,j] \
                        = np.average(np.abs(mu))*100.

    if iplot:
        for i in xrange(len(fns)): ## In each cell of (nx,ny) matric
            ## flow_raw
            fn_fc = fns[i].split('.dat')[0]+'.pck'
            if not(os.path.isfile(fn_fc)):
                raise IOError, 'Could not find the pickled FC object'
            with open(fn_fc,'rb') as pck_fo:
                try:
                    fc_raw = pickle.load(pck_fo)
                except:
                    raise IOError, 'could not load pickle from %s'%fn_fc

            ix, iy = ixys[i]
            irow = iy-1; icol = ix-1
            dat, dicts = reader(fn=fns[i], irst=True)
            sigmas = dicts['sigmas']
            ss     = dicts['ss']
            nbins  = dicts['nbins']
            nfrq   = dicts['nfrq']
            irow   = iy-1; icol = ix-1
            EPS    = dat[0] ## Reduced eps due to nfrq
            _E_    = fc_raw.epsilon_vm  ## Full eps
            ndec=len(_E_[::ss]) ## Reduced DEC acqusition

            ax  = fig.add_subplot(gs[irow,icol])
            for j in xrange(len(sigmas)): ## sDEC
                mu = dat[1+j]
                sigma1 = dat[1+j+len(sigmas)]
                ax.fill_between(
                    EPS, (mu+sigma1)*100,(mu-sigma1)*100,
                    facecolor=fc[j],alpha=0.2)

                ## -----------------------------------------#
                ## Markers where DECs are measured
                if irow==0 and j==0:
                    xs=_E_[::ss]
                    for i_lev in xrange(len(xs)):
                        x=xs[i_lev]
                        ax.plot(x,0,'x',mec='k')
                ## -----------------------------------------#

            if ix==1 and iy==nx: ax_p = ax
            if ix==nx and iy==ny:
                print '(nx,ny):','(,',nx,',',ny,')'
                ax.annotate(s=r'$\mathrm{(N^{DEC}\ ,N^{\psi})}$',
                            horizontalalignment='center',
                            xy=(0.5,-0.2),
                            xycoords='axes fraction')
            tick_bin(ax,5)
            ax.minorticks_on()
            ax.annotate(
                s=r'(%i,%i)'%(ndec,nbins),
                horizontalalignment='center',
                xy=(0.2,0.10),
                xycoords='axes fraction')
            if irow==0:
                ax.annotate(
                    s=r'%i DECs'%ndec,
                    horizontalalignment='center',
                    xy=(0.5,1.2),
                    xycoords='axes fraction')
            if icol==0:
                ax.annotate(
                    s=r'$\mathrm{N^{\psi}=%i}$'%nbins,
                    verticalalignment='center',
                    rotation='vertical',
                    xy=(-0.8,0.5),
                    xycoords="axes fraction")

    precision = np.array(precision)
    bias = np.array(bias)

    print 'npsis'
    print npsis
    print 'nfs'
    print nfs
    print bias.shape
    print precision.shape

    if not(iplot):
        return nfs,npsis,scse,bias,precision

    ms =['^','s','o','>','d']
    ls =['-','--','-.','-','-']
    if iplot:
        for iy in xrange(ny): ## Npsi
            for iz in xrange(nz): ## sCSE
                if iz==0:
                    label=r'$N^{\psi} = %i$'%npsis[iy]
                    ax1.plot(-100,0,ls=ls[iy],color='k',
                             label=label)
                ax1.plot(
                    nfs,bias[:,iy,iz],ls=ls[iy],
                    marker=ms[iz],color='k',mfc='None')
                ax2.plot(
                    nfs,precision[:,iy,iz],ls=ls[iy],
                    marker=ms[iz],color='k',mfc='None')
        for iz in xrange(nz):
            label=r'$s^\mathrm{CSE} = %2.2i \mu$'%(
                sigmas[iz]*1e6)
            ax1.plot(-100,0,ls='None',marker=ms[iz],color='k',
                     label=label,mfc='None')
        #ax1.set_ylim(-1,max(np.array(bias).flatten()))
        ax1.set_xlim(min(nfs),max(nfs))
        ax2.set_xlim(min(nfs),max(nfs))
        ax1.minorticks_on();ax1.grid('on')
        ax2.minorticks_on();ax2.grid('on')
        ax1.set_xlabel(r'$N^\mathrm{DEC}$',dict(fontsize=12))
        ax1.set_ylabel('Bias [%]',dict(fontsize=12))
        ax2.set_xlabel(r'$N^\mathrm{DEC}$',dict(fontsize=12))
        ax2.set_ylabel('Error in precision [%]',
                       dict(fontsize=12))
        mpl_lib.rm_all_lab(fig.axes)
        for i in xrange(len(fig.axes)):
            fig.axes[i].grid(True)
        mpl_lib.ticks_bins_ax_u(fig.axes, 3)
        mpl_lib.ticks_bins_ax_u(fig1.axes,3)

        plt.setp(ax_p.get_xticklabels(),visible=True)
        plt.setp(ax_p.get_yticklabels(),visible=True)
        axes_label.__deco__(iopt=9,ft=13,ax=ax_p)
        ax_p.set_ylabel(
            'Error in \n diffraction stress [%]',
            dict(fontsize=10),multialignment='center')
        ax_p.set_xlabel(
            'Macroscopic \n '+\
            r'Von Mises strain $\mathrm{\bar{E}^{VM}}$',
            dict(fontsize=10))
        ## deco figure
        fig=plt.gcf()
        for i in xrange(len(fig.axes)):
            fig.axes[i].set_ylim(-8,8)
            fig.axes[i].set_xlim(0.,max(EPS))
        fig.savefig(
            'fig_param_%s.pdf'%_fn_pdf_,
            bbox_inches='tight')
        fig1.savefig(
            'fig_trend_%s.pdf'%_fn_pdf_,
            bbox_inches='tight')

    return nfs,npsis,scse,bias,precision


def write_bash_rs(
    f,i,j,noti,que_kinds,ncpu,
    title,**kwargs):
    """
    Thie is to write a bash script prepared
    for a multi-threaded run

    Arguments
    ---------
    f         : file object
    i         : x index
    j         : y index
    noti      : flag to notify the status of job
    que_kinds : type of que
    ncpu

    **kwargs  : key-worded arguments
                to be passed to jac_par_ex1.py
    """
    ## -- block used for CTCMS cluster
    ## ## qsub arguments
    ## f.write('#!/bin/bash\n')
    ## f.write('## RS package parametric study to')
    ## f.write(' evaluate uncertainty\n')
    ## if que_kinds=='fast':
    ##     f.write('#$ -q fast\n')
    ##     f.write('#$ -l short=TRUE\n')
    ## else:
    ##     pass ## wide64/wide2/rack4

    ## f.write('#$ -pe nodal %i\n'%use_cpu)
    ## f.write('#$ -M ynj\n')
    ## if noti: f.write('#$ -m be\n')
    ## f.write('#$ -cwd\n\n')
    ## f.write('# Job command\n')
    ##

    f.write('#!/bin/bash\n')
    #f.write('#PBS -N MonteCarlo-DiffStress-Uncertainty\n')
    f.write('#PBS -N %s\n'%title)
    f.write('#PBS -l select=1:ncpus=%i:mem=%igb,walltime=02:00:00\n'%(ncpu,int(ncpu/2.)))
    if noti: f.write('#PBS -M youngunj@gmail.com -m bea\n')

    f.write('# -------------------------------------------------- #\n')
    f.write('#           Command to be executed below             #\n')
    f.write('# -------------------------------------------------- #\n')

    f.write('\n## change path to working directory\n')
    f.write('cd %s\n\n'%os.getcwd()) ## Work directory

    f.write('python')
    f.write(' jac_par_ex1.py ')
    ## arguments for jac_par_ex1.py
    f.write(" --filename ")
    f.write('rs_par_%1i_%1i.dat'%(i+1,j+1))

    for key in kwargs:
        arg_type = type(kwargs[key]).__name__
        if arg_type=='int':
            f.write(" --%s %i"%(key,kwargs[key]))
        elif arg_type=='list' or arg_type=='ndarray':
            f.write(" --%s '"%key)
            list_args = kwargs[key]
            for element in list_args:
                f.write(" %9.3e"%element)
            f.write("'")
        elif arg_type=='float':
            f.write(" --%s %f"%(key,kwargs[key]))
        else:
            raise IOError, \
                'Unexpected argument type'
    f.write(' --ncpu %i'%ncpu)
    f.write('\n\n')
    pass

def run_211(psi_mx=60,thb=78.2,
            ncpu=8, 
            ensemble=200,
	    title='MC_run211'):
    main_ex01(
        N_psi=[3,5,9,13,30],
        N_F=[99,49,33,19,9],
        ensemble=ensemble,
        sigmas=np.array([50,25,10])*1e-6,
        bragg=thb*np.pi/180., ## {211} Cr
        ird=0.0457,
        bounds=[0,np.sin(psi_mx*np.pi/180)**2],
        nfrq = 9,
        ncpu=ncpu,
        title=title,
	isub=True)

def run_310(psi_mx=60,thb=81.1,
            ncpu=8,
            ensemble=200,
	    title='MC_run211'):
    main_ex01(
        N_psi=[3,5,9,13,30],
        N_F=[99,49,33,19,9],
        ensemble=ensemble,
        sigmas=np.array([50,25,10])*1e-6,
        bragg=thb*np.pi/180.,   ## {310} Co
        ird=0.0457,
        bounds=[0,np.sin(psi_mx*np.pi/180)**2],
        nfrq = 9,
        ncpu=ncpu,
        title=title,
	isub=True)

def copyAndExt(fn):
    """

    Argument
    --------
    fn
    """
    cmd1 = 'cp %s ./'%fn
    os.system(cmd1)
    cmd2 = 'tar -xvf %s'%fn
    os.system(cmd2)

def prep_310():
    """
    copy {310} plane results to the current working directory
    """
    path = '~/repo/evpsc/ipynb/stress_uncertainty/dat/IFsteel/Snapshots/'
    fn = '20150518_BB_IF20k_5deg_310.tar.gz'
    fn = os.path.join(path,fn)
    copyAndExt(fn)

def prep_211():
    """
    copy {211} plane results to the current working directory
    """
    path = '~/repo/evpsc/ipynb/stress_uncertainty/dat/IFsteel/Snapshots/'
    fn = '20150427_BB_IF20k_5deg.tar'
    fn = os.path.join(path,fn)
    copyAndExt(fn)

def ex01(ncpu,ensemble,title):
    prep_211()
    run_211(psi_mx=60,thb=78.2,ncpu=ncpu,ensemble=ensemble,title=title)
def ex02(ncpu,ensemble,title):
    prep_211()
    run_211(psi_mx=35,thb=78.2,ncpu=ncpu,ensemble=ensemble,title=title)
def ex03(ncpu,ensemble,title):
    prep_211()
    run_211(psi_mx=15,thb=78.2,ncpu=ncpu,ensemble=ensemble,title=title)
def ex04(ncpu,ensemble,title):
    prep_310()
    run_310(psi_mx=60,thb=81.1,ncpu=ncpu,ensemble=ensemble,title=title)
def ex05(ncpu,ensemble,title):
    prep_310()
    run_310(psi_mx=35,thb=81.1,ncpu=ncpu,ensemble=ensemble,title=title)
def ex06(ncpu,ensemble,title):
    prep_310()
    run_310(psi_mx=15,thb=81.1,ncpu=ncpu,ensemble=ensemble,title=title)

## command line usage
if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i',type=int,help='Choice of examples')
    parser.add_argument(
        '--ncpu',type=int,
        help='Number of CPU cores',default=8)
    parser.add_argument(
        '--isub',type=int,
        help='Flat to submit the job --  0: no, 1: yes',default=1)
    parser.add_argument(
        '--ensemble',type=int,
        help='the number of ensemble',default=200)
    args=parser.parse_args()

    if   args.isub==0: isub = False
    elif args.isub==1: isub = True

    if args.i==0:
        ## test...
        from RS import rs_ex
        rs_ex.ex_consistency(
            psimx=60., iscatter=False,
            sigma=5e-5, ig_sub=True,ird=0.0457,
            nfrq=9, psi_nbin=5,theta_b=78.2*np.pi/180.,iplot=True)
    elif args.i==-1:
        plot_ex01()
    elif args.i==1:
        ex01(ncpu=args.ncpu,ensemble=args.ensemble,title='MC-ex01')
    elif args.i==2:
        ex02(ncpu=args.ncpu,ensemble=args.ensemble,title='MC-ex02')
    elif args.i==3:
        ex03(ncpu=args.ncpu,ensemble=args.ensemble,title='MC-ex03')
    elif args.i==4:
        ex04(ncpu=args.ncpu,ensemble=args.ensemble,title='MC-ex04')
    elif args.i==5:
        ex05(ncpu=args.ncpu,ensemble=args.ensemble,title='MC-ex05')
    elif args.i==6:
        ex06(ncpu=args.ncpu,ensemble=args.ensemble,title='MC-ex06')
    else:
        raise IOError, 'Not read yet for ex %i'%args.i
