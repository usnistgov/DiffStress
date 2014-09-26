## MS&T 2014 Fall meeting application
## - plot flow stress along various paths

from os import sep, popen
from glob import glob
from RS import rs_ex
from MP.lib import mpl_lib,axes_label
import matplotlib.pyplot as plt
import numpy as np
from MP.mat.mech import find_err
wf = mpl_lib.wide_fig
fl = mpl_lib.fancy_legend

def main_plot_flow_all(hkl='211',sin2psimx=0.5,
                       psi_nbin=1,pmargin=None,iplot_rs=False,
                       path=''):
    """  """
    fns = read(hkl,prefix=path)
    fig = wf(nw=3,nh=1,left=0.2,uw=3.5,
             w0=0,w1=0.3,right=0,iarange=True)
    ax1,ax2,ax3 = fig.axes
    axes_label.__eqv__(ax1,ft=10)
    errors = []
    for i in range(len(fns)):
        fn = fns[i]
        strain_path = fn.split('_')[0]
        untargz(fn)
        model_rs, fwgt, fdsa \
            = rs_ex.ex_consistency(
                sin2psimx=sin2psimx,
                psi_nbin=psi_nbin,
                hkl=hkl,iplot=iplot_rs,
                pmargin=pmargin,path=path)
        fwgt.get_eqv(); fdsa.get_eqv()

        e = find_err(fwgt,fdsa)
        ax3.plot(fwgt.epsilon_vm,e)
        errors.append([fwgt.epsilon_vm,e])

        ax1.plot(fwgt.epsilon_vm,fwgt.sigma_vm,
                 'bx-',label='Weighted Avg', alpha=1.0)
        ax1.plot(fdsa.epsilon_vm,fdsa.sigma_vm,
                 'kx',label='Diff. Stress Analsysis', alpha=1.0)
        ax2.plot(fwgt.sigma[0,0],fwgt.sigma[1,1],'b+')
        ax2.plot(fdsa.sigma[0,0],fdsa.sigma[1,1],'k+')

        ## connector
        npoints = len(fwgt.sigma[0,0])
        wgtx = fwgt.sigma[0,0]; wgty = fwgt.sigma[1,1]
        dsax = fdsa.sigma[0,0]; dsay = fdsa.sigma[1,1]
        for j in range(npoints):
            ax2.plot([wgtx[j],dsax[j]],[wgty[j],dsay[j]],'k-',alpha=0.5)
        plt.show();plt.draw()

    ax1.set_xlim(-0.1,1.1); ax1.set_ylim(-50,700);
    ax2.set_xlim(-10,1000); ax2.set_ylim(-10,1000)
    ax2.set_aspect('equal')
    ax2.set_xlabel(r'$\bar{\Sigma}_{11}$',dict(fontsize=15))
    ax2.set_ylabel(r'$\bar{\Sigma}_{22}$',dict(fontsize=15))
    ax2.locator_params(nbins=3)
    ax2.grid('on'); plt.show()

    fig.savefig('flow_allpaths_%s.pdf'%(hkl))
    plt.close(fig)

    return np.array(errors)[0]

def error_est_conds(psi_nbin=21,sin2psimx=0.5):
    paths=['BB','URD','UTD','PSRD','PSTD']
    for i in range(len(paths)):
        error_est(
            paths[i],psi_nbin=psi_nbin,
            sin2psimx=sin2psimx,iplot_rs=False)

def error_est(
        path='BB',psi_nbin=21,
        sin2psimx=0.5,
        hkls=['200','220','211','310'],
        iplot_rs=False):
    """
    """
    from MP.lib.axes_label import __deco__ as deco
    fig = wf(nw=1,nh=1,left=0.2,uw=3.5,
             w0=0,w1=0.3,right=0,iarange=True)
    ax=fig.axes[0]
    for i in range(len(hkls)):
        err = main_plot_flow_all(
            hkl=hkls[i],sin2psimx=sin2psimx,
            psi_nbin=psi_nbin,pmargin=None,
            iplot_rs=iplot_rs,path=path)
        eps_vm, e = err
        ax.plot(eps_vm,e,label=hkls[i])
    deco(iopt=8,ft=15,ax=ax)
    ax.legend(loc='best',fontsize=9)
    fig.savefig('error_dsa_%s.pdf'%path)
    return fig

def untargz(fn): std = popen('tar -xzf %s'%fn)

def read(hkl='211',prefix=''):
    fns=glob('%s*_%s_*.tar.gz'%(prefix,hkl))
    print '%i files are found'%len(fns)
    return fns

def plot_rsq():
    from MP.lib import mpl_lib
    import numpy as np
    import matplotlib.pyplot as plt
    from pepshkl import reader2 as reader
    import sfig_class
    from MP.mat.mech import FlowCurve as FC
    StressFactor=sfig_class.SF

    flow = FC(); flow.get_model(fn='STR_STR.OUT');flow.get_eqv()
    tdat,usf  = reader(fn='igstrain_fbulk_ph1.out',iopt=1,isort=True)
    vdat,ngrd = return_vf()
    wf = mpl_lib.wide_fig
    nstp,nij,nphi,npsi = tdat.shape[1:]

    Psi = tdat[8,:,:,:,:]
    Phi = tdat[7,:,:,:,:]
    rsq = tdat[9,:,:,:,:]
    sf  = tdat[3,:,:,:,:]

    ## swap sf axes to comply with sfig_class.SF.sf
    ## [nstp,nphi,npsi,nij]

    sf = sf.swapaxes(1,3).swapaxes(1,2)
    rsq=rsq.swapaxes(1,3).swapaxes(1,2)

    psi = Psi[0,0,0,:]
    phi = Phi[0,0,:,0]
    print psi.shape
    print phi.shape
    print sf.shape

    SF=StressFactor()
    SF.add_data(sf=sf,phi=phi,psi=psi)
    SF.flow = flow ## overwrite flow
    SF.add_vf(vdat)
    SF.add_rsq(rsq)
    SF.plot()

    return sf,vdat,rsq,psi


def return_vf():
    """
    Return volume fraction
    """
    from pepshkl import reader4
    tdat,_psis_,_vdat_,_ngrd_ = reader4(
        fn='int_els_ph1.out',ndetector=2,iopt=1)

    ## sorting
    import MP.ssort as sort
    shsort = sort.shellSort
    ss     = sort.ind_swap

    ## sort psi and find ind
    psis, ind = shsort(_psis_)
    nstp, nphi, npsis = np.shape(_vdat_)
    vdat = np.zeros((nstp,nphi,npsis))
    ngrd = np.zeros((nstp,nphi,npsis))

    for istp in range(nstp):
        for iphi in range(nphi):
            dum_v = _vdat_[istp,iphi,:][::]
            vdat[istp,iphi,:] = ss(dum_v[::],ind)[::]

            dum_n = _ngrd_[istp,iphi,:][::]
            ngrd[istp,iphi,:] = ss(dum_n[::],ind)[::]
    return vdat, ngrd


def plot_sf(sff_fn='temp.sff',pmargin=0.1):
    """
    Arguments
    =========
    sff_fn = 'temp.sff'
    pmargin = 0.1
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from pepshkl import reader4
    from RS import sff_converter
    from rs_exp import read_IGSF
    from MP.lib import mpl_lib
    wf = mpl_lib.wide_fig

    sff_converter.main(fn=sff_fn,difile=None,itab=True,
                       ieps0=4,  fn_str='STR_STR.OUT')
    SF, dum = read_IGSF(fn=sff_fn,fn_str='STR_STR.OUT')


    tdat,_psis_,_vdat_,_ngrd_ = reader4(
        fn='int_els_ph1.out',ndetector=2,iopt=1)

    ## sorting
    import MP.ssort as sort
    shsort = sort.shellSort
    ss     = sort.ind_swap

    ## sort psi and find ind
    psis, ind = shsort(_psis_)
    nstp, nphi, npsis = np.shape(_vdat_)
    vdat = np.zeros((nstp,nphi,npsis))
    ngrd = np.zeros((nstp,nphi,npsis))

    for istp in range(nstp):
        for iphi in range(nphi):
            dum_v = _vdat_[istp,iphi,:][::]
            vdat[istp,iphi,:] = ss(dum_v[::],ind)[::]

            dum_n = _ngrd_[istp,iphi,:][::]
            ngrd[istp,iphi,:] = ss(dum_n[::],ind)[::]

    SF.mask_vol() ## shake off all zero values
    SF.add_vf(vdat)
    ## Mask data lacking sufficient volume
    SF.mask_vol(pmargin=pmargin)

    fig=wf(nw=1,nh=1)
    ax=fig.axes[0]
    for istp in range(nstp):
        sign_sin2psi = np.sign(SF.psi)*np.sin(SF.psi*np.pi/180.)**2
        ax.plot(sign_sin2psi,SF.vf[istp,0,:],'--')

    SF.plot(nbin_sin2psi=3,iopt=1)

    return SF,tdat,psis,vdat,ngrd

def plot_sf_psis(
        sff_fn='temp.sff',
        psi_ref=[45.0, 42.4, 39.8, 37.1, 34.3, 31.5,
                 28.5, 25.2, 21.7, 17.5, 12.3, 0]):
    import numpy as np
    import matplotlib.pyplot as plt
    from pepshkl import reader4
    from RS import sff_converter
    from rs_exp import read_IGSF
    from MP.lib import mpl_lib, axes_label
    from RS.rs import find_nearest

    deco = axes_label.__deco__
    wf = mpl_lib.wide_fig

    sff_converter.main(fn=sff_fn,difile=None,itab=True,
                       ieps0=4,  fn_str='STR_STR.OUT')
    SF, dum = read_IGSF(fn=sff_fn,fn_str='STR_STR.OUT')

    tdat,_psis_,_vdat_,_ngrd_ = reader4(
        fn='int_els_ph1.out',ndetector=2,iopt=1)

    ## sorting
    import MP.ssort as sort
    shsort = sort.shellSort
    ss     = sort.ind_swap

    ## sort psi and find ind
    psis, ind = shsort(_psis_)
    nstp, nphi, npsis = np.shape(_vdat_)
    vdat = np.zeros((nstp,nphi,npsis))
    ngrd = np.zeros((nstp,nphi,npsis))

    for istp in range(nstp):
        for iphi in range(nphi):
            dum_v = _vdat_[istp,iphi,:][::]
            vdat[istp,iphi,:] = ss(dum_v[::],ind)[::]

            dum_n = _ngrd_[istp,iphi,:][::]
            ngrd[istp,iphi,:] = ss(dum_n[::],ind)[::]

    ## Find only psi_ref values
    ## sin2psi_ref = np.sin(psi_ref*np.pi/180.)**2
    inds = []
    for i in range(len(psi_ref)):
        inds.append(find_nearest(SF.psi, psi_ref[i]))

    SF.mask_vol() ## shake off all zero values

    fig=wf(nh=nphi,nw=len(psi_ref),iarange=True)
    figv=wf(nh=nphi,nw=len(psi_ref),iarange=True)
    for iphi in range(nphi):
        for ipsi in range(len(psi_ref)):
            #ax = fig.axes[iphi+nphi*ipsi]
            ax = fig.axes[ipsi+len(psi_ref)*iphi]
            axt = figv.axes[ipsi+len(psi_ref)*iphi]
            if iphi==0:
                ax.set_title(
                    r'$\psi= %3.1f^\circ{}$'%\
                    SF.psi[inds[ipsi]]
                )
                axt.set_title(
                    r'$\psi= %3.1f^\circ{}$'%\
                    SF.psi[inds[ipsi]]
                )
            if ipsi==0:
                ax.text(0,0,r'$\phi=%3.1f^\circ{}$'%\
                        SF.phi[iphi],transform=ax.transAxes)
                axt.text(0,0,r'$\phi=%3.1f^\circ{}$'%\
                        SF.phi[iphi],transform=ax.transAxes)

                deco(ax,iopt=6)
                deco(axt,iopt=7)


            x = SF.flow.epsilon_vm
            y = SF.sf[:,iphi,inds[ipsi],0] # f11
            v = vdat[:,iphi,inds[ipsi]]
            ax.plot(x,y*1e12,'b-o')
            axt.plot(x,v,'r-x')

            axt.set_ylim(0.,0.10)

            ax.set_xlim(0.,1.0)
