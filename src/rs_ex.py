"""
A collection of scripts that make use of residual stress module

Func ex_consistency is the core of this module
which will be used for Monte Carlo experiments to estimate
the uncertainty in the stress measurements
"""

import matplotlib as mpl
mpl.use('Agg') ## In case X-window is not available.


import numpy as np
from numpy import pi, sin, cos
from lib import write_args

def ex(ifig=50,
       exp_ref=['exp_dat/Bsteel/EXP_BULGE_JINKIM.txt',
                'exp_dat/Bsteel/uni/avgstr_000.txt'],
       exp_lab=['Exp bulge','Exp uniaxial'],
       mod_ref='STR_STR.OUT'):
    """
    Arguments
    =========
    exp_ref = ['exp_dat'/Bsteel/EXP_BULGE_JINKIM.txt']
    exp_lab = ['Exp bulge', 'Exp uniaxial')
    mod_ref = 'STR_STR.OUT'
    """
    fig=plt.figure(ifig);ax=fig.add_subplot(111)
    myrs = ResidualStress()
    strains = myrs.straine #(20,6)
    # strains = myrs.strainm #(4,6)
    strain_eff = strains.T[0]+strains.T[1]
    strain_eff = strain_eff.T
    stress = []
    for i in xrange(len(strain_eff)):
        dum = myrs.analysis(iopt=0,istp=i)
        stress.append((dum[0] + dum[1])/2.)
    ax.plot(strain_eff,stress,'-x',label='Model SF/IG/ehkl')

    stress = []
    for i in xrange(len(strain_eff)):
        dum = myrs.analysis(iopt=1,istp=i)
        stress.append((dum[0] + dum[1])/2.)
    ax.plot(strain_eff,stress,'-x',label='Exp SF/IG/ehkl')

    stress = []
    for i in xrange(len(strain_eff)):
        dum = myrs.analysis(iopt=2,istp=i)
        stress.append((dum[0] + dum[1])/2.)
    ax.plot(strain_eff,stress,'-x',label='Model SF/IG + Exp ehkl')


    for i in xrange(len(exp_ref)):
        exp_ref = exp_ref[i]
        x,y=np.loadtxt(exp_ref).T
        ax.plot(x,y,'--',label=exp_lab[i])

    dum=np.loadtxt(mod_ref,skiprows=1).T;
    x=dum[2]+dum[3];y=(dum[8]+dum[9])/2.;

    ax.plot(x,y,'--',label='EVPSC biaxial')
    __deco__(ax,iopt=3)

    ax.legend(loc='best',fancybox=True).get_frame().set_alpha(0.5)
    fig.tight_layout()
    plt.show()
    return stress

def proc0(dec_inv_frq,dec_interp):
    from mst_ex import use_intp_sfig, return_vf
    _sf_, _ig_ = use_intp_sfig(dec_inv_frq,iopt=dec_interp,
                               iplot=False,iwgt=False)
    _sf_ = _sf_.sf.transpose(0,3,1,2).copy() * 1e6
    return _sf_, _ig_


def proc_half(model_ngr,_sf_,raw_sfs,model_sfs,
              model_vfs,raw_ehkl):
    filt = model_ngr==0
    _sf_ = _sf_.transpose(1,0,2,3).copy()
    _sf_[:,filt]=np.nan
    _sf_ = _sf_.transpose(1,0,2,3).copy()
    raw_sfs = raw_sfs.transpose(1,0,2,3).copy()
    raw_sfs[:,filt]=np.nan
    raw_sfs = raw_sfs.transpose(1,0,2,3).copy()
    model_sfs = model_sfs.transpose(1,0,2,3).copy()
    model_sfs[:,filt]=np.nan
    model_sfs = model_sfs.transpose(1,0,2,3).copy()
    model_vfs[filt]=np.nan
    model_ngr[filt]=np.nan
    raw_ehkl[filt]=np.nan

    filt = (_sf_==0) | (raw_sfs==0)
    _sf_[filt]=np.nan
    raw_sfs[filt]=np.nan
    model_sfs[filt]=np.nan
    raw_ehkl[filt[:,0,:,:]]=np.nan
    return _sf_,raw_sfs,model_sfs,model_vfs,model_ngr,raw_ehkl


def calc_stress(
        istp,
        _sf_,
        model_rs,
        model_igs,
        model_ehkls,
        tdats,
        ref_psis,wgt):
    """
    Calculate stress - the core part of ex_consistency after
    process the 'raw' data.
    """
    ## masking can be improved by being
    ## specific on phi axis -> require fix in RS.rs
    inds = []
    for ipsi in xrange(len(ref_psis)):
        if not(np.isnan(_sf_[istp,0:2,:,ipsi]).any()):
            inds.append(ipsi)
    model_rs.sf  = _sf_[istp][:,:,inds]
    model_rs.psis = ref_psis[inds]
    model_rs.npsi = len(model_rs.psis)
    model_rs.eps0 = model_igs[istp][:,inds]
    model_rs.ehkl = model_ehkls[istp][:,inds]
    model_rs.tdat = tdats[istp][:,inds]

    #-----------------------------------#
    ## find the sigma ...
    sij_wv = model_rs.dat_model.flow.sigma[:,:,istp]
    stress_wgtavg = np.array(
        [sij_wv[0,0],sij_wv[1,1],sij_wv[2,2],sij_wv[1,2],
         sij_wv[0,2],sij_wv[0,1]])
    ## find the stress
    dsa_sigma = model_rs.find_sigma(
        ivo=[0,1],init_guess=stress_wgtavg,
        weight = wgt) # None
    for i in xrange(6): print '%+7.1f'%(dsa_sigma[i]),
    for i in xrange(6): print '%+7.1f'%(dsa_sigma[i]-stress_wgtavg[i]),
    print ''
    return dsa_sigma

def ex_consistency(
        ifig=50,nxphi=3,exp_ref=[],exp_lab=[],mod_ext=None,
        mod_ref='STR_STR.OUT',sin2psimx=None,iscatter=False,
        sigma=5e-5,psimx=None,psi_nbin=1,ig_sub=True,istep=None,
        hkl=None,iplot=True,iwind=False,wdeg=2,ipsi_opt=1,
        fn_sff=None,pmargin=None,path='',sf_ext=None,ig_ext=None,
        vf_ext=None,iwgt=False,verbose=False,ilog=True,
        dec_inv_frq=1,dec_interp=1,theta_b=None,ird=1.,nfrq=None,
        fnPickle=None):
    """
    Consistency check between 'weighted average' stress and
    the stress obtained following the stress analysis method
    (SF, IG strain)

    ----------------------------------------
    ## Data visualization options
    ifig = 50
    nxphi     : display only first nxphi of results
                along phi axis
    hkl       : hkl of plane, for the sake of labels.
    path      : place holder for strain path

    ----------------------------------------
    ## Data process parameters
    #1. Main options
    sf_ext    : Overwrite stress factor
    ig_ext    : Overwrite IG strain
    vf_ext    : Overwrite grain volume fraction (nstp,nphi,npsi)
    iscatter (False) : introducing a random scattering of 'ehkl'
    sigma     : Level of standard deviation in the Gaussian distribution
                used to 'perturb' the lattice strain.

    #2. tilting angle restriction/treatments
    sin2psimx : limiting the maximum sin2psi values, with which
                stress-fitting is carried out
    psimx     : limiting the maximum psi value
    psi_nbin  : If other than 1, limit the number of data

    #3. IG strain-specific
    ig_sub    : flag whether or not subtract the IG strain

    #4. DEC related
    dec_inv_frq : Inverse Frequency of DEC measures along EVM (nstp)
    dec_interp: Interpolation method for the incomplete DEC

    #5. Misc.
    exp_ref   : experimental reference
    mod_ref   : model's weighted average flow curves are given

    ----------------------------------------
    ## misc. options
    iplot (True) : flag whether or not MPL plot is performed
    ipsi_opt 0: sin2psi
             1: sign(psi) * sin2psi
             2: psi
    istep     : If other than None, analysis is carried out
                only for the given istep
    fnPickle=None
        if given (not None), read 'data' from the given file name

    ----------------------------------------
    ## debugging options
    verbose   : False
    ilog      : False

    ----------------------------------------
    ## Diffraction condition parameters
    theta_b   : None (Bragg's angle) (given in radian)
    ird       : Intensity expected for random distribution
    nfrq        None

    Returns
    -------
    if istep is given, perform the task for that specific step number
                       and returns
        model_rs,
        s11,
        s22,
        dsa_sigma[0],
        dsa_sigma[1],
        raw_psis,
        raw_vfs[istp],
        raw_sfs[istp],
        full_Ei,
        DEC_interp[istp]

    Otherwise, return the results after perform the task for all given steps
    In this case, the returned results are dependent on the
    argument fnPickle.
    if type(fnPickle).__name__=='NoneType':
         return model_rs, flow_weight, flow_dsa, filename_pickle
    else:
         return model_rs, flow_weight, flow_dsa
    """
    np.seterr(all='ignore')
    import time, pickle
    from MP import progress_bar
    from rs import ResidualStress,\
        u_epshkl_geom_inten_vectorize,\
        filter_psi3,psi_reso4
    from mst_ex import  return_vf
    from MP.mat import mech # mech is a module
    uet = progress_bar.update_elapsed_time
    t0 = time.time()
    if ilog:
        fn = 'ex_consistency.log'
        f = open(fn,'w')
        write_args(
            f=f,ihead=True,ifig=ifig,nxphi=nxphi,exp_ref=exp_ref,
            exp_lab=exp_lab,mod_ext=mod_ext,mod_ref=mod_ref,
            sin2psimx=sin2psimx,iscatter=iscatter,sigma=sigma,
            psimx=psimx,psi_nbin=psi_nbin,ig_sub=ig_sub,
            istep=istep,hkl=hkl,iplot=iplot,iwind=iwind,
            wdeg=wdeg,ipsi_opt=ipsi_opt,fn_sff=fn_sff,
            pmargin=pmargin,path=path,sf_ext=sf_ext,
            ig_ext=ig_ext,vf_ext=vf_ext,iwgt=iwgt,
            verbose=verbose,ilog=ilog,nfrq=nfrq,
            ird=ird,theta_b=theta_b,fnPickle=fnPickle)
        f.close()
        print 'log has been saved to ',fn

    #------------------------------------------------------------#
    if not(iwgt): wgt = None # overwrite wgt

    if type(fnPickle).__name__=='NoneType':
        ## i_ip = 1: ioption for the model data
        t0_load=time.time()
        model_rs = ResidualStress(
            mod_ext=mod_ext,
            fnmod_ig='igstrain_fbulk_ph1.out',
            fnmod_sf='igstrain_fbulk_ph1.out',
            i_ip=1)
        uet(time.time()-t0_load,
            'Time spent for loading model_rs in rs_ex.ex_consistency')
        print

        ## Process the sf/eps0/ehkl/wgt and so forth
        ## according to the parameters given
        ivf_ext=True; isf_ext=True; iig_ext=True
        if type(sf_ext)==type(None):isf_ext=False
        if type(ig_ext)==type(None):iig_ext=False
        if type(vf_ext)==type(None):ivf_ext=False

        ## SF/IG/ehkl
        model_ehkls = np.copy(model_rs.dat_model.ehkl)
        if isf_ext and iig_ext:
            model_sfs = sf_ext.copy()
            model_igs = ig_ext.copy()
        elif isf_ext!=iig_ext:
            raise IOError, 'isf_ext should equal to iig_ext'
        ## if isf_ext False and isf_ext False
        else:
            model_sfs = model_rs.dat_model.sf.copy()
            model_igs = model_rs.dat_model.ig.copy()

        ## VF
        if not(ivf_ext): model_vfs, model_ngr = return_vf()
        else: model_vfs = vf_ext.copy()

        ## whether or not vf would be used as weights in
        ## the least-sq estimator

        ## Apply data-analysis
        ## 0. Use of 'interpolated' SF
        ## 1. Limit the range of sin2psi (or psi)
        ## 2. Finite number of tiltings
        ## 3. Assign tdat (subtracting ig strain or not)
        ## 4. Perturb ehkl (common)
        ## 5. Filter based on vf?

        ## Unstaged but the 'raw' arrays (not used for actual calculation)
        raw_psis = model_rs.dat_model.psi.copy()
        raw_vfs  = model_vfs.copy()
        raw_ehkl = np.copy(model_rs.dat_model.ehkl)
        raw_sfs  = model_sfs.copy()

        ## staged arrays for stress estimation (used for actual calcultion)
        model_rs.psis = model_rs.dat_model.psi.copy()
        model_rs.phis = model_rs.dat_model.phi.copy()
        model_rs.npsi = len(model_rs.psis)
        model_rs.nphi = len(model_rs.phis)
        sin2psis_init = np.sin(model_rs.psis)**2

        ## 0. Use interpolated SF
        _sf_, _ig_ = proc0(dec_inv_frq,dec_interp)

        ## 0.5 Mask DEC where volume fraction is depleted...
        ## model_ngr or model_vfs
        _sf_,raw_sfs,model_sfs,model_vfs,model_ngr,raw_ehkl\
            =proc_half(model_ngr,_sf_,raw_sfs,model_sfs,
                       model_vfs,raw_ehkl)

        ## 1. Limit the range of sin2psi (or psi)
        t0_pr1=time.time()
        if type(sin2psimx)!=type(None) or type(psimx)!=type(None):
            if type(sin2psimx)!=type(None): bounds=[0., sin2psimx]
            elif type(psimx)!=type(None): bounds=[0.,np.sin(psimx*np.pi/180.)**2]
            raw_sfs,model_sfs,model_igs,model_vfs,\
                model_ehkls,raw_psis,_sf_,raw_vfs\
                = filter_psi3(sin2psis_init,bounds,
                              raw_sfs,model_sfs,model_igs,model_vfs,
                              model_ehkls,raw_psis,_sf_,raw_vfs)

            ## reduce the psis in model_rs
            model_rs.psis, = filter_psi3(sin2psis_init,bounds,
                                         model_rs.psis.copy())
            model_rs.npsi = len(model_rs.psis)
        uet(time.time()-t0_pr1, 't for pr1');print
        DEC_interp = _sf_.copy()

        ## 2. Finite number of tiltings
        t0_pr2=time.time()
        if psi_nbin!=1:
            model_sfs, model_igs, model_vfs,model_ehkls, _sf_ \
                = psi_reso4(model_rs.psis, psi_nbin,
                            model_sfs,model_igs,
                            model_vfs,model_ehkls,_sf_)
            model_rs.psis, = psi_reso4(model_rs.psis.copy(),psi_nbin,
                                       model_rs.psis.copy())
            model_rs.npsi = len(model_rs.psis)
        uet(time.time()-t0_pr2,'t for pr2');print

        ## 3. Assign tdat
        if ig_sub: model_tdats = model_ehkls - model_igs
        else     : model_tdats = model_ehkls.copy()

        ## 3.5 - save data if required in order to reuse
        ## the processed (analyzed) data for Monte Carlo
        ## calculation
        ## list of items to be saved

        ## model_vfs, model_tdats, model_rs
        ## _sf_, model_igs, model_ehkls

        ## pickle the data
        import MP.lib.temp

        filename_pickle=MP.lib.temp.gen_tempfile(prefix='dsa-data',ext='pck')
        with open(filename_pickle,'wb') as fo:
            pickle.dump(model_vfs,fo)
            pickle.dump(model_tdats,fo)
            pickle.dump(_sf_,fo)
            pickle.dump(model_igs,fo)
            pickle.dump(model_ehkls,fo)
            pickle.dump(model_rs,fo)
            pass
        pass

    elif type(fnPickle).__name__=='str':
        with open(fnPickle,'rb') as fo:
            model_vfs = pickle.load(fo)
            model_tdats=pickle.load(fo)
            _sf_=pickle.load(fo)
            model_igs = pickle.load(fo)
            model_ehkls = pickle.load(fo)
            model_rs=pickle.load(fo)
            pass
        pass

    ## 4. Perturb tdat, i.e., ehkl (performed regardless of fnPickle)
    if iscatter:
        t0_perturb=time.time()
        nstp, nphi, npsi = model_ehkls.shape
        tdats = u_epshkl_geom_inten_vectorize(
            model_vfs,model_tdats,
            ird,sigma,model_rs.psis,theta_b)
        uet(time.time()-t0_perturb,'Time spent for perturbation')
        print
    else: tdats=model_tdats.copy()

    ## end of data process
    #------------------------------------------------------------#

    if mod_ext==None: mod_ref='STR_STR.OUT'
    else            : mod_ref='%s.%s'%(mod_ref.split('.')[0],
                                       mod_ext)

    flow_weight = mech.FlowCurve(name='Model weighted')
    flow_weight.get_model(fn=mod_ref)
    ## calc Von Mises stress/strain
    flow_weight.get_eqv()

    if len(flow_weight.epsilon_vm)<5: lc='k.'
    else:                             lc='k-'

    ## plot all stress factors at individual
    ## deformation levels

    stress = []
    print ('%19s '+'%8s%8s%8s%8s%8s%8s'*2)%(
        '','S11','S22','S33','S23','S13','S12',
        'dS11','dS22','dS33','dS23','dS13','dS12')
    ################################################
    ## *Serial* Loop over the deformation steps
    ref_psis = model_rs.psis.copy()
    nstp = model_rs.dat_model.nstp

    ## Main loop
    ## -------------------------------------------------- ##
    dtCalcstr = 0.
    for istp in xrange(nstp):
        t0_calstr=time.time()
        if type(nfrq).__name__=='int':
            if istp % nfrq !=0: continue
        if iplot==False and type(istep)!=type(None):
            ## exit by conducting only for the given step
            nstp = 1; istp = istep

        print '%13s %2.2i/%2.2i'%('processing:',istp,nstp),
        dsa_sigma = calc_stress(
            istp,_sf_,model_rs,model_igs,
            model_ehkls,tdats,
            ref_psis,wgt)
        stress.append(dsa_sigma)

        if iplot and istp==0:
            from matplotlib import pyplot as plt
            from matplotlib.backends.backend_pdf import PdfPages
            from MP.lib import mpl_lib,axes_label
            plt.ioff()
            wide_fig     = mpl_lib.wide_fig
            fancy_legend = mpl_lib.fancy_legend
            ## Collection of figures at various plastic strains
            fe   = PdfPages('all_ehkl_fits_%s_%s.pdf'%(hkl,path))
            fs   = PdfPages('all_stress_factors_%s_%s.pdf'%(hkl,path))
            f_er = PdfPages('all_Ei-ehkl-e0_%s_%s.pdf'%(hkl,path))
            fig1 = wide_fig(ifig,nw=2,nh=1,left=0.2,uw=3.5,
                            w0=0,w1=0.3,right=0,iarange=True)

            fig_vm = plt.figure(figsize=(3.5,3.5))
            ax_vm = fig_vm.add_subplot(111)

            ## ax1: Equivalent Stress/Strain
            ## ax2: Stress path in the plane stress space (RD/TD)
            ax1 = fig1.axes[0]; ax2 = fig1.axes[1]
            for a in [ax1,ax_vm]:
                a.plot(
                    flow_weight.epsilon_vm,
                    flow_weight.sigma_vm,
                    lc,label=r'$\langle \sigma^c \rangle$',
                    alpha=1.0)
                axes_label.__eqv__(a,ft=10)
        if iplot:
            a=raw_sfs[istp,:2,:,:]
            b=dsa_sigma[:2]
            full_Ei = np.tensordot(a,b,axes=[0,0])
            s11 = model_rs.dat_model.flow.sigma[0,0,istp]
            s22 = model_rs.dat_model.flow.sigma[1,1,istp]

            if type(istep)!=type(None):
                return model_rs, s11,s22, dsa_sigma[0],dsa_sigma[1],\
                    raw_psis.copy(),\
                    raw_vfs[istp].copy(),raw_sfs[istp].copy(),\
                    full_Ei.copy(),\
                    DEC_interp[istp].copy()

            #-----------------------------------#
            if istp==0: ileg=True
            else:       ileg=True #False

            if (istep!=None and istp==istep) or\
               (istep==None and istp==nstp-1):
               fig2,fig3,fig4=__model_fit_plot__(
                   model_rs,ifig=ifig+istp*2+10,
                   istp=istp, nxphi=nxphi,stress_wgt=None,
                   ivo=None,hkl=hkl,ileg=ileg,iwind=iwind,
                   wdeg=wdeg)
            else:
                plt.ioff()
                f1,f2,f3=__model_fit_plot__(
                    model_rs,ifig=ifig+istp*2+10,
                    ## stress_wgt: the mechanical stress.
                    ## will be used as a reference line
                    istp=istp,nxphi=nxphi,stress_wgt=[s11,s22,0,0,0,0],
                    wgt=raw_vfs[istp].copy(),wgt_psi=raw_psis,
                    full_Ei = full_Ei,
                    full_DEC = raw_sfs[istp].copy(),
                    DEC_interp = DEC_interp[istp].copy(),
                    ivo=[0,1],hkl=hkl,ileg=ileg,iwind=False,
                    ipsi_opt=ipsi_opt)
                fs.savefig(f2);fe.savefig(f1);f_er.savefig(f3)
                f1.clf();plt.draw();f2.clf();plt.draw();f3.clf();plt.draw()
                plt.close(f1);plt.close(f2);plt.close(f3);plt.ion()
        dtCalcstr = dtCalcstr + (time.time()-t0_calstr)
    uet(dtCalcstr,'Time spent inside the stress analysis loop');print
    # end of the serial loop over deformation steps
    ############################################################

    if iplot: fe.close(); fs.close(); f_er.close()
    stress   = np.array(stress).T # diffraction stress
    flow_dsa = mech.FlowCurve(name='Diffraction Stress')
    flow_dsa.get_6stress(stress)
    if type(nfrq).__name__=='NoneType':
        flow_dsa.get_33strain(model_rs.dat_model.flow.epsilon)
    elif type(nfrq).__name__=='int':
        flow_dsa.get_33strain(model_rs.dat_model.flow.epsilon[:,:,::nfrq])

    flow_dsa.get_eqv()
    sigma_wgt = flow_weight.sigma
    ## Various plots
    if iplot:
        for a in [ax1, ax_vm]:
            a.plot(flow_dsa.epsilon_vm,flow_dsa.sigma_vm,'k+',
                   label=r'$\sigma^{d}$')
        for i in xrange(len(exp_ref)):
            f = exp_ref[i]; lab = exp_lab[i]
            edat = np.loadtxt(f).T
            ax1.plot(edat[0],edat[1],'-',lw=2,label=lab)

        for a in [ax1,ax_vm]:
            fancy_legend(a,size=10,nscat=1)

        ax2.plot(sigma_wgt[0,0],sigma_wgt[1,1],'k-')
        ax2.plot(flow_dsa.sigma[0,0],flow_dsa.sigma[1,1],'k+')

        ## connector
        npoints = len(sigma_wgt[0,0])
        wgtx = sigma_wgt[0,0];      wgty = sigma_wgt[1,1]
        dsax = flow_dsa.sigma[0,0]; dsay = flow_dsa.sigma[1,1]

        ax2.set_ylim(-100,700); ax2.set_xlim(-100,700)
        ax2.set_aspect('equal')
        ax2.set_xlabel(r'$\bar{\Sigma}_{11}$ [MPa]',dict(fontsize=15))
        ax2.set_ylabel(r'$\bar{\Sigma}_{22}$ [MPa]',dict(fontsize=15))
        ax2.locator_params(nbins=3)
        ax2.set_xticks(np.linspace(300,700,3),dict(fontsize=4))
        ax2.set_yticks(np.linspace(300,700,3),dict(fontsize=4))
        ax2.grid('on'); plt.show()

        ## save figures
        fig1.savefig('flow_%s_%s.pdf'%(hkl,path))
        fig_vm.savefig('all_flow_vm.pdf', bbox_inches='tight')
        fig_vm.savefig('all_flow_vm.ps', bbox_inches='tight')
        # close figures
        try:
            figs=[fig1,fig2,fig3,fig4]
            for f in figs:
                try: plt.close(f)
                except: pass
        except: pass

    uet(time.time()-t0);print
    if type(fnPickle).__name__=='NoneType':
        return model_rs, flow_weight, flow_dsa, filename_pickle
    else:
        return model_rs, flow_weight, flow_dsa

def __model_fit_plot__(
        container,ifig,istp,nxphi=None,hkl=None,
        wgt=None,wgt_psi=None,full_Ei=None,
        full_DEC = None,
        DEC_interp=None,
        stress_wgt=None,ivo=None,fig=None,figs=None,fige=None,
        c1='r',c2='b',m1='--',m2='-',isf=[True,True],
        ileg=True,iwind=False,wdeg=2,
        ipsi_opt=0):
    """
    Plot container's analyzed data
    """
    from MP.lib import mpl_lib, axes_label
    import lib; sin2psi_wind = lib.sin2psi_bounds
    sin2psi_opt = lib.sin2psi_opt
    wide_fig     = mpl_lib.wide_fig
    fancy_legend = mpl_lib.fancy_legend
    rm_lab       = mpl_lib.rm_lab
    tune_x_lim   = mpl_lib.tune_x_lim
    tune_xy_lim  = mpl_lib.tune_xy_lim
    deco         = axes_label.__deco__

    nphi = container.nphi; npsi = container.npsi
    phis = container.phis; psis = container.psis
    if nxphi!=None and nphi>nxphi:
        print 'Probed phis are redundant:', nphi
        print 'Only %i of phis axis are shown'%(nxphi)
        nphi = nxphi

    ivf=True
    #vf   = container.dat_model.vf[istp]
    #ngr  = container.dat_model.ngr[istp]
    if type(wgt)!=type(None): vf   = wgt.copy()
    else: vf =[np.nan]

    if np.all(np.isnan(vf)):
        print 'All volume fractions are nan'
        ivf=False

    iEi = False
    if type(full_Ei)!=type(None): iEi = True

    sf   = container.sf #dat_model.sf[istp]
    tdat = container.tdat; ehkl = container.ehkl
    eps0 = container.eps0; Ei   = container.Ei

    if fig==None: fig = wide_fig(
            ifig,nw=nphi,w0=0.00,ws=0.5,w1=0.0,uw=3.0,
            left=0.15,right=0.25,nh=1,h0=0.2,h1=0,
            down=0.08,up=0.10,iarange=True)
    if figs==None: figs= wide_fig(
            ifig+1,nw=nphi,w0=0.00,ws=0.5,w1=0.0,uw=3.0,
            left=0.15,right=0.25,nh=1,h0=0.2,h1=0,
            down=0.08,up=0.10,iarange=True)
    if fige==None: fige= wide_fig(
            ifig+2,nw=nphi,w0=0.00,ws=0.5,w1=0.0,
            uw=3.0,left=0.12,right=0.10)

    axes  = fig.axes[:nphi]#nphi:nphi*2]
    ax_er = fige.axes[:nphi]; axesf = figs.axes;
    if ivf: axesv = []

    for iphi in xrange(nphi):
        ax = fig.axes[iphi]; axs = figs.axes[iphi]
        ax.set_title( r'$\phi=%3.1f^\circ$'%(phis[iphi]*180/pi))
        axs.set_title(r'$\phi=%3.1f^\circ$'%(phis[iphi]*180/pi))
        ax.locator_params(nbins=4); axs.locator_params(nbins=4)
        if ivf:
            axesv.append(axes[iphi].twinx())
            axesv[iphi].locator_params(nbins=4)

    for iphi in xrange(nphi):
        ax = axes[iphi]; ae = ax_er[iphi]
        if ivf: av = axesv[iphi];
        ## convert psi to user's convenience.
        if iphi==0:
            if iEi: xEi = sin2psi_opt(wgt_psi.copy(),ipsi_opt)
            x  = sin2psi_opt(psis.copy(),ipsi_opt)
            if ivf: xv = sin2psi_opt(wgt_psi.copy(),ipsi_opt)
            if iwind:
                x = sin2psi_opt(psis.copy(),2)
                if ivf: xv = sin2psi_opt(wgt_psi.copy(),2)

        ## E_{i} -- abundant information
        # if ileg: label=r'Fit'
        # else:    label=None
        # ax.plot(x,Ei[iphi]*1e6,'+',mec='k',ms=7,
        #         mfc='None',label=label)
        ## e-e_0
        if type(hkl).__name__=='NoneType' and ileg:
            label=r'$\langle\tilde{\varepsilon}^d \rangle^G$'
        elif type(hkl).__name__!='NoneType' and ileg:
            label=r'$\langle\tilde{\varepsilon}^d \rangle^G$'
        elif ileg!=True: label=None

        y = tdat[iphi]*1e6
        ax.plot(x,y,'k.',label=label)
        if iwind:
            raise IOError, 'iwind is not stable'
            for i in xrange(len(psis)):
                X = psis[i]*180./np.pi; Y = y[i]
                pl, lu, s2l, s2u = sin2psi_wind(
                    w_rad=wdeg,psi0=X)
                ax.plot([s2l,s2u],[Y,Y],'g-')
                ax.plot([s2l,s2u],[Y,Y],'g|')

        if type(hkl)==type(None) and ileg: label=r'$E_{i} - \varepsilon^{\{hkl\}}-\varepsilon^{\{hkl\}}_0$'
        elif type(hkl)!=type(None) and ileg: label='$E_{i} - \varepsilon^{\{%s\}}-\varepsilon^{\{%s\}}_0$'%(hkl,hkl)
        elif ileg!=True: label = None

        if ivf: av.plot(xv,vf[iphi],'k-',lw=2,
                        alpha=0.4)

        ae.plot(x,Ei[iphi]*1e6-y,c2+m2,label=label)
        deco(ax=ax,iopt=0,hkl=hkl,ipsi_opt=ipsi_opt)
        deco(ax=ae,iopt=0,hkl=hkl,ipsi_opt=ipsi_opt)

        ## all_stress_factor_hkl.pdf
        ax=axesf[iphi]
        if type(hkl)==type(None) and ileg:
            lab1=r'$\mathbb{F}_{11}^{I}$ in use'
            lab2=r'$\mathbb{F}_{22}^{I}$ in use'
        elif type(hkl)!=type(None) and ileg:
            lab1=r'$\mathbb{F}^{\{%s\},I}_{11}$ in use'%hkl
            lab2=r'$\mathbb{F}^{\{%s\},I}_{22}$ in use'%hkl
        elif ileg!=True:
            lab1=None; lab2=None

        for i in xrange(len(isf)):
            if i==0:# and ileg:
                lab=lab1
                #st = c1+m1
                # c=c1
                c='k'
                st='r-'
                marker = '.'
            elif i==1:# and ileg:
                lab=lab2
                #st = c2+m2
                # c=c2
                c = 'gray'
                st='b-'
                marker = '+'

            if isf[i]:
                # l, = ax.plot(
                #     xv,sf[i][iphi]*1e6,st,label=lab)

                _l1_=None
                _l2_=None
                if i==0 and iphi==nphi-1:
                    _l1_=r'$\mathbb{F}_{11}$'
                    _l2_=r'$\mathbb{F}^{\ I}_{11}$'
                    ls='-'
                if i==1 and iphi==nphi-1:
                    _l1_=r'$\mathbb{F}_{22}$'
                    _l2_=r'$\mathbb{F}^{\ I}_{22}$'
                    ls='--'

                if type(full_DEC)!=type(None):
                    if (full_DEC[i][iphi]==0).any():
                        raise IOError,'Found zero in full_DEC'
                    ax.plot(xv,full_DEC[i][iphi]*1e6,'-',
                            label=_l1_,color=c)
                if type(DEC_interp)!=type(None):
                    ax.plot(xv,DEC_interp[i][iphi]*1e6,'k--',
                            label=_l2_,color=c)

                for j in xrange(len(sf[i][iphi][:])):
                    if sf[i][iphi][j]!=0:
                        if j==0:
                            ax.plot(
                                x[j],sf[i][iphi][j]*1e6,
                                ls='None',
                                color=c,marker=marker,
                                label=lab)
                        else:
                            ax.plot(
                                x[j],sf[i][iphi][j]*1e6,
                                color=c,marker=marker)
                ## legend
                if iphi==nphi-1 and i==1 and istp==0:
                    fancy_legend(ax,size=11,nscat=1,ncol=1,
                                 bbox_to_anchor=(1.4,1))

        if ivf:
            av.set_ylabel(r'Vol. $f(\phi,\psi)$',
                          dict(fontsize=13))
            # av.tick_params(axis='y',colors='red')
            # av.yaxis.label.set_color('red')

        deco(ax=ax,iopt=1,hkl=hkl,ipsi_opt=ipsi_opt)
        # if iphi==0:fancy_legend(ax,nscat=1)

    if type(stress_wgt)!=type(None):
        container.sigma = np.array(stress_wgt)
        container.calc_Ei(ivo=ivo)

        if ileg:
            label = r'$\mathbb{F}\ : \langle \mathbf{\sigma}^{c} \rangle$'
            lab1 = r'$\mathbb{F}\ : \mathbf{\sigma}^{d}$'
        else:
            label=None
            lab1 =None

        for iphi in xrange(nphi):
            ax=axes[iphi]

            ## Calculate the elastic strain
            ## based on the obtained stress and DEC
            ax.plot(x,container.Ei[iphi]*1e6,'o',
                    mfc='None',mec='black',label=label)
            if iEi: ax.plot(xEi, full_Ei[iphi]*1e6,'k-',label=lab1)

            if iphi==nphi-1:
                if ivf: ax.plot(
                        0,0,'k-',lw=2,alpha=0.4,
                        label='Vol.')## to trick the legend
                fancy_legend(
                    ax,size=11,nscat=1,ncol=1,
                    bbox_to_anchor=(1.4,1))


    for iax in xrange(len(axes)):
        axes[iax].set_ylim(-1500,)

    for iax in xrange(len(axes)):
        axes[iax].set_ylim(-1500,)
        if ivf:
            axesv[iax].set_ylim(0,0.30)

    tune_x_lim(fig.axes,axis='x')
    tune_x_lim(axes,    axis='y')
    tune_xy_lim(ax_er           )
    if ivf: tune_x_lim(axesv,   axis='y')

    tune_x_lim(axesf,   axis='y')

    ## remove redundant axis labels
    for i in xrange(len(axes)-1):
        rm_lab(axes[i+1], axis='y')
        rm_lab(axes[i+1], axis='x')
        rm_lab(axesf[i+1],axis='y')
        rm_lab(axesf[i+1],axis='x')
        rm_lab(ax_er[i+1],axis='y')
        rm_lab(ax_er[i+1],axis='x')
        if ivf: rm_lab(axesv[i],  axis='y')

    return fig, figs, fige

def __model_fit_plot_3d__(container,istp,nxphi=None):
    import phikhi.psikhi2cart as pk
    conv = pk.conv
    convs = pk.convs

    nphi=container.nphi
    if nxphi!=None and nphi>nxphi:
        print 'Probed phis are redundant:', nphi
        print 'Only %i of phis axis are shown'%(nxphi)
        nphi = nxphi

    npsi=container.npsi
    phis=container.phis
    psis=container.psis

    vf = container.dat_model.vf[istp]
    ngr = container.dat_model.ngr[istp]

    sf = container.sf
    tdat = container.tdat
    ehkl = container.ehkl
    eps0 = container.eps0
    Ei   = container.Ei

    from mpl_lib import mpl_lib
    ax3d = mpl_lib.axes3()

    for iphi in xrange(nphi):
        #x = sin(psis)**2
        x = psis
        y = np.ones(len(x))*phis[iphi]
        x, y = convs(k=x,p=y)
        z = Ei[iphi]
        ax3d.plot(x,y,z,'rx')

    for iphi in xrange(nphi):
        x = psis
        y = np.ones(len(x)) * phis[iphi]
        x, y = convs(k=x,p=y)
        z = tdat[iphi]
        ax3d.plot(x,y,z,'b+')

    #ax3d.set_xlabel(r'$\sin^2{\psi}$')
    ax3d.set_xlabel(r'$\psi$')
    ax3d.set_ylabel(r'$\phi$')
    ax3d.set_xlim(-1,1)
    ax3d.set_ylim(-1,1)



def sf_scan(fn='sf_ph1.out',npair=1):
    import MP.read_blocks as rb
    from MP.mat import voigt
    from MP.lib import mpl_lib

    wide_fig = mpl_lib.wide_fig
    figs=wide_fig(nw=2)
    ijv = voigt.ijv
    read=rb.read_sf_scan
    sig,eps=read(fn=fn)
    print len(sig),len(eps)
    nstp = len(sig)
    if np.mod(nstp,npair)!=0:
        raise IOError, 'Cannot be paired with given npair...'

    for i in xrange(nstp/npair):
        for j in xrange(npair):
            s = sig[i*npair+j]
            e = eps[i*npair+j]
            for k in xrange(2):
                figs.axes[k].plot(s[k], e[10])
                figs.axes[k].plot(s[k], e[1])
                figs.axes[k].plot(s[k], e[100])
                figs.axes[k].plot(s[k], e[200])


def stress_plot(ext=['311','211','220','200'],ifig=1,nphi=3,istp=2,sin2psimx=None):
    """
    Plot results on various hkl
    """
    from MP.lib import mpl_lib
    from MP.lib import axes_label
    wide_fig = mpl_lib.wide_fig
    fig  = wide_fig(ifig,nw=nphi,w0=0.00,ws=0.5,w1=0.0,uw=3.0,
                    left=0.15,right=0.10,
                    nh=1,h0=0.2,h1=0,down=0.08,up=0.10,
                    iarange=True)
    figs = wide_fig(ifig+1,nw=nphi,w0=0.00,ws=0.5,w1=0.0,uw=3.0,
                    left=0.12,right=0.10)
    fige = wide_fig(ifig+2,nw=nphi,w0=0.00,ws=0.5,w1=0.0,uw=3.0,
                    left=0.12,right=0.10)
    hkl=ext


    c1=['r','g','b','k']
    c2=['m','y','b','r']

    m1=['x','+','-','--']
    m2=['+','+','+','+']

    diffs=[]
    f_w=[]
    f_d=[]

    for i in xrange(len(ext)):
        # weighted flow, flow based on diffraction analysis
        rst = ex_consistency(iplot=False,mod_ext=ext[i],
                             sin2psimx=sin2psimx)
        diff, flow_w, flow_d = rst[0],rst[1],rst[2]
        diffs.append(diffs)
        f_w.append(flow_w)
        f_d.append(flow_d)

        print ext[i]
        f1,f2,f3=__model_fit_plot__(
            container=diff,ifig=ifig,istp=istp,nxphi=3,hkl=ext[i],
            stress_wgt=None,ivo=None,fig=fig,figs=figs,fige=fige,
            c1=c1[i],c2=c2[i],
            m1=m1[i],m2=m2[i],
            isf=[True,False]
            )

    return diffs, f_w, f_d
