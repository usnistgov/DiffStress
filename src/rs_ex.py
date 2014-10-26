"""
A collection of scripts that make use of residual stress module
"""

import numpy as np
from numpy import pi, sin, cos

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
    for i in range(len(strain_eff)):
        dum = myrs.analysis(iopt=0,istp=i)
        stress.append((dum[0] + dum[1])/2.)
    ax.plot(strain_eff,stress,'-x',label='Model SF/IG/ehkl')

    stress = []
    for i in range(len(strain_eff)):
        dum = myrs.analysis(iopt=1,istp=i)
        stress.append((dum[0] + dum[1])/2.)
    ax.plot(strain_eff,stress,'-x',label='Exp SF/IG/ehkl')

    stress = []
    for i in range(len(strain_eff)):
        dum = myrs.analysis(iopt=2,istp=i)
        stress.append((dum[0] + dum[1])/2.)
    ax.plot(strain_eff,stress,'-x',label='Model SF/IG + Exp ehkl')


    for i in range(len(exp_ref)):
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

def ex_consistency(
        ifig=50,nxphi=3,
        exp_ref=[],exp_lab=[],mod_ext=None,
        mod_ref='STR_STR.OUT',sin2psimx=None,
        iscatter=False,sigma=5e-5,psimx=None,psi_nbin=1,
        ig_sub=True,istep=None,hkl=None,iplot=True,
        iwind=False,wdeg=2,ipsi_opt=1,fn_sff=None,
        pmargin=None,path='',
        sf_ext=None,ig_ext=None,iwgt=False):
    """
    Consistency check between 'weighted average' stress and
    the stress obtained following the stress analysis method
    (SF, IG strain)

    ifig = 50
    nxphi     : display only first nxphi of results
                along phi axis
    exp_ref   : experimental reference
    mode_ref  : model's weighted average flow curves are given
    sin2psimx : limiting the maximum sin2psi values, with which
                stress-fitting is carried out
    iscatter (False) : introducing a random scattering of 'ehkl'
    psimx     : limiting the maximum psi value
    psi_nbin  : If other than 1, limit the number of data
    ig_sub    : flag whether or not subtract the IG strain
    istep     : If other than None, analysis is carried out
                for a specific istep
    hkl       : hkl of plane, for the sake of labels.
    iplot (True) : flag whether or not MPL plot is performed
    ipsi_opt 0: sin2psi
             1: sing(psi) * sin2psi
             2: psi
    pmargin   : portional margin of volume that should exceed to
                contribute to the ehkl/SF/IG in model_rs
    path      : place holder for strain path
    sf_ext    : Overwrite stress factor
    ig_ext    : Overwrite IG strain
    iwgt      : Whether or not accounting for 'weight'
    """
    from rs import ResidualStress,u_epshkl,filter_psi,\
        filter_psi2,psi_reso, psi_reso2, psi_reso3

    from MP.mat import mech # mech is a module
    FlowCurve = mech.FlowCurve

    if iplot:
        from matplotlib import pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
        from MP.lib import mpl_lib,axes_label
        wide_fig     = mpl_lib.wide_fig
        fancy_legend = mpl_lib.fancy_legend

        ## Collection of figures at various plastic strains
        fe   = PdfPages('all_ehkl_fits_%s_%s.pdf'%(hkl,path))
        fs   = PdfPages('all_stress_factors_%s_%s.pdf'%(hkl,path))
        f_er = PdfPages('all_Ei-ehkl-e0_%s_%s.pdf'%(hkl,path))

        # fig1 is 'flow_**.pdf'
        fig1 = wide_fig(ifig,nw=2,nh=1,left=0.2,uw=3.5,
                        w0=0,w1=0.3,right=0,iarange=True)

        ## ax1: Equivalent Stress/Strain
        ## ax2: Stress path in the plane stress space (RD/TD)
        ax1 = fig1.axes[0]; ax2 = fig1.axes[1]
        ax2.set_axis_bgcolor('0.95')

    ## i_ip = 1: model case
    model_rs = ResidualStress(
        mod_ext=mod_ext,
        fnmod_ig='igstrain_fbulk_ph1.out',
        fnmod_sf='igstrain_fbulk_ph1.out',
        i_ip=1)

    ## masking array element based on diffraction volume
    model_rs.dat_model.mask_vol()
    if pmargin!=None:
        model_rs.dat_model.mask_vol_margin(pmargin)

    if mod_ext==None: mod_ref='STR_STR.OUT'
    else:             mod_ref='%s.%s'%(
            mod_ref.split('.')[0],mod_ext)

    flow_weight = FlowCurve(name='Model weighted')
    flow_weight.get_model(fn=mod_ref)
    flow_weight.get_eqv() ## calc Von Mises stress/strain

    if len(flow_weight.epsilon_vm)<5: lc='bx'
    else:                             lc='bx-'

    if iplot:
        ax1.plot(flow_weight.epsilon_vm,flow_weight.sigma_vm,
                 lc,label='Weighted avg',alpha=1.0)
        axes_label.__eqv__(ax1,ft=10)

    ## plot all stress factors at individual deformation levels
    stress = []
    print '%8s%8s%8s%8s%8s%8s'%(
        'S11','S22','S33','S23','S13','S12')

    for istp in range(model_rs.dat_model.nstp):
        """
        Dimensions of data arrays for:
        ==============================
        sf     (nstp, k, nphi, npsi)
        ig     (nstp, nphi, npsi)
        ehkl   (nstp, nphi, npsi)
        strain (nstp, 6)
        vf     (nstp, nphi, npsi)
        """
        model_rs.sf   = model_rs.dat_model.sf[istp]
        model_rs.eps0 = model_rs.dat_model.ig[istp]
        model_rs.ehkl = model_rs.dat_model.ehkl[istp]

        ## whether or not intergranular strain is subtracted.
        if ig_sub: model_rs.tdat = model_rs.ehkl - model_rs.eps0
        else:      model_rs.tdat = model_rs.ehkl
        tdat_ref = model_rs.tdat[::]

        ## Inducing 'scatters'
        if iscatter:
            tdat_scatter = []
            for iphi in range(len(tdat_ref)):
                dum = u_epshkl(tdat_ref[iphi],sigma=sigma)
                tdat_scatter.append(dum)
            tdat_scatter = np.array(tdat_scatter)
            model_rs.tdat = tdat_scatter

        model_rs.phis = model_rs.dat_model.phi
        model_rs.psis = model_rs.dat_model.psi
        model_rs.nphi = len(model_rs.phis)
        model_rs.npsi = len(model_rs.psis)
        wgt           = model_rs.dat_model.vf[istp][::]

        if sin2psimx!=None or psimx!=None:
            filter_psi(model_rs,sin2psimx=sin2psimx,psimx=psimx)
            wgt = filter_psi2(
                wgt,sin2psi=np.sin(model_rs.psis)**2,
                bounds =[0., sin2psimx])
            if sf_ext!=None:
                model_rs.sf = sf_ext[istp]
            elif ig_ext!=None:
                model_rs.ig = ig_ext[istp]

        if psi_nbin!=1:
            wgt = psi_reso3(wgt,psi=model_rs.psis,ntot=psi_nbin)
            psi_reso2(model_rs,ntot=psi_nbin)

        #-----------------------------------#
        ## find the sigma ...
        s11 = model_rs.dat_model.flow.sigma[0,0][istp]
        s22 = model_rs.dat_model.flow.sigma[1,1][istp]

        if iwgt: pass
        else: wgt = None # overwrite wgt

        dsa_sigma = model_rs.find_sigma(
            ivo=[0,1],
            init_guess=[0,0,0,0,0,0],
            #init_guess=[s11,s22,0,0,0,0],
            weight = wgt # None
            )

        for i in range(6): print '%+7.1f'%(dsa_sigma[i]),
        print ''
        stress.append(dsa_sigma)
        #-----------------------------------#
        if istp==0: ileg=True
        else:       ileg=False

        if (istep!=None and istp==istep) or\
           (istep==None and istp==model_rs.dat_model.nstp-1)\
           and iplot:
            fig2,fig3,fig4=__model_fit_plot__(
                model_rs,ifig=ifig+istp*2+10,
                istp=istp, nxphi=nxphi, stress_wgt=None,
                ivo=None,hkl=hkl,ileg=ileg,iwind=iwind,
                wdeg=wdeg)
        elif iplot:
            plt.ioff()
            f1,f2,f3=__model_fit_plot__(
                model_rs,ifig=ifig+istp*2+10,istp=istp,
                nxphi=nxphi,stress_wgt=[s11,s22,0,0,0,0],
                ivo=[0,1],hkl=hkl,ileg=ileg,iwind=False,
                ipsi_opt=ipsi_opt)
            fs.savefig(f2);fe.savefig(f1);f_er.savefig(f3)
            plt.close(f1);plt.close(f2);plt.close(f3)
            plt.ion()

    if iplot: fe.close(); fs.close(); f_er.close()

    stress   = np.array(stress).T # diffraction stress
    flow_dsa = FlowCurve(name='Diffraction Stress')
    flow_dsa.get_6stress(stress)
    flow_dsa.get_33strain(model_rs.dat_model.flow.epsilon)
    flow_dsa.get_eqv()
    if iplot:
        ax1.plot(flow_dsa.epsilon_vm,flow_dsa.sigma_vm,'k+',
                 label='Stress Analysis')
        for i in range(len(exp_ref)):
            f = exp_ref[i]; lab = exp_lab[i]
            edat = np.loadtxt(f).T
            ax1.plot(edat[0],edat[1],'-',lw=2,label=lab)
            ## ax1.set_ylim(0.,800)
        fancy_legend(ax1,size=10)

    sigma_wgt = flow_weight.sigma

    if iplot:
        ax2.plot(sigma_wgt[0,0],sigma_wgt[1,1],'bx')
        ax2.plot(flow_dsa.sigma[0,0],flow_dsa.sigma[1,1],'k+')

        ## connector
        npoints = len(sigma_wgt[0,0])
        wgtx = sigma_wgt[0,0];      wgty = sigma_wgt[1,1]
        dsax = flow_dsa.sigma[0,0]; dsay = flow_dsa.sigma[1,1]
        for i in range(npoints):
            ax2.plot([wgtx[i],dsax[i]],[wgty[i],dsay[i]],'k-',alpha=0.2)

        ax2.set_ylim(-100,700); ax2.set_xlim(-100,700)
        ax2.set_aspect('equal')
        ax2.set_xlabel(r'$\bar{\Sigma}_{11}$',dict(fontsize=15))
        ax2.set_ylabel(r'$\bar{\Sigma}_{22}$',dict(fontsize=15))
        ax2.locator_params(nbins=3)
        ax2.set_xticks(np.linspace(300,700,3),dict(fontsize=4))
        ax2.set_yticks(np.linspace(300,700,3),dict(fontsize=4))
        ax2.grid('on'); plt.show()

        ## save figures
        fig1.savefig('flow_%s_%s.pdf'%(hkl,path))
        ## fig2.savefig('ehkl_%s_fit_%s.pdf'%(hkl,path))
        ## fig3.savefig('sf_%s_%s.pdf'%(hkl,path))
        ## fig4.savefig('ehkl_fit_err_%s_%s.pdf'%(hkl,path))
        # close figures
        plt.close(fig1); plt.close(fig2); plt.close(fig3); plt.close(fig4)

    return model_rs, flow_weight, flow_dsa

def __model_fit_plot__(container,ifig,istp,nxphi=None,hkl=None,
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

    vf   = container.dat_model.vf[istp]
    ngr  = container.dat_model.ngr[istp]
    sf   = container.dat_model.sf[istp]
    tdat = container.tdat; ehkl = container.ehkl
    eps0 = container.eps0; Ei   = container.Ei

    if fig==None: fig = wide_fig(
            ifig,nw=nphi,w0=0.00,ws=0.5,w1=0.0,uw=3.0,
            left=0.15,right=0.10,nh=1,h0=0.2,h1=0,
            down=0.08,up=0.10,iarange=True)
    if figs==None: figs= wide_fig(
            ifig+1,nw=nphi,w0=0.00,ws=0.5,w1=0.0,
            uw=3.0,left=0.12,right=0.10)
    if fige==None: fige= wide_fig(
            ifig+2,nw=nphi,w0=0.00,ws=0.5,w1=0.0,
            uw=3.0,left=0.12,right=0.10)

    axes  = fig.axes[:nphi]#nphi:nphi*2]
    ax_er = fige.axes[:nphi]; axesf = figs.axes; axesv = []

    for iphi in range(nphi):
        ax = fig.axes[iphi]; axs = figs.axes[iphi]
        axesv.append(axes[iphi].twinx())
        ax.set_title( r'$\phi=%3.1f^\circ$'%(phis[iphi]*180/pi))
        axs.set_title(r'$\phi=%3.1f^\circ$'%(phis[iphi]*180/pi))
        ax.locator_params(nbins=4); axs.locator_params(nbins=4)
        axesv[iphi].locator_params(nbins=4)

    for iphi in range(nphi):
        ax = axes[iphi]; av = axesv[iphi]; ae = ax_er[iphi]

        x  = sin2psi_opt(psis,ipsi_opt)
        xv = sin2psi_opt(container.dat_model.psi,ipsi_opt)
        if iwind:
            x = sin2psi_opt(psis,2)
            xv = sin2psi_opt(container.dat_model.psi,2)

        ## E_{i}
        if ileg: label=r'$E_{i}$ (fitting)'
        else:    label=None
        ax.plot(x,Ei[iphi]*1e6,'o',mec='g',mfc='None',label=label)

        ## e-e_0
        if hkl==None and ileg:
            label=r'$\mathrm{\varepsilon^{hkl}-\varepsilon^{hkl}_0}$'
        elif hkl!=None and ileg:
            label=r'$\varepsilon^{%s}-\varepsilon^{%s}_0$'%(hkl,hkl)
        elif ileg!=True: label=None

        y = tdat[iphi]*1e6
        ax.plot(x,y,'bx',label=label)
        if iwind:
            xerr = []
            for i in range(len(psis)):
                X = psis[i]*180./np.pi; Y = y[i]
                pl, lu, s2l, s2u = sin2psi_wind(
                    w_rad=wdeg,psi0=X)
                ax.plot([s2l,s2u],[Y,Y],'b-')
                ax.plot([s2l,s2u],[Y,Y],'b|')

        if hkl==None and ileg: label=r'$E_{i} - \varepsilon^{hkl}-\varepsilon^{hkl}_0$'
        elif hkl!=None and ileg: label=r'$E_{i} - \varepsilon^{%s}-\varepsilon^{%s}_0$'%(hkl,hkl)
        elif ileg!=True: label = None
        av.plot(xv,vf[iphi],'r-')
        ae.plot(x,Ei[iphi]*1e6-y,c2+m2,label=label)

        deco(ax=ax,iopt=0,hkl=hkl,ipsi_opt=ipsi_opt)
        deco(ax=ae,iopt=0,hkl=hkl,ipsi_opt=ipsi_opt)
        if iphi==0 and ileg:
            ax.legend(loc='upper right',fontsize=9,fancybox=True).\
                get_frame().set_alpha(0.5)
            # av.legend(loc='lower right',fontsize=9,fancybox=True).\
            #     get_frame().set_alpha(0.5)

        ## all_stress_factor_hkl.pdf
        ax=axesf[iphi]
        if hkl==None and ileg:
            lab1=r'$F_{11}$'; lab2=r'$F_{22}$'
        elif hkl!=None and ileg:
            lab1=r'$F^{%s}_{11}$'%hkl; lab2=r'$F^{%s}_{22}$'%hkl
        elif ileg!=True:
            lab1=None; lab2=None

        for i in range(len(isf)):
            if i==0:# and ileg:
                lab=lab1
                #st = c1+m1
                c=c1
                st='r-'
            elif i==1:# and ileg:
                lab=lab2
                #st = c2+m2
                c=c2
                st='b-'

            if isf[i]:
                l, = ax.plot(
                    xv,sf[i][iphi]*1e6,st,label=lab)
                ax.plot(xv,sf[i][iphi]*1e6,color=c,marker='.')

        av.set_ylabel(r'Vol. $f(\phi,\psi)$',dict(fontsize=13))
        av.tick_params(axis='y',colors='red')
        av.yaxis.label.set_color('red')
        deco(ax=ax,iopt=1,hkl=hkl,ipsi_opt=ipsi_opt)
        if iphi==0:fancy_legend(ax)

    if stress_wgt!=None:
        container.sigma = np.array(stress_wgt)
        container.calc_Ei(ivo=ivo)

        # label = 'E_{i} \mathrm{(fitting) with}'
        # for i in range(len(ivo)):
        #     label = '%s \sigma_{%1i}: %3.1f'%(label,i+1,stress_wgt[i])
        # label=r'$%s$'%label
        if ileg: label = r'$E_{i}$ with given stress'
        else: label=None

        for iphi in range(nphi):
            ax=axes[iphi]
            ##x=sin(psis)**2
            ax.plot(x,container.Ei[iphi]*1e6,'k--',label=label)
            if iphi==0: fancy_legend(ax,size=10)

    tune_x_lim(fig.axes,axis='x')
    tune_x_lim(axes,    axis='y')
    tune_xy_lim(ax_er           )
    tune_x_lim(axesv,   axis='y')
    tune_x_lim(axesf,   axis='y')

    ## remove redundant axis labels
    for i in range(len(axes)-1):
        rm_lab(axes[i+1], axis='y')
        rm_lab(axes[i+1], axis='x')
        rm_lab(axesf[i+1],axis='y')
        rm_lab(axesf[i+1],axis='x')
        rm_lab(ax_er[i+1],axis='y')
        rm_lab(ax_er[i+1],axis='x')
        rm_lab(axesv[i],  axis='y')

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

    for iphi in range(nphi):
        #x = sin(psis)**2
        x = psis
        y = np.ones(len(x))*phis[iphi]
        x, y = convs(k=x,p=y)
        z = Ei[iphi]
        ax3d.plot(x,y,z,'rx')

    for iphi in range(nphi):
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

    for i in range(nstp/npair):
        for j in range(npair):
            s = sig[i*npair+j]
            e = eps[i*npair+j]
            for k in range(2):
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

    for i in range(len(ext)):
        # weighted flow, flow based on diffraction analysis
        diff, flow_w, flow_d = ex_consistency(iplot=False,mod_ext=ext[i],
                                              sin2psimx=sin2psimx)
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
