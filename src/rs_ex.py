"""
A collection of scripts that make use of residual stress module
"""

import numpy as np
pi = np.pi
sin = np.sin
cos = np.cos

def ex(ifig=50,
       exp_ref=['exp_dat/Bsteel/EXP_BULGE_JINKIM.txt',
                'exp_dat/Bsteel/uni/avgstr_000.txt'],
       exp_lab=['Exp bulge','Exp uniaxiai'],
       mod_ref='STR_STR.OUT'):
    """
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
        ifig=50,
        nxphi=3,
        #exp_ref=['exp_dat/Bsteel/bulge/EXP_BULGE_JINKIM.txt',
        #'exp_dat/Bsteel/uni/avgstr_000.txt'],
        #exp_lab=['Exp bulge','Exp uniaxiai'],
        #exp_ref=['exp_dat/Bsteel/uni/avgstr_000.txt'],
        #exp_lab=['Exp uniaxiai'],
        exp_ref=[], exp_lab=[],

        mod_ref='STR_STR.OUT',
        sin2psimx=None,
        iscatter=False,
        psimx=None,
        psi_nbin=1,ig_sub=True,
        istep=None):
    """
    Consistency check between 'weighted average' stress and the stress obtained
    following the stress analysis method (SF, IG strain)

    ifig = 50
    nxphi   : display only first nxphi of results along phi axis
    exp_ref : experimental reference
    mode_ref: model's weighted average flow curves are given
    """
    from rs import ResidualStress
    from rs import u_epshkl
    from rs import filter_psi
    from rs import psi_reso, psi_reso2
    from matplotlib import pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    fe=PdfPages('all_ehkl_fits.pdf')
    fs=PdfPages('all_stress_factors.pdf')
    from MP.lib import mpl_lib # mpl_lib is a module
    from MP.lib import axes_label
    from MP.mat import mech # mech is a module

    wide_fig = mpl_lib.wide_fig
    fancy_legend = mpl_lib.fancy_legend

    FlowCurve = mech.FlowCurve

    fig1 = wide_fig(ifig,nw=3,nh=1,left=0.2,
                    uw=3.5,w0=0,w1=0.3,right=0,iarange=True)

    ax1 = fig1.axes[0]
    ax2 = fig1.axes[1]
    ax3 = fig1.axes[2]

    ax2.set_axis_bgcolor('0.95')
    model_rs = ResidualStress(i_ip=1)

    ## masking array element based on diffraction volume
    model_rs.dat_model.mask_vol()

    flow_weight = FlowCurve(name='Model weighted')
    flow_weight.get_model(fn=mod_ref)
    #flow_weight = model_rs.dat_model.flow
    flow_weight.get_eqv()

    ax1.plot(flow_weight.epsilon_vm,flow_weight.sigma_vm,
             'b--',label='Average',alpha=1.0)
    axes_label.__eqv__(ax1,ft=10)

    # stress3x3 = model_rs.dat_model.flow.sigma
    # strain3x3 = model_rs.dat_model.flow.epsilon
    # x_weight = strain3x3[0,0] + strain3x3[1,1]
    # y_weight = (stress3x3[0,0]+stress3x3[1,1])/2.
    # ax1.plot(x_weight,y_weight,
    #          'bx',label='Average',alpha=1.0)

    # ax1.set_xlabel(r'$\bar{E}^{\mathrm{eq}}$',dict(fontsize=12))
    # ax1.set_ylabel(r'$\bar{\Sigma}^{\mathrm{eq}}}$',dict(fontsize=12))

    ## plot all stress factors at individual deformation levels
    stress = []
    print '%10s%11s%11s%11s%11s%11s'%('S11','S22','S33','S23','S13','S12')
    for istp in range(model_rs.dat_model.nstp):
        # istp is the computational step
        # at which the SF probing was
        # performed.
        # for istp in range(1):
        #     istp = model_rs.dat_model.nstp-1
        """
        sf (nstp, k, nphi,npsi)
        ig (nstp, nphi, npsi)
        ehkl (nstp, nphi,npsi)
        strain (nstp,6)
        vf (nstp,nphi,npsi)
        """
        model_rs.sf   = model_rs.dat_model.sf[istp]
        model_rs.eps0 = model_rs.dat_model.ig[istp]
        model_rs.ehkl = model_rs.dat_model.ehkl[istp]
        if ig_sub:
            model_rs.tdat = model_rs.ehkl - model_rs.eps0
        else: model_rs.tdat = model_rs.ehkl

        tdat_ref = model_rs.tdat[::]

        if iscatter:
            tdat_scatter = []
            for iphi in range(len(tdat_ref)):
                dum = u_epshkl(tdat_ref[iphi],sigma=5e-5)
                tdat_scatter.append(dum)
            tdat_scatter = np.array(tdat_scatter)
            model_rs.tdat=tdat_scatter

        model_rs.phis = model_rs.dat_model.phi
        model_rs.psis = model_rs.dat_model.psi
        model_rs.nphi = len(model_rs.phis)
        model_rs.npsi = len(model_rs.psis)

        if sin2psimx!=None or psimx!=None:
            filter_psi(model_rs,sin2psimx=sin2psimx,psimx=psimx)
        if psi_nbin!=1:
            #psi_reso(model_rs,nbin=psi_nbin)
            psi_reso2(model_rs,ntot=psi_nbin)

        #-----------------------------------#
        ## coeff.
        # model_rs.coeff() # assign sf to cffs
        ## find the sigma ...
        dsa_sigma = model_rs.find_sigma(ivo=[0,1])
        for i in range(6): print '%+10.1f'%(dsa_sigma[i]),
        print ''
        stress.append(dsa_sigma)
        #-----------------------------------#

        plt.ioff()
        f1,f2=__model_fit_plot__(
            model_rs,ifig=ifig+istp*2+10,
            istp=istp, nxphi=nxphi)
        fe.savefig(f2);fs.savefig(f1)
        plt.close(f1);plt.close(f2)
        plt.ion()

        if (istep!=None and istp==istep) or\
           (istep==None and istp==model_rs.dat_model.nstp-1):
            fig2,fig3=__model_fit_plot__(
                model_rs,ifig=ifig+istp*2+10,
                istp=istp, nxphi=nxphi)

    fe.close(); fs.close()

    stress=np.array(stress).T # diffraction stress
    flow_dsa = FlowCurve(name='Diffraction Stress')
    flow_dsa.get_6stress(stress)
    flow_dsa.get_33strain(model_rs.dat_model.flow.epsilon)
    flow_dsa.get_eqv()

    ax1.plot(flow_dsa.epsilon_vm,flow_dsa.sigma_vm,'k+',
             label='Stress Analysis')

    # stress3x3 = flow.sigma
    # x_dsa = strain3x3[0,0] + strain3x3[1,1]
    # y_dsa = (stress3x3[0,0] + stress3x3[1,1])/2.

    # ax1.plot(x_dsa,y_dsa,'k+',
    #          label='Stress Analysis')

    # ax3.plot(x_dsa[1:],np.abs((y_dsa[1:]-y_weight[1:])/y_dsa[1:]),'o',
    #          label='error')
    ax3.set_ylim(0.,)
    fancy_legend(ax3,loc='best')

    #ax1.legend(loc='best',fancybox=True).get_frame().set_alpha(0.5)

    for i in range(len(exp_ref)):
        f = exp_ref[i]
        lab = exp_lab[i]
        edat=np.loadtxt(f).T
        ax1.plot(edat[0],edat[1],'-',lw=2,label=lab)
        ax1.set_ylim(0.,800)

    fancy_legend(ax1,size=10)

    #sigma_wgt = model_rs.dat_model.flow.sigma
    sigma_wgt = flow_weight.sigma
    ax2.plot(sigma_wgt[0,0],sigma_wgt[1,1],'b--',label='weight')
    ax2.plot(flow_dsa.sigma[0,0],flow_dsa.sigma[1,1],'k+',label='DSA')

    ax2.set_ylim(-100,700); ax2.set_xlim(-100,700)
    ax2.set_aspect('equal')
    ax2.set_xlabel(r'$\bar{\Sigma}_{11}$',dict(fontsize=15))
    ax2.set_ylabel(r'$\bar{\Sigma}_{22}$',dict(fontsize=15))
    ax2.locator_params(nbins=3)
    ax2.set_xticks(np.linspace(300,700,3),dict(fontsize=4))
    ax2.set_yticks(np.linspace(300,700,3),dict(fontsize=4))
    ax2.grid('on')
    plt.show()

    ## save figures
    fig1.savefig('flow.pdf')
    fig2.savefig('ehkl_fit.pdf')
    fig3.savefig('sf.pdf')

    return model_rs

def __model_fit_plot__(container,ifig,istp,nxphi=None):
    from MP.lib import mpl_lib
    wide_fig = mpl_lib.wide_fig
    fancy_legend=mpl_lib.fancy_legend
    rm_lab = mpl_lib.rm_lab
    tune_x_lim = mpl_lib.tune_x_lim
    from MP.lib import axes_label
    __deco__ = axes_label.__deco__

    # from mpl_lib.mpl_lib import wide_fig
    # from mpl_lib.mpl_lib import fancy_legend
    # from mpl_lib.mpl_lib import rm_lab
    # from mpl_lib.mpl_lib import tune_x_lim
    # from mpl_lib.axes_label import __deco__

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

    sf = container.dat_model.sf[istp]
    tdat = container.tdat
    ehkl = container.ehkl
    eps0 = container.eps0
    Ei   = container.Ei

    fig = wide_fig(ifig,nw=nphi,
                   w0=0.00,ws=0.5,w1=0.0,uw=3.0,
                   left=0.15,right=0.10,
                   nh=1,h0=0.2,h1=0,down=0.08,up=0.10,
                   iarange=True)

    figs= wide_fig(ifig+1,nw=nphi,
                   w0=0.00,ws=0.5,w1=0.0,uw=3.0,
                   left=0.12,right=0.10)

    #axeig = fig.axes[:nphi]
    axes= fig.axes[:nphi]#nphi:nphi*2]
    #ax_er= fig.axes[nphi*2:nphi*3]

    axesf = figs.axes
    axesv = []

    for iphi in range(nphi):
        ax=fig.axes[iphi]
        axesv.append(axes[iphi].twinx())
        axs=figs.axes[iphi]
        #axg=axeig[iphi]

        ax.set_title(r'$\phi=%3.1f^\circ$'%(phis[iphi]*180/pi))
        axs.set_title(r'$\phi=%3.1f^\circ$'%(phis[iphi]*180/pi))
        ax.locator_params(nbins=4)
        axs.locator_params(nbins=4)
        #axg.locator_params(nbins=4)
        axesv[iphi].locator_params(nbins=4)

    for iphi in range(nphi):
        ax=axes[iphi]
        av=axesv[iphi]
        #ag=axeig[iphi]
        #ae=ax_er[iphi]

        x=sin(psis)**2
        xv = sin(container.dat_model.psi)**2

        ax.plot(x,Ei[iphi]*1e6,'r--',
                label=r'$E_{i}$ (fitting)')

        # y = container.dat_model.ehkl[istp][iphi]-\
        #     container.dat_model.ig[istp][iphi]
        # y = y*1e6
        y = tdat[iphi]*1e6
        #y = tdat_ref[iphi]*1e6
        ax.plot(x,y,'bx',
                label=r'$\varepsilon^{hkl}-'\
                    '\varepsilon^{hkl}_0$')


        av.plot(xv,vf[iphi],'g-',label='Vol. Fraction')
        # ag.plot(
        #     x,eps0[iphi]*1e6,'+',color='k',
        #     label=r'$\varepsilon_{\mathrm{IG}}^{hkl}$')

        # ae.plot(
        #     x,np.abs(((eps0[iphi]-tdat[iphi]))/tdat[iphi]),
        #     'r--',label='Error')

        __deco__(ax=ax,iopt=0)
        #__deco__(ax=ag,iopt=2)

        if iphi==0:
            ax.legend(loc='upper right',fontsize=9,fancybox=True).\
                get_frame().set_alpha(0.5)
            av.legend(loc='lower right',fontsize=9,fancybox=True).\
                get_frame().set_alpha(0.5)
        ax=axesf[iphi]
        l, = ax.plot(
            xv,sf[0][iphi]*1e6,'b-o',
            label=r'$F_{11}$')
        l, = ax.plot(
            xv,sf[1][iphi]*1e6,'r-x',
            label=r'$F_{22}$')
        av.set_ylabel(r'$f(\phi,\psi)$',dict(fontsize=15))
        __deco__(ax=ax,iopt=1)
        if iphi==0:
            fancy_legend(ax)

    tune_x_lim(fig.axes,axis='x')
    tune_x_lim(axes,axis='y')
    tune_x_lim(axesv,axis='y')
    tune_x_lim(axesf,axis='y')
    # tune_x_lim(axeig,axis='y')
    # tune_x_lim(ax_er,axis='y')
    # tune_x_lim(ax_er,axis='x')

    ## remove redundant axis labels
    for i in range(len(axes)-1):
        rm_lab(axes[i+1],axis='y')
        rm_lab(axes[i+1],axis='x')
        # rm_lab(ax_er[i+1],axis='y')
        # rm_lab(ax_er[i+1],axis='x')
        # rm_lab(axeig[i+1],axis='y')
        rm_lab(axesf[i+1],axis='y')
        rm_lab(axesf[i+1],axis='x')
        rm_lab(axesv[i],axis='y')

    #for i in range(len(axes)):
        # rm_lab(axeig[i],axis='x')


    #fig.tight_layout()
    #figs.tight_layout()
    return fig, figs

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
            
