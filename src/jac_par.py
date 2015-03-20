"""
Applications for JAC paper
"""
from RS.rs_ex import ex_consistency as main
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from MP.lib import mpl_lib, axes_label
from RS import lib
import numpy as np
sin2psi_opt = lib.sin2psi_opt
tune_xy_lim  = mpl_lib.tune_xy_lim
tune_x_lim  = mpl_lib.tune_x_lim
fancy_legend = mpl_lib.fancy_legend
deco         = axes_label.__deco__
sin2psi_opt = lib.sin2psi_opt

##
sigmas=[5e-5,1e-4,2e-4] # best to worst
nbins=[29,15,11,7] # best to worst

gs=gridspec.GridSpec(
    len(sigmas),len(nbins),
    wspace=0,hspace=0,left=0.25,right=0.8,top=0.8)

def tabular_figs02():
    fig=plt.figure()
    axes=[]
    for i in range(len(sigmas)):
        axes.append([])
        for j in range(len(nbins)):
            if sigmas[i]==0: iscatter=False
            else:            iscatter= True

            rst = main(
                ## main variables
                sigma=sigmas[i], psi_nbin=nbins[j],

                ## minor
                sin2psimx=0.5, iscatter=iscatter,
                iplot=False,dec_inv_frq=1,
                dec_interp=1)
            model_rs, flow_weight, flow_dsa = rst
            ax=fig.add_subplot(gs[i,j])
            ax.locator_params(nbins=4)
            axes[i].append(ax)

            # von mises
            lab1=r'$(\langle \sigma^c \rangle)^{VM}$'
            lab2=r'$(\sigma^\mathrm{RS})^{VM}$'
            ax.plot(flow_weight.epsilon_vm,
                    flow_weight.sigma_vm,'k-',label=lab1)
            ax.plot(flow_dsa.epsilon_vm,
                    flow_dsa.sigma_vm,'kx',label=lab2)

            if i==len(sigmas)-1 and j==0:
                #VM Deco
                axes_label.__vm__(ax=ax,ft=10)
            else:
                mpl_lib.rm_lab(ax,axis='x')
                mpl_lib.rm_lab(ax,axis='y')


    for iax in range(len(fig.axes)):
        fig.axes[iax].set_ylim(0,)
        fig.axes[iax].set_xlim(0,)
    tune_xy_lim(fig.axes)
    tune_x_lim(fig.axes,axis='y')



    ## annotations
    for j in range(len(nbins)):
        if j==0: s=r'$N^\psi$=%i'%nbins[j]
        else:
            s = '%i'%nbins[j]
        axes[0][j].annotate(
            s=s,horizontalalignment='center',
            size=14,xy=(0.5,1.2),
            xycoords='axes fraction')

    for i in range(len(sigmas)):
        if i==0:
            s=r'$s^\mathrm{CSE}$=%i'%(
                sigmas[i]*(10**6))
        elif i==len(sigmas)-1:
              s=r'%i $\mu$ strain'%(
                  sigmas[i]*(10**6))
        else: s='%i'%(sigmas[i]*(10**6))

        axes[i][-1].annotate(
            s=s,
            verticalalignment='center',
            horizontalalignment='center',
            rotation=270,
            size=14,xy=(1.20,0.5),
            xycoords='axes fraction')

    fancy_legend(
        fig.axes[0],size=11,nscat=1,ncol=1,
        bbox_to_anchor=(-0.4,1))

    # vm
    fig.axes[len(sigmas)*len(nbins)-len(nbins)].\
        set_xticks(np.arange(0.,1.01,0.5))
    fig.savefig('tab2.pdf')
    fig.savefig('tab2.png')
    fig.clf()

def tabular_figs01():
    """
    Ehkl vs sin2psi plot

    y (sigmas)
    x (nbins)
    """
    fig=plt.figure()
    axes=[]
    for i in range(len(sigmas)):
        axes.append([])
        for j in range(len(nbins)):
            ax=fig.add_subplot(gs[i,j])
            axes[i].append(ax)
            ax.locator_params(nbins=4);
            if sigmas[i]==0: iscatter=False
            else:            iscatter= True
            rst = main(
                ## main variables
                sigma=sigmas[i], psi_nbin=nbins[j],

                ## minor
                sin2psimx=0.5, iscatter=iscatter,
                istep=2,iplot=False,dec_inv_frq=1,
                dec_interp=1)

            model_rs, sig11, sig22, dsa11,dsa22,\
                raw_psis,raw_vfs, raw_sfs,\
                full_Ei,dec_intp  = rst

            # raise IOError
            x = sin2psi_opt(model_rs.psis, 1)
            ax.plot(
                x, model_rs.tdat[0]*1e6,'k.',
                label=\
                r'$\tilde{ \langle\varepsilon^e \rangle}^G$')
            ax.plot(
                x, model_rs.Ei[0]*1e6,'kx',
                label=\
                r'$\mathbb{F}_{ij} \sigma^\mathrm{RS}_{ij}$')

            x = sin2psi_opt(raw_psis, 1)
            # ax.plot(x, full_Ei[0]*1e6,'k-') ## continuous Ei
            ## ax.plot(raw_psis, dec_intp[0][0]*1e6,'k--')

            ## True Ei? = F_ij * <sigma>^c_{ij}

            model_rs.psis = raw_psis.copy()
            model_rs.npsi = len(model_rs.psis)
            model_rs.cffs = raw_sfs.copy()
            ## plug the weight average stress
            model_rs.sigma=[sig11, sig22, 0, 0, 0, 0]
            model_rs.calc_Ei()
            ax.plot(
                x,model_rs.Ei[0]*1e6,'k-',
                label=r'$\mathbb{F}_{ij} \langle\sigma\rangle^c_{ij}$'
            )

            if i==len(sigmas)-1 and j==0:
                deco(ax=ax,iopt=0,hkl=None,ipsi_opt=1)
            else:
                mpl_lib.rm_lab(ax,axis='x')
                mpl_lib.rm_lab(ax,axis='y')

    for j in range(len(nbins)):
        if j==0: s=r'$N^\psi$=%i'%nbins[j]
        else:
            s = '%i'%nbins[j]
        axes[0][j].annotate(
            s=s,horizontalalignment='center',
            size=14,xy=(0.5,1.2),
            xycoords='axes fraction')

    for i in range(len(sigmas)):
        if i==0:
            s=r'$s^\mathrm{CSE}$=%i'%(
                sigmas[i]*(10**6))
        elif i==len(sigmas)-1:
              s=r'%i $\mu$ strain'%(
                  sigmas[i]*(10**6))
        else: s='%i'%(sigmas[i]*(10**6))

        axes[i][-1].annotate(
            s=s,
            verticalalignment='center',
            horizontalalignment='center',
            rotation=270,
            size=14,xy=(1.20,0.5),
            xycoords='axes fraction')

    for iax in range(len(fig.axes)):
        fig.axes[iax].set_ylim(-1500,)
    tune_xy_lim(fig.axes)
    tune_x_lim(fig.axes,axis='y')

    # plt.annotate(s=r'$\varepsilon^\mathrm{CSE}$ Counting Statistics Error',
    #              size=20,rotation=90,
    #              verticalalignment='center',
    #              horizontalalignment='center',
    #              xy=(1.01,0.5),
    #              xycoords='figure fraction')

    # plt.annotate(s=r'Number of $\psi$',size=20,
    #              horizontalalignment='center',
    #              xy=(0.6,0.93),
    #              xycoords='figure fraction')


    fig.axes[len(sigmas)*len(nbins)-len(nbins)].\
        set_xticks(np.arange(-0.5,0.501,0.5))

    fancy_legend(
        fig.axes[0],size=11,nscat=1,ncol=1,
        bbox_to_anchor=(-0.4,1))
    fig.savefig('tab.pdf')
    fig.savefig('tab.png')
    fig.clf()


def DEC_intp(ss=[1,2,4],intps=[0,3,4],inds=[79,90,120]):
    """
    DEC interpolation comparison

    intps=0: NN (piecewise linear interpolation)
         =1: Assign nearest data
         =2: Cubic
         =3: Quadratic
         =4: Linear fit
         =5: Poly2
         =6: poly3
         =7: Place-holder for power law fit
         =8: zero
         =9: slinear
    inds=[70,40]
    """

    from MP.mat import mech # mech is a module
    FlowCurve = mech.FlowCurve
    from mst_ex import use_intp_sfig, return_vf
    from rs import ResidualStress
    model_rs = ResidualStress(
        mod_ext=None,fnmod_ig='igstrain_fbulk_ph1.out',
        fnmod_sf='igstrain_fbulk_ph1.out',i_ip=1)
    psi = model_rs.dat_model.psi
    sin2psi = sin2psi_opt(psi,1)

    DEC_raw = model_rs.dat_model.sf.copy()
    fc = FlowCurve(name='Weighted Average')
    fc.get_model(fn='STR_STR.OUT')
    fc.get_eqv()

    evm = fc.epsilon_vm.copy()
    fig = plt.figure()

    gs=gridspec.GridSpec(
        len(ss),len(intps),
        wspace=0,hspace=0,left=0.25,right=0.8,top=0.8)

    axes=[]
    ms=['o','x','+','d','>','t']
    for i in range(len(ss)):
        axes.append([])
        for j in range(len(intps)):
            ## interpolated DECs
            print i,j
            print ss[i],intps[j]
            # if i==1 and j==0:
            #     raise IOError
            _sf_, _ig_ = use_intp_sfig(
                ss[i],iopt=intps[j],
                iplot=False,iwgt=False)
            dum = _sf_.sf.swapaxes(-1,-2).swapaxes(-2,-3)
            DEC_intp = dum.copy()

            ax = fig.add_subplot(gs[i,j])
            ax.locator_params(nbins=4)
            axes[i].append(ax)

            for k in range(len(inds)):
                ind = inds[k]
                val_sin2psi = sin2psi[ind]
                val_psi     = psi[ind] * 180./np.pi
                lab1, lab2 = None, None
                if k==0:
                    lab1=r'Actual $\mathbb{F}_{ij}$'
                    lab2=r'Interpolated $\mathbb{F}^{\ I}_{ij}$'
                y_raw = DEC_raw[:,0,0,ind]*1e6
                ax.plot(evm, y_raw,'k-',
                        label=lab1)
                y_intp =  DEC_intp[:,0,0,ind]*1e12
                ax.plot(evm,y_intp,'k--',
                        label=lab2)
                ax.fill_between(evm,y_raw,y_intp,
                                facecolor='gray',alpha=0.5)
                ax.plot(
                    evm[::ss[i]],
                    DEC_raw[:,0,0,ind][::ss[i]]*1e6,
                    ms[k],mfc='None',mec='black',
                    label=r'$\psi=%2.0f^\circ{}$'%val_psi)

            if i==len(ss)-1 and j==0:
                pass
            else:
                mpl_lib.rm_lab(ax,axis='x')
                mpl_lib.rm_lab(ax,axis='y')

    deco(ax=axes[-1][0],ft=10,iopt=6,ij=[1,1])
    axes[-1][0].grid('off')
    tune_xy_lim(fig.axes)
    tune_x_lim(fig.axes,axis='y')


    ## annotations
    for j in range(len(intps)):
        if j==0: s = 'Piecewise'
        if j==1: s = 'Quadratic'
        if j==2: s = 'Linear fit'
        if j==3: s = 'Power law fit'

        axes[0][j].annotate(
            s=s, horizontalalignment='center',
            size=10,xy=(0.5,1.2),
            xycoords='axes fraction')

    for i in range(len(ss)):
        s = r'$f^{\ \mathbb{F}}=^1/_{%i}$'%ss[i]
        axes[i][-1].annotate(
            s=s,
            verticalalignment='center',
            horizontalalignment='center',
            rotation=270,
            size=10,xy=(1.20,0.5),
            xycoords='axes fraction')


    fancy_legend(axes[0][0],size=10,nscat=1,ncol=1,
                 bbox_to_anchor=(-0.1,1))
    fig.savefig('dec_intp.pdf')
    fig.savefig('dec_intp.png')
    return model_rs
