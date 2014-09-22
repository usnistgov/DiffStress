## MS&T 2014 Fall meeting application
## - plot flow stress along various paths

from os import sep, popen
from glob import glob
from RS import rs_ex
from MP.lib import mpl_lib,axes_label
import matplotlib.pyplot as plt
wf = mpl_lib.wide_fig
fl = mpl_lib.fancy_legend

def main_plot_flow_all(hkl='211',sin2psimx=0.5,psi_nbin=12,pmargin=None):
    fns = read(hkl)
    fig = wf(nw=2,nh=1,left=0.2,uw=3.5,
             w0=0,w1=0.3,right=0,iarange=True)
    ax1,ax2 = fig.axes
    axes_label.__eqv__(ax1,ft=10)

    for i in range(len(fns)):
        fn = fns[i]
        strain_path = fn.split('_')[0]
        untargz(fn)
        model_rs, flow_wgt, flow_dsa \
            = rs_ex.ex_consistency(
                sin2psimx=sin2psimx,
                psi_nbin=psi_nbin,
                hkl=hkl,iplot=False,pmargin=pmargin)

        flow_wgt.get_eqv(); flow_dsa.get_eqv()

        ax1.plot(flow_wgt.epsilon_vm,flow_wgt.sigma_vm,
                 'bx-',label='Weighted Avg', alpha=1.0)
        ax1.plot(flow_dsa.epsilon_vm,flow_dsa.sigma_vm,
                 'kx',label='Diff. Stress Analsysis', alpha=1.0)
        ax2.plot(flow_wgt.sigma[0,0],flow_wgt.sigma[1,1],'b+')
        ax2.plot(flow_dsa.sigma[0,0],flow_dsa.sigma[1,1],'k+')

        ## connector
        npoints = len(flow_wgt.sigma[0,0])
        wgtx = flow_wgt.sigma[0,0]; wgty = flow_wgt.sigma[1,1]
        dsax = flow_dsa.sigma[0,0]; dsay = flow_dsa.sigma[1,1]
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

def untargz(fn): std = popen('tar -xzf %s'%fn)

def read(hkl='211'):
    fns=glob('*_%s_*.tar.gz'%hkl)
    print '%i files are found'%len(fns)
    return fns



def plot_sf(sff_fn='temp.sff'):
    from RS import sff_converter, sfig_class
    from rs_exp import read_IGSF
    sff_converter.main(fn=sff_fn,difile=None,itab=True,
                       ieps0=4,fn_str='STR_STR.OUT')
    SF,dum = read_IGSF(fn=sff_fn,fn_str='STR_STR.OUT')

    ## mask data lacking sufficient volume

    SF.plot(nbin_sin2psi=3)
