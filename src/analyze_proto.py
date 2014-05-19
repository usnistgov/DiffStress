##
import numpy as np
from RS import rs_exp
from MP.lib import mpl_lib
wide_fig=mpl_lib.wide_fig
fancy_legend = mpl_lib.fancy_legend

def main_reader(path='../dat/23JUL12', fref='Bsteel_BB_00.txt',
                fn_sf ='YJ_Bsteel_BB.sff',icheck=False):
    import copy
    EXP, SF, IG = rs_exp.read_main(path,fref,fn_sf,icheck)

    ## interpolate based on experimental phis, psis, IG

    # EXP.plot(istps=[0,10,20])
    # SF.plot()
    # IG.plot()

    print
    print '#-----------------------------------------------------#'
    print ' Interpolate SF and IG for matching with D-spacings'

    SF_orig = copy.deepcopy(SF)
    IG_orig = copy.deepcopy(IG)
    SF_orig.flow.get_vm_strain()
    IG_orig.flow.get_vm_strain()

    SF.interp_strain(EXP.flow.epsilon_vm)
    IG.interp_strain(EXP.flow.epsilon_vm)

    SF.interp_psi(EXP.psis)
    IG.interp_psi(EXP.psis)

    if SF.phi!=EXP.phis:
        print '  ** Phi angles of SF are different **'
    if IG.phi!=EXP.phis:
        print '  ** Phi angles of IG are different **'
    print '#-----------------------------------------------------#'

    # SF.plot()
    # IG.plot()

    SF.determine_phis(phi_new=EXP.phi)
    IG.determine_phis(phi_new=EXP.phi)

    SF.plot()
    IG.plot()

    #EXP.plot(istps=np.arange(EXP.flow.nstp)[::7])

    return EXP, SF, IG, SF_orig, IG_orig

class StressAnalysis:
    def __init__(self,
                 path='/Users/yj/repo/rs_pack/dat/23Jul12',
                 fref='Bsteel_BB_00.txt',
                 fn_sf='YJ_Bsteel_BB.sff'):
        self.EXP,self.SF,self.IG,self.SF_orig,self.IG_orig\
            = main_reader(path,fref,fn_sf,icheck=False)

        self.EXP.get_ehkl() # get ehkl based on d_avg
        self.nstp = self.EXP.flow.nstp

    def get_visualize(self,sigma=[0,0],istp=0,uni='MPa',d0=None):
        """
        Let's visualize the procedure. (sigma in Pascal)
        """
        sigma = np.array(sigma) / 1e6
        if d0!=None:
            self.EXP.assign_d0(d0)
        import matplotlib.pyplot as plt
        plt.close('all')
        nphi = self.EXP.nphi
        npsi = self.EXP.npsi
        figs = wide_fig(nw=nphi)
        # calculate Ei based on the current stress & SF
        Ei = np.zeros((nphi, npsi))
        for iphi in range(nphi):
            for ipsi in range(npsi):
                for k in range(len(sigma)):
                    Ei[iphi,ipsi]\
                        = Ei[iphi,ipsi] + \
                        self.SF.sf[istp,iphi,ipsi,k] * \
                        sigma[k]
        ehkl = self.EXP.ehkl[istp].copy()
        ig   = self.IG.ig[istp].copy()
        for iphi in range(nphi):
            #x = np.sin(self.EXP.psi*np.pi**180)**2
            x = self.EXP.psi[::]
            l, = figs.axes[iphi].plot(x,Ei[iphi,:],'-rx',
                                      label='Fitting')
            # figs.axes[iphi].plot(x,ehkl[iphi,:],'-ko',
            #                      label=r'$\varepsilon^{hkl}$')
            figs.axes[iphi].plot(
                x,ehkl[iphi,:]-ig[iphi,:],'-ko',
                label=r'$\varepsilon^{hkl}-\varepsilon^{IG}$')
        fancy_legend(figs.axes[iphi])
        if d0!=None:
            self.EXP.get_ehkl()

    def calc_Ei(self,ivo=None):
        """
        Given the 'assumed' macro stress (sigma),
        esitmate the elastic strains...
        """
        nphi = self.EXP.nphi
        npsi = self.EXP.npsi

        self.Ei = np.zeros((nphi,npsi)) # initialize self.Ei
        sf = self.SF.sf[self.istp]

        # for k in range(6):
        #     self.Ei[:,:] = self.Ei[:,:] + \
        #         sf[:,:,k] * self.sigma[k]

        for iphi in range(nphi):
            for ipsi in range(npsi):
                for k in range(6):
                    if ivo==None or (ivo!=None and k in ivo):
                        self.Ei[iphi,ipsi] = \
                         self.Ei[iphi,ipsi] + \
                         sf[iphi,ipsi,k] * self.sigma[k]

    def f_least_Ei_d0(self,array=[1.701,0,0,0,0,0],ivo=None):
        d0 = array[0]
        self.EXP.assign_d0(d0)
        return self.f_least_Ei(stress=array[1:],ivo=ivo)

    def f_least_Ei_fixed_d0(self,array=[0,0,0,0,0,0],
                            ivo=None,d0=None):
        self.EXP.assign_d0(d0)
        return self.f_least_Ei(stress=array,ivo=ivo)

    def f_least_Ei(self,stress=[0,0,0,0,0,0],ivo=None):
        self.sigma=np.array(stress)
        nphi = self.EXP.nphi
        npsi = self.EXP.npsi
        self.calc_Ei(ivo=ivo)

        f_array=[]
        for iphi in range(nphi):
            for ipsi in range(npsi):
                d = self.EXP.ehkl[self.istp,iphi,ipsi] -\
                    self.IG.ig[self.istp,iphi,ipsi] -\
                    self.Ei[iphi,ipsi]
                if np.isnan(d): d = 0
                f_array.append(d)
        return np.array(f_array)

    def residuals(self,a,x,y,f):
        return y-f(x,a)

    def find_sigma(self,ivo=None,istp=0,iplot=False):
        from scipy import optimize
        fmin = optimize.leastsq
        self.istp=istp

        # print '#----------------------#'
        # print ' epsilon_{VM} %4.2f'%self.EXP.flow.epsilon_vm[istp]
        # print ' ivo:', ivo
        # --core

        dat = fmin(self.f_least_Ei_d0,
                   [1.1710,100,100,0,0,0,0],ivo,
                   full_output=True)
        ## --core
        stress = dat[0][1:]
        d0 = dat[0][0]

        # dat = fmin(self.f_least_Ei_fixed_d0,
        #            [100,100,0,0,0,0],args=(ivo,d0),
        #            full_output=True)
        # stress=dat[0]

        cov_x, infodict, mesg, ier = dat[1:]
        if not(ier in [1,2,3,4]):
            raise IOError, "solution was not found"

        if iplot:
            Ei = self.Ei.copy()
            ehkl = self.EXP.ehkl[self.istp,:,:]
            sf = self.SF.sf[self.istp,:,:,:]

            nphi = self.EXP.nphi
            x = self.EXP.psi
            figs = wide_fig(nw=nphi)
            for iphi in range(nphi):
                l, = figs.axes[iphi].plot(x,Ei[iphi,:],'--')
                figs.axes[iphi].plot(x,ehkl[iphi,:],'-x',
                                     color=l.get_color())

        stress = np.array(stress) / 1e6

        return stress, d0

def main(path='/Users/yj/repo/rs_pack/dat/23Jul12',
         fref='Bsteel_BB_00.txt',fn_sf='YJ_Bsteel_BB.sff',
         fexp=None,
         iso_SF=False,
         ishow=False,
         ind_plot=False):
    """
    """
    if fexp==None:
        fexp='/Users/yj/repo/evpsc-dev/'\
            'exp_dat/Bsteel/bulge/EXP_BULGE_JINKIM.txt',
    from MP.mat import mech
    from MP.lib import axes_label
    from MP.lib import mpl_lib
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    eqv = axes_label.__eqv__
    tune_xy_lim = mpl_lib.tune_xy_lim
    deco = axes_label.__deco__
    rm_inner =mpl_lib.rm_inner
    rm_lab = mpl_lib.rm_lab
    ticks_bin_u = mpl_lib.ticks_bins_ax_u

    plt.ioff()

    RS_graphs = PdfPages('RS_Graphs.pdf')

    mystress = StressAnalysis(path=path,fref=fref,fn_sf=fn_sf)
    mystress.nstp
    mystress.EXP.plot_all()
    RS_graphs.savefig(plt.gcf())

    # isotropic stress factor?
    if iso_SF: mystress.SF.get_iso_sf(E=204e9,nu=0.3)

    # calc stress
    eps_vm = mystress.EXP.flow.epsilon_vm
    dknot = []
    s11 = []; s22 = []
    Eis = []; eps = []; igs = []

    for istp in range(mystress.nstp):
        stress, d0 = mystress.find_sigma(
            ivo=[0,1],istp=istp,iplot=False)

        dknot.append(d0)
        s11.append(stress[0])
        s22.append(stress[1])

        Eis.append(mystress.Ei) # fitted Ei
        ehkl = mystress.EXP.ehkl[istp] # experimental eps_hkl
        ig   = mystress.IG.ig[istp]    # experimental eps_ig
        eps.append(ehkl.copy())
        igs.append(ig.copy())

    # macro flow object
    mystress.flow = mech.FlowCurve()
    mystress.flow.get_stress(s11,0,0)
    mystress.flow.get_stress(s22,1,1)
    mystress.flow.set_zero_sigma_ij(i=2,j=2)
    mystress.flow.set_zero_shear_stress()
    mystress.flow.get_vm_stress()

#------------------------------------------------------------#
    # plots that illustrate the fitted curves...
    Eis = np.array(Eis)
    eps = np.array(eps)
    igs = np.array(igs)
    figs = wide_fig(nw=mystress.EXP.nphi,w0=0,w1=0,left=0.2,
                    right=0.15)

    from MP.lib import mpl_lib
    cmap, c = mpl_lib.norm_cmap(
        mx = eps_vm[mystress.EXP.flow.nstp-1],
        mn = eps_vm[0]) # c.to_rbga()

    for iphi in range(mystress.EXP.nphi):
        x = mystress.EXP.psi
        ax = figs.axes[iphi]
        ax.set_title(
            r'$\phi: %3.1f^\circ{}$'%mystress.EXP.phi[iphi])
        for istp in range(mystress.EXP.flow.nstp):
            y_fit = Eis[istp,iphi]
            y1=eps[istp,iphi] # e^{hkl}
            y0=igs[istp,iphi] # e^{ig}
            y_exp = y1-y0
            ax.plot(x,y_fit,'--',color=c.to_rgba(eps_vm[istp]))
            ax.plot(x,y_exp,'x',color=c.to_rgba(eps_vm[istp]))

    deco(ax=figs.axes[0],iopt=5)
    rm_inner(figs.axes)
    ticks_bin_u(figs.axes)
    tune_xy_lim(figs.axes)

    b = figs.axes[-1].get_position()
    axcb = figs.add_axes([0.88,b.y0,0.03,b.y1-b.y0])
    mpl_lib.add_cb(axcb,cmap=cmap,filled=True,
                   ylab='Equivalent Strain')

    RS_graphs.savefig(plt.gcf())

#------------------------------------------------------------#
    # Save fitted curves into individual plots to pdf files
    if ind_plot:
        from matplotlib.backends.backend_pdf import PdfPages
        pdf_pages = PdfPages('Fit_results.pdf')
        x = mystress.EXP.psi
        for istp in range(mystress.EXP.flow.nstp):
            figx = wide_fig(nw=mystress.EXP.nphi+1,
                            w0=0,w1=0,left=0.2,
                            right=0.15)
            for iphi in range(mystress.EXP.nphi):
                x = mystress.EXP.psi
                ax = figx.axes[iphi]
                ax.set_title(
                    r'$\phi: %3.1f^\circ{}$'%\
                    mystress.EXP.phi[iphi])

                y_fit = Eis[istp,iphi]
                y1=eps[istp,iphi] # e^{hkl}
                y0=igs[istp,iphi] # e^{ig}
                y_exp = y1-y0
                ax.plot(x,y_fit,'r-',label='Fitting')
                ax.plot(x,y_exp,'bx',
                        label=r'$\varepsilon^{hkl}'\
                        '-\varepsilon^{IG}$')

            ax = figx.axes[-1]
            eqv(ax,ft=8,zero_xy=False)
            ax = ax.twinx()
            eqv(ax,ft=8,zero_xy=False)
            ax.plot(eps_vm,mystress.flow.sigma_vm,'b-')
            ax.plot(eps_vm[istp],
                    mystress.flow.sigma_vm[istp],'ro')

            for i in range(len(figx.axes)-2):
                figx.axes[i].set_ylim(-0.0008,0.0004)

            fancy_legend(figx.axes[0])
            rm_inner(figx.axes[:4])
            ticks_bin_u(figx.axes[:4])

            #tune_xy_lim(figx.axes)
            ax.set_xlim(0.,1)
            rm_lab(ax=figx.axes[-2],axis='y')
            deco(figx.axes[0],iopt=5)

            pdf_pages.savefig(figx)
            plt.close(figx)
        pdf_pages.close()
        plt.ion()
#------------------------------------------------------------#
    ## flow stress curve plotting

    figs = wide_fig(nw=2,w1=0.2)
    figs.axes[0].plot(eps_vm, mystress.flow.sigma_vm,'x')
    figs.axes[1].plot(eps_vm, dknot,'o')
    eqv(figs.axes[0],ft=8,zero_xy=True)
    eqv(figs.axes[1],ft=8,zero_xy=False)
    figs.axes[1].set_ylabel(r'$\mathrm{d}_o$',dict(fontsize=12))
    RS_graphs.savefig(plt.gcf())

#------------------------------------------------------------#
    ## flow stress in plane stress space
    figs = wide_fig(nw=2)
    figs.axes[0].plot(
        mystress.flow.sigma[0,0],
        mystress.flow.sigma[1,1],'-x')
    figs.axes[1].plot(
        mystress.EXP.flow.epsilon[0,0],
        mystress.EXP.flow.epsilon[1,1],'-x')
    axes_label.__plane__(ax=figs.axes[0],ft=10,iopt=0)
    axes_label.__plane__(ax=figs.axes[1],ft=10,iopt=1)
    RS_graphs.savefig(plt.gcf())
    RS_graphs.close()
    if not(ishow): plt.close('all')
    plt.ion()
    return mystress
