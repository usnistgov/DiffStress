"""
Analyze experimental data from proto
"""

## dependents
import numpy as np
import os
from os import sep
from RS import rs_exp
from MP.lib import mpl_lib
wide_fig=mpl_lib.wide_fig
fancy_legend = mpl_lib.fancy_legend

def main_reader(path='../dat/23JUL12', fref='Bsteel_BB_00.txt',
                fn_sf ='YJ_Bsteel_BB.sff',icheck=False,isym=False,
                fc=None,fn_str=None):
    """
    Arguments
    =========
    path (path that contains the data file)
    fref (reference data that gives strain and date-named for X-ray)
    fn_sf (stress factor file name)
    icheck
    isym
    fc     (flow curve object)
    fn_str (stress-strain data file)
    """
    if type(fc)==type(None) and type(fn_str)==type(None):
        print '---------------------------------------------'
        print 'Strain information where SF/IG were measured'
        print 'requires either fc or fn_str specified'
        print 'If not the strain column in fn_sf is used,'
        print 'subsequent analysis is performed by assuming'
        print 'that the sample is in equibiaxial strain'
        print '---------------------------------------------\n'
        # raw_input('press enter to proceed>>')

        # # dat = open('%s%s%s'%(path,sep,fn_sf)).read()
        # fn = os.path.join(path,fn_sf)
        # if os.path.isfile(fn):
        #     dat = open(fn).read()
        if os.path.isfile(fn_sf):
            dat = open(fn_sf).read()
        else:
            print 'Could not find fn_sf: '+fn_sf
            raise IOError, 'tip: use the fully pathed file name'

        ## Find proper line-breaker
        lbs=['\r','\n']; lns = []
        for i in xrange(len(lbs)):
            lns.append(len(dat.split(lbs[i])))
        lns = np.array(lns)
        lb = lbs[lns.argmax()]
        ## --

        dat_line = dat.split(lb)
        e1 = np.array(map(float,dat_line[4].split()[1:]))
        e2 = e1[::]
        e3 = -e1-e2
        from MP.mat.mech import FlowCurve as FC
        flow = FC()
        flow.get_strain(e1,0,0)
        flow.get_strain(e2,1,1)
        flow.get_strain(e3,2,2)
        flow.set_zero_shear_strain()
        flow.get_vm_strain()
        fc = flow
    elif type(fc)!=type(None):
        print 'fc was given'
        pass

    import copy
    EXP, SF, IG = rs_exp.read_main(
        path=path,fref=fref,fn_sf=fn_sf,fc=fc,fn_str=fn_str,
        icheck=icheck,isym=isym)


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
    SF.interp_psi(EXP.psis)
    ## in case that IG strain is available
    IG.interp_strain(EXP.flow.epsilon_vm)
    IG.interp_psi(EXP.psis)

    ## if not?



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
                 fn_sf='YJ_Bsteel_BB.sff',isym=False,
                 fc_sf=None,fn_str=None):
        """
        fc_sf
        """
        import copy
        self.EXP,self.SF,self.IG,self.SF_orig,self.IG_orig\
            = main_reader(
                path=path,fref=fref,fn_sf=fn_sf,
                fc=fc_sf,fn_str=fn_str,icheck=False,
                isym=isym)

        self.nstp = self.EXP.flow.nstp

        ############################################################
        ## True below if-block in order to enable
        ## 'flattening' d-spacing based on initial scatter in dspc
        if False:
            d_bar = np.zeros((self.EXP.nphi))
            d_diff = np.zeros((self.EXP.nphi,self.EXP.npsi))
            for iphi in xrange(self.EXP.nphi):
                ds = []
                for ipsi in xrange(self.EXP.npsi):
                    d = self.EXP.P_scan[0].protophi[iphi].ppscans[ipsi].dspc
                    ds.append(d)
                m = np.array(ds).mean()
                e = np.array(ds).std()
                print 'STD (d) at phi %4.0f = %f'%(self.EXP.phis[iphi],e)
                print 'Mean(d) at phi %4.0f = %f'%(self.EXP.phis[iphi],m)
                d_bar[iphi] = m
                for ipsi in xrange(self.EXP.npsi):
                    d = self.EXP.P_scan[0].protophi[iphi].ppscans[ipsi].dspc
                    d_diff[iphi,ipsi] = d - d_bar[iphi]

            # ## should I flatten the error????
            for istp in xrange(self.nstp):
                for iphi in xrange(self.EXP.nphi):
                    for ipsi in xrange(self.EXP.npsi):
                        self.EXP.P_scan[istp].protophi[iphi].ppscans[ipsi].dspc\
                            = self.EXP.P_scan[istp].protophi[iphi].ppscans[ipsi].dspc \
                            - d_diff[iphi,ipsi]
            pass
        ############################################################
        self.EXP.get_ehkl() # get ehkl based on d_avg

    def put_psi_offset(self,offset=0.0):
        self.EXP.put_psi_offset(offset=offset)
        self.SF.interp_psi(self.EXP.psis)
        self.IG.interp_psi(self.EXP.psis)

    def apply_sym(self):
        pass

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
        for iphi in xrange(nphi):
            for ipsi in xrange(npsi):
                for k in xrange(len(sigma)):
                    Ei[iphi,ipsi]\
                        = Ei[iphi,ipsi] + \
                        self.SF.sf[istp,iphi,ipsi,k] * \
                        sigma[k]
        ehkl = self.EXP.ehkl[istp].copy()
        ig   = self.IG.ig[istp].copy()
        for iphi in xrange(nphi):
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
        estimate the elastic strains...
        """
        nphi = self.EXP.nphi
        npsi = self.EXP.npsi

        self.Ei = np.zeros((nphi,npsi)) # initialize self.Ei
        sf = self.SF.sf[self.istp]

        # for k in xrange(6):
        #     self.Ei[:,:] = self.Ei[:,:] + \
        #         sf[:,:,k] * self.sigma[k]

        for iphi in xrange(nphi):
            for ipsi in xrange(npsi):
                for k in xrange(6):
                    if ivo==None or (ivo!=None and k in ivo):
                        self.Ei[iphi,ipsi] = \
                         self.Ei[iphi,ipsi] + \
                         sf[iphi,ipsi,k] * self.sigma[k]

    def f_least_Ei_d0(self,array=[1.701,0,0,0,0,0],ivo=None,
                      wgt=None):
        """
        Call the least square objective function
        with allowing d0 to vary
        """
        d0 = array[0]
        self.EXP.assign_d0(d0)
        return self.f_least_Ei(stress=array[1:],ivo=ivo,wgt=wgt)

    def f_least_Ei_fixed_d0(self,array=[0,0,0,0,0,0],
                            ivo=None,d0=None,wgt=None):
        """
        Call the least square objective function
        without allowing d0 to vary
        """
        self.EXP.assign_d0(d0)
        return self.f_least_Ei(stress=array,ivo=ivo,wgt=wgt)

    def f_least_Ei(self,stress=[0,0,0,0,0,0],ivo=None,
                   wgt=None):
        """
        The least square objective function

        Arguments
        =========
        stress = [0,0,0,0,0,0]
        ivo    = None
        wgt    = None
        """
        self.sigma=np.array(stress)
        nphi = self.EXP.nphi
        npsi = self.EXP.npsi
        self.calc_Ei(ivo=ivo)

        f_array=[]

        for iphi in xrange(nphi):
            for ipsi in xrange(npsi):
                d = self.EXP.ehkl[self.istp,iphi,ipsi] -\
                    self.IG.ig[self.istp,iphi,ipsi] -\
                    self.Ei[iphi,ipsi]
                if np.isnan(d): d = 0
                if type(wgt)!=type(None):
                    d = d * wgt[iphi,ipsi]
                f_array.append(d)
        return np.array(f_array)

    def residuals(self,a,x,y,f): return y-f(x,a)

    def find_sigma(self,ivo=None,istp=0,iplot=False,
                   iwgt=True,ifix_d0=False,d0=1.701):
        """
        Find stress.

        Arguments
        =========
        ivo = None
        istp=0
        iplot=False
        iwgt=True
        ifix_d0=False
        """
        from scipy import optimize
        fmin = optimize.leastsq
        self.istp=istp

        # print '#----------------------#'
        # print ' epsilon_{VM} %4.2f'%self.EXP.flow.epsilon_vm[istp]
        # print ' ivo:', ivo
        # --core

        if iwgt==False: wgt=None
        elif iwgt:
            wgt = np.zeros((self.EXP.nphi,self.EXP.npsi))
            for iphi in xrange(self.EXP.nphi):
                for ipsi in xrange(self.EXP.npsi):
                    wgt[iphi,ipsi]=self.EXP.P_scan[
                        istp].protophi[iphi].ppscans[ipsi].ints

        if ifix_d0:
            dat = fmin(self.f_least_Ei_fixed_d0,
                       [100,100,0,0,0,0],args=(ivo,d0,wgt),
                       full_output=True)
            stress=dat[0]
        elif not(ifix_d0):
            dat = fmin(self.f_least_Ei_d0,
                       [1.1710,100,100,0,0,0,0],args=(ivo,wgt),
                       full_output=True)
            ## --core
            stress = dat[0][1:]
            d0 = dat[0][0]

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
            for iphi in xrange(nphi):
                l, = figs.axes[iphi].plot(x,Ei[iphi,:],'--')
                figs.axes[iphi].plot(x,ehkl[iphi,:],'-x',
                                     color=l.get_color())

        stress = np.array(stress) / 1e6

        return stress, d0

def main(path='../dat/BB/',
         fref='Bsteel_fref_DIC.txt',fn_sf='YJ_BB_10times.sff',
         fexp=None,iso_SF=False,ishow=False,ind_plot=False,
         psi_offset=0.0,psi_sym=False,fc=None,fn_str=None,ifix_d0=False,
         d0_ref=1.17025,
         iwgt=False):
    """
    Arguments
    =========
    path
    fref
    fn_sf
    fexp
    iso_SF
    ishow
    ind_plot
    psi_offset
    psi_sym
    fc
    fn_str
    ifix_d0
    d0_ref=1.17025
    iwgt=False
    """
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

    if type(fexp)==type(None):
        fexp='/Users/yj/repo/evpsc-dev/'\
            'exp_dat/Bsteel/bulge/EXP_BULGE_JINKIM.txt',

    RS_graphs = PdfPages('RS_Graphs.pdf')
    mystress = StressAnalysis(path=path,fref=fref,
                              fn_sf=fn_sf,isym=psi_sym,
                              fc_sf=fc,fn_str=fn_str)

    if psi_offset!=0: mystress.put_psi_offset(psi_offset)
    if psi_sym:
        print 'psi symmetry has been applied.'
        mystress.apply_sym()

    mystress.nstp
    mystress.EXP.plot_all()
    RS_graphs.savefig(plt.gcf())

    # isotropic stress factor?
    if iso_SF: mystress.SF.get_iso_sf(E=204e9,nu=0.3)

    # calc stress
    eps_vm = mystress.EXP.flow.epsilon_vm
    dknot = []
    s11 = []; s22 = [];# s12 = []
    Eis = []; eps = []; igs = []

    d_ehkl = np.zeros((mystress.EXP.nphi))
    for istp in xrange(mystress.nstp):
        stress, d0 = mystress.find_sigma(
            ivo=[0,1],istp=istp,iplot=False,iwgt=iwgt,
            ifix_d0=ifix_d0,d0=d0_ref)

        dknot.append(d0)
        s11.append(stress[0])
        s22.append(stress[1])

        Eis.append(mystress.Ei) # fitted Ei
        ehkl = mystress.EXP.ehkl[istp] # experimental eps_hkl
        ig   = mystress.IG.ig[istp]    # experimental eps_ig
        eps.append(ehkl.copy())
        igs.append(ig.copy())

        if istp==0:
            for iphi in xrange(mystress.EXP.nphi):
                d_ehkl[iphi] = np.array(ehkl[iphi][::]).std()

    print '-----------------------------------'
    print 'Standard deviation in d_ehkl at istp=0\n'
    print 'phi:',
    for iphi in xrange(mystress.EXP.nphi):
        print '%7.0f '%mystress.EXP.phis[iphi],
    print 'avg'
    print 'std:',
    for iphi in xrange(mystress.EXP.nphi):
        print '%7.1e '%d_ehkl[iphi],
    print '%7.1e '%d_ehkl.mean()

    # macro flow object
    mystress.flow = mech.FlowCurve()
    mystress.flow.get_stress(s11,0,0)
    mystress.flow.get_stress(s22,1,1)
    mystress.flow.set_zero_sigma_ij(i=2,j=2)
    # mystress.flow.get_stress(s12,0,1)
    # mystress.flow.get_stress(s12,1,0)


    mystress.flow.set_zero_sigma_ij(i=0,j=1)
    mystress.flow.set_zero_sigma_ij(i=1,j=0)
    mystress.flow.set_zero_sigma_ij(i=2,j=1)
    mystress.flow.set_zero_sigma_ij(i=1,j=2)
    mystress.flow.set_zero_sigma_ij(i=0,j=2)
    mystress.flow.set_zero_sigma_ij(i=2,j=0)
    mystress.flow.get_vm_stress()

#------------------------------------------------------------#
    # plots that illustrate the fitted curves...
    Eis = np.array(Eis)
    eps = np.array(eps)
    igs = np.array(igs)
    figs = wide_fig(nw=mystress.EXP.nphi,w0=0,w1=0,left=0.2,
                    right=0.15)

    from MP.lib import mpl_lib
    # mx = eps_vm[mystress.EXP.flow.nstp-1],
    # mn = eps_vm[0])
    mx = 1.0; mn = 0.0
    norm = mpl.colors.Normalize(vmin=mn, vmax=mx)
    cmap, c = mpl_lib.norm_cmap(
        mx=mx,mn=mn) # c.to_rbga()

    for iphi in xrange(mystress.EXP.nphi):
        x = mystress.EXP.psi
        ax = figs.axes[iphi]
        ax.set_title(
            r'$\phi: %3.1f^\circ{}$'%mystress.EXP.phi[iphi])
        for istp in xrange(mystress.EXP.flow.nstp):
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
                   ylab='Equivalent Strain',norm=norm)

    RS_graphs.savefig(plt.gcf())

#------------------------------------------------------------#
    # Save fitted curves into individual plots to pdf files
    if ind_plot:
        from matplotlib.backends.backend_pdf import PdfPages
        pdf_pages = PdfPages('Fit_results.pdf')
        x = mystress.EXP.psi
        for istp in xrange(mystress.EXP.flow.nstp):
            figx = wide_fig(nw=mystress.EXP.nphi+1,
                            w0=0,w1=0,left=0.2,
                            right=0.15)
            for iphi in xrange(mystress.EXP.nphi):
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

            for i in xrange(len(figx.axes)-2):
                figx.axes[i].set_ylim(-0.0010,0.0006)

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
    figs.axes[0].plot(mystress.SF_orig.flow.epsilon_vm,
                      np.zeros((mystress.SF_orig.flow.nstp)),
                      'o',ms=8,mfc='None',mec='r')


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
