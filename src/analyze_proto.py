##
import numpy as np
from RS import rs_exp
from MP.lib import mpl_lib
wide_fig=mpl_lib.wide_fig
fancy_legend = mpl_lib.fancy_legend

def main_reader(path='../dat/23JUL12', fref='Bsteel_BB_00.txt',
         fn_sf ='YJ_Bsteel_BB.sff',icheck=True):
    EXP, SF, IG = rs_exp.read_main(path,fref,fn_sf,icheck)

    ## interpolate based on experimental phis, psis, IG

    # EXP.plot(istps=[0,10,20])
    # SF.plot()
    # IG.plot()

    print
    print '#-----------------------------------------------------#'
    print ' Interpolate SF and IG for matching with D-spacings'

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

    return EXP, SF, IG

class StressAnalysis:
    def __init__(self,
                 path='/Users/yj/repo/rs_pack/dat/23Jul12',
                 fref='Bsteel_BB_00.txt',fn_sf='YJ_Bsteel_BB.sff'):
        self.EXP,self.SF,self.IG=main_reader(path,fref,fn_sf,icheck=True)
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
        figs = wide_fig(ifig=3031,nw=nphi)
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

    def find_sigma(self,ivo=None,istp=0,iplot=False):
        from scipy import optimize
        fmin = optimize.leastsq
        self.istp=istp

        print '#----------------------#'
        print ' epsilon_{VM} %4.2f'%self.EXP.flow.epsilon_vm[istp]
        print ' ivo:', ivo
        # --core
        # dat = fmin(self.f_least_Ei,[100,100,0,0,0,0],ivo,
        #            full_output=True)
        dat = fmin(self.f_least_Ei_d0,[1.704,100,100,0,0,0,0],ivo,
                   full_output=True)

        # --core
        stress = dat[0][1:]
        d0 = dat[0][0]

        #print 'stress:', stress
        print 'd0:', d0

        print '#----------------------#'
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
                figs.axes[iphi].plot(x,ehkl[iphi,:],'-x',color=l.get_color())

        stress = np.array(stress) / 1e6
        return stress, d0

def main(path='/Users/yj/repo/rs_pack/dat/23Jul12',
         fref='Bsteel_BB_00.txt',fn_sf='YJ_Bsteel_BB.sff'):
    from MP.mat import mech
    from MP.lib import axes_label
    eqv = axes_label.__eqv__

    figs = wide_fig(nw=2)
    mystress = StressAnalysis(path=path,fref=fref,fn_sf=fn_sf)
    mystress.nstp

    mystress.EXP.plot_all()

    eps_vm = mystress.EXP.flow.epsilon_vm
    dknot = []
    s11 = []
    s22 = []
    for istp in range(mystress.nstp):
        stress, d0 = mystress.find_sigma(ivo=[0,1],
                            istp=istp,iplot=False)
        dknot.append(d0)
        s11.append(stress[0])
        s22.append(stress[1])

    flow = mech.FlowCurve()
    flow.get_stress(s11,0,0)
    flow.get_stress(s22,1,1)

    flow.set_zero_sigma_ij(i=2,j=2)
    flow.set_zero_shear_stress()
    flow.get_vm_stress()

    figs.axes[0].plot(eps_vm, flow.sigma_vm,'x')
    figs.axes[1].plot(eps_vm, dknot,'o')

    eqv(figs.axes[0],ft=8)

    return mystress
