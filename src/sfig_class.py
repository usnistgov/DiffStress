# class stress factor
# class intergranular strain
from MP.mat import mech
import numpy as np
fc = mech.FlowCurve
sin = np.sin
cos = np.cos
from MP.lib import mpl_lib
wide_fig = mpl_lib.wide_fig
from RS import rs
interp = rs.interpolate

class SF:
    def add_data(self,sf=None,phi=None,psi=None):
        """
        sf[nstp,nphi,npsi,6]
        """
        self.flow = fc()
        if sf!=None:
            self.sf  = sf
            self.sf.shape
            # if len(self.sf.shape)!=4:
            #     self.sf = np.array([self.sf])
            self.nstp = len(self.sf[0])

            self.nij = self.sf.shape[-1]
            if self.nij!=6:
                print 'Stress factor is not fully given in 6D'
        if phi!=None:
            self.phi = phi
            self.nphi = len(self.phi)
        if psi!=None:
            self.psi = psi
            self.npsi = len(self.psi)
    def add_flow(self,eps,i,j):
        self.flow.get_strain(eps,i,j)

    def get_iso_sf(self,E,nu):
        self.add_iso_sf(E=E,nu=nu)
        sf_new = np.zeros(self.sf.shape)
        sf_new[:] = self.sf_iso.copy()
        self.sf_old = self.sf.copy()
        self.sf = sf_new.copy()

    def restore_sf(self):
        self.sf = self.sf_old.copy()

    def add_iso_sf(self,E,nu):
        E  = 204e9
        nu = 0.3
        phi = self.phi.copy()
        psi = self.psi.copy()
        self.sf_iso = np.zeros((len(phi),len(psi),6))
        for iphi in range(len(phi)):
            for ipsi in range(len(psi)):
                self.sf_iso[iphi,ipsi] = calc_iso_sf(
                    phi=phi[iphi],psi=psi[ipsi],
                    nu=nu,Y=E,iopt=1).copy()

    def interp_strain(self,epsilon_vm):
        self.sf_old = self.sf.copy()
        self.flow.get_vm_strain()
        self.flow.epsilon_vm #
        self.sf_new = np.zeros(
            (len(epsilon_vm),self.nphi,self.npsi,self.nij))
        for iphi in range(self.nphi):
            for ipsi in range(self.npsi):
                for k in range(self.nij):
                    y = self.sf_old[:,iphi,ipsi,k]
                    self.sf_new[:,iphi,ipsi,k] = interp(
                        epsilon_vm,self.flow.epsilon_vm,y)

        self.sf = self.sf_new
        # Overwrite the flow? self.nstp also needs change
        self.flow = fc()
        self.flow.epsilon_vm = np.array(epsilon_vm)
        self.nstp = len(self.flow.epsilon_vm)

    def interp_psi(self,psi=None,iopt=1):
        """
        iopt=0: interpolate along psi
        iopt=1: interpolate along +- sin2psi
        """
        psi = np.array(psi)
        from numpy import sign
        if iopt==0:
            x = psi
            x0 = self.psi
        if iopt==1:
            x = sign(psi) * sin(psi*np.pi/180.)**2
            x0 = sign(self.psi) * sin(self.psi*np.pi/180.)**2

        self.sf_old = self.sf.copy()
        self.psi_old = self.psi.copy()
        self.npsi # also need change

        self.sf_new = np.zeros(
            (self.nstp,self.nphi,len(x),self.nij))
        for istp in range(self.nstp):
            for iphi in range(self.nphi):
                for k in range(self.nij):
                    y = self.sf_old[istp,iphi,:,k]
                    self.sf_new[istp,iphi,:,k] = interp(
                        x,x0,y)

        self.sf = self.sf_new.copy()
        self.psi = psi.copy()
        self.npsi = len(self.psi)

    def change_phi(self,phi):
        """
        """
        self.sf_old = self.sf.copy()
        self.phi_old = self.phi.copy()
        self.nphi # also need change
        pass

    def determine_phis(self,phi_new):
        deg = 180./np.pi
        nphi = self.nphi
        npsi = self.npsi
        nstp = self.flow.nstp
        phi_old = self.phi.tolist()
        nphi_new = len(phi_new)
        ind = []

        for i in range(len(phi_new)):
            try:
                j = phi_old.index(phi_new[i])
            except ValueError:
                #print i
                #print 'Mirror?'
                #print 'Is this okay to set -phi = phi?'
                #print 'Warning!!!> Go ahead only',
                #print ' if you know what you are getting to'
                j = phi_old.index(abs(phi_new[i]))
                ind.append(j)
            else: ind.append(j)

        sf_new = np.zeros((
            self.nstp,nphi_new,self.npsi,self.nij))

        for i in range(len(ind)):
            sf_new[:,i,:,:] = self.sf[:,ind[i],:,:]
        self.sf = sf_new.copy()
        self.phi = np.array(phi_new)
        self.nphi = len(self.phi)
        
    def plot(self):
        from MP.lib import axes_label
        from MP.lib import mpl_lib
        import matplotlib as mpl
        import matplotlib.cm as cm

        figs = wide_fig(nw=self.nphi,w0=0,w1=0,
                        left=0.2,right=0.15)

        mx = max(self.flow.epsilon_vm)
        mn = min(self.flow.epsilon_vm)

        cmap, c = mpl_lib.norm_cmap(mx=mx,mn=mn)
        colors=[]
        self.flow.nstp = len(self.flow.epsilon_vm)
        for i in range(self.nphi):
            for j in range(self.flow.nstp):
                eps = self.flow.epsilon_vm[j] 
                cl = c.to_rgba(eps)
                if i==0: colors.append(cl)

                l, = figs.axes[i].plot(
                    sin(self.psi*np.pi/180.)**2,
                    self.sf[j,i,:,0],'-',color=cl)
                figs.axes[i].plot(
                    sin(self.psi*np.pi/180.)**2,
                    self.sf[j,i,:,1],'--',color=cl)

                if j==0:
                    figs.axes[i].set_title(
                        r'$\phi: %3.1f^\circ{}$'%self.phi[i])

        
        deco = axes_label.__deco__
        rm_inner =mpl_lib.rm_inner
        ticks_bin_u = mpl_lib.ticks_bins_ax_u

        deco(figs.axes[0],iopt=1)
        mpl_lib.tune_xy_lim(figs.axes)
        rm_inner(figs.axes)
        ticks_bin_u(figs.axes,n=4)

        # color bar
        b = figs.axes[-1].get_position()
        axcb = figs.add_axes([0.88,b.y0,0.03,b.y1-b.y0])
        mpl_lib.add_cb(ax=axcb,filled=False,
                       levels=self.flow.epsilon_vm,
                       colors=colors,ylab='Equivalent strain')

class IG:
    def add_data(self,ig,phi,psi):
        self.flow = fc()
        #self.ig.shape(nstp,nphi,npsi)
        self.ig = ig
        self.phi=phi
        self.nphi = len(self.phi)
        self.psi=psi
        self.npsi = len(self.psi)
        self.nstp = self.ig.shape[0]
    def add_flow(self,eps,i,j):
        self.flow.get_strain(eps,i,j)
    def interp_strain(self,epsilon_vm):
        self.ig_old = self.ig.copy()
        self.flow.get_vm_strain()
        self.flow.epsilon_vm #
        self.ig_new = np.zeros(
            (len(epsilon_vm),self.nphi,self.npsi))
        for iphi in range(self.nphi):
            for ipsi in range(self.npsi):
                y = self.ig_old[:,iphi,ipsi]
                self.ig_new[:,iphi,ipsi] = interp(
                    epsilon_vm,self.flow.epsilon_vm,y)

        self.ig = self.ig_new
        # Overwrite the flow? self.nstp also needs change
        self.flow = fc()
        self.flow.epsilon_vm = np.array(epsilon_vm)
        self.nstp = len(self.flow.epsilon_vm)
    def interp_psi(self,psi=None,iopt=1):
        """
        iopt=0: interpolate along psi
        iopt=1: interpolate along +- sin2psi
        """
        psi = np.array(psi)
        from numpy import sign
        if iopt==0:
            x = psi
            x0 = self.psi
        if iopt==1:
            x = sign(psi) * sin(psi*np.pi/180.)**2
            x0 = sign(self.psi) * sin(self.psi*np.pi/180.)**2

        self.ig_old = self.ig.copy()
        self.psi_old = self.psi.copy()
        self.npsi # also need change

        self.ig_new = np.zeros(
            (self.nstp,self.nphi,len(x)))
        for istp in range(self.nstp):
            for iphi in range(self.nphi):
                y = self.ig_old[istp,iphi,:]
                self.ig_new[istp,iphi,:] = interp(
                    x,x0,y)

        self.ig = self.ig_new.copy()
        self.psi = psi.copy()
        self.npsi = len(self.psi)

    def determine_phis(self,phi_new):
        deg = 180./np.pi
        nphi = self.nphi
        npsi = self.npsi
        nstp = self.flow.nstp
        phi_old = self.phi.tolist()
        nphi_new = len(phi_new)
        ind = []

        for i in range(len(phi_new)):
            try:
                j = phi_old.index(phi_new[i])
            except ValueError:
                # print i
                # print 'Mirror?'
                # print 'Is this okay to set -phi = phi?'
                # print 'Warning!!!> Go ahead only',
                # print ' if you know what you are getting to'
                j = phi_old.index(abs(phi_new[i]))

                ind.append(j)
            else: ind.append(j)

        ig_new = np.zeros((
            self.nstp,nphi_new,self.npsi))

        for i in range(len(ind)):
            ig_new[:,i,:] = self.ig[:,ind[i],:]
        self.ig = ig_new.copy()
        self.phi = np.array(phi_new)
        self.nphi = len(self.phi)
        
    def plot(self):
        from MP.lib import axes_label
        from MP.lib import mpl_lib
        import matplotlib as mpl
        import matplotlib.cm as cm

        figs = wide_fig(nw=self.nphi,w0=0,w1=0,
                        left=0.2,right=0.15)

        mx = max(self.flow.epsilon_vm)
        mn = min(self.flow.epsilon_vm)

        cmap, c = mpl_lib.norm_cmap(mx=mx,mn=mn)
        colors=[]
        for i in range(self.nphi):
            for j in range(self.nstp):
                eps = self.flow.epsilon_vm[j] 
                cl = c.to_rgba(eps)
                if i==0: colors.append(cl)

                figs.axes[i].plot(
                    sin(self.psi*np.pi/180.)**2,
                    self.ig[j,i,:],'-x',color=cl)

                if j==0:
                    figs.axes[i].set_title(
                        r'$\phi: %3.1f^\circ{}$'%self.phi[i])

        deco = axes_label.__deco__
        rm_inner =mpl_lib.rm_inner
        ticks_bin_u = mpl_lib.ticks_bins_ax_u

        deco(figs.axes[0],iopt=2)
        mpl_lib.tune_xy_lim(figs.axes)
        rm_inner(figs.axes)
        ticks_bin_u(figs.axes,n=4)

        # color bar
        b = figs.axes[-1].get_position()
        axcb = figs.add_axes([0.88,b.y0,0.03,b.y1-b.y0])
        mpl_lib.add_cb(ax=axcb,filled=False,
                       levels=self.flow.epsilon_vm,
                       colors=colors,ylab='Equivalent strain')

def calc_iso_sf(phi,psi,nu,Y,iopt=0):
    """
    Calculate the isotropic elastic stress factor
    iopt=0: return sf[3,3]
    else:   return sf[6]
    """
    phi = phi * np.pi / 180.
    psi = psi * np.pi / 180.

    sin = np.sin
    cos = np.cos

    sf = np.zeros((3,3))
    sf6 = np.zeros((6,))
    nuy = (1. + nu) / Y

    S1 = -nu/Y
    S2 = 2*(1+nu)/Y

    # sf[0,0] = nuy * (sin(psi)**2) * (cos(phi)**2) - nu / Y
    # sf[1,1] = nuy * (sin(psi)**2) * (sin(phi)**2) - nu / Y
    # sf[2,2] = nuy * (cos(psi)**2)                 - nu / Y
    # sf[0,1] = nuy * (sin(psi)**2) * sin(2*phi)
    # sf[0,2] = nuy * sin(2*psi)    * cos(phi)
    # sf[1,2] = nuy * sin(2*psi)    * sin(phi)

    # 2011/JAC/Gnaeupel-Herold, Creuziger, and Iadicola
    sf[0,0] = S1  + 0.5 * S2 * (cos(phi)**2) * (sin(psi)**2)
    sf[1,1] = S1  + 0.5 * S2 * (sin(phi)**2) * (sin(psi)**2)
    sf[2,2] = S1  + 0.5 * S2 *                 (cos(psi)**2)
    sf[0,1] = 0.5 * 0.5 * S2 *  sin(2*phi)   * (sin(psi)**2)
    sf[0,2] = 0.5 * 0.5 * S2 *  cos(phi)     *  sin(2*psi)
    sf[1,2] = 0.5 * 0.5 * S2 *  sin(phi)     *  sin(2*psi)

    # Apply symmetry
    sf[1,0] = sf[0,1]
    sf[2,0] = sf[0,2]
    sf[1,2] = sf[2,1]

    if iopt==0:
        return sf

    from MP.mat import voigt
    vij = voigt.vij
    for i in range(3):
        for j in range(3):
            k = vij[i,j]
            sf6[k] = sf[i,j]

    return sf6
