from numpy import sin,pi,sign
import numpy as np
def sin2psi_bounds(w_rad=2.,psi0=0):
    """
    Given the 'radius' of the angular window
    Calculate psi0+psi and psi0-psi, and thus,
    sin2(psi0+psi), sin2(psi0-psi)

    All given in 'degree'
    """
    psi_l = psi0-w_rad; psi_u = psi0+w_rad

    sin2   = sin(psi0*pi/180.)**2
    sin2_l = sin(psi_l*pi/180.)**2
    sin2_u = sin(psi_u*pi/180.)**2

    # print sin2_l, sin2, sin2_u
    # print
    # print sin2 - sin2_l
    # print sin2_u - sin2

    #print psi_l, psi_u, sin2_l, sin2_u
    return psi_l, psi_u, sin2_l, sin2_u

def ex1(mx=0.5,w_rad=5.):
    import matplotlib.pyplot as plt
    ax = plt.gca()

    npsi = 10
    sin2psi = np.linspace(0,mx,npsi)

    psis = np.arcsin(np.sqrt(sin2psi)) * 180./np.pi
    xs = sin2psi[::]
    for i in range(npsi):
      x = xs[i]
      psi = psis[i]
      psil, psiu, sin2l, sin2u = sin2psi_bounds(w_rad=w_rad,psi0=psi)
      
      l, = ax.plot([x,x],[sin2l,sin2u])
      ax.plot(x,sin2psi[i],'x',color=l.get_color())
      ax.plot(x,(sin2l+sin2u)/2.,'+',color=l.get_color())

def sin2psi_opt(psi,iopt):
    """
    psi is radian

    iopt=0 (sin2psi)
    iopt=1 (+-sin2psi)
    iopt=2 (psi)
    """
    signs = sign(psi)
    sin2psi = sin(psi)**2
    if iopt==0: return sin2psi
    if iopt==1:
        for i in range(len(sin2psi)):
            sin2psi[i] = sin2psi[i] * signs[i]
        return sin2psi
    if iopt==2:
        return psi[::]
    
