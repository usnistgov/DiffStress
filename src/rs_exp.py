"""
experimental residual stress analysis
"""
def check(EXP=None):
    from MP.lib import mpl_lib
    from phikhi import draw_circ
    draw_circ = draw_circ.main
    wide_fig = mpl_lib.wide_fig    

    figs = wide_fig(nw=2)

    if EXP==None:EXP=main()
    EXP.list()
    d_avg = []

    for istp in range(len(EXP.P_scan)):
        d_avg.append(EXP.P_scan[istp].d_avg)

    figs.axes[0].plot(-EXP.flow.epsilon[2,2],d_avg,'x')

    pf = figs.axes[1]
    istp = -1
    PF_scans = EXP.P_scan[istp]
    counts = []
    X = []; Y = []
    for iphi in range(len(PF_scans.protophi)):
        pp = PF_scans.protophi[iphi].ppscans
        for i in range(len(pp)):
            counts.append(pp[i].ints)
            X.append(pp[i].x)
            Y.append(pp[i].y)

    mxc = max(counts)
    mnc = min(counts)
    k = 0
    for iphi in range(len(PF_scans.protophi)):
        pp = PF_scans.protophi[iphi].ppscans
        for i in range(len(pp)):
            alpha = (counts[k]-mnc)/(mxc-mnc)
            pf.plot(X[k],Y[k],'ko',ms=4,
                    alpha=alpha)
            k = k + 1
    pf.set_xlim(-1,1)
    pf.set_ylim(-1,1)
    draw_circ(pf)
    pf.set_aspect('equal')
    pf.set_axis_off()

def read_main(path='../dat/23JUL12',fref='Bsteel_BB_00.txt',
              fn_sf='YJ_Bsteel_BB.sff',
              icheck=False):
    import read_proto, os
    ExpProto = read_proto.ProtoExp(path,fref)
    if icheck:check(ExpProto)

    if path==None: path='.'
    SF, IG = read_IGSF(fn='%s%s%s'%(path,os.sep,fn_sf))
    return ExpProto, SF, IG

def read_IGSF(fn='YJ_Bsteel_BB.sff'):
    """
    Read from the stress factor file
    """
    import numpy as np
    import RS; reload(RS)
    from RS import sff_plot
    from RS import sfig_class

    StressFactor = sfig_class.SF
    IGstrain = sfig_class.IG
    read_sff = sff_plot.reader

    SF = StressFactor()  # class
    IG = IGstrain()      # class

    if fn!=None:
        # sf[nstp,nphi,npsi,6]
        sf, ig, phi, psi, eps_macro = read_sff(fn)
        sf = sf * 1e-12 # Pascal
        eps_macro = np.array(eps_macro)
        SF.add_data(sf=sf,phi=phi,psi=psi)
        IG.add_data(ig=ig,phi=phi,psi=psi) 

        # ## assign the flow
        SF.add_flow(eps_macro,0,0)
        SF.add_flow(eps_macro,1,1)
        SF.add_flow(-eps_macro*2,2,2)
        #SF.flow.set_zero_epsilon_ij(2,2)
        SF.flow.set_zero_shear_strain()

        IG.add_flow(eps_macro,0,0)
        IG.add_flow(eps_macro,1,1)
        IG.add_flow(-eps_macro*2,2,2)
        #IG.flow.set_zero_epsilon_ij(2,2)
        IG.flow.set_zero_shear_strain()

    return SF, IG
