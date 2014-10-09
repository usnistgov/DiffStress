"""
Experimental residual stress analysis
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
              fn_sf='YJ_Bsteel_BB.sff',fc=None,fn_str=None,
              icheck=False,isym=False):
    """
    Make ProtoExp, SF, IG and return them.
    """
    import read_proto, os
    ExpProto = read_proto.ProtoExp(path,fref,isym)
    if icheck:check(ExpProto)
    if path==None: path='.'
    SF, IG = read_IGSF(fn='%s%s%s'%(path,os.sep,fn_sf),
                       fc=fc,fn_str=fn_str)
    return ExpProto, SF, IG

def read_IGSF(fn='YJ_Bsteel_BB.sff', fc=None, fn_str=None):
    """
    Read from the stress factor file
    ## Discard the 'strain' recorded in the *.sff file
       if either fc or fn_str is given.

    Arguments
    =========
    fn = 'YJ_Bsteel_BB.sff' ; file name of sff file
    fc = None (should be a 'mech.FlowCurve' object)
    fn_str = 'STR_STR.OUT'
    """
    import numpy as np
    import RS; reload(RS)
    from RS import sff_plot, sfig_class

    if fc!=None and fn_str!=None: raise IOError,\
       'Only one of fc and fn_str should be given'

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

        # ## assign the flow -- only biaxial ...
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

    elif fn==None:
        print 'SFF file name is missing.',
        print ' Empty SF, IG objects are returned'
    if fc!=None:
        print "Warning: overwrite the SF.flow / IG.flow"
        SF.flow = fc; IG.flow = fc
    elif fn_str!=None:
        print "Warning: overwirte the SF.flow / IG.flow"
        from MP.mat import mech
        fc=mech.FlowCurve(name='Flow for SF/IG')
        fc.get_model(fn=fn_str)
        fc.get_eqv()
        SF.flow = fc; IG.flow = fc
    return SF, IG
