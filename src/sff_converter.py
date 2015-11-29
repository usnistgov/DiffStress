### sff file converter
import numpy as np
from nist_diff import intepsphiout, fij
def main(fn='temp.sff', difile=None, iph=1, factor=1e6, itab=False,
         ieps0=4,irev=False,fn_str=None, flow=None):
    """
    Arguments
    =========
    fn       = 'temp.sff'
    difile   = 'None'
    iph      = 1
    factor   = 1e6
    itab     = False (use '\t' as the delimiter)
    # ieps0 = 0,1,2,3 (option for ntepsphiout)
    #          0: int_eps_ph - at unloads (eps^hkl)
    #          1: int_els_ph - at loads   (eps^hkl)
    #          2: igstrain_ph - at loads  (Eps-eps^hkl)
    #          3: igstrain_ph - at unloads  (Eps-eps^hkl)
    ieps0 is deprecated. The option for IG and SF is now hardwired
             # 4: igstrain_bulk - at loads (eps^hkl - Fij<Sij>)
             # 5: igstrain_bulk - at loads (eps^hkl - Fij/Fij^bulk E)
             IG strain is calculated from igstrain_loads_avg.out

    irev    = False (if True, eps0 = -eps0) - This argument is for
             the one time. (Mistake in diffwrite subroutine...)

    * Note that stress factor file for Thomas's PF software
    takes 1/TPa unit.
    Either VPSC or EVPSC has itself no imposition of unit.
    However, single crystal file usually takes on MPa unit
    for both single crystal moduli and hardening parameters.
    So, it is most likely that the stress factor should be
    converted from 1/MPa to 1/TPa.  10^6 MPa = 1 TPa.

    There are various methods to calculate SF/IG strains in EVPSC.
    1. igstrain_loads_avg.out: this file is not
       one of many EVPSC-generated-output files. This file is
       calculated based on ...
    2. igstrain_bix_ph1.out (diff prior to *unloading*)
    3. igstrain_fbulk_ph1.out
    4. int_els_ph1.out (prior to unloading)
    5. int_eps_ph1.out (post unloading)

    """
    if not(itab): raise IOError,'use itab==True'

    if difile==None:
        print 'difile is none: find it from EVPSC.IN'
        dl     = open('EVPSC.IN','r').readlines()
        nph    = int(dl[1].split()[0])
        n      = 12*(iph-1) + 3 + 11
        difile = dl[n][:-1]
    print difile

    # condition
    difl, nphi, phis, npsi, dum1, dum2 = condition(fn=difile)

    # ------------------------------------------------------- #
    fn_ig_avg='igstrain_loads_avg.out' ## IG strain from avg.
    #fn_ig_avg='igstrain_fbulk_ph1.out'
    from pepshkl import reader3 as reader
    import os.path

    if flow!=None and fn_str!=None:
        raise IOError, 'Either flow or fn_str should be given'
    elif flow!=None: flow.get_eqv()
    elif fn_str!=None:
        from MP.mat.mech import FlowCurve
        flow = FlowCurve()
        flow.get_model(fn=fn_str)
        flow.get_eqv()

    neps    = flow.nstp
    strains = flow.epsilon_vm[::]

    if not(os.path.isfile(fn_ig_avg)):
        from pepshkl import ex_igb_bix_t
        print 'Run pepshkl.ex_igb_cf or pepshkl.ex_igb_bix_t'\
            ' to obtain igstrain_loads_avg.out'
        dummy_rst = ex_igb_bix_t(fn='igstrain_bix_ph1.out',
                                 fnout=fn_ig_avg,
                                 flow=flow)
    dat,psis,fij = reader(fn_ig_avg,isort=True)
    # ------------------------------------------------------- #

    psis = psis[0]
    print 'difile:',difile
    print 'npsi', npsi
    npsi = len(psis)
    print 'npsi', npsi

    # (step, nphi, psi)

    if ieps0<4: raise IOError, 'Use if ieps0<4 is deprecated'

    sff = open(fn, 'w') # stress factor file
    sff.write('Header')
    for i in xrange(neps*10-2): sff.write('\t')
    sff.write('\r\n')
    sff.write('#strainLevels \t %i'%neps)
    for i in xrange(neps*10-3): sff.write('\t')
    sff.write('\r\n')
    sff.write('#phivalues \t%i'%nphi)
    for i in xrange(neps*10-3): sff.write('\t')
    sff.write('\r\n')
    sff.write('#psivalues \t%i'%npsi)
    for i in xrange(neps*10-3): sff.write('\t')
    sff.write('\r\n')
    sff.write(('%s'%('e^{eff} \t%5.3f'%(strains[0]))).ljust(14))
    for i in xrange(8): sff.write('\t')

    if neps>1:
        sff.write(('\t%5.3f'%(strains[1])).rjust(115))
        for i in xrange(9): sff.write('\t')
        for i in xrange(neps-2):
            sff.write(('\t%5.3f'%(strains[i+2])).rjust(118))
            for j in xrange(9):
                if j==8 and i==neps-3: pass
                else: sff.write('\t')

    sff.write('\r\n')
    for i in xrange(neps):
        stemp = '%10s%9s ' ; ftemp = '%9s'%'%6.2f'
        ftemp = ftemp + '%8s '%'%6.2f'
        for j in xrange(7):
            stemp = stemp + '%13s '
            ftemp = ftemp + '% +12.6e '
        sff.write(stemp%('Phi\t','Psi\t','F11\t','F22\t',
                  'F33\t','F23\t','F13\t','F12\t','eps0'))
        if i==neps-1: pass
        else: sff.write('\t')
        if i==neps-1: pass
        else: sff.write('\t')


    print 'nphi,neps,npsi:', nphi,neps,npsi

    sff.write('\r\n') # line breaker
    for i in xrange(nphi):
        for j in xrange(npsi):
            for k in xrange(neps):
                ph = phis[i]
                ps = psis[j]
                sf = fij[k,i,j,:]*factor
                eps0 = dat[k,i,j]

                if irev: eps = - eps
                sff.write('%9.2f\t  %9.2f\t %+12.4e\t'\
                              ' %+12.4e\t %+12.4e\t'\
                              ' %+12.4e\t %+12.4e\t'\
                              ' %+12.4e\t %+12.4e'%(
                        ph,ps,sf[0],sf[1],sf[2],sf[3],sf[4],sf[5],eps0))
                if k==neps-1: pass
                else: sff.write('\t')
                if k==neps-1: pass
                else: sff.write('\t')
            sff.write('\r\n')
        pass
    sff.close()
    pass # end of main


# print """
# definition main
# """
# print main.__doc__

def condition(fn=None):
    """
    return conditions of the last run based on the
    dif file given as the argument fn

    Argument
    ========
    fn = None
    """
    if fn==None:
        # print 'dif filename must be given'
        EVPSC_IN = open('EVPSC.IN','r')
        lines    = EVPSC_IN.readlines()
        EVPSC_IN.close()
        nph      = int(lines[1].split()[0])
        if nph>1: raise IOError, 'Only single phase is of consideration.'
        fn       = lines[14].split()[0]

    difl = open(fn, 'r').readlines()
    nphi = int(difl[5].split()[0])
    phis = map(float,difl[6].split()[:nphi])
    nbeta= int(difl[7].split()[0])
    neps = int(difl[11].split()[0])
    eps = map(float, difl[12].split()[:neps])
    eps = eps[::-1] #; eps.append(0)
    eps = np.array(eps[::-1])
    return difl, nphi, phis, nbeta, neps, eps

