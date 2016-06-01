## a module that aims to quickly load igstrain_fbulk_ph1.out
import numpy as np

def reader_ehkl(fn='int_els_ph1.out',isort=True):
    """
    Arguments
    ---------
    fn='int_els_ph1.out'
    isort=True


    orig         new
    0 istep       0
    1 phi         1
    2 beta
    3 psi1        2
    4 psi2
    5 eps1        3
    6 eps2
    7 sig1        4
    8 sig2
    9 n1          5
    10 n2
    11 v1         6
    12 v2
    13 sig11      7
    14 sig22      8
    15 sig33      9
    16 sig23     10
    17 sig13     11
    18 sig12     12
    19 eps11     13
    20 eps22     14
    21 eps33     15
    22 eps23     16
    23 eps13     17
    24 eps12     18
    """
    import time
    t0 = time.time()

    with open(fn,'r') as fo:
        fo.readline()
        npb = int(fo.readline().split()[0])

    dat = np.loadtxt(fn,skiprows=2) ## nrow, ncol

    steps = np.unique(dat[:,0])

    nrow, ncol = dat.shape
    nstp = nrow/npb
    dat_a = np.zeros((nstp,npb,ncol))
    for istp in xrange(nstp):
        i0 = istp*npb
        i1 = (istp+1)*npb
        dat_a[istp,:,:] = dat[i0:i1,:]

    phis  = np.unique(dat_a[0,:,1])
    betas = np.unique(dat_a[0,:,2])
    psis1 = np.unique(dat_a[0,:,3])
    psis2 = np.unique(dat_a[0,:,4])
    psis = np.append(psis1,psis2)
    nphis = len(phis)
    nbetas = len(betas)
    npsis =  nbetas*2

    dat_b=np.zeros((nstp,nphis,nbetas,ncol))
    for istp in xrange(nstp):
        for iphi in xrange(nphis):
            i0 = iphi*nbetas
            i1 = i0 + nbetas
            dat_b[istp,iphi,:,:] = dat_a[istp,i0:i1,:]


    dat_c = np.zeros((nstp,nphis,npsis,ncol-6))

    for ibet in xrange(nbetas):
        ## isteps
        icol_in_b=0
        icol_in_c=0
        dat_c[:,:,ibet*2, icol_in_c]   = dat_b[:,:,ibet,icol_in_b]
        dat_c[:,:,ibet*2+1,icol_in_c] = dat_b[:,:,ibet,icol_in_b]

        ## phis
        icol_in_b=1
        icol_in_c=1
        dat_c[:,:,ibet*2,icol_in_c]   = dat_b[:,:,ibet,icol_in_b]
        dat_c[:,:,ibet*2+1,icol_in_c] = dat_b[:,:,ibet,icol_in_b]

        ## psis
        icol_in_b0=3
        icol_in_c=2
        dat_c[:,:,ibet*2,icol_in_c]   = dat_b[:,:,ibet,icol_in_b0]
        dat_c[:,:,ibet*2+1,icol_in_c] = dat_b[:,:,ibet,icol_in_b0+1]

        ## ehkl
        icol_in_b0=5
        icol_in_c=3
        dat_c[:,:,ibet*2,icol_in_c]   = dat_b[:,:,ibet,icol_in_b0]
        dat_c[:,:,ibet*2+1,icol_in_c] = dat_b[:,:,ibet,icol_in_b0+1]

        ## sig
        icol_in_b0=7
        icol_in_c=4
        dat_c[:,:,ibet*2,icol_in_c]   = dat_b[:,:,ibet,icol_in_b0]
        dat_c[:,:,ibet*2+1,icol_in_c] = dat_b[:,:,ibet,icol_in_b0+1]

        ## ngr
        icol_in_b0=9
        icol_in_c=5
        dat_c[:,:,ibet*2,icol_in_c]   = dat_b[:,:,ibet,icol_in_b0]
        dat_c[:,:,ibet*2+1,icol_in_c] = dat_b[:,:,ibet,icol_in_b0+1]

        ## vgr
        icol_in_b0=11
        icol_in_c=6
        dat_c[:,:,ibet*2,icol_in_c]   = dat_b[:,:,ibet,icol_in_b0]
        dat_c[:,:,ibet*2+1,icol_in_c] = dat_b[:,:,ibet,icol_in_b0+1]

        ## sigij
        ## not needed...

    if isort:
        psis = dat_c[0,0,:,2]
        ind  = np.argsort(psis)
        dat_c = dat_c[:,:,ind,:] ## sorted
        psis = psis[ind]

    print 'elpased time:',time.time()-t0
    return psis, phis, dat_c[:,:,:,3], dat_c[:,:,:,5],\
        dat_c[:,:,:,6], steps

def reader_ig(fn='igstrain_fbulk_ph1.out',isort=True):
    """
    Read igstrain_fbulk_ph1.out and return the data
    in ndarray in the shape of (nstp, nfij, nphi,npsis,ncol)
    col0: phi
    col1: psi
    col2: sin2psi
    col3: ehkl
    col4: E
    col5: ehkl-E
    col6: f(hkl)
    col7: F(bulk)
    col8: ehkl-f/F*E
    col9: Sij
    col10: R^2

    Arguments
    ---------
    fn='igstrain_fbulk_ph1.out'
    isort=True

    Returns
    -------
    data in the dimension of (ncol, nstp,nFij,nphis,npsis)
    each column is decrbied above.
    """
    import time
    t0 = time.time()
    with open(fn,'r') as fo:
        completeStrings= fo.read()
    blocks = completeStrings.split('--\n')
    blocks = blocks[1:]
    lines0=blocks[0].split('\n')[:-1]
    nb = len(blocks)
    ncol = len(lines0[0].split())
    nrow = len(lines0)
    print 'nb:',len(blocks)
    print 'ncol,nrow',ncol,nrow
    dat = np.zeros((nb,nrow,ncol))
    Fij = []
    for i in xrange(nb): ##  each Fij
        lines = blocks[i].split('\n')
        # lines = blocks[i].split('\n')[:-1]
        for j in xrange(nrow):
            dat[i,j,:] = map(float,lines[j].split())
        if i<2:
            Fij.append(map(int,lines[j+1].split('F')[-1].split(',')))

    print 'Fij:',Fij

    ## restructre the data...
    #dat = np.zeros((nb,nrow,ncol))
    phis = dat[0,:,0]; psis = dat[0,:,1]
    phis = np.unique(phis); psis = np.unique(psis)
    nFij,nphis,npsis= len(Fij),len(phis),len(psis)
    nstp = nb/nFij
    print 'nstp:',nstp

    #
    d1_a=np.zeros((nstp,nFij,nrow,ncol))
    d1_b=np.zeros((nstp,nFij,nphis,npsis,ncol))
    for istp in xrange(nstp):
        i0=istp*nFij
        i1=(istp+1)*nFij
        d1_a[istp,:,:,:]=dat[i0:i1,:,:]
        for ij in xrange(nFij):
            for iphi in xrange(nphis):
                #each fij-nphi block
                i0=iphi*npsis
                i1=i0+npsis
                d1_b[istp,ij,iphi,:,:]=d1_a[istp,ij,i0:i1,:] # nstp, nfij, nphis * npsi,ncols

    if isort:
        psis = d1_b[0,0,0,:,1]
        inds = np.argsort(psis)
        d1_b = d1_b[:,:,:,inds,:]
    print time.time()-t0, 'seconds'
    # return d1_b.swapaxes(0,-1) ## ncol,nstp,fij,nphi,npsi,ncol
    return d1_b.transpose(4,0,1,2,3)

if __name__=='__main__':
    ## test
    # reader()
    reader_ehkl()
