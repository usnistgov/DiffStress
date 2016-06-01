## a module that aims to quickly load igstrain_fbulk_ph1.out
import numpy as np
def reader(fn='igstrain_fbulk_ph1.out',isort=True):
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
    reader()

