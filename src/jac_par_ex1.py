"""
Applications for JAC  - one
"""
from RS import mst_ex

sigmas       = [1e-5, 2e-5, 5e-5, 1e-4]
ss           = 3
bounds       = [0.,0.5]
nbins        = 10
iwgt         = False
nsample      = 4
intp_opt     = 0
iplot        = False
DEC_freq_sym = True
NCPU         = 0

if __name__=='__main__':
    mst_ex.influence_of_cnts_stats()

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--sigmas', type=str,help='a list of standard deviation for eps(hkl) noise')
    parser.add_argument('--ss', type=int,help='inverse frequency in DEC')
    parser.add_argument('--bounds', type=str,help='sin2psi bounds')
    parser.add_argument('--nbins', type=int,help='Number of tilting angles (psi)')
    parser.add_argument('--nsample', type=int,help='Number of ensembles')
    parser.add_argument('--ncpu', type=int,help='Number of CPUS for multiprocessing.pool')

    args=parser.parse_args()

    sigmas  = args.sigmas
    ss      = args.ss
    bounds  = args.bounds
    nbins   = args.nbins
    nsample = args.nsample
    ncpus   = args.ncpu

    sigmas=map(float,sigmas)
    bounds=map(float,bounds)

    x,M,s = mst_ex.influence_of_cnts_stats(
        sigmas,ss,bounds,nbins,iwgt,nsample,intp_opt,iplot,DEC_freq_sym,NCPU)


    f = open('jac_ex01_dat.txt','w')
    f.write('strain ')
    for i in range(len(s)):
        f.write('%7.1e'%s[i])
        f.write('\n')
    
    f.close()
