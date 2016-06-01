import numpy as np
from RS import rs_ex
# import os
# home =os.getcwd()
# os.chdir('/home/younguj/repo/DiffStress/src/')
# import rs_ex
# os.chdir(home)

def main():
    rs_ex.ex_consistency(
        psi_nbin=15, sigma=5e-5, psimx=35, ig_sub=True,
        iplot=False,ird=0.0457, theta_b=78.2*np.pi/180.,iscatter=False)

if __name__=='__main__':
    main()
