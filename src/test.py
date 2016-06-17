##
import numpy as np
def test1():
    """
    Test rs_ex.ex_consistency to schematically illustrate the
    spread of noise as a function of (psi,phi)
    """
    import rs_ex
    rst = rs_ex.ex_consistency(iscatter=True,sigma=25e-6, psimx=60.,
                               iplot=False,iwgt=False,ird=0.182, 
                               theta_b=78.2*np.pi/180.,fnPickle=None)
    print rst[-1]
