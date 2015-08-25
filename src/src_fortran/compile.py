## compilation for <cs.f>
import os, shutil

def main():
    cmd = "f2py -c -m cs cs.f mat_lib.f VOIGT.f --f77flags='-O3 -fdefault-real-8 -fdefault-double-8 -fbounds-check -finit-local-zero -finit-integer=zero -finit-real=zero'  --f90flags='-O3 -fdefault-real-8 -fdefault-double-8 -fbounds-check -finit-local-zero -finit-integer=zero -finit-real=zero' >> /tmp/dum_f2py_compile"
    iflag=os.system(cmd)
    if iflag!=0: raise IOError, 'Failed compilation'



if __name__=='__main__':
    main()
