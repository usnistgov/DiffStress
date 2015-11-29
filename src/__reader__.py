"""
Reads X-Ray file and returns data in an organized way 
for further analysis
"""

def main(fn='24JUL12_0022Data.txt', mode='tr'):

    if mode=='tr':
        s, e = __readtr__(fn)
        return s, e
    elif mode=='lin':
        detector1, detector2 = __readlin__(fn)
        return detector1, detector2
    raise IOError, 'Unexpected flow.'

def __readlin__(fn='24JUL12_0004Data1Phi-90.txt'):
    """ Reads psi sine square vs d-spacing file """
    import numpy as np
    cnts = open(fn, 'r').readlines()
    phi = float(cnts[3].split('=')[1]) # phi
    tmp = cnts[4].split('=')[1] 
    s = float(tmp.split()[0])   # stress
    d = float(tmp.split()[2])   # deviation
    s1 = cnts[12:22]
    s2 = cnts[24:34]

    beta, sin2psi, dspacing, th2, eps, intensity\
        = _detector_(s1)
    det1 = dict(beta=beta, sin2psi=sin2psi, \
                    dspacing=dspacing, th2=th2, \
                    eps=eps, intensity=intensity)

    beta, sin2psi, dspacing, th2, eps, intensity \
        = _detector_(s2)
    det2 = dict(beta=beta, sin2psi=sin2psi, \
                    dspacing=dspacing, th2=th2,\
                    eps=eps, intensity=intensity)

    return det1, det2

def _detector_(s):
    """ Reads the detector block and returns proper
    data """
    beta = s[0].split()[2:]
    npoints = len(beta)
    psi = map(float,s[1].split(':')[1].split())
    sin2psi = map(float, s[2].split(':')[1].split())
    dspacing = map(float, s[3].split(':')[1].split())
    th2 = map(float, s[4].split(':')[1].split())
    strainE3 = map(float, s[5].split(':')[1].split())
    fwhm = map(float, s[6].split(':')[1].split())
    fwhmfit = map(float, s[7].split(':')[1].split())
    breadth = map(float, s[8].split(':')[1].split())
    intensity = map(float, s[9].split(':')[1].split())

    print 'beta',
    print beta
    print 'psi'
    print psi
    print 'sin2psi'
    print sin2psi
    print 'dspacing'
    print dspacing
    print 'th2'
    print th2

    return beta, sin2psi, dspacing, th2, strainE3, \
        intensity
    
        
def __readtr__(fn):
    import numpy as np
    cnts = open(fn, 'r').readlines()
    #print cnts[0].split('(')[1].split(')')[0]
    s11, s12, s13, dum, e11, e12, e13 = cnts[6].split()
    s22, s23, dum, e22, e23 = cnts[7].split()
    s33, dum, e33 = cnts[8].split()

    s11, s12, s13 = map(float, [s11, s12, s13])
    s22, s23      = map(float, [s22, s23])
    s33           = float(s33)

    e11, e12, e13 = map(float, [e11, e12, e13])
    e22, e23      = map(float, [e22, e23])
    e33           = float(e33)

    # print '%f %f %f '%(s11, s12, s13)
    # print '%f %f '%(s22, s23)
    # print '%f'%s33
    # print '%f %f %f '%(e11, e12, e13)
    # print '%f %f '%(e22, e23)
    # print '%f'%e33

    ## To symmetry tensor.
    s = [[ s11, s12, s13],
         [ s12, s22, s23],
         [ s13, s23, s33]]
    e = [[ e11, e12, e13],
         [ e12, e22, e23],
         [ e13, e23, e33]]

    return np.array(s), np.array(e)
    
def __logreader__(fn='log', loc=None):
    """ log file reader """
    import os
    # if location is not given, fn = fn
    # otherwise, the full-path filename is given.
    if loc==None: fn = fn
    else: fn = os.path.join(loc, fn)

    cnt = open(fn, 'r').readlines()
    cnt = cnt[2:]
    filenames = []
    ln_fn     = []
    images    = []

    for i in xrange(len(cnt)):

        if len(cnt[i])<3: pass
        else:
            dum = cnt[i].split()[0]
            dum = dum + 'Data.txt'
            filenames.append(dum)
            dum = cnt[i].split()[0]

            fn_phis = []
            fn_phis.append(dum + 'Data1Phi-90.txt')
            fn_phis.append(dum + 'Data2Phi0.txt')
            fn_phis.append(dum + 'Data3Phi45.txt')
            fn_phis.append(dum + 'Data4Phi135.txt')
            
            ln_fn.append([])
            for j in xrange(len(fn_phis)):
                ln_fn[i].append(fn_phis[j])

            images.append([])

            images[i].append(int(cnt[i].split(',')[0].split()[1]))
            if len(cnt[i].split(','))==2:
                images[i].append(int(cnt[i].split(',')[1]))
            elif len(cnt[i].split(','))==1: pass
            else: raise IOError, 'Unexpected case...'

    return filenames, ln_fn, images
