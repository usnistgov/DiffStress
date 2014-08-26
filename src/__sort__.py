"""
Sorting the files in the directory based on the data labels
of the file names.
"""
from os import sep
def __findalltr__(path='.'):
    """ Finds all triaxial data files and returned the file names """
    import glob
    fn = glob.glob('%s%s*.txt'%(path,sep))
    fntr = [] # triaxial files only.
    for i in range(len(fn)):
        cfn = fn[i]
        cfn = cfn.split('.txt')[0]
#        print cfn[-4:]
        if cfn[-4:]=='Data':
            fntr.append(fn[i])
    return fntr
    
def __alltr_ids__(path='.'):
    """ Returns all triaxial data file's ids..."""
    fntr = __findalltr__(path=path)
    ids = []
    for i in range(len(fntr)):
        cfn = fntr[i]
        iid = int(cfn.split('_')[-1].split('Data')[0])
        ids.append(iid)
    return ids
    
def __findlin__(n):
    """
    Finds a fixed-phi scan
    """
    import glob
    fn = glob.glob('*_%sData?Phi*.txt'%str(n).zfill(4))
    return fn
    
def main(mode='tr',i=None,path='.'):
    fntr = __findalltr__(path=path)
    ids = __alltr_ids__(path=path)
    if mode=='tr':
        "Triaxial results"
        return fntr, ids
    elif mode=='lin':
        if i==None or type(i).__name__!='int':
            print i
            raise IOError, 'The argument i should'\
                ' be an integer.'
        elif not(i in ids):
            raise IOError, \
                'The given i does not exist.'
        import glob
        fnlins = glob.glob('*%sData?Phi*.txt'\
                           %str(i).zfill(4))
        return fnlins

        
