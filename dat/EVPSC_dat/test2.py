import numpy as np
from RS import rs
import time
t0=time.time()
rs.reader(fn='int_els_ph1.out',isort=True)
print 'elapsed time:', time.time()-t0
