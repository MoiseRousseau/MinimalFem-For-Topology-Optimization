import numpy as np
from scipy.sparse import coo_matrix
import struct

def read_jacobian(prefix):
    filename = prefix + "_jacobian.bin"
    return read_matrix(filename)


def read_sensitivity(prefix):
    filename = prefix + "_sensitivity.bin"
    return read_matrix(filename)


def read_matrix(filename):
    f = open(filename, 'rb')
    rows = struct.unpack('q', f.read(8))[0]
    cols = struct.unpack('q', f.read(8))[0]
    nnzs = struct.unpack('q', f.read(8))[0]
    outS = struct.unpack('q', f.read(8))[0]
    innS = struct.unpack('q', f.read(8))[0]
      
    val = np.zeros(nnzs, 'f8')
    for count in range(nnzs):
      	val[count] = struct.unpack('d', f.read(8))[0]
    i = np.zeros(outS, 'i8')
    j = np.zeros(nnzs, 'i8')
    for count in range(outS):
      	i[count] = struct.unpack('I', f.read(4))[0]
    for count in range(nnzs):
      	j[count] = struct.unpack('I', f.read(4))[0]
    ii = np.zeros(nnzs, 'i8')
    count = 0
    for k in range(outS-1):
      	ii[i[count]:i[count+1]] = count
      	count += 1
    ii[i[count]:] = count
    coo_mat = coo_matrix((val, (j, ii)), shape=(rows,cols))
    
    return coo_mat
