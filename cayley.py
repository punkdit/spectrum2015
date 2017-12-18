#!/usr/bin/env python

import sys, os

import math 

import numpy
from numpy import log2 
from scipy import sparse
from scipy.sparse.linalg import eigs, eigsh

from solve import shortstr, shortstrx, parse, eq2, dot2, zeros2, array2, identity2
from models import build_model
#from action import mulclose, mulclose_fast

from argv import argv


def mulclose(gen, verbose=False, maxsize=None):
    els = set(gen)
    bdy = list(els)
    changed = True 
    while bdy:
        if verbose:
            print "mulclose:", len(els)
        _bdy = []
        for A in gen:
            for B in bdy:
                C = A*B  
                if C not in els: 
                    els.add(C)
                    _bdy.append(C)
                    if maxsize and len(els)>=maxsize:
                        return list(els)
        bdy = _bdy
    return els


class Op(object):
    def __init__(self, v, sign=+1):
        self.v = array2(v)
        self.sign = sign
        n = len(self.v)
        assert n%2==0
        self.n = n//2 # qubits

    def __str__(self):
        return "Op(%s, %s)"%(self.v, self.sign)

    def __hash__(self):
        return hash((self.v.tostring(), self.sign))

    def __eq__(self, other):
        return eq2(self.v, other.v) and self.sign == other.sign

    def __ne__(self, other):
        return not eq2(self.v, other.v) or self.sign != other.sign

    def __mul__(self, other):
        assert self.n == other.n
        n = self.n
        v = self.v
        u = other.v
        vz = v[n:]
        ux = u[:n]
        sign = -1 if (vz*ux).sum()%2 else +1
        op = Op((u+v)%2, sign*self.sign*other.sign)
        return op

    def __neg__(self):
        return Op(self.v, -self.sign)



def test():

    I = Op([0, 0], +1)
    X = Op([1, 0], +1)
    Z = Op([0, 1], +1)

    assert I*I == I
    assert X*I == X
    assert X*X == I
    assert Z*Z == I
    assert X*Z == -(Z*X)
    assert X*Z == (-Z)*X

    G = mulclose([X, Z])
    assert len(G)==8

    II = Op([0, 0, 0, 0], +1)
    XI = Op([1, 0, 0, 0], +1)
    ZI = Op([0, 1, 0, 0], +1)
    IX = Op([0, 0, 1, 0], +1)
    IZ = Op([0, 0, 0, 1], +1)
    
    G = mulclose([XI, ZI, IX, IZ])
    assert len(G)==32

test()

def cayley():

    model = build_model()

    print model
    print

    gen = []
    n = model.n
    pad = array2([0]*n)
    for v in model.Gx:
        v = numpy.concatenate((v, pad))
        op = Op(v, +1)
        gen.append(op)
    
    for v in model.Gz:
        v = numpy.concatenate((pad, v))
        op = Op(v, +1)
        gen.append(op)
    
    size = 2**(2*len(model.Rx)+1) * 2**(len(model.Hx)+len(model.Hz))
    print "size:", size
    G = mulclose(gen, maxsize=size, verbose=True)
    #G = mulclose_fast(gen)
    assert len(G) == size

    print "G:", len(G)

    lookup = {}
    for i, op in enumerate(G):
        lookup[op] = i

    print "building H..."
    N = len(G)
    offset = N+1 # make H positive definite
    H = {}
    for i in range(N):
        H[i, i] = offset
    for i in range(N):
        for g in gen:
            j = lookup[g * G[i]]
            assert i != j
            H[i, j] = H.get((i, j), 0) + 1

    del lookup
    print "building sparse..."
    keys = H.keys()
    keys.sort()
    data = [] 
    rows = [] 
    cols = [] 
    for idx in keys:
        #H1[idx] = H[idx]
        data.append(H[idx])
        rows.append(idx[0])
        cols.append(idx[1])
    del H
    del keys

    H1 = sparse.coo_matrix((data, (rows, cols)), (N, N))
    H1 = sparse.csr_matrix(H1, dtype=numpy.float64)

    print "eigsh..."
    k = argv.get("k", 40)
    vals, vecs = sparse.linalg.eigsh(H1, k=min(H1.shape[0]-1, k), which="LM")
    vals -= offset

    print ' '.join("%.6f"%x for x in vals)


if __name__ == "__main__":
    test()

    cayley()



