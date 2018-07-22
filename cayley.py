#!/usr/bin/env python

from __future__ import print_function

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
            print("mulclose:", len(els))
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

    def is_identity(self):
        return self.v.sum() == 0 and self.sign==1

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

    def mkcsr(self, G, lookup):
        A = {}
        N = len(G)
        for i in range(N):
            j = lookup[G[i]*self] # symmetries act on the "right"
            #j = lookup[self*G[i]] # symmetries act on the ???
            A[i, j] = A.get((i, j), 0) + 1
        A = mkcsr(A, N)
        return A



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


def mkcsr(H, N):
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
    return H1


EPSILON = 1e-8
def is_close(A, B):
    x = (A-B).power(2).sum()
    return x < EPSILON


def mkops(Ax, flip=False):
    n = Ax.shape[1]
    pad = array2([0]*n)
    xgen = []
    for v in Ax:
        if flip:
            v = numpy.concatenate((pad, v)) # Z ops
        else:
            v = numpy.concatenate((v, pad)) # X ops
        op = Op(v, +1)
        xgen.append(op)
    return xgen

mkxops = mkops
mkzops = lambda Ax : mkops(Ax, True)

def cayley():

    model = build_model()

    print(model)
    print("Gx:")
    print(shortstr(model.Gx))
    print("Gz:")
    print(shortstr(model.Gz))

    xgen = mkxops(model.Gx)
    zgen = mkzops(model.Gz)
    gen = xgen + zgen

    n = model.n
    pad = array2([0]*n)
    v = numpy.concatenate((pad, pad))
    I = Op(v, +1)
    
    v = numpy.concatenate((pad, pad))
    nI = Op(v, -1)
    
    size = 2**(2*len(model.Rx)+1) * 2**(len(model.Hx)+len(model.Hz))
    print("size:", size)
    G = mulclose(gen, maxsize=size, verbose=False)
    #G = mulclose_fast(gen)
    assert len(G) == size

    print("G:", len(G))

    lookup = {}
    for i, op in enumerate(G):
        lookup[op] = i

    print("building H...")
    N = len(G)
    H = {}
    #for i in range(N):
    #    H[i, i] = offset
    for i in range(N):
        for g in gen:
            j = lookup[g * G[i]] # the "left" cayley graph
            assert i != j
            H[i, j] = H.get((i, j), 0) + 1
    H = mkcsr(H, N)
    I = I.mkcsr(G, lookup)

    if argv.proj:
#        for X in xgen:
#            X = X.mkcsr(G, lookup)
#            assert is_close(X.dot(X), I)
#        
#            P = 0.5 * (I + X)
#            assert is_close(X.dot(H), H.dot(X))
#            assert is_close(P, P.transpose())
#    
#            #H = P.dot(H.dot(P))
#
#        for Z in zgen:
#            Z = Z.mkcsr(G, lookup)
#            assert is_close(Z.dot(Z), I)
#        
#            P = 0.5 * (I + Z)
#            assert is_close(Z.dot(H), H.dot(Z))
#            assert is_close(P, P.transpose())
#    
#            #H = P.dot(H.dot(P))

        Rx = [X.mkcsr(G, lookup) for X in mkxops(model.Rx)]
        Rz = [Z.mkcsr(G, lookup) for Z in mkzops(model.Rz)]
        sxgen = mkxops(model.Hx)
        szgen = mkzops(model.Hz)
        stab = mulclose(sxgen + szgen)
        assert len(stab) == 2**(len(sxgen+szgen))
        stab = [S.mkcsr(G, lookup) for S in stab]
        print("stab:", len(stab))

        #nI = nI.mkcsr(G, lookup)
        #assert is_close(nI.dot(nI), I)
    
        #P = 0.5 * (I - nI)
        #assert is_close(nI.dot(H), H.dot(nI))
        #assert is_close(P, P.transpose())

        r = len(Rx)
        P = mkcsr({}, N) # zero
        count = 0

        for S in stab:
            P += S # character is +1 (up to multiple)
            count += 1
            for i in range(r):
                X = Rx[i]
                Z = Rz[i]
                #assert is_close(X.dot(Z), -Z.dot(X)) # nope
    
                A = X.dot(Z.dot(X.dot(Z))) # -I
                assert is_close(S.dot(X), X.dot(S))
                assert is_close(S.dot(Z), Z.dot(S))
                A = S.dot(A)
                P += -A # character is -1 (up to multiple)
                count += 1
                break

        P *= 1./count
        assert is_close(P, P.transpose())
        assert is_close(P.dot(P), P)
        H = P.dot(H.dot(P))

    assert is_close(H, H.transpose())
    #return

    print("building sparse...")

    offset = N+1 # make H positive definite
    H = H + offset*I # spectral shift operator

    del lookup

    if argv.dense:
        H = H.toarray()
        from numpy.linalg import eigh
        vals, vecs = eigh(H)

    else:
    
        print("eigsh...")
        k = argv.get("k", 40)
        which = argv.get("which", "LM")
        vals, vecs = sparse.linalg.eigsh(H, k=min(H.shape[0]-1, k), which=which)

    vals -= offset

    print(' '.join("%.6f"%x for x in vals))


if __name__ == "__main__":
    test()

    cayley()



