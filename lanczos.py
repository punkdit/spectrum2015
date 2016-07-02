#!/usr/bin/env python

import sys

import numpy
import scipy.sparse.linalg as la

scalar = numpy.float64


class Op(la.LinearOperator):

    def __init__(self, n):
        la.LinearOperator.__init__(self, 
            (n, n),
            dtype=scalar, 
            matvec=self.matvec)

    def matvec(self, v):
        n = len(v)
        w = numpy.zeros(n)
        for i in range(n):
            w[i] = (i%2)*v[i]
        return w





class XOp(Op):
    def __init__(self, n, idxs):
        "apply X operator to specified qubits "
        if type(idxs) in (int, long):
            idxs = [idxs]
        for idx in idxs:
            assert 0<=idx<n
        perm = []
        assert n<=30
        ns = range(n)
        mask = [1 if i in idxs else 0 for i in ns]
        #print "XOp:", mask
        for i in range(2**n):
            #print i,
            bits = []
            for flip in mask:
                bit = (1 - i%2) if flip else i%2
                bits.append(bit)
                i >>= 1
            #print bits,
            j = 0
            for bit in reversed(bits):
                j += bit
                j <<= 1
            j >>= 1
            #print j
            perm.append(j)
        self.perm = perm
        Op.__init__(self, 2**n)

    def matvec(self, v, w=None):
        assert len(v) == len(self.perm)
        if w is None:
            w = numpy.zeros(len(v))
        #for i, j in enumerate(self.perm):
        #    w[j] += v[i]
        v = v[self.perm]
        w += v
        return w


class ZOp(Op):
    def __init__(self, n, idxs):
        "apply Z operator to specified qubits "
        if type(idxs) in (int, long):
            idxs = [idxs]
        for idx in idxs:
            assert 0<=idx<n
        phases = []
        assert n<=30
        ns = range(n)
        mask = [1 if i in idxs else 0 for i in ns]
        #print "XOp:", mask
        for i in range(2**n):
            #print i,
            phase = 1
            bits = []
            for flip in mask:
                if flip and i%2:
                    phase *= -1
                i >>= 1
            phases.append(phase)
        #print "ZOp:", idxs, phases
        self.phases = numpy.array(phases)
        Op.__init__(self, 2**n)

    def matvec(self, v, w=None):
        assert len(v) == len(self.phases)
        if w is None:
            w = numpy.zeros(len(v))
        #for i, phase in enumerate(self.phases):
        #    w[i] += phase * v[i]
        w += self.phases * v
        return w


class SumOp(Op):
    def __init__(self, n, ops):
        self.ops = ops
        self.n = 2**n
        self.count = 0
        Op.__init__(self, self.n)

    def matvec(self, v, w=None):
        #print "SumOp.matvec"
        assert len(v)==self.n
        if w is None:
            w = numpy.zeros(len(v))
        for op in self.ops:
            op.matvec(v, w)
        self.count += 1
        sys.stdout.write('.');sys.stdout.flush()
        return w


def build_compass(l):

    print "build_compass..."
    n = l**2

    keys = [(i, j) for i in range(l) for j in range(l)]
    coords = {}  
    for i, j in keys:
        for di in range(-l, l+1):
          for dj in range(-l, l+1):
            coords[i+di, j+dj] = keys.index(((i+di)%l, (j+dj)%l))

    m = n 
    ops = []

    idx = 0 
    for i in range(l):
        print "build_compass: i=", i
        for j in range(l):
            ops.append(XOp(n, [coords[i, j], coords[i, j+1]]))
            ops.append(ZOp(n, [coords[i, j], coords[i+1, j]]))
            idx += 1
            print idx

    assert idx == m
    assert len(ops) == 2*n

    op = SumOp(n, ops)
    print "done"
    return op


def test():

    if 0:
        op = XOp(3, [])
        op = XOp(3, [2])
        op = XOp(3, [1,2])
    
        op = ZOp(3, [])
        op = ZOp(3, [2])
        op = ZOp(3, [1,2])
    
        X0 = XOp(3, 0)
        X1 = XOp(3, 1)
        X2 = XOp(3, 2)
    
        Z0 = ZOp(3, 0)
        Z1 = ZOp(3, 1)
        Z2 = ZOp(3, 2)
    
        A = SumOp(3, [X0, X1, X2, Z0, Z1, Z2])
    
        vals, vecs = la.eigsh(A, k=2, v0=None, which='LA')
        #print vals

    l = argv.get("l", 3)
    A = build_compass(l)
    vals, vecs = la.eigsh(A, k=4, v0=None, which='LA', maxiter=1000) #, tol=1e-8)
    print

    print vals
    print "iters:", A.count


    return

    n = 128
    A = Op((n, n))

    vals, vecs = la.eigsh(A, k=2, v0=None)

    print vals
    print vecs





from argv import Argv 
argv = Argv()

if __name__ == "__main__":

    test()




