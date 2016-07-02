#!/usr/bin/env python

import sys
from heapq import heappush, heappop, heapify

import numpy
import scipy.sparse.linalg as la

EPSILON = 1e-8
scalar = numpy.float64

def write(s):
    sys.stdout.write(str(s)+' ')
    sys.stdout.flush()


class Op(la.LinearOperator):

    def __init__(self, n):
        self.n = n
        la.LinearOperator.__init__(self, 
            (n, n),
            dtype=scalar, 
            matvec=self.matvec)

    def __len__(self):
        return self.n

    def matvec(self, v):
        n = len(v)
        w = numpy.zeros(n)
        for i in range(n):
            w[i] = (i%2)*v[i]
        return w





class XOp(Op):
    def __init__(self, n, idxs):
        " X operator to specified qubits "
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
        " Z operator to specified qubits "
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


from smap import SMap
def strop(l, idx):
    keys = [(i, j) for i in range(l) for j in range(l)]
    s = SMap()
    bits = []
    for i in range(l):
      for j in range(l):
        s[i, j] = 'X' if idx%2 else '.'
        idx //= 2
    return str(s)
    

def build_compass(l):

    print "build_compass...",
    n = l**2

    keys = [(i, j) for i in range(l) for j in range(l)]
    coords = {}  
    for i, j in keys:
        for di in range(-l, l+1):
          for dj in range(-l, l+1):
            coords[i+di, j+dj] = keys.index(((i+di)%l, (j+dj)%l))

    m = n 
    xops = []
    zops = []

    idx = 0 
    for i in range(l):
        for j in range(l):
            op = XOp(n, [coords[i, j], coords[i, j+1]])
            xops.append(op)

            op = ZOp(n, [coords[i, j], coords[i+1, j]])
            zops.append(op)

            idx += 1
            write(idx)

    assert idx == m
    assert len(xops) == len(zops) == n

    gauge = SumOp(n, xops+zops)
    print "done"

    stabs = []
    for i in range(l-1):
        idxs = []
        for j in range(l):
            idxs.append(coords[j, i])
            idxs.append(coords[j, i+1])
        op = XOp(n, idxs)
        stabs.append(op)

    return gauge, xops, zops, stabs




def test():

    l = argv.get("l", 3)
    A, xops, zops, stabs = build_compass(l)

    n = A.n
    v0 = None
    v0 = numpy.zeros(len(A))
    v0[0] = 1
    
    
    vals, vecs = la.eigsh(A, k=2, v0=v0, which='LA', maxiter=1000) #, tol=1e-8)
    print

    print vals
    print "iters:", A.count

    v0 = vecs[:, -1]
    v0 = numpy.abs(v0)

    count = 0
    idxs = []
    values = []
    for i in range(n):
        if abs(v0[i]) > EPSILON:
            #write(i)
            idxs.append(i)
            values.append(v0[i])
    print
    print "nnz:", len(idxs)

    propagate = SumOp(l**2, xops)

    syndrome = numpy.zeros(n)
    for op in zops:
        syndrome += op.phases

    idxs = [i for i in range(n) if v0[i]>EPSILON]
    idxs.sort(key = lambda i : -syndrome[i])

    value = v0[idxs[0]]
    best = []
    for idx in idxs:
        if v0[idx] >= value-EPSILON:
            best.append(idx)
        else:
            break

    print "best:", best
    marked = set(best)
    for i in idxs:
        for xop in xops:
            j = xop.perm[i]
            if j in marked:
                continue
            marked.add(j)
            continue # <---------
            if syndrome[j] > syndrome[i]:
                #print "*", i, '->', j, 'and', syndrome[i], syndrome[j]
                print strop(l, i), syndrome[i]
                print "->"
                print strop(l, j), syndrome[j]
                print v0[i], "->", v0[j]
                print

    for i in idxs:
        for xop in xops:
            j = xop.perm[i]
            if syndrome[j] > syndrome[i] and v0[j]<v0[i]:
                #print "*", i, '->', j, 'and', syndrome[i], syndrome[j]
                print strop(l, i), syndrome[i]
                print "->"
                print strop(l, j), syndrome[j]
                print v0[i], "->", v0[j]
                print



from argv import Argv 
argv = Argv()

if __name__ == "__main__":

    test()




