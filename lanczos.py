#!/usr/bin/env python

import sys
from heapq import heappush, heappop, heapify

import numpy
import scipy.sparse.linalg as la

from code import texstr

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

    def __rmul__(self, alpha):
        return RMulOp(self, alpha)

    def matvec(self, v):
        n = len(v)
        w = numpy.zeros(n)
        for i in range(n):
            w[i] = (i%2)*v[i]
        return w



class XOp(Op):
    def __init__(self, n, idxs, alpha=1.0):
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
        self.alpha = alpha
        Op.__init__(self, 2**n)

    def matvec(self, v, w=None):
        assert len(v) == len(self.perm)
        if w is None:
            w = numpy.zeros(len(v))
        #for i, j in enumerate(self.perm):
        #    w[j] += v[i]
        v = v[self.perm]
        w += self.alpha * v
        return w


class ZOp(Op):
    def __init__(self, n, idxs, alpha=1.0):
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
        self.alpha = alpha
        Op.__init__(self, 2**n)

    def matvec(self, v, w=None):
        assert len(v) == len(self.phases)
        if w is None:
            w = numpy.zeros(len(v))
        #for i, phase in enumerate(self.phases):
        #    w[i] += phase * v[i]
        w += self.alpha * self.phases * v
        return w


def mkop(tp, s):
    n = len(s)
    idxs = []
    for i, x in enumerate(s):
        if s[i]=='1':
            idxs.append(i)
    return tp(n, idxs)


class SumOp(Op):
    def __init__(self, n, ops):
        self.ops = ops
        self.n = 2**n
        self.count = 0
        Op.__init__(self, self.n)

    def matvec(self, v, w=None, verbose=True):
        #print "SumOp.matvec"
        assert len(v)==self.n
        if w is None:
            w = numpy.zeros(len(v))
        for op in self.ops:
            op.matvec(v, w)
        self.count += 1
        if verbose:
            sys.stdout.write('.');sys.stdout.flush()
        return w


class RMulOp(Op):
    def __init__(self, op, alpha=1.0):
        self.op = op
        self.alpha = alpha
        Op.__init__(self, op.n)

    def matvec(self, v, w=None, verbose=True):
        assert len(v)==self.n
        if w is None:
            w = numpy.zeros(len(v))
        w += self.alpha * self.op.matvec(v, verbose)
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
    

class Model(object):
    pass

class CompassModel(Model):
    def __init__(self, l):
    
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
    
        xstabs = []
        for i in range(l-1):
            idxs = []
            for j in range(l):
                idxs.append(coords[j, i])
                idxs.append(coords[j, i+1])
            op = XOp(n, idxs)
            xstabs.append(op)
    
        self.n = n # qubits
        self.A = gauge
        self.xops = xops
        self.zops = zops
        self.xstabs = xstabs


gcolor_gauge = """
1111...........
11..11.........
1.1.1.1........
..11..11.......
.1.1.1.1.......
....1111.......
11......11.....
1.1.....1.1....
........1111...
..11......11...
.1.1.....1.1...
1...1...1...1..
........11..11.
.1...1...1...1.
....11......11.
........1.1.1.1
..1...1...1...1
....1.1.....1.1
""".strip().split()

gcolor_stab = """
11111111.......
1111....1111...
11..11..11..11.
1.1.1.1.1.1.1.1
""".strip().split()

class GColorModel(Model):
  def __init__(self):

    xops = []
    zops = []

    for op in gcolor_gauge:
        xops.append(mkop(XOp, op))
        zops.append(mkop(ZOp, op))

    xstabs = []
    for op in gcolor_stab:
        xstabs.append(mkop(XOp, op))
        #stabs.append(mkop(ZOp, op))

    n = len(op)
    self.A = SumOp(n, xops+zops)
    self.n = n
    self.xops = xops
    self.zops = zops

    self.xstabs = xstabs

  def get_syndrome(self, v):

    syndrome = []
    for op in self.xstabs:
        v1 = op.matvec(v)
        r = numpy.dot(v, v1) # real values
        syndrome.append(r)
    return syndrome



def test():

    if argv.compass:
        l = argv.get("l", 3)
        model = CompassModel(l)

    elif argv.gcolor:
        model = GColorModel()

    else:
        return


    dim = model.A.n

    v0 = None

    if argv.stabdist:
        v0 = numpy.zeros(len(model.A))
        v0[0] = 1
    
    k = argv.get("k", 2)

    A = model.A
    alpha = 0.0001
    perturb = [alpha * op for op in model.xstabs]
    A = SumOp(model.n, A.ops+perturb)
    
    vals, vecs = la.eigsh(A, k=k, v0=v0, which='LA', maxiter=None) #, tol=1e-8)
    print

    # vals go from smallest to highest
    print vals
    print "iterations:", model.A.count

    for i in range(k):
        print i, "syndrome", model.get_syndrome(vecs[:, i])

    v0 = vecs[:, -1]

    v0 = numpy.abs(v0)

    count = 0
    idxs = []
    values = []
    for i in range(dim):
        if abs(v0[i]) > EPSILON:
            #write(i)
            idxs.append(i)
            values.append(v0[i])
    print "nnz:", len(idxs)


    return

    propagate = SumOp(model.n, model.xops)

    syndrome = numpy.zeros(dim)
    for op in model.zops:
        syndrome += op.phases

    idxs = [i for i in range(dim) if v0[i]>EPSILON]
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
        for xop in model.xops:
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
        for xop in model.xops:
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




