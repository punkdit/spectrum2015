#!/usr/bin/env python

import sys

import numpy
import scipy.sparse.linalg as la


EPSILON = 1e-8
scalar = numpy.float64

def write(s):
    sys.stdout.write(str(s)+' ')
    sys.stdout.flush()


class Op(la.LinearOperator):

    verbose = False
    def __init__(self, n):
        self.n = n
        la.LinearOperator.__init__(self, 
            (n, n),
            dtype=scalar, 
            matvec=self.matvec)

    def __len__(self):
        return self.n

    def __mul__(self, other):
        return MulOp(self, other)

    def __add__(self, other):
        return SumOp(self.n, [self, other])

    def __sub__(self, other):
        return SumOp(self.n, [self, -1.0*other])

    def __rmul__(self, alpha):
        return RMulOp(self, alpha)

    def matvec(self, v, w=None, verbose=False):
        if self.verbose:
            write('.')
        assert len(v) == self.n
        if w is None:
            return v
        w += v
        return w

    def todense(self):
        n = self.n
        A = numpy.zeros((n, n))
        v = numpy.zeros(n)
        for i in range(n):
            v[i] = 1
            A[:, i] = self.matvec(v)
            v[i] = 0
        return A


class IdOp(Op):
    pass

class XOp(Op):
    def __init__(self, n, idxs, alpha=1.0):
        " X operator to specified qubits "

        if type(idxs) in (int, long):
            idxs = [idxs]
        for idx in idxs:
            assert 0<=idx<n
        assert len(set(idxs))==len(idxs), idxs
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
        self.idxs = list(idxs)
        self.alpha = alpha
        Op.__init__(self, 2**n)

    def __str__(self):
        return "XOp(%d, %s, %s)"%(self.n, self.idxs, self.alpha)

    def matvec(self, v, w=None, verbose=False):
        if self.verbose:
            write('.')
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
        assert len(set(idxs))==len(idxs), idxs
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
        self.idxs = list(idxs)
        self.alpha = alpha
        Op.__init__(self, 2**n)

    def __str__(self):
        return "ZOp(%d, %s, %s)"%(self.n, self.idxs, self.alpha)

    def matvec(self, v, w=None, verbose=False):
        if self.verbose:
            write('.')
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
        self.n = n
        self.count = 0
        Op.__init__(self, self.n)

    def matvec(self, v, w=None, verbose=True):
        #print "SumOp.matvec"
        if self.verbose:
            write('.')
        assert len(v)==self.n
        if w is None:
            w = numpy.zeros(len(v))
        for op in self.ops:
            op.matvec(v, w)
        self.count += 1
        return w


class RMulOp(Op):
    def __init__(self, op, alpha=1.0):
        self.op = op
        self.alpha = alpha
        Op.__init__(self, op.n)

    def matvec(self, v, w=None, verbose=True):
        if self.verbose:
            write('.')
        assert len(v)==self.n
        if w is None:
            w = numpy.zeros(len(v))
        w += self.alpha * self.op.matvec(v, verbose=verbose)
        return w


class MulOp(Op):
    def __init__(self, a, b):
        # first do b then a !!
        self.ops = a, b
        assert a.n == b.n, (a.n, b.n)
        Op.__init__(self, a.n)

    def matvec(self, v, w=None, verbose=True):
        if self.verbose:
            write('.')
        assert len(v)==self.n
        if w is None:
            w = numpy.zeros(len(v))
        a, b = self.ops
        v = b.matvec(v, verbose=verbose)
        v = a.matvec(v, verbose=verbose)
        w += v
        return w


def test():

    if argv.threads:
        os.environ['OMP_NUM_THREADS'] = str(argv.threads)

    n = argv.get("n", 1024)

    k = argv.get("k", 2)

    A = IdOp(n)

    A.verbose = True
    which = argv.get("which", 'LA')
    vals, vecs = la.eigsh(A, k=k, which=which, maxiter=None) #, tol=1e-8)

    print 
    print vals
    

from argv import Argv 
argv = Argv()


if __name__ == "__main__":

    test()




