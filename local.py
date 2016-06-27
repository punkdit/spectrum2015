#!/usr/bin/env python

import sys 
from math import sqrt
from random import random, choice
from operator import mul, add 
#import cPickle as pickle
import pickle

import numpy

from qupy.codes.smap import SMap

from qupy.abstract import AbstractQu
from qupy.dense import Qu, Gate, genidx
from qupy.sparse import SparseQu, EPSILON
from qupy.ldpc.solve import parse, dot2, shortstr, solve, span, enum2


def cross(itemss):
    if len(itemss)==0:
        yield ()
    else:
        for head in itemss[0]:
            for tail in cross(itemss[1:]):
                yield (head,)+tail

class FuzzyDict(object):
    def __init__(self, epsilon=EPSILON):
        self._items = []
        self.epsilon = epsilon

    def __getitem__(self, x):
        epsilon = self.epsilon
        for key, value in self._items:
            if abs(key - x) < epsilon:
                return value
        raise KeyError

    def get(self, x, default=None):
        epsilon = self.epsilon
        for key, value in self._items:
            if abs(key - x) < epsilon:
                return value
        return default

    def __setitem__(self, x, y):
        epsilon = self.epsilon
        items = self._items
        idx = 0
        while idx < len(items):
            key, value = items[idx]
            if abs(key - x) < epsilon:
                items[idx] = (x, y)
                return
            idx += 1
        items.append((x, y))

    def __str__(self):
        return "{%s}"%(', '.join("%s: %s"%item for item in self._items))



class IdentityQu(AbstractQu):

    # XX cache these objects ...
    # def __new__..

    def __init__(self, shape, valence):
        d = shape[0]
        if shape != (d, d):
            raise ValueError, "square only"
        AbstractQu.__init__(self, shape, valence)

    def __eq__(self, other):
        if isinstance(other, IdentityQu):
            return self.space==other.space
        elif self.space==other.space:
            d = self.shape[0]
            for i in range(d):
                for j in range(d):
                    r = other[i, j]
                    if (r==1) != (i==j): # EPSILON? symbolic?
                        return False
        return False

    def __or__(self, other):
        self.match_valence(other)
        return other
    __ror__ = __or__

    def clone(self):
        return self

    def tosparse(self):
        box = SparseQu(self.shape, self.valence)
        d = self.shape[0]
        assert self.shape == (d, d)
        for i in range(d):
            box[i, i] = 1.
        return box

    def tofp(self):
        d = self.shape[0]
        assert self.shape == (d, d)
        return Gate.identity(d)

    def eigs(self):
        basis = Qu.basis(self.shape, self.valence)
        items = [(v, 1.) for v in basis]
        return items


class LazyTensorQu(AbstractQu):

    def __init__(self, ops):
        shape = reduce(add, [op.shape for op in ops])
        valence = ''.join(op.valence for op in ops)
        AbstractQu.__init__(self, shape, valence)
        self.ops = ops

    def __str__(self):
        return "LazyTensorQu(%s)"%(', '.join(op.name or str(op) for op in self.ops))

    def __mul__(A, B):
        return LazyTensorQu(A.ops + B.ops)

    def __rmul__(A, r):
        op = A.clone()
        op[0] *= r
        return op

    def __neg__(A):
        return (-1.)*A

    def __div__(A, r):
        assert float(r)==r, "scalar only"
        return (1./r)*A

    def _match(A, B):
        if [op.space for op in A.ops] != [op.space for op in B.ops]:
            raise ValueError

    def __or__(A, B):
        if isinstance(B, SparseQu):
            return A.__or_sparse(B)
        if isinstance(B, LazyTensorQu):
            A._match(B)
            ops = [op1|op2 for (op1, op2) in zip(A.ops, B.ops)]
            return LazyTensorQu(ops)
        else:
            raise ValueError, type(B)

    def __or_sparse(A, B): 
        d = B.shape[0]
        v = SparseQu(B.shape, B.valence)
        enum_A_items = list(enumerate(A.items()))
        mulovers = list(genidx((d,)*len(enum_A_items)))
        data = v.data
        for basis, phase in B.items():
            for mulover in mulovers:
                basis1 = list(basis)
                for i, (idx, op) in enum_A_items:
                    basis1[idx] = mulover[i]
                r = phase
                for i, (idx, op) in enum_A_items:
                    r *= op.v[mulover[i], basis[idx]]
                basis1 = tuple(basis1)
                r0 = data.get(basis1, None)
                if r0 is None:
                    if r==0:
                        pass
                    elif type(r)!=float or abs(r) > EPSILON:
                        data[basis1] = r 
                else:
                    r += r0
                    if r==0:
                        del data[basis1]
                    elif type(r)!=float or abs(r) > EPSILON:
                        data[basis1] = r 
                    else:
                        del data[basis1]
        return v

#    def __add__(A, B):
#        A._match(B) # should only need spaces to match
#        idxs = []
#        for idx in range(len(A.ops)):
#            aop, bop = A.ops[idx], B.ops[idx]
#            if not isinstance(aop, IdentityQu) or \
#                not isinstance(bop, IdentityQu):
#                # XXX TOO HARD :-( XXX

    def __eq__(A, B):
        if A.space != B.space:
            return False
        if A.ops == B.ops:
            # sufficient but not necessary
            return True
        else:
            assert 0, "not implemented; use .todense()?"

    def clone(self):
        return LazyTensorQu([op.clone() for op in self.ops])

    def insert(self, idx, op):
        ops = [op.clone() for op in self.ops]
        ops.insert(idx, op)
        return LazyTensorQu(ops)

    @classmethod
    def promote(self, n, d, op):
        if isinstance(op, LazyTensorQu):
            return op
        if isinstance(op, dict):
            ops = [op.get(i, IdentityQu((d, d), 'ud')) for i in range(n)]
            op = LazyTensorQu(ops)
            return op
        assert 0, type(op)

    def todense(self):
        return reduce(mul, self.ops)

    def tofp(self):
        ops = [op.tofp() for op in self.ops]
        return LazyTensorQu(ops)

    def eigs(self):
        es = [op.eigs() for op in self.ops]
        items = []
        for item in cross(es):
            val = 1.
            vec = None
            for _val, _vec in item:
                val *= _val
                vec = _vec if vec is None else vec*_vec
            items.append((val, vec))
        return items

    # 
    # ---------- this stuff makes us look like a dict -------
    #

    def keys(self):
        keys = [i for i, op in enumerate(self.ops)
            if not isinstance(op, IdentityQu)]
        return keys

    def values(self):
        keys = self.keys()
        return [self.ops[key] for key in keys]

    def items(self):
        return zip(self.keys(), self.values())

    def get(self, i, default=None):
        op = self.ops[i]
        if isinstance(op, IdentityQu):
            return default
        return op

    def __getitem__(self, i):
        op = self.ops[i]
        if isinstance(op, IdentityQu):
            raise KeyError
        return op

    def __setitem__(self, i, box):
        assert self.ops[i].shape == box.shape
        self.ops[i] = box

    # 
    # ---------- end of stuff that makes us look like a dict -------
    #


class LocalQu(AbstractQu):
    """ Sum of LazyTensorQu's 
    """

    def __init__(self, d, n, ops=[]):
        self.d = d
        self.n = n # == rank/2
        shape = (d, d) * n
        valence = 'ud'*n
        AbstractQu.__init__(self, shape, valence)

        self.ops = []
        for op in ops:
            op = LazyTensorQu.promote(self.n, self.d, op)
            self.ops.append(op)

    def clone(self):
        box = LocalQu(self.d, self.n)
        for op in self.ops:
            op = op.clone()
            box.ops.append(op)
        return box

#    def compact(self):
#        items = []
#        for opmap in self.ops:
#            keys = opmap.keys()
#            keys.sort()
#            item = tuple((key, opmap[key]) for key in keys)
#            items.append(item)
#        uniq = list(set(items))
#        print (len(items), len(uniq))
#        #for item in uniq:
#        #    print items.count(item),
#        #print
#        assert 0, "not implemented"

    some_descs = [
        (Gate.I, 'I'),
        (Gate.X, 'X'),
        (Gate.Z, 'Z'),
        (Gate.Y, 'Y'),
        (Qu((2, 2), 'ud', [[1, 0], [0, 0]]), '|0><0|'),
        (Qu((2, 2), 'ud', [[0, 0], [0, 1]]), '|1><1|'),
        (Qu((2, 2), 'ud', [[0, 1], [0, 0]]), '|0><1|'),
        (Qu((2, 2), 'ud', [[0, 0], [1, 0]]), '|1><0|'),
    ]

    # XX this should be a Qu method XX
    @classmethod
    def get_desc(cls, B):
        for A, desc in cls.some_descs:
            if B.space==A.space and A.is_close(B):
                return desc
        return str(B)

    def __str__(self):
        descs = []
        for op in self.ops:
            desc = []
            for key, A in op.items():
                desc.append("%s:%s" % (key, LocalQu.get_desc(A)))
            descs.append(','.join(desc))
        return "LocalQu(%d, %d, %s)"%(
            self.d, self.n, ' + '.join(descs))

    def add(self, *items):
        d, n = self.d, self.n
        ops = [None] * self.n
        idx = 0
        while idx+1 < len(items):
            i = items[idx]
            box = items[idx+1]
            assert box.shape == (d, d)
            assert box.valence == 'ud'
            assert ops[i] is None
            ops[i] = box.clone()
            idx += 2
        for i, op in enumerate(ops):
            if op is None:
                op = IdentityQu((d, d), 'ud')
                ops[i] = op
        op = LazyTensorQu(ops)
        self.ops.append(op)

    def __add__(self, other):
        assert self.space == other.space
        return LocalQu(self.d, self.n, self.ops + other.ops)

    def __mul__(self, r):
        self = self.clone()
        for op in self.ops:
            if not op:
                continue
            key = op.keys()[0]
            #print op[key]
            op[key] = r*op[key]
            #print op[key]
            #print
        return self

    __rmul__ = __mul__

    def __neg__(self):
        return -1. * self

    def __sub__(self, other):
        return self + (-1. * other)

    def __div__(self, r):
        return (1./r) * self

    def __or__(self, v):
        if (v.shape == (self.d,)*self.n and
            v.valence == 'u'*self.n):
            assert isinstance(v, SparseQu)
            u = SparseQu(v.shape, v.valence)
            for op in self.ops:
                u0 = op | v
                u += u0
        else:
            if not isinstance(v, LocalQu):
                raise ValueError, type(v)
            A, B = self, v
            ops = [op1|op2 for op1 in A.ops for op2 in B.ops]
            u = LocalQu(self.d, self.n, ops)
        return u

    def __eq__(A, B):
        if A.space != B.space:
            return False
        if A.ops == B.ops:
            return True
        # not necessarily False... :-(
        #return False
        raise ValueError

    def __ne__(A, B):
        return not A.__eq__(B)

    def tosparse(self):
        A = SparseQu(self.shape, self.valence)
        n, d = self.n, self.d
        I = SparseQu.identity(d)
        for op in self.ops:
            A0 = 1.
            for i in range(n):
                A1 = op.get(i, I)
                A1 = A1.tosparse()
                A0 = A0*A1 # tensor
            A += A0
        return A

    def todense(self):
        A = self.tosparse()
        A = A.todense()
        return A

    def tofp(self):
        ops = [op.tofp() for op in self.ops]
        return LocalQu(self.d, self.n, ops)



def build_double(G):
    n = len(G)
    Lps, Lns, Tps, Tns = {}, {}, {}, {}
    zero = lambda n : Qu((n, n), 'ud')
    for g in G:
        Lp, Ln, Tp, Tn = zero(n), zero(n), zero(n), zero(n)
        for h in G:
            Lp[G.coord(g*h, h)] = 1
            Ln[G.coord(h*~g, h)] = 1
        Tp[G.coord(g, g)] = 1
        Tn[G.coord(~g, ~g)] = 1
        Lps[g] = Lp
        Lns[g] = Ln
        Tps[g] = Tp
        Tns[g] = Tn
    return Lps, Lns, Tps, Tns


class Lattice(object):
    def __init__(self, G, l):

        double = build_double(G)
        self.Lps, self.Lns, self.Tps, self.Tns = double
        self.G = G
        self.d = d = len(self.G)
        self.l = l
        self.n = n = 2*self.l**2

        self.ijs = ijs = [(i, j) for i in range(l) for j in range(l)]
        self.keys = [(i, j, k) for i,j in ijs for k in (0, 1)]

        G = self.G
        As = {}
        Bs = {}
        all_ops = []
        for i, j in ijs:
            if i==l-1 and j==l-1:
                continue
            for g in G:
                A = self.make_star(i, j, g)
                As[i, j, g] = A
            A = reduce(add, [As[i, j, g] for g in G])
            A = (1./len(G))*A
            all_ops.append(A)
            B = self.make_plaq(i, j)
            Bs[i, j] = B
            all_ops.append(B)
        self.As = As
        self.Bs = Bs
        self.all_ops = all_ops

    def make_star(self, i, j, g):
        A = self.mkop({
            (i, j, 0) : self.Lns[g], 
            (i, j, 1) : self.Lns[g], 
            (i-1, j, 0) : self.Lps[g], 
            (i, j-1, 1) : self.Lps[g],
        })
        return A

    def make_plaq(self, i, j, g=None):
        G = self.G
        if g is None:
            g = G.id
        Tps, Tns = self.Tps, self.Tns
        #if (i+j)%2==1:
        #    Tps, Tns = Tns, Tps
        op = LocalQu(self.d, self.n)
        for h1 in G:
         for h2 in G:
          for h3 in G:
            h4 = ~(h3*h2*h1)
            assert h4*h3*h2*h1 == G.id
            op0 = self.mkop({
                (i, j, 0) : Tns[h1], 
                (i, j, 1) : Tps[h2], 
                (i, j+1, 0) : Tps[h3], 
                (i+1, j, 1) : Tns[h4],
            })
            op = op + op0
        #print "make_plaq", (i, j), len(op.ops),
        #print op.compact() # useless..
        if len(self.G)==2 and self.l<3:
            op = op.tosparse() # much much faster than a LocalQu..
            #print op.nnz # << d**n
        return op

    def get_index(self, i, j, k):
        idx = self.keys.index((i%self.l, j%self.l, k))
        return idx

    def mkop(self, coords):
        if type(coords) is dict:
            coords = coords.items()
        op = {}
        for key, value in coords:
            i, j, k = key
            idx = self.keys.index((i%self.l, j%self.l, k))
            value = value.todense()
            # or tosparse ?
            op[idx] = value
        op = LocalQu(self.d, self.n, [op])
        return op


class Heisenberg(object):
    def __init__(self, l, **kw):
        self.d = d = 2
        self.n = n = l**2
        self.keys = keys = [(i, j) for i in range(l) for j in range(l)]
        coords = {}
        for i, j in keys:
            coords[i, j] = keys.index((i, j))
            coords[i+1, j] = keys.index(((i+1)%l, j))
            coords[i, j+1] = keys.index((i, (j+1)%l))
            coords[i+1, j+1] = keys.index(((i+1)%l, (j+1)%l))

        X, Y, Z = Gate.X, Gate.Y, Gate.Z
        op = LocalQu(d, n)
        jx, jy, jz = 1., 1., 1.
        h = -0.1
        for i, j in keys:
            for lop in (jx*X, jy*Y, jz*Z):
                op.add(coords[i, j], lop, coords[i+1, j], lop)
                op.add(coords[i, j], lop, coords[i, j+1], lop)
            op.add(coords[i, j], h*Z)

        self.H = op # negative H

        
class Heisenberg3(object):
    def __init__(self, l, **kw):
        self.d = d = 2
        keys = [
            (i, j, k) 
            for i in range(l)
            for j in range(l)
            for k in range(l)]
        keys2 = ijs = [(i, j) for i in range(l) for j in range(l)]

        # create handy lookup table
        coords = {}
        for i, j, k in keys:
            coords[i, j, k] = keys.index((i, j, k%l))
            coords[i+1, j, k] = keys.index(((i+1)%l, j, k%l))
            coords[i, j+1, k] = keys.index((i, (j+1)%l, k%l))
            coords[i+1, j+1, k] = keys.index(((i+1)%l, (j+1)%l, k%l))
            coords[i, j, k+1] = keys.index((i, j, (k+1)%l))
            coords[i+1, j, k+1] = keys.index(((i+1)%l, j, (k+1)%l))
            coords[i, j+1, k+1] = keys.index((i, (j+1)%l, (k+1)%l))
            coords[i+1, j+1, k+1] = keys.index(((i+1)%l, (j+1)%l, (k+1)%l))

        keys3 = list(keys)
        for i, j in ijs:
            keys.append((i, j))
        keys = keys3+keys2
        self.n = n = len(keys)

        for i, j in ijs:
            coords[i, j] = keys.index((i, j))
            coords[i+1, j] = keys.index(((i+1)%l, j))
            coords[i, j+1] = keys.index((i, (j+1)%l))
            coords[i+1, j+1] = keys.index(((i+1)%l, (j+1)%l))

        X, Y, Z = Gate.X, Gate.Y, Gate.Z
        op = LocalQu(d, n)
        jx, jy, jz = 1., 1., 1. # note: actually j* = sqrt(J)
        h = -0.1
        for i, j, k in keys3:
            for lop in (jx*X, jy*Y, jz*Z):
                op.add(coords[i, j, k], lop, coords[i+1, j, k], lop)
                op.add(coords[i, j, k], lop, coords[i, j+1, k], lop)
                op.add(coords[i, j, k], lop, coords[i, j, k+1], lop)
            op.add(coords[i, j, k], h*Z)

        A = -0.1 # A<<J
        for i, j in ijs:
            op.add(
                coords[i+1, j], Z, coords[i, j], Y,
                coords[i, j+1], Z, coords[i+1, j+1], Y) # W_p
            op.add(
                coords[i+1, j], A*Z, coords[i, j], Y,
                coords[i, j+1], Z, coords[i+1, j+1], Y,
                coords[i, j, 0], X) # H_int

        self.H = op # negative H


class RandomBond(object):
    def __init__(self, l, **kw):
        self.d = d = 2
        self.n = n = l**2
        self.keys = keys = [(i, j) for i in range(l) for j in range(l)]
        coords = {}
        for i, j in keys:
            coords[i, j] = keys.index((i, j))
            coords[i+1, j] = keys.index(((i+1)%l, j))
            coords[i, j+1] = keys.index((i, (j+1)%l))
            coords[i+1, j+1] = keys.index(((i+1)%l, (j+1)%l))

        X, Y, Z = Gate.X, Gate.Y, Gate.Z
        op = LocalQu(d, n)
        jx, jy, jz = 1., 1., 1.
        h = -0.1
        for i, j in keys:
            lop = choice((jx*X, jy*Y, jz*Z))
            op.add(coords[i, j], lop, coords[i+1, j], lop)
            lop = choice((jx*X, jy*Y, jz*Z))
            op.add(coords[i, j], lop, coords[i, j+1], lop)
            #op.add(coords[i, j], h*Z)

        self.H = op # negative H

        
class Ising(object):
    def __init__(self, l, **kw):
        self.d = d = 2
        self.n = n = l**2
        self.keys = keys = [(i, j) for i in range(l) for j in range(l)]
        coords = {}
        for i, j in keys:
            coords[i, j] = keys.index((i, j))
            coords[i+1, j] = keys.index(((i+1)%l, j))
            coords[i, j+1] = keys.index((i, (j+1)%l))
            coords[i+1, j+1] = keys.index(((i+1)%l, (j+1)%l))

        X, Y, Z = Gate.X, Gate.Y, Gate.Z
        op = LocalQu(d, n)
        jx, jy, jz = 1., 1., 1.
        h = 1.0
        for i, j in keys:
            lop = Z
            op.add(coords[i, j], lop, coords[i+1, j], lop)
            op.add(coords[i, j], lop, coords[i, j+1], lop)
            op.add(coords[i, j], h*X)

        self.H = op # negative H

        
class Toric(object):
    def __init__(self, l, Jz=1., Jx=1., **kw):
        self.d = d = 2
        self.n = n = 2*l**2
        self.keys = keys = [(i, j, k)
            for i in range(l) for j in range(l) for k in range(2)]
        coords = {}
        for i, j, k in keys:
            coords[i, j, k] = keys.index((i, j, k))
            coords[i+1, j, k] = keys.index(((i+1)%l, j, k))
            coords[i, j+1, k] = keys.index((i, (j+1)%l, k))
            coords[i+1, j+1, k] = keys.index(((i+1)%l, (j+1)%l, k))
            coords[i-1, j, k] = keys.index(((i-1)%l, j, k))
            coords[i, j-1, k] = keys.index((i, (j-1)%l, k))

        X, Y, Z = Gate.X, Gate.Y, Gate.Z
        op = LocalQu(d, n)
        for i in range(l):
          for j in range(l):
            op.add(
                coords[i, j, 0], Jz*Z,
                coords[i, j, 1], Z,
                coords[i, j+1, 0], Z,
                coords[i+1, j, 1], Z,
            )
            op.add(
                coords[i, j, 0], Jx*X,
                coords[i, j, 1], X,
                coords[i-1, j, 0], X,
                coords[i, j-1, 1], X,
            )

        self.H = op # negative H

        
class Compass(object):
    "2d compass model"
    def __init__(self, l, Jz=1., Jx=1., periodic=True, **kw):
        self.d = d = 2
        self.n = n = l**2
        self.keys = keys = [(i, j) for i in range(l) for j in range(l)]
        coords = {}
        for i, j in keys:
            coords[i, j] = keys.index((i, j))
            coords[i+1, j] = keys.index(((i+1)%l, j))
            coords[i-1, j] = keys.index(((i-1)%l, j))
            coords[i, j-1] = keys.index((i, (j-1)%l))
            coords[i, j+1] = keys.index((i, (j+1)%l))
            coords[i+1, j+1] = keys.index(((i+1)%l, (j+1)%l))

        X, Y, Z = Gate.X, Gate.Y, Gate.Z
        op = LocalQu(d, n)
        for i, j in keys:
            if i==l-1 and not periodic:
                continue
            if j==l-1 and not periodic:
                continue
            #print coords[i, j], "--", coords[i+1, j]
            op.add(coords[i, j], Jz*Z, coords[i+1, j], Z)
            op.add(coords[i, j], Jx*X, coords[i, j+1], X)
        self.H = op # negative H


class Gapped(object):
    def __init__(self, n):
        self.n = n
        self.d = 2
        X, Y, Z = Gate.X, Gate.Y, Gate.Z
        op = LocalQu(self.d, n)

        for i in range(n):
            op.add(i, X)
            op.add(i, Z)
        self.H = op # negative H


class CompassBlock(object):
    def __init__(self, l, **kw):
        self.d = d = 2
        self.n = n = l**2
        self.keys = keys = [(i, j) for i in range(l) for j in range(l)]
        coords = {}
        for i, j in keys:
            coords[i, j] = keys.index((i, j))
            coords[i+1, j] = keys.index(((i+1)%l, j))
            coords[i-1, j] = keys.index(((i-1)%l, j))
            coords[i, j-1] = keys.index((i, (j-1)%l))
            coords[i, j+1] = keys.index((i, (j+1)%l))
            coords[i+1, j+1] = keys.index(((i+1)%l, (j+1)%l))

        X, Y, Z = Gate.X, Gate.Y, Gate.Z
        op = LocalQu(d, n)
        for i, j in keys:
            op.add(coords[i, j], X, coords[i, j+1], X)
            op.add(coords[i, j], Z, coords[i+1, j], Z)

        for i in range(l):
            op.add(coords[i, 0], X)
            op.add(coords[0, i], Z)

        for i in range(l-1):
            ops = []
            for j in range(l):
                ops.extend([coords[j, i], X, coords[j, i+1], X])
            op.add(*ops)
            ops = []
            for j in range(l):
                ops.extend([coords[i, j], Z, coords[i+1, j], Z])
            op.add(*ops)

        ops = []
        for i in range(l):
            ops.extend([coords[i, 0], X])
        op.add(*ops)

        ops = []
        for i in range(l):
            ops.extend([coords[0, i], Z])
        op.add(*ops)

        self.H = op # negative H


class Model(object):
    def __init__(self, d, n):
        self.d = d
        self.n = n
        I = Gate.I
        In2 = LocalQu(d, n)
        In2.add(0, 0.5 * I)
        self.In2 = In2

    def build_op(self, v, tp):
        op = LocalQu(self.d, self.n)
        items = []
        for j in range(self.n):
            if v[j]:
                items.append(j)
                items.append(tp)
        op.add(*items)
        return op

    def measure(self, g, v):
        g = 0.5 * g
        pos = self.In2 + g
        neg = self.In2 - g

        pv = pos|v
        nv = neg|v
        pr = pv.dot(pv)
        nr = nv.dot(nv)
        assert abs((pr+nr)-1.)<EPSILON, abs((pr+nr)-1.)
        #print pr, nr, v.nnz
        if random() <= pr:
            eigv = 1
            v = pv
        else:
            eigv = -1
            v = nv
        v = v.normalized()
        return eigv, v


class Gauge(Model):
    def __init__(self, code, Jz=1., Jx=1., **kw):
        Model.__init__(self, 2, code.n)

        X, Y, Z = Gate.X.clone(), Gate.Y, Gate.Z.clone()
        #Z.v = Z.v.real
        #X.v = X.v.real
        H = LocalQu(self.d, self.n)
        Hx = LocalQu(self.d, self.n)
        Hz = LocalQu(self.d, self.n)
        for i in range(code.gz):
            g = code.Gz[i]
            items = []
            for j in range(self.n):
                if g[j]:
                    items.append(j)
                    items.append(Z)
            items[1] = Jz*Z
            H.add(*items)
            Hz.add(*items)
        for i in range(code.gx):
            g = code.Gx[i]
            items = []
            for j in range(self.n):
                if g[j]:
                    items.append(j)
                    items.append(X)
            items[1] = Jx*X
            H.add(*items)
            Hx.add(*items)
        self.H = H # negative H
        self.Hxs = Hx
        self.Hzs = Hz

        if code.Lx is not None:
            self.Lx = self.build_op(code.Lx[0], X)
        if code.Lz is not None:
            self.Lz = self.build_op(code.Lz[0], Z)

        self.Hx = [self.build_op(code.Hx[i], X) for i in range(code.mx)]
        self.Hz = [self.build_op(code.Hz[i], Z) for i in range(code.mz)]
        self.Gx = [self.build_op(code.Gx[i], X) for i in range(code.gx)]
        self.Gz = [self.build_op(code.Gz[i], Z) for i in range(code.gz)]

        op = LocalQu(self.d, self.n)
        for h in self.Hx + self.Hz:
            op += h
        self.H_stab = op

        self.mx = code.mx
        self.mz = code.mz
        self.gx = code.gx
        self.gz = code.gz
        self.k = 1


def get_model():

    l = argv.get('l', 1)
    n = argv.get('n', 16)

    Jz = argv.get('Jz', 1./2**0.5)
    Jx = argv.get('Jx', 1./2**0.5)
    print "Jz:", Jz, "Jx:", Jx

    if argv.heisenberg:
        model = Heisenberg(l)
    elif argv.heisenberg3:
        assert l%2, "l needs to be odd i think"
        model = Heisenberg3(l)
    elif argv.random:
        model = RandomBond(l)
    elif argv.ising:
        model = Ising(l)
    elif argv.compass:
        model = Compass(l, Jz, Jx)
    elif argv.compassblock:
        model = CompassBlock(l)
    elif argv.gapped:
        model = Gapped(n)
    elif argv.toric:
        model = Toric(l, Jz, Jx)
    elif argv.gcolor:
        from qupy.ldpc.gcolor import Lattice
        lattice = Lattice(l)
        code = lattice.build_code()
        model = Gauge(code, Jz, Jx)
    else:
        model = None

    return model


def test_model():
    """
        This is the power method. We look for the largest
        eigenvalue(s) of -H.
    """

    l = argv.get('l', 1)
    model = get_model()

    n, d = model.n, model.d
    H = model.H

    print "n =", n, "2**n =", 2**n

    K = argv.get('K', 1)

    sigma = argv.get('sigma', 10.)
    max_count = argv.get('count', None)
    threshold = argv.get('threshold')

    normalize = argv.get('normalize', True)
    orthogonal = argv.get('orthogonal', False)

    if argv.exact:
        H = H.todense()
        print "shape:", H.shape, H.valence
        H = H.get_flatop().do(H)
        #eigs = H.eigs()
        import numpy
        #eigvals, v = numpy.linalg.eig(H.v)
        eigvals = numpy.linalg.eigvals(H.v)
        eigvals = [float(value.real) for value in eigvals]
        eigvals.sort(reverse=True)
        print eigvals#[:50]
        #for val, vec in eigs[:4]:
        #    print val
        return

    vs = [None]*K

    stem = argv.get('load')
    if stem:
        for j in range(K):
            name = '%s-%d.txt'%(stem, j)
            try:
                v = SparseQu.load(name)
                vs[j] = v
                if normalize:
                    v /= v.norm()
            except IOError, e:
                print e
                continue

    if argv.killstab:
      for i in range(K):
        v = vs[i]
        keys = list(v.keys())
        for key in keys:
            r = v[key]
            idx = ''.join(str(x) for x in key)
            assert abs(r.imag)<EPSILON
            u = parse(idx)
            u = u.transpose()
            syndrome = dot2(code.Hz, u)
            if syndrome.sum():
                del v[key]

    if argv.stabilize:
        v = vs[-1]
        for _ in range(model.mx):
          for h in model.Hx:
            v = v+(h|v)
        #v /= v.norm()
          print "norm:", v.norm()
        vs[-1] = v

    if argv.show is not None:

        idx = argv.show
        assert idx is not True, "please set show=idx"
        v = vs[idx]
        keys = list(v.keys())
        keys.sort(key = lambda idx : -abs(v[idx]))
        for key in keys:
            r = v[key]
            if abs(r) < threshold:
                break
            idx = ''.join(str(x) for x in key)
            assert abs(r.imag)<EPSILON
            u = parse(idx)
            u = u.transpose()
            #w = solve(code.Gx.transpose(), u)
            syndrome = dot2(code.Gz, u)
            #syndrome = dot2(code.Hz, u)
            info = shortstr(syndrome.transpose())
            #info = shortstr(w.transpose()) if w is not None else '-'
            print idx, info, sum(key), r.real

            #u = SparseQu((d,)*n, 'u'*n)
            #u[key] = 1.
            #w = H|u
            #print w.shortstr()
            #print u.dot(w)
            #print
        return

    if argv.measure is not None:
        idx = argv.measure
        assert idx is not True, "please set measure=idx"
        v = vs[idx]
        v /= v.norm()

        syndrome = []
        assert model.gx == model.gz
        for _ in range(0):
         for i in range(model.gx):
          for G in [model.Gx, model.Gz]:

            eigv, v = model.measure(G[i], v)
            #print v.nnz

            u = model.H | v
            E = v.dot(u)
            print "<E>", E.real

        #print syndrome
        #print v.nnz
        #print v.str()

        for S in [model.Hx, model.Hz]:
          for i in range(model.mx):
            #eigv, v = model.measure(S[i], v)
            r = v - (S[i] | v)
            print r.norm()

        vs[idx] = v
        return

    HH = None
    if argv.simultaneous:
        HH = -model.H_stab

    if argv.H_stab:
        H = -model.H_stab

    if argv.expect:
      for idx in range(K):
        v = vs[idx]
        if v is None:
            continue

        v /= v.norm()

        u = H|v
        E = v.dot(u)
        print "<E> =", E.real
      return

    if argv.stabdist:

        v = SparseQu((d,)*n, 'u'*n)
        #S = model.Hx

        #k = 2**(code.mx-1)
        #values = [-1] * k + [1] * k
        #for u in span(code.Hx):
        #    key = tuple(u)
        #    print key
        #    v[key] = values.pop(choice(range(len(values))))
        #assert not values

        values = [choice([-1, +1]) for i in range(code.mx)]
        for u in enum2(code.mx):
            r = 1
            for i, ii in enumerate(u):
                if ii:
                    r *= values[i]
            S = dot2(u, code.Hx)
            S = tuple(S)
            print S, r
            v[S] = 1.*r

        #return

        v /= v.norm()

        vs[K-1] = v

        print v.str()

    if argv.compass and 1:
        for i in range(min(2, K)):
            v = SparseQu((d,)*n, 'u'*n)
            for key in genidx((2,)*l):
                if sum(key)%2==i:
                    #print key*l
                    v[key*l] = 1. # set each stabilizer*logop to 1.
            v /= v.norm()
            vs[i] = v
            print "[%d] lnnz: %d"%(i, v.nnz)

    if argv.logop:
        v = vs[0]
        v = model.Lx | v
        v /= v.norm()
        #v = model.Lx | v
        #v /= v.norm()
        vs[1] = v
        print "logop:", v.norm(), v.nnz

    for i in range(K):
        if vs[i] is not None:
            continue
        v = SparseQu((d,)*n, 'u'*n)
        nnz = argv.get('nnz', 10)
        #v = SparseQu.random((d,)*n, 'u'*n, nnz)
        v = SparseQu((d,)*n, 'u'*n)
        for j in range(nnz):
            idx = ''.join([choice('01') for _ in range(n)])
            v[idx] = choice([-1., 1.])
        if nnz==0:
            v['0'*n] = 1
#            v['1'*n] = 1
        v = v / v.norm()
        if argv.todense:
            v = v.todense()
        print "setting vs[%d] nnz=%d"%(i, v.nnz)
        vs[i] = v

    dot = lambda u, v : u.dot(v)
    norm2 = lambda r : r.real**2 + r.imag**2

    if argv.todense:
        H = H.todense()
        #dot = lambda u, v : (~u|v)

    i = 0
    while max_count is None or i<max_count:
    #for i in range(max_count):

      i += 1

      for j in range(K):

        if argv.lastonly and j<K-int(argv.lastonly):
            continue

        v = vs[j]

        u = (H | v)

        #eigv = u.norm()
        eigv = dot(v, u).real

        u = u + sigma*v

        if normalize:
            r = u.norm()
            u /= r

        if not orthogonal:
          for jj in range(j):
            r = u.dot(vs[jj])
            u = u - r*vs[jj]

        if j and normalize:
            u /= u.norm()

        delta = (v-u).norm()

        print i, "[%d] delta ="%(j,), delta, norm2(dot(u, v)),
        #print "<E> = %.11f, <E>/n = %f" % (eigv, eigv/n),
        print "<E> = %.11f," % (eigv,),

        v = u

        if threshold:
            v = v.threshold(threshold)
            v /= v.norm()

        print "nnz =", v.nnz

        vs[j] = v

        if delta<1e-6 and v.nnz==1:
            break

        stem = argv.get('save')
        if stem:
            name = '%s-%d.txt'%(stem, j)
            try:
                vs[j].save(name)
            except KeyboardInterrupt:
                vs[j].save(name) # Really save!!!
                break

        if argv.species:
            counts = {}
            values = list(v.values())
            values.sort()
            r0 = values[0].real
            counts[r0] = 1
            i = 1
            while i < len(values):
                r1 = values[i].real
                if abs(r1-r0) / (r0+r1) < 1e-4:
                    counts[r0] += 1
                else:
                    r0 = r1
                    counts[r0] = 1
                i += 1
            keys = counts.keys()
            keys.sort()
            print "species:", len(keys)
            for r0 in keys:
                print "%3d"%counts[r0], r0

      if HH is not None:
        H, HH = HH, H

      if K>2:
        print

    print v

    return # <------------ return

    eigs = []
    for j in range(K):
        eig = dot(vs[j], H|vs[j]).real
        eigs.append(eig)
        print "eigv =", eig
        for k in range(j):
            print norm2(dot(vs[j], vs[k]))
        #if v.nnz <= 16:
        #    print v
    for j in range(1, K):
        print "gap", eigs[j] - eigs[j-1]


def test_lazy_tensor():

    X, I = Gate.X, Gate.I

    n = 8
    op = [I] * n
    op[3] = X
    op = LazyTensorQu(op)

    #print op

    assert (op|op)==(op|op)


def test_local():

    X, I = Gate.X, Gate.I
    d, n = 2, 3 # three qubits
    op = LocalQu(d, n)
    op.add(0, X, 2, X) # flip first and last qubits

    assert (2.*op.todense()) == (2.*op).todense()

    v = SparseQu((d,)*n, 'u'*n, {(0, 0, 0):1.})
    w = SparseQu((d,)*n, 'u'*n, {(1, 0, 1):1.})
    assert op|v == w

    v = SparseQu((d,)*n, 'u'*n, {(0, 0, 0):1., (0, 1, 1):1.})
    w = SparseQu((d,)*n, 'u'*n, {(1, 0, 1):1., (1, 1, 0):1.})
    assert op|v == w

    assert (op|op)|v == v

    # --------------------------------------------------------------

    # three qutrits
    d = 3
    n = 3
    op = LocalQu(d, n)

    # flip 1<->2
    X_12 = Qu((d,d), 'ud')
    X_12[0, 0] = 1.
    X_12[1, 2] = 1. # flip
    X_12[2, 1] = 1. # flip

    # put the flip on third site
    op.add(2, X_12)

    assert (2.22*op.todense()) == (2.22*op).todense()

    v = SparseQu((d,)*n, 'u'*n, {(0, 0, 1): 3.0})
    v2 = SparseQu((d,)*n, 'u'*n, {(0, 0, 2): 3.0}) 
    v1 = op | v
    assert v1==v2

    v = SparseQu((d,)*n, 'u'*n,
        {(2, 0, 2): 3.0, (0, 0, 0): 1.0, (1, 0, 1): 2.0})
    v1 = op | v

    v2 = SparseQu((d,)*n, 'u'*n, 
        {(1, 0, 2): (2+0j), (0, 0, 0): (1+0j), (2, 0, 1): (3+0j)})

    assert v1==v2

    op = op.tosparse()
    v1 = op | v
    assert v1==v2

    # --------------------------------------------------------------

    c0 = Compass(3, -1., 0.)
    c1 = Compass(3, 0., -1./1000)
    c2 = Compass(3, -1., -1./1000)

    assert (c0.H.todense() + c1.H.todense()) == c2.H.todense()

    #print c1.H.todense().norm()
    #print c1.H.todense().nnz


def test_toric_ground():

    from qupy.symbolic import groups
    if argv.cyclic:
        d = argv.get('d', 2)
        G = groups.Group.cyclic(d)
    elif argv.sym:
        G = groups.Group.symmetric(3)
    else:
        return


    l = argv.get('l', 2)
    lattice = Lattice(G, l)
    d = lattice.d
    n = lattice.n

    As = lattice.As
    Bs = lattice.Bs

    if argv.get('test'):
        op = reduce(add, [As[0, 0, g] for g in G])
        r = 1./len(G)
        if len(G)==2:
            assert (r*op.todense()) == (r*op).todense()
        op = r*op
        sparse_ops = [op]
        sparse_ops += Bs.values()
    
        # check projection
        for op in sparse_ops:
            v = SparseQu.random((d,)*n, 'u'*n, 64)
            u1 = op|v
            u2 = op|(op|v)
            assert u1==u2
    
        # check commutative
        sn = len(sparse_ops)
        for i in range(sn):
            for j in range(i+1, sn):
                print (i, j)
                v = SparseQu.random((d,)*n, 'u'*n, 16)
                opi = sparse_ops[i]
                opj = sparse_ops[j]
                u1 = opj|(opi|v)
                u2 = opi|(opj|v)
                if u1!=u2:
                    print "FAIL:", (u1-u2).norm()

    all_ops = lattice.all_ops

    v = SparseQu((d,)*n, 'u'*n)
    print "n =", n
    v['0'*n] = 1.
    v /= v.norm()

    if G.n == 2 and l==2:
        sparse_ops = [A.tosparse() for A in all_ops]

        for i in range(len(all_ops)):
            A = all_ops[i]
            w1 = A|v
            w2 = sparse_ops[i]|v
    
            assert w1 == w2

    max_count = argv.get('count', 100)

    name = argv.get('load')
    if name:
        v = SparseQu.load(name)

    print "project:"
    for op in all_ops:
        u = op|v
        u1 = op|u
        assert u==u1
        r = u.norm()
        print r, (u-v).norm(), u.nnz
        if r<1e-6:
            print "nullspace"
            break
        v = u/r

    name = argv.get('save')
    if name:
        v.save(name)

    print "result:"
    for op in all_ops:
        u = op|v
        r = u.norm()
        u /= r
        print abs(u.dot(v))

    #print u
    print "nonzeros:"
    keys = v.keys()
    keys.sort()
    for key in keys:
        print '\t', key


    #I, X, Y, Z = Gate.I, Gate.X, Gate.Y, Gate.Z

    I = Gate.identity(d)
    Z = Gate.Z
    X = lattice.Lps[G[1]]
    assert X != I
    get_index = lattice.get_index

    #print I.eigs()

    v = v.todense()

    dist = FuzzyDict()

    sp = 0.
    op = [I]*n
    #for i in range(l):
        #op[get_index(i, 0, 1)] = X # (0.5, 0.5)
        #op[get_index(0, i, 1)] = Z # (1., 0.)

    op[get_index(0, 0, 0)] = X
    op[get_index(0, 0, 1)] = X
    #op[get_index(1, 0, 0)] = X
    #op[get_index(1, 0, 1)] = X

    op = LazyTensorQu(op)
    for eval, evec in op.eigs():
        r = ~evec|v
        r = complex(r)
        p = (r.real)**2 + (r.imag)**2
        sp += p
        r = sqrt(p)
        dist[eval] = dist.get(eval, 0.) + p
        #print r,
        if r:
            print r, eval
    print sp
    print dist


if __name__ == "__main__":

    from qupy.tool.argv import Argv
    argv = Argv()

    test_lazy_tensor()
    test_local()
    test_toric_ground()
    test_model()




