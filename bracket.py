#!/usr/bin/env python

import sys, os
from fractions import gcd

import numpy

import networkx as nx

from solve import get_reductor, array2, row_reduce, dot2, shortstr, zeros2, shortstrx
from solve import u_inverse
from isomorph import write

import models
from models import genidx

zero = 0

#def primes(n):
#    primfac = []
#    d = 2
#    while d*d <= abs(n):
#        while (n % d) == 0:
#            primfac.append(d)
#            n //= d
#        d += 1
#    if n > 1:
#       primfac.append(n)
#    return primfac
#
#
#assert primes(2) == [2]
#assert primes(-2) == [2]
#assert primes(3) == [3]
#assert primes(4) == [2,2]
#assert primes(6) == [2,3]


def genpow(idxs, n):
    if n==0:
        yield ()
    else:
        for idx in idxs:
            for tail in genpow(idxs, n-1):
                yield (idx,)+tail



class Space(object):
    "Finite dimensional vector space"

#    _cache = {}
#    def __new__(cls, idxs):
#        if type(idxs) in (int, long):
#            ob = cls._cache.get(idxs)
#            if ob is not None:
#                return ob
#        ob = object.__new__(cls, idxs)
#        if type(idxs) in (int, long):
#            cls._cache[idxs] = ob
#        return ob

    def __init__(self, idxs, factors=None):
        if type(idxs) in (int, long):
            idxs = range(idxs)
        self.idxs = list(idxs)
        self.elems = set(idxs)
        if factors is None:
            factors = [self]
        self.factors = factors

    def __str__(self):
        return "Space(%s)"%(self.idxs,)
    __repr__ = __str__

    # Expose tensor product structure in __len__ and __getitem__
    def __len__(self):
        return len(self.factors)

    def __getitem__(self, i):
        return self.factors[i]

    def __add__(self, other):
        idxs = [(0,idx) for idx in self.idxs] + [(1,idx) for idx in other.idxs]
        return Space(idxs)

    def __mul__(self, other):
        return TensorSpace(self, other)

    def dual(self):
        return self # ???

    def __eq__(self, other):
        return self.idxs == other.idxs

    def __ne__(self, other):
        return self.idxs != other.idxs

    def __pow__(self, n):
        if n==0:
            return self.zero
        idxs = [idx for idx in genpow(self.idxs, n)]
        space = Space(idxs)
        return space

    def identity(self):
        elems = {}
        for i in self.idxs:
            elems[i, i] = 1
        return Operator(elems, self*self.dual())

    def basis(self, i=None):
        if i is None:
            return [Vector({idx : 1}, self) for idx in self.idxs]
        return Vector({self[i] : 1}, self)

    def contains(self, idx):
        return idx in self.elems


Space.zero = Space(0)


class TensorSpace(Space):
    def __init__(self, a, b):
        idxs = []
        for idx in a.idxs:
            for jdx in b.idxs:
                idxs.append((idx, jdx))
        Space.__init__(self, idxs, [a, b])


class Hom(object):
    def __init__(self, src, tgt, f):
        """ f is a dict mapping src.idxs to tgt.idxs""" 
        self.src = src
        self.tgt = tgt
        self.f = f


class Vector(object):
    def __init__(self, elems={}, space=None, name=""):
        assert space is None or isinstance(space, Space), repr(space)
        self.elems = dict(elems) # map (i, j) -> value
        for key in elems.keys():
            assert space.contains(key), "%s not in %s" % (key, space)
        self.support = set(elems.keys())
        self.space = space
        self.name = name

    def check(self):
        for idx in self.elems:
            if not self.space.contains(idx):
                print "%s not in %s" % (idx, self.space)
                assert 0

    def __str__(self):
        if self.name:
            return self.name
        elems = self.elems
        keys = elems.keys()
        keys.sort()
        s = ["%s:%s"%(key, elems[key]) for key in keys]
        #return self.name or "%s(%s)"%(self.__class__.__name__, self.elems)
        return "%s(%s)"%(self.__class__.__name__, ' '.join(s))
    __repr__ = __str__

    def __len__(self):
        return len(self.elems)

    def __eq__(self, other):
        assert self.space == other.space or not other or not self, (self.space, other.space)
        keys = set(self.elems.keys()+other.elems.keys())
        for key in keys:
            if self.elems.get(key, zero) != other.elems.get(key, zero):
                return False
        return True

    def __ne__(self, other):
        return not (self==other)

    def __hash__(self):
        elems = self.elems
        keys = elems.keys()
        keys.sort()
        key = [(k, elems[k]) for k in keys if elems[k]!=zero]
        return hash(tuple(key))

    def __add__(self, other):
        assert self.space == other.space or not other or not self
        keys = set(self.elems.keys()+other.elems.keys())
        elems = {}
        for key in keys:
            value = self.elems.get(key, zero) + other.elems.get(key, zero)
            if value != zero:
                elems[key] = value
        return self.__class__(elems, self.space)

    def __sub__(self, other):
        assert self.space == other.space or not other or not self
        keys = set(self.elems.keys()+other.elems.keys())
        elems = {}
        for key in keys:
            value = self.elems.get(key, zero) - other.elems.get(key, zero)
            if value != zero:
                elems[key] = value
        return self.__class__(elems, self.space)

    def __rmul__(self, r):
        elems = {}
        for key, value in self.elems.items():
            value = r*value
            elems[key] = value
        return self.__class__(elems, self.space)

    def __neg__(self):
        return (-1)*self

    def reduce(self):
        values = [abs(v) for v in self.elems.values()]
        values = list(set(values))
        if not values:
            return self
        factor = reduce(gcd, values)
        elems = {}
        for key, value in self.elems.items():
            elems[key] = value // factor
        return self.__class__(elems, self.space)

    def dot(self, other):
        assert self.space == other.space or not self or not other
        keys = self.support.intersection(other.support)
        value = 0
        #for key, v in self.elems.items():
        for key in keys:
            value += self.elems[key] * other.elems[key]
        return value

    def norm(self):
        value = 0
        for key, v in self.elems.items():
            value += v * v
        return value**0.5

    def eigval(A, B):
        "return val such that B == val*A, or None"
        assert len(A)
        #if A.support != B.support:
        #    return None
        keys = iter(A.support)
        key = keys.next()
        while A.elems[key] == 0:
            key = keys.next()
        a, b = A.elems[key], B.elems.get(key, 0)
        assert b % a == 0
        val = b//a
        if val*A != B:
            return None
        return val

Vector.zero = Vector({}, Space.zero)


class Operator(Vector):
    def __init__(self, elems={}, space=None, name=""):
        assert space is None or isinstance(space, Space), repr(space)
        #_space = space*space.dual() # as a vector, i live here
        assert len(space)==2, space
        Vector.__init__(self, elems, space, name)
        src = space[1].dual()
        tgt = space[0]
        rows = {} # map i -> list of cols [j..]
        cols = {} # map j -> list of rows [i..]
        for (i, j) in elems.keys():
            assert space.contains((i, j))
            assert src.contains(j)
            assert tgt.contains(i)
            rows.setdefault(i, []).append(j)
            cols.setdefault(j, []).append(i)
        self.rows = rows # output space
        self.cols = cols # input space
        self.tgt = tgt # output
        self.src = src # input

    def __call__(self, v):
        if not isinstance(v, Vector):
            raise TypeError
        assert v.space == self.src or not v
        elems = {}
        for i in self.rows.keys():
            elems[i] = 0
        for (i, j), a in self.elems.items():
            elems[i] += a*v.elems.get(j, 0)
        for key in list(elems.keys()):
            if elems[key] == 0:
                del elems[key]
        return Vector(elems, self.tgt)

    def __mul__(self, other):
        if not isinstance(other, Operator):
            raise TypeError
        assert self.src == other.tgt
        elems = {}
        # elems[i,j] = sum_k self[i,k] * other[k,j]
        rows = self.rows
        #cols = other.cols
        src, tgt = self.src, self.tgt
        for i, ks in rows.items():
            for k in ks:
                a = self.elems[i, k]
                for j in other.rows.get(k, []):
                    value = elems.get((i, j), 0) + a * other.elems[k, j]
                    if value:
                        elems[i, j] = value
                    else:
                        del elems[i, j]
        return Operator(elems, self.tgt * other.src.dual(), "(%s)(%s)"%(self.name, other.name) if self.name and other.name else "")

    def tensor(self, other):
        elems = {}
        _elems = other.elems
        for a, v in self.elems.items():
            for b, u in _elems.items():
                elems[((a[0], b[0]), (a[1], b[1]))] = v*u
        space = self.space * other.space
        #print elems
        return Operator(elems, space)

    def directsum(self, other):
        elems = {}
        for (i, j), v in self.elems.items():
            key = (0, i), (0, j)
            elems[key] = v
        for (i, j), v in other.elems.items():
            key = (1, i), (1, j)
            elems[key] = v
        space = (self.space[0] + other.space[0]) * (self.space[1] + other.space[1])
        return Operator(elems, space)

    def lietensor(self, other):
        IA = self.space[0].identity()
        IB = other.space[0].identity()
        assert self*IA == self, (IA, self)
        assert other*IB == other
        return self.tensor(IB) + IA.tensor(other)

        elems = {}
        _elems = other.elems

        for v in self.cols.keys():
          for w in other.cols.keys():
            vw = v + w # tuple!
            for Xv in self.cols[v]:
                #  key = [ output   , input  ]
                key = (Xv+w, vw)
                value = elems.get(key, 0) + self.elems[Xv, v]
                elems[key] = value
                if value == 0:
                    del elems[key]

            for Xw in other.cols[w]:
                #  key = [ output   , input  ]
                key = (v+Xw, vw)
                value = elems.get(key, 0) + other.elems[Xw, w]
                elems[key] = value
                if value == 0:
                    del elems[key]
        return Operator(elems)

    def dual(self):
        elems = {}
        for (i, j), v in self.elems.items():
            elems[j, i] = -v
        space = self.space[1] * self.space[0].dual()
        return Operator(elems, space)

    def bracket(self, other):
        return self*other - other*self

    @classmethod
    def xop(cls, u):
        n = len(u)
        elems = {}
        for idx in genidx((2,)*n):
            jdx = list(idx)
            for i, ii in enumerate(u):
                if ii:
                    jdx[i] = 1-jdx[i]
            jdx = tuple(jdx)
            elems[idx, jdx] = 1
        space = Space(2)**n
        return cls(elems, space*space.dual(), "xop%s"%(u,))
    
    @classmethod
    def eop(cls, u, v):
        "make operator: |u><v|"
        n = len(u)
        assert n==len(v)
        space = Space(2)**n
        u, v = tuple(u), tuple(v)
        assert space.contains(u)
        assert space.contains(v)
        elems = {(u, v):1}
        return cls(elems, space*space.dual(), "|%s><%s|"%(u, v))
    
    @classmethod
    def uop(cls, u):
        "upper operator"
        n = len(u)
        elems = {}
        for idx in genidx((2,)*n):
            jdx = list(idx)
            for i, ii in enumerate(u):
                if ii:
                    jdx[i] = 1-jdx[i]
            jdx = tuple(jdx)
            if idx<jdx:
                elems[idx, jdx] = 1 # UPPER operator
        space = Space(2)**n
        return cls(elems, space*space.dual())
    
    @classmethod
    def lop(cls, u):
        "lower operator"
        n = len(u)
        elems = {}
        for idx in genidx((2,)*n):
            jdx = list(idx)
            for i, ii in enumerate(u):
                if ii:
                    jdx[i] = 1-jdx[i]
            jdx = tuple(jdx)
            if idx>jdx:
                elems[idx, jdx] = 1 # UPPER operator
        space = Space(2)**n
        return cls(elems, space*space.dual())
    
    @classmethod    
    def zop(cls, u):
        n = len(u)
        elems = {}
        for idx in genidx((2,)*n):
            z = 1
            for i, ii in enumerate(u):
                if ii and idx[i]:
                    z *= -1
            elems[idx, idx] = z
        space = Space(2)**n
        return cls(elems, space*space.dual(), "zop%s"%(u,))

Operator.zero = Operator({}, Space.zero*Space.zero.dual())


class Rep(object):
    def __init__(self, ops, mh):
        "first mh elems of ops form basis for Cartan subalgebra"
        assert 0<=mh<len(ops)
        self.ops = list(ops)
        self.mh = mh 
        self.hops = ops[:mh]
        self.eops = ops[mh:]
        space = ops[0].space[0]
        for op in ops:
            assert op.space[0] == space
            op.check()
        self.space = space

    def __len__(self):
        return len(self.ops)

    def __getitem__(self, i):
        return self.ops[i]

    def check(self):
        for A in self.hops:
            for B in self.hops:
                assert A.bracket(B) == Operator.zero
            for C in self.eops:
                if A.bracket(C) != Operator.zero:
                    break
            else:
                assert 0, "Cartan subalgebra is not maximal"

    def tensor(self, other):
        assert len(self.ops)==len(other.ops)
        n = len(self.ops)
        ops = []
        for i in range(n):
            ops.append(self.ops[i].lietensor(other.ops[i]))
        return Rep(ops, self.mh)
    __mul__ = tensor

    def directsum(self, other):
        assert len(self.ops)==len(other.ops)
        n = len(self.ops)
        ops = []
        for i in range(len(self.ops)):
            ops.append(self.ops[i].directsum(other.ops[i]))
        return Rep(ops, self.mh)
    __add__ = directsum

    def dual(self):
        ops = [op.dual() for op in self.ops]
        return Rep(ops, self.mh)

    def vector(self, op):
        v = []
        for _op in self.ops:
            v.append(op.dot(_op))
        return v

    def isomorphic(self, other, verbose=True):
        assert len(self)==len(other)
        n = len(self)
        for i in range(n):
          for j in range(n):
            if i==j:
                continue
            op = self[i].bracket(self[j])
            v = self.vector(op)
            op = other[i].bracket(other[j])
            _v = other.vector(op)
            if verbose:
                print "isomorphic:"
                print v
                print _v
            if v != _v and not verbose:
                return False
        return True



class Sl(Rep):
    def __init__(self, n):
        ops = []
        assert n>1
        space = Space(n)
        print space
        space2 = space*space.dual()
        for i in range(n-1):
            op = Operator({(i, i) : 1, (i+1, i+1) : -1}, space2)
            op.check()
            ops.append(op)
        mh = len(ops)
    
        for i in range(n):
            for j in range(n):
                if i==j:
                    continue
                op = Operator({(i, j) : 1}, space2)
                ops.append(op)
                op.check()
    
        Rep.__init__(self, ops, mh)
        assert len(ops) == n**2-1


def test_sl():

    n = argv.get("n", 3)

    sl = Sl(n)

    #I = Operator.identity(n)
    for A in sl:
      for B in sl:
        lhs = A.lietensor(A).bracket(B.lietensor(B)) 
        rhs = A.bracket(B).lietensor(A.bracket(B))
        assert lhs==rhs
    #print sl.isomorphic(sl + sl, verbose=True)

    rep = sl + sl
    #rep = sl * sl

    for op in rep.hops:
        print op
        assert op.space[0] == rep.space
    print rep.space
    for e in rep.space.basis():
        print e,
        for H in rep.hops:
            #print H, H(e), ',',
            print e.eigval(H(e)),
        print

    for X in rep.eops:
        for H in rep.hops:
            print X.eigval(H.bracket(X)),
        print


def test_1():

    #Rx, Rz = models.build_reduced()
    #print shortstr(Rx)

    #Rx = [Operator.xop(u) for u in Rx]
    #Ru = [Operator.uop(u) for u in Rx]
    #Rl = [Operator.lop(u) for u in Rx]
    #Rz = [Operator.zop(u) for u in Rz]

    Ru = [Operator.uop([1,0]), Operator.uop([0,1])]
    Rl = [Operator.lop([1,0]), Operator.lop([0,1])]

    ZZ = Operator.zop([1,1])
    ZI = Operator.zop([1,0])
    UI, IU = Ru
    LI, IL = Rl

    I = Operator({(0,0):1, (1,1):1})
    U = Operator({(0, 1):1})
    L = Operator({(1, 0):1})

    assert UI == U.tensor(I)

    B = ZI.bracket(UI)
    print UI
    print B
    print B.dot(UI)

    B = ZZ.bracket(UI - IU)
    print UI+IU
    print B
    print B.dot(UI+IU)

    return

    for A in Ru+Rl:
        for z in Rz:
            B = z.bracket(A)
            r = B.dot(A)
            #assert B.support == A.support
            #print B.norm(),
            print len(B.support.intersection(A.support)),
        print

    #search(Ru + Rl + Rz)


def test():

    space = Space(2)
    space2 = space * space.dual()
    U = Operator({(0, 1):1}, space2)
    L = Operator({(1, 0):1}, space2)
    X = Operator({(1, 0):1, (0, 1):1}, space2)
    Z = Operator({(0, 0):1, (1, 1):-1}, space2)

    assert U.space[0] == space
    assert (U+L).space[0] == space
    assert X == U+L
    ZX = Z*X

    a = Vector({0 : 1}, space)
    b = Vector({1 : 1}, space)
    assert U(a) == Vector.zero
    assert U(b) == a
    assert X(a) == b
    assert Z(a) == a
    assert Z(b) == -b, (Z(b), -b)

    assert Z.bracket(X) == 2*ZX
    assert ZX.bracket(X) == 2*Z

    assert ZX + X == 2*U
    assert ZX - X == -2*L

    assert Z.bracket(U) == 2*U
    assert Z.bracket(L) == -2*L

    I = space.identity()
    XX = X.tensor(X)
    ZZ = Z.tensor(Z)
    ZI = Z.tensor(I)

    assert XX.bracket(ZZ) == Operator.zero
    assert XX.bracket(ZI) != Operator.zero
    assert XX.bracket(ZI) == 2*XX*ZI

#    assert Operator.xop([0]).space[0] == space
#    assert Operator.xop([0]) == I
#    assert I.tensor(I).space[0] == (space*space)
#    assert I.tensor(I) == (space*space).identity()

    xop = Operator.xop
    zop = Operator.zop

    assert xop([0,0]) == I.tensor(I)
    assert xop([0,1]) == I.tensor(X)
    assert xop([1,0]) == X.tensor(I)
    assert xop([1,1]) == X.tensor(X)

    assert zop([0,0]) == I.tensor(I)
    assert zop([1,0]) == Z.tensor(I)
    assert zop([0,1]) == I.tensor(Z)
    assert zop([1,1]) == Z.tensor(Z)

#    assert xop([0,1,1]) == I.tensor(X.tensor(X)) # FAILS ASSOC
#    assert zop([1,0,1]) == Z.tensor(I.tensor(Z)) # FAILS ASSOC
    assert xop([0,1,1]) * xop([1,1,0]) == xop([1,0,1])
    assert xop([0,1,1]) * zop([1,0,0]) == zop([1,0,0]) * xop([0,1,1])
    assert xop([0,1,1]) * zop([1,1,0]) == -zop([1,1,0]) * xop([0,1,1])

    UIU = Operator.uop([1,0,1])
    LIL = Operator.lop([1,0,1])
    assert UIU + LIL == Operator.xop([1,0,1])

    assert X.lietensor(X) == X.tensor(I) + I.tensor(X)
    assert X.lietensor(Z) == X.tensor(I) + I.tensor(Z)
    assert Z.lietensor(X) == Z.tensor(I) + I.tensor(X)

    #print Space(5).identity()

    print "OK"


def test_pauli():
    """
        We build sl_{2**n) using pauli operators
    """

    n = argv.get("n", 3)
    space = Space(2)**n
    xop, zop, eop = Operator.xop, Operator.zop, Operator.eop

    ops = []
    pos, neg = [], []
    for u in genidx((2,)*n):
      for v in genidx((2,)*n):
        if u==v:
            continue
        op = eop(u, v) # |u><v|
        #print op
        ops.append(op)
        if v < u:
            pos.append(op)
        else:
            neg.append(op)

    h = [zop(u) for u in genidx((2,)*n) if sum(u)]

    assert len(h) + len(ops) == (2**n)*(2**n)-1

    basis = []
    for u in genidx((2,)*n):
      for v in genidx((2,)*n):
        if sum(u)==sum(v)==0:
            continue
        op = xop(u) * zop(v)
        basis.append(op)
    assert len(basis) == (2**n)*(2**n)-1

    for op1 in ops:
        #print op1
        op = None
        for op2 in basis:
            r = op1.dot(op2)
            if r:
                #print r, op2
                if op is None:
                    op = r*op2
                else:
                    op = op + r*op2
            #print "op = ", op
        #print
        assert op == (2**n)*op1

    v = Vector({(0,)*n : 1}, space)
    for op in neg:
        assert op(v) == Vector.zero
    #for op in pos:
    #    print op(v)


class LISet(object):
    "linear independent set, row-reduced"
    def __init__(self, space, ops):
        self.space = space
        self.ops = []
        self.idxs = {} # map idx -> op
        for op in ops:
            assert op.space == space
            self.add(op)

    def __len__(self):
        return len(self.idxs)

    def __getitem__(self, i):
        return self.ops[i]

    def add(self, op):
        if op == Operator.zero:
            return
        #print "add", op
        for idx in self.space.idxs:
            a = op.elems.get(idx)
            if a is None or a==0:
                continue
            op1 = self.idxs.get(idx)
            if op1 is None:
                if a < 0:
                    op = -op
                self.idxs[idx] = op
                self.ops.append(op) # slightly redundant
                return True # <-------- return
            b = op1.elems[idx]
            assert b > 0
            if a < 0:
                op = -op
                a = -a
            c = a*b / gcd(a, b)
            assert c%a == 0
            assert c%b == 0
            if c>a:
                op = (c//a)*op
            else:
                assert c==a
            assert op.elems[idx] == c, (op.elems[idx] , c)
            if c>b:
                op1 = (c//b)*op1
            else:
                assert c==b
            assert op1.elems[idx] == c, (op1.elems[idx] , c)
            op = op - op1 # row reduce op
            assert op.elems.get(idx, 0) == 0, op.elems.get(idx,0)
            if op == Operator.zero:
                break
        return False


def closure(ops):
    "lie bracket closure"
    op = ops[0]
    space = op.space
    li = LISet(space, ops)

    while 1:
        sz = len(li)
        print sz
        items = list(li)
        for A in items:
          for B in items:
            if A==B:
                continue
            C = A.bracket(B)
            if li.add(C):
                write('.')
        if len(li)==sz:
            break
    print len(li)
    return li


def main():

    xop, zop, eop = Operator.xop, Operator.zop, Operator.eop

    ops = [xop([1,0]), xop([0,1]), zop([1,0]), zop([0,1])]

    closure(ops)


def test_model():

    Gx, Gz, Hx, Hz = models.build()

    Px = get_reductor(Hx) # projector onto complement of rowspan of Hx
    Pz = get_reductor(Hz)
    Rz = [dot2(Pz, g) for g in Gz] 
    Rz = array2(Rz)
    Rz = row_reduce(Rz, truncate=True)

    Rx = [dot2(Px, g) for g in Gx] 
    Rx = array2(Rx)
    Rx = row_reduce(Rx, truncate=True)

    n = Rx.shape[1]

    print shortstrx(Rx, Rz)

    r, n = Rx.shape
    N = 2**r
    gz = len(Gz)

    Qx = u_inverse(Rx)
    RR = dot2(Gz, Rx.transpose())
    PxtQx = dot2(Px.transpose(), Qx)

    if N<=1024:
        H = numpy.zeros((N, N))
    else:
        H = None
    A = {}
    U = []

    basis = []
    lookup = {}
    for i, v in enumerate(genidx((2,)*r)):
        v = array2(v)
        lookup[v.tostring()] = i
        basis.append(v)

    space = Space(len(basis))
    space2 = space*space.dual()

    ops = []

    # zops
    for gz in Gz:
        elems = {}
        for i, v in enumerate(basis):
            bit = dot2(gz, Rx.transpose(), v)
            elems[(i, i)] = 1 - 2*bit # send 0,1 -> +1,-1
        ops.append(Operator(elems, space2))

    # xops
    for gx in Gx:
        elems = {}
        for i, v in enumerate(basis):
            u = (v+dot2(gx, PxtQx))%2
            j = lookup[u.tostring()]
            elems[j, i] = 1
        for (i, j) in elems:
            assert elems[i, j] == elems[j, i]
        ops.append(Operator(elems, space2))

    print "ops:", len(ops)
    print ops[0]

    closure(ops)



from argv import Argv

argv = Argv()

if __name__ == "__main__":

    test()

    if argv.profile:
        import cProfile as profile
        profile.run("main()")

    else:
        fn = argv.next()
        fn = eval(fn)
        fn()


