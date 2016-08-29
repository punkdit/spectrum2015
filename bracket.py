#!/usr/bin/env python

import sys, os
from fractions import gcd

#import numpy

import networkx as nx
from solve import get_reductor, array2, row_reduce, dot2, shortstr, zeros2

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





class Vector(object):
    def __init__(self, elems={}):
        self.elems = dict(elems) # map (i, j) -> value
        self.support = set(elems.keys())

    def __str__(self):
        return "%s(%s)"%(self.__class__.__name__, self.elems)
    __repr__ = __str__

    def __len__(self):
        return len(self.elems)

    def __eq__(self, other):
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
        keys = set(self.elems.keys()+other.elems.keys())
        elems = {}
        for key in keys:
            value = self.elems.get(key, zero) + other.elems.get(key, zero)
            if value != zero:
                elems[key] = value
        return self.__class__(elems)

    def __sub__(self, other):
        keys = set(self.elems.keys()+other.elems.keys())
        elems = {}
        for key in keys:
            value = self.elems.get(key, zero) - other.elems.get(key, zero)
            if value != zero:
                elems[key] = value
        return self.__class__(elems)

    def __rmul__(self, r):
        elems = {}
        for key, value in self.elems.items():
            value = r*value
            elems[key] = value
        return self.__class__(elems)

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
        return self.__class__(elems)

    def dot(self, other):
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



class Operator(Vector):
    def __init__(self, elems={}):
        Vector.__init__(self, elems)
        rows = {} # map i -> list of cols [j..]
        cols = {} # map j -> list of rows [i..]
        for (i, j) in elems.keys():
            rows.setdefault(i, []).append(j)
            cols.setdefault(j, []).append(i)
        self.rows = rows # output space
        self.cols = cols # input space

    def __call__(self, v):
        if not isinstance(v, Vector):
            raise TypeError
        elems = {}
        for i in self.rows.keys():
            elems[i] = 0
        for (i, j), a in self.elems.items():
            elems[i] += a*v.elems.get(j, 0)
        for key in list(elems.keys()):
            if elems[key] == 0:
                del elems[key]
        return Vector(elems)

    def __mul__(self, other):
        if not isinstance(other, Operator):
            raise TypeError
        elems = {}
        # elems[i,j] = sum_k self[i,k] * other[k,j]
        rows = self.rows
        #cols = other.cols
        for i, ks in rows.items():
            for k in ks:
                a = self.elems[i, k]
                for j in other.rows.get(k, []):
                    value = elems.get((i, j), 0) + a * other.elems[k, j]
                    if value:
                        elems[i, j] = value
                    else:
                        del elems[i, j]
        return Operator(elems)

    def _promote(self):
        # use tuples for row and col indexes 
        keys = self.elems.keys()
        if not keys:
            return self
        (i, j) = keys[0]
        assert type(i)==type(j)
        if type(i) is tuple:
            return self
        elems = {}
        for (i, j), value in self.elems.items():
            i = (i,)
            j = (j,)
            elems[i, j] = value
        return Operator(elems)

    def tensor(self, other):
        self = self._promote()
        other = other._promote()
        elems = {}
        _elems = other.elems
        for a, v in self.elems.items():
            for b, u in _elems.items():
                elems[(a[0]+b[0], a[1]+b[1])] = v*u
        return Operator(elems)

    def lietensor(self, other):
        self = self._promote()
        other = other._promote()
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
        return Operator(elems)

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
        return cls(elems)
    
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
        return cls(elems)
    
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
        return cls(elems)
    
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
        return cls(elems)


def shrink(ops, A):
    changed = True
    while changed:
        changed = False
        for B in ops:
            r = A.dot(B)
            if r==0:
                continue
            if r>0:
                A1 = A - B
            else:
                A1 = A + B
            if A1.norm() < A.norm():
                A = A1
                changed = True
    return A


def test():
    U = Operator({(0, 1):1})
    L = Operator({(1, 0):1})
    X = Operator({(1, 0):1, (0, 1):1})
    Z = Operator({(0, 0):1, (1, 1):-1})

    assert X == U+L
    ZX = Z*X

    a = Vector({0 : 1})
    b = Vector({1 : 1})
    ZERO = Vector()
    assert U(a) == ZERO
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

    I = Operator({(0, 0):1, (1, 1):1})
    XX = X.tensor(X)
    ZZ = Z.tensor(Z)
    ZI = Z.tensor(I)

    ZERO = Operator()
    assert XX.bracket(ZZ) == ZERO
    assert XX.bracket(ZI) != ZERO
    assert XX.bracket(ZI) == 2*XX*ZI

    assert Operator.xop([0,0]) == I.tensor(I)
    assert Operator.xop([0,1]) == I.tensor(X)
    assert Operator.xop([1,0]) == X.tensor(I)
    assert Operator.xop([1,1]) == X.tensor(X)

    assert Operator.zop([0,0]) == I.tensor(I)
    assert Operator.zop([1,0]) == Z.tensor(I)
    assert Operator.zop([0,1]) == I.tensor(Z)
    assert Operator.zop([1,1]) == Z.tensor(Z)

    assert Operator.xop([0,1,1]) == I.tensor(X.tensor(X))
    assert Operator.zop([1,0,1]) == Z.tensor(I.tensor(Z))

    UIU = Operator.uop([1,0,1])
    LIL = Operator.lop([1,0,1])
    assert UIU + LIL == Operator.xop([1,0,1])

    assert (5*XX).reduce() == XX
    assert (5*ZI).reduce() == ZI

    assert shrink([XX, ZI], XX + ZI) == ZERO
    assert shrink([XX, ], XX + ZI) == ZI
    assert shrink([XX, ], 2*XX - ZI) == -ZI


    assert X.lietensor(X) == X.tensor(I) + I.tensor(X)
    assert X.lietensor(Z) == X.tensor(I) + I.tensor(Z)
    assert Z.lietensor(X) == Z.tensor(I) + I.tensor(X)

    print "OK"


class Rep(object):
    def __init__(self, ops):
        self.ops = list(ops)


class sl(Rep):
    def __init__(self, n):
        ops = []
        Rep.__init__(self, ops)


def main():

    n = argv.get("n", 2)

    h = []
    assert n>1
    for i in range(n-1):
        op = Operator({(i, i) : 1, (i+1, i+1) : -1})
        h.append(op)

    basis = [Vector({i:1}) for i in range(n)]
    for e in basis:
        print e,
        for H in h:
            #print H, H(e), ',',
            print e.eigval(H(e)),
        print

    if argv.adjoint:
        ops = []
        for i in range(n):
            for j in range(n):
                if i==j:
                    continue
                ops.append(Operator({(i, j) : 1}))

        for X in ops:
            for H in h:
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


def search(ops):
    found = set(ops)

    pairs = [(A, B) for A in ops for B in ops]
    while 1:
        new = []
        for A, B in pairs:
            C = A.bracket(B)
            C = C.reduce()
            if C not in found and -C not in found:
                new.append(C)
                #found.add(C)
                #print C
        print "new:", len(new),;sys.stdout.flush()
        _new = []
        for A in new:
            A = shrink(ops, A)
            #print "shrink:", len(A)
            if A not in found and -A not in found:
                _new.append(A)
                found.add(A)
            #else:
            #    print "*",
        new = _new
        pairs = [(A, B) for A in ops for B in new]
        print "new:", len(new), "pairs:", len(pairs)
        if not new:
            break
        ops.extend(new)

    print "ops:", len(ops)


from argv import Argv

argv = Argv()

if __name__ == "__main__":

    test()
    main()




