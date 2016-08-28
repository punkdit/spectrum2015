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





class Operator(object):
    def __init__(self, elems={}):
        self.elems = dict(elems) # map (i, j) -> value
        rows = {} # map i -> list of cols [j..]
        cols = {} # map j -> list of rows [i..]
        for (i, j) in elems.keys():
            rows.setdefault(i, []).append(j)
            cols.setdefault(j, []).append(i)
        self.rows = rows
        self.cols = cols

    def __str__(self):
        return "Operator(%s)"%(self.elems)
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
        return Operator(elems)

    def __sub__(self, other):
        keys = set(self.elems.keys()+other.elems.keys())
        elems = {}
        for key in keys:
            value = self.elems.get(key, zero) - other.elems.get(key, zero)
            if value != zero:
                elems[key] = value
        return Operator(elems)

    def __rmul__(self, r):
        elems = {}
        for key, value in self.elems.items():
            value = r*value
            elems[key] = value
        return Operator(elems)

    def __neg__(self):
        return (-1)*self

    def __mul__(self, other):
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

    def promote(self):
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
        self = self.promote()
        other = other.promote()
        elems = {}
        _elems = other.elems
        for a, v in self.elems.items():
            for b, u in _elems.items():
                elems[(a[0]+b[0], a[1]+b[1])] = v*u
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

    def reduce(self):
        values = [abs(v) for v in self.elems.values()]
        values = list(set(values))
        if not values:
            return self
        factor = reduce(gcd, values)
        elems = {}
        for key, value in self.elems.items():
            elems[key] = value // factor
        return Operator(elems)

    def dot(self, other):
        value = 0
        for key, v in self.elems.items():
            value += v * other.elems.get(key, 0)
        return value

    def norm(self):
        value = 0
        for key, v in self.elems.items():
            value += v * v
        return value**0.5


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
    assert X == U+L
    Z = Operator({(0, 0):1, (1, 1):-1})
    ZX = Z*X

    assert Z.bracket(X) == 2*ZX
    assert ZX.bracket(X) == 2*Z

    assert ZX + X == 2*U
    assert ZX - X == -2*L

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

    print "OK"


def main():

    Rx, Rz = models.build_reduced()

    print shortstr(Rx)

    #Rx = [Operator.xop(u) for u in Rx]
    Ru = [Operator.uop(u) for u in Rx]
    Rl = [Operator.lop(u) for u in Rx]
    Rz = [Operator.zop(u) for u in Rz]

    ops = Ru + Rl + Rz
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


from argv import Argv

argv = Argv()

if __name__ == "__main__":

    test()
    main()




