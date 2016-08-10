#!/usr/bin/env python

import sys

import numpy
from numpy import linalg as la
from numpy import dot, kron, eye

from code import texstr
from lanczos import show_eigs, write
from isomorph import from_ham, search

from argv import Argv
argv = Argv()

"""
See: Fulton chapter 11.

This is representations of sl_2.
"""

EPSILON = 1e-10

I = numpy.array([[1, 0], [0, 1.]])
H = numpy.array([[1, 0], [0, -1.]]) # diagonal
X = numpy.array([[0, 1], [0, 0.]]) # raising
Y = numpy.array([[0, 0], [1, 0]]) # lowering


def is_square(A):
    return len(A.shape)==2 and A.shape[0] == A.shape[1]

def bracket(A, B):
    assert is_square(A) and is_square(B)
    return dot(A, B) - dot(B, A)

def directsum(A, B):
    assert is_square(A) and is_square(B)
    a, b = A.shape[0], B.shape[0]
    n = a+b
    C = numpy.zeros((n, n))
    C[:a,:a] = A
    C[a:,a:] = B
    return C

def directrmul(A, r):
    assert r>=1
    assert is_square(A)
    a = A.shape[0]
    n = r*a
    C = numpy.zeros((n, n))
    for i in range(r):
        C[i*a:(i+1)*a, i*a:(i+1)*a] = A
    return C

def tensor(A, B):
    assert is_square(A) and is_square(B)
    a, b = A.shape[0], B.shape[0]
    C = kron(A, eye(b)) + kron(eye(a), B)
    return C


class Repr(object):
    def __init__(self, H, X, Y):
        assert H.shape == X.shape == Y.shape
        self.H = H
        self.X = X
        self.Y = Y
        self.n = len(H)

    def __str__(self):
        vals = [self.H[i, i] for i in range(self.n)]
        uniq = list(set(vals))
        uniq.sort(reverse=True)
        s = ', '.join("%d(%d)"%(val, vals.count(val)) for val in uniq)
        return "Repr(%s)"%s
            

    def check(self):
        H, X, Y = self.H, self.X, self.Y
        assert abs(H.trace()) < EPSILON
        assert abs(X.trace()) < EPSILON
        assert abs(Y.trace()) < EPSILON
        assert numpy.allclose(bracket(H, X), 2*X)
        assert numpy.allclose(bracket(H, Y), -2*Y)
        assert numpy.allclose(bracket(X, Y), H)

    @classmethod
    def standard(cls, n):
        assert n>=1
        H = numpy.zeros((n, n))
        X = numpy.zeros((n, n))
        Y = numpy.zeros((n, n))
        k = n-1
        for i in range(n):
            H[i, i] = k
            k -= 2
        for i in range(n-1):
            X[i, i+1] = i+1
            Y[i+1, i] = n-1-i

        return cls(H, X, Y)

    def __add__(self, other):
        H = directsum(self.H, other.H)
        X = directsum(self.X, other.X)
        Y = directsum(self.Y, other.Y)
        return Repr(H, X, Y)

    def __rmul__(self, r):
        H = directrmul(self.H, r)
        X = directrmul(self.X, r)
        Y = directrmul(self.Y, r)
        return Repr(H, X, Y)

    def __mul__(self, other):
        "tensor product"
        H = tensor(self.H, other.H)
        X = tensor(self.X, other.X)
        Y = tensor(self.Y, other.Y)
        return Repr(H, X, Y)

    def __pow__(self, r):
        assert r>=1
        rep = self
        for i in range(r-1):
            rep = self*rep
        return rep


def build_orbigraph(A):
    import networkx as nx

    n = len(A)
    graph = nx.Graph()
    for i in range(n):
        graph.add_node(i)

    bag0 = from_ham(A)
    bag1 = from_ham(A)

    count = 0
    fs = set()
    for fn in search(bag0, bag1):
        f = [None]*n
        #print fn
        for i, j in fn.items():
            if i>=n:
                continue
            assert i<n and j<n
            f[i] = j
            graph.add_edge(i, j)
        f = tuple(f)
        if f in fs:
            #write('/')
            continue # <---- continue
        fs.add(f)
        #write('.')
        count += 1
    print
    print "isomorphisms:", count

    equs = nx.connected_components(graph)
    m = len(equs)

    print "components:", m

    P = numpy.zeros((n, m))
    Q = numpy.zeros((m, n))
    for i, equ in enumerate(equs):
        for j in equ:
            P[j, i] = 1
        Q[i, j] = 1

    #print shortstr(P)
    #print shortstr(Q)

    A = numpy.dot(Q, numpy.dot(A, P))
    return A


def main():

    n = argv.get("n", 7)

    rep = Repr.standard(n)
    rep.check()

    rep = Repr.standard(3) + Repr.standard(4)
    rep.check()
    #print rep.H

    rep = Repr.standard(3) * Repr.standard(4)
    rep.check()

    sym = lambda n : Repr.standard(n+1)

    s2 = sym(1)
    rep = s2 ** 6
    print rep

#    rep = sym(6)+5*sym(4)+9*sym(2)+5*sym(0)

    A = rep.H + rep.X + rep.Y
    print A
    vals, vecs = la.eig(A)

    r2 = 2**0.5
    show_eigs(vals/r2)

    A = build_orbigraph(A)
    print A
    vals, vecs = la.eig(A)
    show_eigs(vals/r2)

    print "OK"



if __name__=="__main__":
    main()

    


