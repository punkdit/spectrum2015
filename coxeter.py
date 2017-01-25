#!/usr/bin/env python

"""
Represent coxeter group as reflection group of a root system.

See: Humphreys, p63-66
"""


import sys, os


from action import Perm, Group, mulclose
from argv import argv



def factorial(n):
    r = 1
    for i in range(1, n+1):
        r *= i
    return r
 

def choose(m, n):
    return factorial(m) // factorial(n)


class Weyl(object):
    def __init__(self, roots, gen):
        self.roots = roots # list of tuples
        self.gen = gen # list of Weyl group generators (Perm's)
        self.identity = Perm.identity(roots)
        for g in gen:
            assert g*g == self.identity

    @classmethod
    def build_A(cls, n):
    
        roots = []
        lookup = {}
        for i in range(n+1):
          for j in range(n+1):
    
            if i==j:
                continue
            root = [0]*(n+1)
            root[i] = 1
            root[j] = -1
    
            root = tuple(root)
            lookup[root] = len(roots)
            roots.append(root)
    
        #assert len(pos_roots) == choose(n+1, 2) # number of positive roots
        assert len(lookup) == len(roots)
    
        gen = []
        for i in range(n):
            perm = []
            for idx, root in enumerate(roots):
                _root = list(root)
                _root[i], _root[i+1] = _root[i+1], _root[i] # swap i, i+1 coords
                jdx = lookup[tuple(_root)]
                perm.append(jdx)
            perm = Perm(perm, roots)
            gen.append(perm)
        return cls(roots, gen)

    @classmethod
    def build_B(cls, n, k=1):
    
        roots = []
        lookup = {}

        # long roots
        for i in range(n):
          for j in range(i+1, n):
    
            for a in [-1, 1]:
             for b in [-1, 1]:
                root = [0]*n
                root[i] = a
                root[j] = b
                root = tuple(root)
                lookup[root] = len(roots)
                roots.append(root)

        # short roots
        for i in range(n):
          for a in [-1, 1]:
            root = [0]*n
            root[i] = a*k
            root = tuple(root)
            lookup[root] = len(roots)
            roots.append(root)
        assert len(lookup) == len(roots)
    
        gen = []
        for i in range(n-1):
            perm = []
            for idx, root in enumerate(roots):
                _root = list(root)
                _root[i], _root[i+1] = _root[i+1], _root[i] # swap i, i+1 coords
                jdx = lookup[tuple(_root)]
                perm.append(jdx)
            perm = Perm(perm, roots)
            gen.append(perm)

        perm = []
        for idx, root in enumerate(roots):
            _root = list(root)
            _root[n-1] = -_root[n-1]
            jdx = lookup[tuple(_root)]
            perm.append(jdx)
        perm = Perm(perm, roots)
        gen.append(perm)

        return cls(roots, gen)

    @classmethod
    def build_C(cls, n):
        return cls.build_B(n, k=2)

    @classmethod
    def build_D(cls, n):

        assert n>=2

        roots = []
        lookup = {}

        for i in range(n):
          for j in range(i+1, n):
    
            for a in [-1, 1]:
             for b in [-1, 1]:
                root = [0]*n
                root[i] = a
                root[j] = b
                root = tuple(root)
                lookup[root] = len(roots)
                roots.append(root)

        assert len(lookup) == len(roots)
        assert len(roots) == 2*n*(n-1)

        gen = []
        for i in range(n-1):
            perm = []
            for idx, root in enumerate(roots):
                _root = list(root)
                _root[i], _root[i+1] = _root[i+1], _root[i] # swap i, i+1 coords
                jdx = lookup[tuple(_root)]
                perm.append(jdx)
            perm = Perm(perm, roots)
            gen.append(perm)

        perm = []
        for idx, root in enumerate(roots):
            _root = list(root)
            _root[n-1], _root[n-2] = -_root[n-2], -_root[n-1] # swap & negate last two coords
            jdx = lookup[tuple(_root)]
            perm.append(jdx)
        perm = Perm(perm, roots)
        gen.append(perm)

        return cls(roots, gen)



def test(n):


    G = Weyl.build_A(n)

#    print "%d roots" % len(G.roots)
#    print G.roots
#    for g in G.gen:
#        print [g(root) for root in G.roots]

    gen = G.gen
    for i in range(n):
      for j in range(n):
        gi = gen[i]
        gj = gen[j]
        if i==j:
            assert gi*gj == G.identity
        elif abs(i-j)==1:
            assert gi*gj != G.identity
            assert (gi*gj)**2 != G.identity
            assert (gi*gj)**3 == G.identity
        else:
            assert gi*gj != G.identity
            assert (gi*gj)**2 == G.identity

    if n < 5:
        assert len(mulclose(G.gen)) == factorial(n+1)

    # ---------------------------------------------------------

    G = Weyl.build_B(n)

    gen = G.gen
    for g in gen:
        assert g*g == G.identity

    for i in range(n-1):
      for j in range(i+1, n-1):
        gi = gen[i]
        gj = gen[j]
        if abs(i-j)==1:
            assert gi*gj != G.identity
            assert (gi*gj)**2 != G.identity
            assert (gi*gj)**3 == G.identity
        else:
            assert gi*gj != G.identity
            assert (gi*gj)**2 == G.identity
        if i < n-1:
            gj = gen[n-1]
            assert gi*gj != G.identity
            assert (gi*gj)**2 == G.identity

    if n>2:
        gi = gen[n-2]
        gj = gen[n-1]
        assert (gi*gj) != G.identity
        assert (gi*gj)**2 != G.identity
        assert (gi*gj)**3 != G.identity
        assert (gi*gj)**4 == G.identity

    if n < 5:
        assert len(mulclose(G.gen)) == (2**n)*factorial(n)

    # ---------------------------------------------------------
    G = Weyl.build_C(n)

    # ---------------------------------------------------------

    if n<3:
        return # <---------------------- return

    G = Weyl.build_D(n)

    gen = G.gen
    for i in range(n-1):
        gi = gen[i]
        for j in range(i+1, n-1):
            gj = gen[j]
            if abs(i-j)==1:
                assert gi*gj != G.identity
                assert (gi*gj)**2 != G.identity
                assert (gi*gj)**3 == G.identity
            else:
                assert gi*gj != G.identity
                assert (gi*gj)**2 == G.identity

        gj = gen[n-1]
        if i < n-3 or i==n-2:
            assert gi*gj != G.identity
            assert (gi*gj)**2 == G.identity
        elif i==n-3:
            assert gi*gj != G.identity
            assert (gi*gj)**2 != G.identity
            assert (gi*gj)**3 == G.identity

    if n < 5:
        assert len(mulclose(G.gen)) == (2**(n-1))*factorial(n)


def main():

    n = argv.n
    if n:
        test(n)

    else:
        for n in range(2, 6):
            test(n)


if __name__ == "__main__":

    main()






