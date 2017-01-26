#!/usr/bin/env python

"""
Represent coxeter group as reflection group of a root system.

See: 
    Humphreys, p63-66
    https://en.wikipedia.org/wiki/Root_system
"""


import sys, os
from fractions import Fraction
from operator import mul


from action import Perm, Group, mulclose
from argv import argv


def cross(itemss):
    if len(itemss)==0:
        yield ()
    else:
        for head in itemss[0]:
            for tail in cross(itemss[1:]):
                yield (head,)+tail


def factorial(n):
    r = 1
    for i in range(1, n+1):
        r *= i
    return r
 

def choose(m, n):
    return factorial(m) // factorial(n)


class Weyl(object):
    def __init__(self, roots, gen, simple=None):
        self.roots = roots # list of tuples
        self.gen = gen # list of Weyl group generators (Perm's)
        self.identity = Perm.identity(roots)
        self.simple = simple
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

    @classmethod
    def build_simple(cls, roots, simple):
        "use generators from reflections of simple roots"
        for root in simple:
            assert root in roots
        n = len(roots[0])
        idxs = range(n)
        lookup = dict((root, i) for (i, root) in enumerate(roots))
        gen = []
        for alpha in simple:
            #print "alpha:", alpha
            r0 = sum(alpha[i]*alpha[i] for i in idxs)
            perm = []
            for root in roots:
                #print "    root:", root
                r = sum(alpha[i]*root[i] for i in idxs)
                _root = tuple(root[i] - 2*Fraction(r, r0)*alpha[i] for i in idxs)
                #print "    _root:", _root
                perm.append(lookup[_root])
            perm = Perm(perm, roots)
            gen.append(perm)
        return cls(roots, gen, simple)

    @classmethod
    def build_E8(cls):
        D8 = cls.build_D(8)
        half = Fraction(1, 2)
        roots = list(D8.roots)
        for signs in cross([(-1, 1)]*8):
            r = reduce(mul, signs)
            if r != 1:
                continue
            root = tuple(sign*half for sign in signs)
            roots.append(root)

        assert len(roots) == 240

        simple = []
        for i in range(6):
            root = [0]*8
            root[i] = 1
            root[i+1] = -1
            simple.append(tuple(root))

        root = [0]*8
        root[i] = 1
        root[i+1] = 1
        simple.append(tuple(root))

        root = [-half]*8
        simple.append(tuple(root))

        return cls.build_simple(roots, simple)

    @classmethod
    def build_E7(cls):
        E8 = cls.build_E8()
        idxs = range(8)
        # delete one root:
        root0 = E8.simple[0]
        roots = []
        for root in E8.roots:
            if sum(root0[i]*root[i] for i in idxs)==0:
                roots.append(root)
        assert len(roots)==126
        simple = [roots[0]] + [root for root in E8.simple if root in roots]
        return cls.build_simple(roots, simple)

    @classmethod
    def build_E6(cls):
        E8 = cls.build_E8()
        idxs = range(8)
        # delete two roots:
        root0 = E8.simple[0]
        root1 = E8.simple[1]
        roots = []
        for root in E8.roots:
            if sum(root0[i]*root[i] for i in idxs)==0 and \
                sum(root1[i]*root[i] for i in idxs)==0:
                roots.append(root)
        assert len(roots)==72
        simple = [(0, 0, 0, 1, 0, 0, 0, -1)]
        simple += [root for root in E8.simple if root in roots]
        return cls.build_simple(roots, simple)

    @classmethod
    def build_F4(cls):
        roots = []
        idxs = range(4)
        for root in cross([(-1,0,1)]*4):
            d = sum(root[i]**2 for i in idxs)
            if d==1 or d==2:
                roots.append(root)
        half = Fraction(1, 2)
        for root in cross([(-half,0,half)]*4):
            d = sum(root[i]**2 for i in idxs)
            if d==1 or d==2:
                roots.append(root)
        assert len(roots)==48
        simple = [
            (1, -1, 0, 0),
            (0, 1, -1, 0),
            (0, 0, 1, 0),
            (-half, -half, -half, -half)]
        return cls.build_simple(roots, simple)

    @classmethod
    def build_G2(cls):
        roots = []
        for root in cross([(-2, -1, 0, 1, 2)]*3):
            if sum(root) != 0:
                continue
            d = sum(root[i]**2 for i in range(3))
            if d==2 or d==6:
                roots.append(root)
        print roots
        assert len(roots)==12
        simple = [(1, -1, 0), (-1, 2, -1)]
        return cls.build_simple(roots, simple)

    def matrix(self):
        "coxeter matrix"
        m = {}
        gen = self.gen
        n = len(gen)
        for i in range(n):
          for j in range(i+1, n):
            gi = gen[i]
            gj = gen[j]
            g = gigj = gi*gj
            k = 1
            while 1:
                if g == self.identity:
                    if k > 2:
                        m[i, j] = k
                    break
                g = gigj * g
                k += 1
                assert k<10
        return m




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

    #else:
    #    for n in range(2, 6):
    #        test(n)

    if argv.E8:
        E8 = Weyl.build_E8()
        print E8.matrix()

    if argv.E7:
        E7 = Weyl.build_E7()
        print E7.matrix()

    if argv.E6:
        E6 = Weyl.build_E6()
        print E6.matrix()
        #items = mulclose(E6.gen, verbose=True)
        #print "|E6|=", len(items)

    if argv.F4:
        F4 = Weyl.build_F4()
        print F4.matrix()

    if argv.G2:
        G2 = Weyl.build_G2()
        print G2.matrix()
        assert len(mulclose(G2.gen)) == 12 # dihedral group D_6


if __name__ == "__main__":

    if argv.profile:
        import cProfile as profile
        profile.run("main()")

    else:

        main()






