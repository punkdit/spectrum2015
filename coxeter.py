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


def rdot(a, b):
    assert len(a)==len(b)
    r = sum(ai*bi for (ai,bi) in zip(a, b))
    return r


def rnorm2(a):
    r = sum(ai*ai for ai in a)
    return r


def rscale(v, a):
    assert Fraction(v)==v
    assert type(a) is tuple
    a = tuple(v*ai for ai in a)
    return a


class Weyl(object):
    def __init__(self, roots, gen, simple=None):
        self.roots = roots # list of tuples
        self.gen = gen # list of Weyl group generators (Perm's)
        self.identity = Perm.identity(roots)
        self.simple = simple
        for g in gen:
            assert g*g == self.identity

        #if simple is not None:
        #    self.check_simple()

    def check_simple(self):
        """
        A set of simple roots for a root system \Phi is
        a set of roots that form a basis for the Euclidean
        space spanned by \Phi with the special property that each
        root has components with respect to this basis that are
        either all nonnegative or all nonpositive.
        
        https://en.wikipedia.org/wiki/E8_(mathematics)#Simple_roots
        """

        simple = self.simple
        print "simple:", simple
        roots = self.roots
        n = len(simple)
        for root in roots:
            if root in simple or rscale(-1, root) in simple:
                continue
            print "root:", root
        print "OK"

    @classmethod
    def build_A(cls, n):
    
        roots = []
        simple = []
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
            if j==i+1:
                simple.append(root)
    
        #assert len(pos_roots) == choose(n+1, 2) # number of positive roots
        assert len(lookup) == len(roots)
        assert len(roots) == n*(n+1)
    
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
        return cls(roots, gen, simple)

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
        assert len(roots) == 2*n**2
    
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
        assert len(roots)==12
        simple = [(1, -1, 0), (-1, 2, -1)]
        return cls.build_simple(roots, simple)

    def matrix(self, desc=None):
        "coxeter matrix"
        m = {}
        gen = self.gen
        n = len(gen)
        for i in range(n):
          for j in range(i+1, n):
            key = (i, j)
            if desc:
                key = (desc[i], desc[j])
            gi = gen[i]
            gj = gen[j]
            g = gigj = gi*gj
            k = 1
            while 1:
                if g == self.identity:
                    if k > 2:
                        m[key] = k
                    break
                g = gigj * g
                k += 1
                assert k<10
        return m



def mulclose_pri(els, verbose=False, maxsize=None):
    els = set(els)
    changed = True
    while changed:
        if verbose:
            print "mulclose:", len(els)
        changed = False
        _els = list(els)
        pairs = [(g, h) for g in _els for h in _els]
        pairs.sort(key = lambda (g,h) : len(g.word+h.word))

        for A, B in pairs:
            C = A*B 
            if C not in els:
                els.add(C)
                if maxsize and len(els)>=maxsize:
                    return list(els)
                changed = True
    return els 





def test_monoid(G):
    "build the Coxeter-Bruhat monoid, represented as a monoid of functions."
    "this representation is not faithful."

    roots = G.roots
    print "roots:", len(roots)
#    print roots
#    print

    names = 'ABCDEF'
    gen = G.gen
    for i, g in enumerate(gen):
        g.word = names[i]
        print "%s:"%g.word,
        print g.str()
    print
    n = len(gen)
    identity = Perm(dict((r, r) for r in roots), roots, '')
    weyl = mulclose_pri([identity]+gen)
    weyl = list(weyl)
    weyl.sort(key = lambda g : (len(g.word), g.word))
    print "weyl:",
    for w in weyl:
        print w.word,
    print
    print "weyl:", len(weyl)

    r0 = roots[0]
    bdy = set([r0])
    seen = set()
    identity = dict((r, r) for r in roots)
    perms = [dict(identity) for g in gen]
    while bdy:

        _bdy = set()
        seen.update(bdy)
        for r0 in bdy:
          for i in range(n):
            g = gen[i]
            perm = perms[i]
            r1 = g*r0
            assert perm.get(r0) == r0
            if r1 not in seen:
                perm[r0] = r1
                _bdy.add(r1)

        bdy = _bdy

    gen = [Perm(perms[i], roots, gen[i].word) for i in range(len(perms))]
    identity = Perm(identity, roots)
    for g in gen:
        print g.str()
        assert g*g == g
    
    monoid = mulclose_pri([identity]+gen)
    print "monoid:", len(monoid)
    #monoid = mulclose(monoid)
    #print "monoid:", len(monoid)

    desc = "ABCDEFG"[:n]

    def translate(word):
        "monoid to weyl"
        g = identity
        for c in word:
            g = g * G.gen[desc.index(c)]
        return g

    monoid = list(monoid)
    monoid.sort(key = lambda g : (len(g.word), g.word))
    for g in monoid:
        tgt = list(set(g.perm.values()))
        g = translate(g.word)
        print "%6s"%g.word, len(tgt)
    print

    return

    m = G.matrix(desc)

    from bruhat import BruhatMonoid
    monoid = BruhatMonoid(desc, m, bruhat=True, build=True)
    #for w in monoid.words:
    #    print w,
    #print

    def translate(word):
        "translate to function monoid"
        g = identity
        for c in word:
            g = g * gen[desc.index(c)]
        return g

    lookup = {}
    for w0 in monoid.words:
        g = translate(w0)
        w1 = lookup.get(g)
        if w1 is not None:
            print w0, "=", w1
        else:
            lookup[g] = w0
            print w0

    for w0 in monoid.words:
      for w1 in monoid.words:
        w2 = monoid.mul[w0, w1]
        if translate(w0)*translate(w1) == translate(w2):
            pass
        else:
            print "%r*%r = %r" % (w0, w1, w2),
            print " ****************** FAIL"



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

    elif argv.test:
        for n in range(2, 6):
            test(n)

    G = None

    if argv.E_8:
        G = Weyl.build_E8()
        assert G.matrix() == {
            (0, 1):3, (1, 2):3, (2, 3):3, (3, 4):3, (4, 5):3, (4, 6):3, (6, 7):3}

    if argv.E_7:
        G = Weyl.build_E7()
        assert G.matrix() == {
            (0, 6): 3, (1, 2): 3, (2, 3): 3, (3, 4): 3, (3, 5): 3, (5, 6): 3}

    if argv.E_6:
        G = Weyl.build_E6()
        assert G.matrix() == {
            (0, 1): 3, (1, 2): 3, (2, 3): 3, (2, 4): 3, (4, 5): 3}
        #items = mulclose(G.gen, verbose=True)
        #print "|E6|=", len(items)

    if argv.F_4:
        G = Weyl.build_F4()
        assert G.matrix() == {(0, 1): 3, (1, 2): 4, (2, 3): 3}
        #items = mulclose(G.gen, verbose=True)
        #print "|F4|=", len(items) # == 1152

    if argv.G_2:
        G = Weyl.build_G2()
        assert G.matrix() == {(0, 1): 6}
        assert len(mulclose(G.gen)) == 12 # dihedral group D_6

    for arg in argv:
        if len(arg) != 3 or arg[1]!="_":
            continue
        try:
            n = int(arg[2:])
        except:
            continue

        print "building %s"%arg
        if arg.startswith("A"):
            G = Weyl.build_A(n)
        if arg.startswith("B"):
            G = Weyl.build_B(n)
        if arg.startswith("C"):
            G = Weyl.build_C(n)
        if arg.startswith("D"):
            G = Weyl.build_D(n)

    if G is None:
        return

    if argv.monoid:
        test_monoid(G)


if __name__ == "__main__":

    if argv.profile:
        import cProfile as profile
        profile.run("main()")

    else:

        main()






