#!/usr/bin/env python

import sys

from util import factorial, all_subsets
from argv import argv


def mulclose(els, verbose=False, maxsize=None):
    els = set(els)
    changed = True
    while changed:
        if verbose:
            print "mulclose:", len(els)
        changed = False
        _els = list(els)
        for A in _els:
            for B in _els:
                C = A*B 
                if C not in els:
                    els.add(C)
                    if maxsize and len(els)>=maxsize:
                        return list(els)
                    changed = True
    return els

def identity(items):
    return dict((i, i) for i in items)

class Set(object):
    def __init__(self, els):
        self.els = els

    def __str__(self):
        return str(self.els)


class Perm(object):

    """
    A permutation of a list of items.
    """
    def __init__(self, perm, items, word=''):
        #if isinstance(perm, list):
        #    perm = tuple(perm)
        if perm and isinstance(perm, (list, tuple)) and isinstance(perm[0], (int, long)):
            perm = list(items[i] for i in perm)
        if not isinstance(perm, dict):
            perm = dict((perm[i], items[i]) for i in range(len(perm)))
        self.perm = perm # map item -> item
        self.items = set(items)
        assert len(perm) == len(items), (perm, items)
        self.n = len(perm)
        self._str = None
        if len(items)>10:
            items = set(items)
        for key, value in perm.items():
            assert key in items, repr(key)
            assert value in items, repr(value)
        self.word = word

    @classmethod
    def identity(cls, items, *args, **kw):
        n = len(items)
        perm = dict((item, item) for item in items)
        return Perm(perm, items, *args, **kw)

    def restrict(self, items, *args, **kw):
        perm = dict((i, self.perm[i]) for i in items)
        return Perm(perm, items, *args, **kw)

    def fixes(self, items):
        items = set(items)
        for item in items:
            item = self(item)
            if item not in items:
                return False
        return True

    def _X_str__(self):
        #return str(dict((i, self.perm[i]) for i in range(self.n)))
        #return str(dict((i, self.perm[i]) for i in range(self.n)))
        if self._str:
            return self._str
        perm = self.perm
        keys = perm.keys()
        keys.sort()
        items = ["%s:%s"%(key, perm[key]) for key in keys]
        s = "{%s}"%(', '.join(items))
        self._str = s
        return s

    def __str__(self):
        remain = set(self.items)
        s = []
        while remain:
            item = iter(remain).next()
            orbit = [item]
            item1 = self*item
            while item1 != item:
                orbit.append(item1)
                item1 = self*item1
            s.append("(%s)"%(' '.join(str(item) for item in orbit)))
            assert orbit
            n = len(remain)
            for item in orbit:
                remain.remove(item)
            assert len(remain) < n
        #return "Perm(%s)"%(''.join(s))
        return ''.join(s)
    __repr__ = __str__

    def str(self):
        perm = self.perm
        items = self.items
        s = []
        for i, item in enumerate(items):
            j = items.index(perm[item])
            s.append("%d:%d"%(i, j))
        s = "{%s}"%(', '.join(s))
        return s

    def __hash__(self):
        return hash(str(self))

    def __mul__(self, other):
        if not isinstance(other, Perm):
            item = self.perm[other]
            return item
        assert self.items == other.items
        perm = {}
        #for item in self.items:
            #perm[item] = other.perm[self.perm[item]]
        for item in other.items:
            perm[item] = self.perm[other.perm[item]]
        return Perm(perm, self.items, self.word+other.word)
    __call__ = __mul__

    def __pow__(self, n):
        assert int(n)==n
        if n==0:
            return Perm.identity(self.items)
        if n<0:
            self = self.__invert__()
            n = -n
        g = self
        for i in range(n-1):
            g = self*g
        return g

    def __invert__(self):
        perm = {}
        for item in self.items:
            perm[self.perm[item]] = item
        return Perm(perm, self.items)

    def __eq__(self, other):
        assert self.items == other.items
        return self.perm == other.perm

    def __ne__(self, other):
        assert self.items == other.items
        return self.perm != other.perm

    def __getitem__(self, idx):
        return self.perm[idx]

    def fixed(self):
        items = []
        for item in self.items:
            if self.perm[item] == item:
                items.append(item)
        return items


#class Species(object):
#    def __call__(self, group, items):
#        pass

#class PointedSpecies(Species):

class Item(object):
    def __init__(self, item, name=None):
        self.item = item
        if name is None:
            name = str(item)
        self.name = name # a canonical name

    def __hash__(self):
        return hash(self.name)

    def __str__(self):
        return self.name
    def __repr__(self):
        return "I(%s)"%(self.name)
    def __eq__(self, other):
        return self.name == other.name
    def __ne__(self, other):
        return self.name != other.name


class TupleItem(Item):
    def __init__(self, items):
        Item.__init__(self, items)


class SetItem(Item):
    "a hashable unordered set"
    def __init__(self, items):
        items = list(items)
        #for i in range(len(items)):
        #    item = items[i]
        #    if isinstance(item, Item):
        #        pass
        #    elif isinstance(items[i], (int, long, str)):
        #        items[i] = Item(item)
        #    else:
        #        assert 0, repr(item)
        items.sort(key = lambda item : str(item))
        #items = tuple(items)
        Item.__init__(self, items)

    def __iter__(self):
        return iter(self.item)

    def __len__(self):
        return len(self.item)


def disjoint_union(items, _items):
    items = [(0, item) for item in items]
    _items = [(1, item) for item in _items]
    return items + _items
    

def all_functions(source, target):
    m, n = len(source), len(target)
    source = list(source)
    target = list(target)
    assert n**m < 1e8, "too big"
    if m==0:
        yield {}
    elif n==0:
        return # no functions here
    elif m==1:
        for i in range(n):
            yield dict([(source[0], target[i])])
    else:
        for func in all_functions(source[1:], target):
            for i in range(n):
                _func = dict(func)
                _func[source[0]] = target[i]
                yield _func

assert len(list(all_functions('ab', 'abc'))) == 3**2
assert len(list(all_functions('abc', 'a'))) == 1
assert len(list(all_functions('a', 'abc'))) == 3


def __choose(items, k):
    "choose k elements"
    n = len(items)
    assert 0<=k<=n
    if k==0:
        yield [], items, [] # left, right, chosen
    elif k==n:
        yield items, [], [] # left, right, chosen
    elif k==1:
        for i in range(n):
            yield items[:i], items[i+1:], [items[i]]
    else:
        for left, right, chosen in __choose(items, k-1):
            n = len(right)
            for i in range(n):
                yield left + right[:i], right[i+1:], chosen+[right[i]]

def _choose(items, *ks):
    "yield a tuple "
    if len(ks)==0:
        yield items,
    elif len(ks)==1:
        k = ks[0]
        for left, right, chosen in __choose(items, k):
            yield items, chosen
    else:
        k = ks[0]
        for flag in _choose(items, *ks[:-1]):
            chosen = flag[-1]
            for chosen, _chosen in _choose(chosen, ks[-1]):
                yield flag + (_chosen,)


def choose(items, *ks):
    "choose k elements"
    items = list(items)
    _items = []
    #for left, right, chosen in _choose(items, k):
    #    _items.append((SetItem(left+right), SetItem(chosen)))
    for flag in _choose(items, *ks):
        flag = tuple(SetItem(item) for item in flag)
        _items.append(flag)
    return _items

items4 = list('abcd')
assert len(choose(items4, 0)) == 1
assert len(choose(items4, 1)) == 4
assert len(choose(items4, 2)) == 4*3//2
assert len(choose(items4, 3)) == 4
assert len(choose(items4, 4)) == 1


class Action(object):
    """
    A list of Perm's, closed under multiplication.
    """

    def __init__(self, perms, items, check=False):
        perms = list(perms)
        self.perms = perms
        self.items = set(items)
        for perm in perms:
            assert perm.items == self.items, (perm.items, items)

    @classmethod
    def generate(cls, perms, *args, **kw):
        items = perms[0].items
        perms = list(mulclose(perms, *args))
        return cls(perms, items, **kw)

    @property
    def identity(self):
        p = Perm.identity(self.items)
        return p

    @classmethod
    def trivial(cls, items_or_n=1, check=False):
        "the trivial action on items fixes every item"
        if type(items_or_n) in (int, long):
            items = range(items_or_n)
        else:
            items = list(items_or_n)
        perm = Perm.identity(items)
        G = Action([perm], items, check=check)
        return G

    @classmethod
    def symmetric(cls, items_or_n, check=False):
        if type(items_or_n) in (int, long):
            items = range(items_or_n)
        else:
            items = list(items_or_n)
        perms = []
        n = len(items)
        for i in range(n-1):
            perm = dict((item, item) for item in items)
            perm[items[i]] = items[i+1]
            perm[items[i+1]] = items[i]
            perms.append(perm)
        perms = [Perm(perm, items) for perm in perms]
        G = Action.generate(perms, check=check)
        assert len(G) == factorial(n)
        return G

    @classmethod
    def cyclic(cls, items_or_n, check=False):
        if type(items_or_n) in (int, long):
            items = range(items_or_n)
        else:
            items = list(items_or_n)
        perms = []
        n = len(items)
        perms = [dict((items[i], items[(i+k)%n]) for i in range(n))
            for k in range(n)]
        assert len(perms) == n
        perms = [Perm(perm, items) for perm in perms]
        G = Action(perms, items, check=check)
        return G

    @classmethod
    def dihedral(cls, items_or_n, check=False):
        if type(items_or_n) in (int, long):
            items = range(items_or_n)
        else:
            items = list(items_or_n)
        perms = []
        n = len(items)
        perms = [
            dict((items[i], items[(i+1)%n]) for i in range(n)),
            dict((items[i], items[(-i)%n]) for i in range(n))]
        perms = [Perm(perm, items) for perm in perms]
        G = Action.generate(perms, check=check)
        assert len(G) == 2*n
        return G

    def __str__(self):
        return "Action(%s, %s)"%(self.perms, self.items)

    def __len__(self):
        return len(self.perms)

    def __getitem__(self, idx):
        return self.perms[idx]
    
    def __add__(self, other):
        items = disjoint_union(self.items, other.items)
        perms = []
        for perm in self:
            _perm = {}
            for item in self.items:
                _perm[0, item] = 0, perm[item]
            for item in other.items:
                _perm[1, item] = 1, item # identity
            _perm = Perm(_perm, items)
            perms.append(_perm)
        for perm in other:
            _perm = {}
            for item in self.items:
                _perm[0, item] = 0, item # identity
            for item in other.items:
                _perm[1, item] = 1, perm[item]
            _perm = Perm(_perm, items)
            perms.append(_perm)
        perms = list(mulclose(perms))
        return Action(perms, items)

    def __mul__(self, other):
        "direct product of groups"
        items = [(i, j) for i in self.items for j in other.items]
        perms = []
        for g in self:
          for h in other:
            perm = {}
            for i, j in items:
                perm[i, j] = g[i], h[j]
            perm = Perm(perm, items)
            perms.append(perm)
        group = Action(perms, items)
        return group

    def square(self, other=None):
        if other is None:
            other = self
        assert self.isomorphic(other) # bit of a hack...
        items = [(i, j) for i in self.items for j in other.items]
        perms = []
        send_perms = {}
        send_items = dict((i, (i, i)) for i in self.items) # the diagonal
        diagonal = [(i, i) for i in self.items]
        for i in range(len(self)):
            g = self[i]
            h = other[i]
            perm = {}
            for i, j in items:
                perm[i, j] = g[i], h[j]
            perm = Perm(perm, items)
            perms.append(perm)
            #send_perms[g] = perm.restrict(diagonal)
            send_perms[g] = perm
        action = Action(perms, items)
        #print len(send_perms), len(self)
        #print self
        #print action
        #print send_perms
        hom = Hom(self, action, send_perms, send_items)
        hom.check()
        return hom

    def fixed(self):
        fixed = {}
        for el in self.perms:
            fixed[el] = el.fixed()
        return fixed

    def stab(self, item):
        "stabilizer subgroup"
        perms = []
        for g in self.perms:
            if g[item]==item:
                perms.append(g)
        return Action(perms, self.items)

    def orbit(self, item):
        items = set(g*item for g in self.perms)
        return items

    def orbits(self):
        #print "orbits"
        #print self.perms
        #print self.items
        remain = set(self.items)
        orbits = []
        while remain:
            #print "remain:", remain
            item = iter(remain).next()
            orbit = set(g*item for g in self.perms)
            #print "orbit:", orbit
            for item in orbit:
                remain.remove(item)
            orbits.append(orbit)
        return orbits

    def components(self):
        orbits = self.orbits()
        actions = [Action([perm.restrict(orbit) for perm in self.perms], orbit) for orbit in orbits]
        return actions

    def check(group):
        #print "check"
        orbits = group.orbits()
        assert sum(len(orbit) for orbit in orbits) == len(group.items)

        # orbit stabilizer theorem
        n = sum(len(group.stab(item)) for item in group.items)
        assert n == len(group) * len(orbits)

        # Cauchy-Frobenius lemma
        assert n == sum(len(g.fixed()) for g in group)

    def all_subactions(self):
        "All subgroups, acting on the same items. Return as list of homs into self."
        # do the brute-force stupid algorithm
        n = len(self.perms)
        permss = []
        homs = []
        send_items = identity(self.items)
        for idxs in all_subsets(n):
            if not idxs:
                continue
            perms = [self.perms[idx] for idx in idxs]
            perms = set(mulclose(perms))
            if perms in permss:
                continue
            permss.append(perms)
            send_perms = identity(perms)
            action = Action(perms, self.items)
            hom = Hom(action, self, send_perms, send_items)
            homs.append(hom)
        #actions = [Action(perms, self.items) for perms in permss]
        #return actions
        return homs

#    def choice(group, k):
#        "choose k elements"
#    
#        items = group.items
#        _items = choose(items, k)
#        #_group = set()
#        _group = []
#        for g in group:
#            perm = g.perm
#            _perm = {}
#            for left, right in _items:
#                _left = SetItem(tuple(perm[item] for item in left))
#                _right = SetItem(tuple(perm[item] for item in right))
#                _perm[left, right] = _left, _right
#            _g = Perm(_perm, _items)
#            #_group.add(_g)
#            _group.append(_g)
#        return Action(_group, _items)

    def choice(group, *ks):
        "choose k elements"
    
        items = group.items
        _items = choose(items, *ks)
        _group = []
        for g in group:
            perm = g.perm
            _perm = {}
            for flag in _items:
                _flag = tuple(SetItem(tuple(perm[item] for item in __items)) for __items in flag)
                _perm[flag] = _flag
            _g = Perm(_perm, _items)
            _group.append(_g)
        return Action(_group, _items)

    def is_hom(self, items, func):
        #print "is_hom", items, func
        # If we use a set here this could make a "subgroup" of self
        #group = set()
        group = []
        for g in self:
            perm = {} # try to send g to this guy
            for i in self.items:
                gi = g[i]
                f_i = func[i]
                f_gi = func[gi]
                j = perm.get(f_i)
                if j is None:
                    perm[f_i] = f_gi
                elif j != f_gi:
                    return # fail
            #group.add(Perm(perm, items))
            group.append(Perm(perm, items))
        group = Action(group, items)
        return group

    def all_homs(self, items, surjective=True):
        for func in all_functions(self.items, items):
            if surjective:
                if len(set(func.values())) < len(items):
                    continue
            group = self.is_hom(items, func)
            if group is not None:
                group.check()
                yield group, func

    def isomorphic(self, other):
        assert isinstance(other, Action)
        if len(self)!=len(other):
            return False
        n = len(self)
        for i in range(n):
          for j in range(n):
            g = self.perms[i] * self.perms[j]
            g = self.perms.index(g) # slow
            h = other.perms[i] * other.perms[j]
            h = other.perms.index(h) # slow
            if g!=h:
                return False
        return True

    def is_hom_iso(self, other, func):
        for i in range(len(self)):
            g = self[i]
            h = other[i]
            for i in self.items:
                gi = g[i]
                fi = func[i]
                fgi = func[gi]
                hfi = h(fi)
                if hfi != fgi:
                    return # fail
        return True

    def all_homs_iso(self, other, surjective=True):
        assert self.isomorphic(other)
        for func in all_functions(self.items, other.items):
            if surjective:
                if len(set(func.values())) < len(other.items):
                    continue
            if self.is_hom_iso(other, func):
                yield func


class Hom(object):
    """
        A map of Action's is a map of the perms and a 
        map of the inderlying items, such that these two maps commute.
    """
    def __init__(self, G, H, send_perms, send_items, check=False):
        assert isinstance(G, Action)
        assert isinstance(H, Action)
        self.G = G
        self.H = H
        self.send_perms = dict(send_perms) # map G.perms to H.perms
        self.send_items = dict(send_items) # map G.items to H.items
        if check:
            self.check()

    @property
    def src(self):
        return self.G

    @property
    def tgt(self):
        return self.H

    def attrs(self):
        return (self.G, self.H, self.send_perms, self.send_items)

    # Equality on-the-nose:
    def __eq__(self, other):
        return (self.G==other.G and self.H==other.H and self.send_perms==other.send_perms 
            and self.send_items==other.send_items)

    def __ne__(self, other):
        return (self.G!=other.G or self.H!=other.H or self.send_perms!=other.send_perms 
            or self.send_items!=other.send_items)

    @classmethod
    def identity(cls, G, check=False):
        send_perms = dict((g, g) for g in G)
        send_items = dict((item, item) for item in G.items)
        return cls(G, G, send_perms, send_items, check=check)

    def check(self):
        G, H, send_perms, send_items = self.G, self.H, self.send_perms, self.send_items
        assert len(send_perms)==len(G.perms)
        assert len(send_items)==len(G.items)
        for item in G.items:
            assert item in send_items
            assert send_items[item] in H.items
        for perm in G.perms:
            assert perm in send_perms
            assert send_perms[perm] in H.perms

        # G and H should really be the same (isomorphic) as abstract groups. 
        # Here we merely check that we have a homomorphism of groups.
        for g1 in G.perms:
          h1 = send_perms[g1]
          for g2 in G.perms:
            h2 = send_perms[g2]
            assert send_perms[g1*g2] == h1*h2
        #assert len(G) == len(H)

        for g in G.perms:
            # g is a perm on G.items
            h = send_perms[g]
            # h is a perm on H.items, check it agrees with send_items of g.
            perm = dict((send_items[item], send_items[g[item]]) for item in G.items)
            items = set(send_items[item] for item in G.items)
            h1 = Perm(perm, items)
            assert h.restrict(items) == h1
            #for item in G.items:

    def __mul__(other, self):
        "composition of Hom's"
        assert isinstance(other, Hom)
        assert isinstance(self, Hom)
        assert self.tgt is other.src
        send_perms = {}
        for perm in self.G.perms:
            send_perms[perm] = other.send_perms[self.send_perms[perm]]
        send_items = {}
        for item in self.G.items:
            send_items[item] = other.send_items[self.send_items[item]]
        hom = Hom(self.G, other.H, send_perms, send_items)
        return hom


def test_hom():

    n = argv.get("n", 3)
    items = range(n)

    G = Action.trivial(items, check=True)
    assert len(G.components()) == len(items)

    G = Action.cyclic(items, check=True)
    assert len(G.components()) == 1

    hom = G.square()
    hom.check()
    assert len(hom.src) == len(G)
    for G1 in hom.src.components():
        assert len(G1)==n

    F = Action.dihedral(items, check=True)
    assert len(F.components()) == 1
    homFF = F.square()
    homFF.check()
    for G1 in homFF.src.components():
        print len(G1),
    print

    return

    homs = F.all_subactions()
    for A in homs:
      for B in homs:
        print len(hom.G)
        A1 = homFF * A
        B1 = homFF * B

    print
#    for T1 in Ts:
#      for T2 in Ts:
#        print "%s * %s = %s" % (len(T1), len(T2), [len(T3) for T3 in (T1*T2).components()])

    G = Action.symmetric(items, check=True)
    assert len(G.components()) == 1

    f = Hom.identity(G, check=True)
    ff = f*f
    ff.check()
    assert ff == f

    print "OK"



def is_hom(hom):
    for g in hom.keys():
      for h in hom.keys():
        k = g*h
        if k in hom:
            if hom[k] != hom[g]*hom[h]:
                return False
    return True


def close_hom(hom):
    hom = dict(hom)
    for g in hom.keys():
      for h in hom.keys():
        k = g*h
        if k in hom:
            if hom[k] != hom[g]*hom[h]:
                return None
        else:
            hom[k] = hom[g]*hom[h]
    return hom


def find_all_homs(G, H, hom=None, remain=None):
    """
        Iterate through all homomorphisms from Action G to Action H.
        Yield each homomorphism as a dict : Perm -> Perm.
    """
    if hom is None:
        GI = G.identity
        hom = {GI : H.identity}
        #remain = [p for p in G if p != GI]

    for g in G:
        if g in hom:
            continue

        for h in H:
            hom[g] = h

            hom1 = close_hom(hom)
            if hom1 is None:
                pass

            elif len(hom1) == len(G):
                yield hom1

            else:
                for _hom in find_all_homs(G, H, hom1):
                    yield _hom


#class Action(object):
#    """
#        A group action is a homomorphism from a Action to another Action.
#    """
#    def __init__(self, group, hom):
#        self.group = group
#        self.hom = dict(hom) # dict : Perm -> Perm



def test():

    items = list('abc')

    e = Perm.identity(items)

    g = Perm((1, 0, 2), items)
    assert g*g == e
    assert ~g == g

    h = Perm((1, 2, 0), items)
    assert h*h*h == e
    assert h*h != e
    assert h*h != g
    assert not (h*h == e)

    S3 = Action.generate([g, h])
    S3.check()
    assert len(S3) == 6

    assert len(list(S3.all_homs_iso(S3))) == 1

    #for g in S3:
    #    print g, g.fixed()

    items4 = list('abcd')
    g = Perm((1, 2, 3, 0), items4)
    h = Perm((1, 0, 2, 3), items4)
    S4 = Action.generate([g, h])
    assert len(S4)==24

    S4_22 = S4.choice(2)
    S4_22.check()

    # Pauli group
    X = Perm((3, 1, 2, 0), items4)
    Z = Perm((2, 3, 0, 1), items4)
    I = Perm((0, 1, 2, 3), items4)
    w = Perm((3, 2, 1, 0), items4)
    assert X*X==I
    assert Z*Z==I
    assert Z*X != X*Z
    assert Z*X*Z*X == w
    P1 = Action.generate([X, Z])
    assert len(P1)==8

    if 0:
        import numpy
        group = P1.square(P1)
        print "orbits:", len(group.orbits())
        for orbit in group.orbits():
            print orbit
            A = numpy.zeros((4, 4))
            for (i, j) in orbit:
                i, j = items4.index(i), items4.index(j)
                A[i, j] = 1
            print A

    # Pauli group
    items = "+00 -00 +01 -01 +10 -10 +11 -11".split()
    #          0   1   2   3   4   5   6   7
    II = Perm((0, 1, 2, 3, 4, 5, 6, 7), items)
    XI = Perm((4, 5, 6, 7, 0, 1, 2, 3), items)
    IX = Perm((2, 3, 0, 1, 6, 7, 4, 5), items)
    ZI = Perm((0, 1, 2, 3, 5, 4, 7, 6), items)
    IZ = Perm((0, 1, 3, 2, 4, 5, 7, 6), items)
    w2 = Perm((1, 0, 3, 2, 5, 4, 7, 6), items)
    assert XI*XI==II
    assert ZI*ZI==II
    assert IX*IX==II
    assert IZ*IZ==II
    assert ZI*XI != XI*ZI
    assert ZI*XI == w2*XI*ZI
    assert ZI*XI*ZI*XI == w2
    assert IZ*IX != IX*IZ
    assert IZ*IX == w2*IX*IZ
    assert IZ*IX*IZ*IX == w2

    assert ZI*IX == IX*ZI

    P2 = Action.generate([XI, ZI, IX, IZ])
    assert len(P2)==32

    #homs = [f for f in find_all_homs(P1, P2, {I:II, w:w2})] # 152 solutions
    #assert len(homs)==152

    hom = {II:I} # 963 solutions
    hom = {II:I, w2:I} # 963 solutions
    hom = {II:I, w2:w, XI:X} # 0 solutions
    hom = {II:I, XI:X} # 64 solutions
    hom = {II:I, XI:X, IX:X} # 16 solutions
    hom = {II:I, w2:w} # 0 solutions
    count = 0
    for f in find_all_homs(P2, P1, hom):
        count += 1
    assert count==0, count

    if 1:
        #G1 = S4.choice(2)
        G1 = S4.choice(3, 2, 1)
        G2 = S4.choice(2)
        group = G1.square(G2)
        print "SxS:", len(group.items)
        print "orbits:", len(group.orbits())
        print

        group = S4.choice(2).square()
        print "##\n##"
        print "SxS:", len(group.items)
        print "orbits:", len(group.orbits())
    
        group = S4.choice(2, 1).square()
        print "##\n#\n#"
        print len(group)
        print "SxS:", len(group.items)
        print "orbits:", len(group.orbits())
    
        group = S4.choice(3, 2, 1).square()
        print "#\n#\n#\n#"
        print len(group)
        print "SxS:", len(group.items)
        print "orbits:", len(group.orbits())

    S4_211 = S4.choice(2, 1)
    assert len(S4_211.items)==12
    S4_211.check()

    assert S4.isomorphic(S4_22)

    assert len(list(S4_22.all_homs_iso(S4_22))) == 2
    assert len(list(S4.all_homs_iso(S4))) == 1

    #print len(list(S4.all_homs(list('ab'))))

    Z4 = Action.generate([Perm((1, 2, 3, 0), items4)])
    Z4.check()
    assert len(Z4)==4

    #print len(list(Z4.all_homs(list('abcd'))))

    assert len(list(Z4.all_homs(list('ab')))) == 2
    assert len(list(Z4.all_homs_iso(Z4))) == 4

    Z4_Z4 = Z4 + Z4
    assert len(Z4_Z4)==16
    Z4_Z4.check()
    assert len(Z4_Z4.orbits()) == 2

    group = Z4 * Z4
    group.check()
    #print len(group), len(group.orbits())

    #for g in Z4:
    #    print g, g.fixed()

    Z22 = Action.generate([Perm((1, 0, 3, 2), items4), Perm((2, 3, 0, 1), items4)])
    assert len(Z22)==4

    assert len(list(Z22.all_homs_iso(Z22))) == 4
    assert len(list(Z22.all_homs(list('ab')))) == 6
    #for group, func in Z22.all_homs(list('ab')):
    #    print func

    if 0:
        group = Z22.choice(2)
        print "Z22.choice(2): orbits", [len(orbit) for orbit in group.orbits()]
        group = Z22.choice(2, 1)
        print "Z22.choice(2, 1): orbits", [len(orbit) for orbit in group.orbits()]
        group = Z22.choice(3, 2, 1)
        print "Z22.choice(3, 2, 1): orbits", [len(orbit) for orbit in group.orbits()]
    
        group = Z4.choice(2)
        print "Z4.choice(2): orbits", [len(orbit) for orbit in group.orbits()]
        group = Z4.choice(2, 1)
        print "Z4.choice(2, 1): orbits", [len(orbit) for orbit in group.orbits()]
        group = Z4.choice(3, 2, 1)
        print "Z4.choice(3, 2, 1): orbits", [len(orbit) for orbit in group.orbits()]

#    print "fixed:",
#    for g in group:
#        #print g,
#        print len(g.fixed()),
#        #for item in g.fixed():
#        #    print "\t", item
#    print
#
#    print "orbits:"
##    for item in group.items:
##        print item, len(group.orbit(item)), ":", #repr(group.orbit(item)),
##        #for item in group.orbit(item):
##        #    print item,
##        print
#
#    print len(group.orbits())

    print "OK"



if __name__ == "__main__":

    test_hom()

    if argv.test:
        test()


