#!/usr/bin/env python

import sys


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


class Set(object):
    def __init__(self, els):
        self.els = els

    def __str__(self):
        return str(self.els)


class Perm(object):

    """
    A permutation of a list of items.
    """
    def __init__(self, perm, items):
        #if isinstance(perm, list):
        #    perm = tuple(perm)
        if perm and isinstance(perm, (list, tuple)) and isinstance(perm[0], (int, long)):
            perm = list(items[i] for i in perm)
        if not isinstance(perm, dict):
            perm = dict((perm[i], items[i]) for i in range(len(perm)))
        self.perm = perm # map item -> item
        self.items = items
        assert len(perm) == len(items), (perm, items)
        self.n = len(perm)
        self._str = None
        if len(items)>10:
            items = set(items)
        for key, value in perm.items():
            assert key in items, repr(key)
            assert value in items, repr(value)

    @classmethod
    def identity(cls, items):
        n = len(items)
        perm = dict((item, item) for item in items)
        return Perm(perm, items)

    def __str__(self):
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

    def __hash__(self):
        return hash(str(self))

    def __mul__(self, other):
        if not isinstance(other, Perm):
            item = self.perm[other]
            return item
        assert self.items is other.items
        perm = {}
        for item in self.items:
            perm[item] = other.perm[self.perm[item]]
        return Perm(perm, self.items)
    __call__ = __mul__

    def __invert__(self):
        perm = {}
        for item in self.items:
            perm[self.perm[item]] = item
        return Perm(perm, self.items)

    def __eq__(self, other):
        assert self.items is other.items
        return self.perm == other.perm

    def __ne__(self, other):
        assert self.items is other.items
        return self.perm != other.perm

    def __getitem__(self, idx):
        return self.perm[idx]

    def fixed(self):
        items = []
        for item in self.items:
            if self.perm[item] == item:
                items.append(item)
        return items


class Species(object):
    def __call__(self, group, items):
        pass


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


class Group(object):
    """
    This is really a group *action*, since each of our
    elements is presented as an action on a set of "items".
    """

    def __init__(self, els, items):
        els = list(els)
        self.els = els
        for el in els:
            assert el.items is items
        self.items = items

    def __len__(self):
        return len(self.els)

    def __getitem__(self, idx):
        return self.els[idx]
    
    @classmethod
    def generate(cls, els, *args):
        items = els[0].items
        els = list(mulclose(els, *args))
        return cls(els, items)

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
        return Group(perms, items)

    def __mul__(self, other):
        items = [(i, j) for i in self.items for j in other.items]
        perms = []
        for g in self:
          for h in other:
            perm = {}
            for i, j in items:
                perm[i, j] = g[i], h[j]
            perm = Perm(perm, items)
            perms.append(perm)
        group = Group(perms, items)
        return group

    def fixed(self):
        fixed = {}
        for el in self.els:
            fixed[el] = el.fixed()
        return fixed

    def stab(self, item):
        "stabilizer subgroup"
        els = []
        for g in self.els:
            if g[item]==item:
                els.append(g)
        return Group(els, self.items)

    def orbit(self, item):
        items = set()
        changed = True
        while changed:
            changed = False
            for g in self.els:
                _item = g*item
                if _item not in items:
                    changed = True
                    items.add(_item)
        return SetItem(items)

    def orbits(self):
        items = set(self.orbit(item) for item in self.items)
        return list(items)

    def check(group):
        orbits = group.orbits()
        assert sum(len(orbit) for orbit in orbits) == len(group.items)

        # orbit stabilizer theorem
        n = sum(len(group.stab(item)) for item in group.items)
        assert n == len(group) * len(orbits)

        # Cauchy-Frobenius lemma
        assert n == sum(len(g.fixed()) for g in group)

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
#        return Group(_group, _items)

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
        return Group(_group, _items)

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
        group = Group(group, items)
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
        assert isinstance(other, Group)
        if len(self)!=len(other):
            return False
        n = len(self)
        for i in range(n):
          for j in range(n):
            g = self.els[i] * self.els[j]
            g = self.els.index(g) # slow
            h = other.els[i] * other.els[j]
            h = other.els.index(h) # slow
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

    S3 = Group.generate([g, h])
    S3.check()
    assert len(S3) == 6

    assert len(list(S3.all_homs_iso(S3))) == 1

    #for g in S3:
    #    print g, g.fixed()

    items4 = list('abcd')
    g = Perm((1, 2, 3, 0), items4)
    h = Perm((1, 0, 2, 3), items4)
    S4 = Group.generate([g, h])
    assert len(S4)==24

    S4_22 = S4.choice(2)
    S4_22.check()

    S4_211 = S4.choice(2, 1)
    assert len(S4_211.items)==12
    S4_211.check()

    assert S4.isomorphic(S4_22)

    assert len(list(S4_22.all_homs_iso(S4_22))) == 2
    assert len(list(S4.all_homs_iso(S4))) == 1

    #print len(list(S4.all_homs(list('ab'))))

    Z4 = Group.generate([Perm((1, 2, 3, 0), items4)])
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

    Z22 = Group.generate([Perm((1, 0, 3, 2), items4), Perm((2, 3, 0, 1), items4)])
    assert len(Z22)==4

    assert len(list(Z22.all_homs_iso(Z22))) == 4
    assert len(list(Z22.all_homs(list('ab')))) == 6
    #for group, func in Z22.all_homs(list('ab')):
    #    print func

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



if __name__ == "__main__":

    test()
        
    print "OK"


