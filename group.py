#!/usr/bin/env python

"""
failed experiment..

"""

import sys, os


class Element(object):
    def __init__(self, parent, name, data=None):
        self.parent = parent
        self.name = name
        self.data = data
        self.parent.add(self)

    #def __hash__(self):
    #    return hash(self.name)

    def __str__(self):
        return self.name

    def __repr__(self):
        return "%s:%s" % (self.name, self.parent.name)

    def __mul__(self, other):
        parent = self.parent
        result = parent.mul(self, other)
        return result

    def __invert__(self):
        parent = self.parent
        result = parent.invert(self)
        return result


class Category(object):

    def __init__(self, name='cat'):
        self.name = name
        self.els = set([])

    def add(self, el):
        assert el not in self.els
        self.els.add(el)

    def __len__(self):
        return len(self.els)

    #def __getitem__(self, idx):
    #    return self.els[idx]
    def __iter__(self):
        return iter(self.els)

    def __str__(self):
        return self.name

    def mul(self, g, h):
        pass

    def invert(self, g):
        pass


class AbstractGroup(Category):

    def __init__(self, name='group'):
        Category.__init__(self, name)
        self.els = set([])
        self.ident = None
        self.mtable = {} # map (g, h) -> gh
        self.itable = {} # map g -> ~g

    def check(self):
        els = self.els
        ident = self.ident
        for g in els:
            assert ~g in els
            for h in els:
                assert g*h in els
            assert ident * g is g
            assert g * ident is g
            assert (g*~g) is ident
            assert (~g*g) is ident
        
        for g in els:
         for h in els:
          for k in els:
            assert (g*h)*k is g*(h*k)
    
    def mul(self, g, h):
        gh = self.mtable[g, h]
        return gh

    def invert(self, g):
        ig = self.itable[g]
        return ig

    def __mul__(self, other):
        return ProductGroup(self, other)

    def dump(self):
        els = list(self.els)
        els.sort(key = lambda el : str(el))
        w = max([len(str(el)) for el in els])
        print " "*w + " | ",
        for h in els:
            print str(h).rjust(w),
        print
        print "-"*w + "-+-" + "-"*((w+1)*len(els))
        for g in els:
            print str(g).rjust(w) + " | ",
            for h in els:
                gh = g*h
                print str(gh).rjust(w),
            print


class ProductGroup(AbstractGroup):
    def __init__(self, G, H):
        name = "(%s*%s)"%(G, H)
        AbstractGroup.__init__(self, name)
        els = []
        pairs = {}
        for g in G:
          for h in H:
            el = Element(self, "(%s, %s)"%(g.name, h.name), data=(g, h))
            pairs[g, h] = el
        self.ident = pairs[G.ident, H.ident]
        itable = self.itable
        mtable = self.mtable
        for g in G:
          for h in H:
            gh = pairs[g, h]
            igh = pairs[~g, ~h]
            itable[gh] = igh
            for g1 in G:
              for h1 in H:
                gh1 = pairs[g1, h1]
                ghgh1 = pairs[g*g1, h*h1]
                mtable[gh, gh1] = ghgh1


class CyclicGroup(AbstractGroup):
    def __init__(self, size):
        name = "Z_%d"%(size)
        AbstractGroup.__init__(self, name)
        els = []
        for i in range(size):
            el = Element(self, str(i))
            els.append(el)
        self.ident = els[0]
        for i, g in enumerate(els):
            self.itable[g] = els[(size-i)%size]
            for j, h in enumerate(els):
                self.mtable[g, h] = els[(i+j)%size]
        self.check()
        

def test():

    z2 = CyclicGroup(2)
    z5 = CyclicGroup(5)
    assert len(z5) == 5

    z2z5 = z2*z5
    assert len(z2z5) == 10

    z2z5.dump()

    print "OK"



if __name__ == "__main__":
    test()






