#!/usr/bin/env python

import sys
from heapq import heappush, heappop, heapify

import numpy
#import scipy.sparse.linalg as la

from code import texstr

EPSILON = 1e-8
scalar = numpy.float64

def write(s):
    sys.stdout.write(str(s)+' ')
    sys.stdout.flush()


gcolor_gauge = """
1111...........
11..11.........
1.1.1.1........
..11..11.......
.1.1.1.1.......
....1111.......
11......11.....
1.1.....1.1....
........1111...
..11......11...
.1.1.....1.1...
1...1...1...1..
........11..11.
.1...1...1...1.
....11......11.
........1.1.1.1
..1...1...1...1
....1.1.....1.1
"""

gcolor_stab = """
11111111.......
1111....1111...
11..11..11..11.
1.1.1.1.1.1.1.1
"""

cube_ham = """
6111....
14..11..
1.4.1.1.
1..4.11.
.11.2..1
.1.1.2.1
..11..21
....1110
"""


def parse(s):
    s = s.replace('.', '0') 
    lines = s.split()
    lines = [l.strip() for l in lines if l.strip()]
    rows = [list(int(c) for c in l) for l in lines]
    if rows:
        n = len(rows[0])
        for row in rows:
            assert len(row)==n, "rows have varying lengths"
    a = numpy.array(rows, dtype=numpy.int32)
    return a


class Point(object):
    def __init__(self, desc, idx, nbd=None, colour=""):
        self.desc = desc
        self.colour = colour
        self.idx = idx
        if nbd is None:
            nbd = []
        self.nbd = nbd

    def get_desc(self, depth=0, source=None):
        assert self.nbd
        assert depth>=0
        desc = self.desc+str(self.colour)
        if depth==0:
            return desc
        if source is None:
            source = []
        else:
            assert self not in source
        descs = [a.get_desc(depth-1, source+[self]) for a in self.nbd if a not in source]
        descs.sort()
        desc = "%s[%s]"%(desc, ' '.join(descs))
        return desc

    def __str__(self):
        return "Point(%s: %s)"%(self.desc, descs)


class Bag(object):
    def __init__(self, m, n, points):
        self.m = m
        self.n = n
        self.points = points

    def __len__(self):
        return len(self.points)

    def __getitem__(self, idx):
        return self.points[idx]

    @classmethod
    def build(cls, Gx):
    
        m, n = Gx.shape
        points = []
        for i in range(m):
            g = Gx[i]
            assert g.sum()==4
            weights = []
            for j in numpy.where(g)[0]:
                weights.append(Gx[:, j].sum())
            weights.sort()
            desc = ''.join(str(w) for w in weights)
            a = Point(desc, i)
            points.append(a)
        #print [a.desc for a in points]
        
        for i in range(m):
            g = Gx[i]
            a = points[i]
            for j in numpy.where(g)[0]:
                for i1 in numpy.where(Gx[:, j])[0]:
                    if i1 != i:
                        a.nbd.append(points[i1])

        return cls(m, n, points)

    def map(self, fn):
        points = [None]*len(self)
        for p in self.points:
            p = Point(p.desc, fn[p.idx])
            points[p.idx] = p
        for p in self.points:
            for p1 in p.nbd:
                points[fn[p.idx]].nbd.append(points[fn[p1.idx]]) # whoops.. tricky
        
        return self.__class__(self.m, self.n, points)

    def get_desc(self, depth=0):
        return [v.get_desc(depth) for v in self.points]

    def get_orbits(self, depth=0):
        orbits = {}
        for p in self.points:
            key = p.get_desc(depth)
            orbit = orbits.setdefault(key, [])
            orbit.append(p)
        return orbits


class Tanner(Bag):
    @classmethod
    def build(cls, Gx):
        # This is the Tanner graph
        m, n = Gx.shape
        checks = [Point('c', i) for i in range(m)]
        bits = [Point('b', i+m) for i in range(n)]
        for i in range(m):
            for j in range(n):
                if Gx[i, j]==0:
                    continue
                checks[i].nbd.append(bits[j])
                bits[j].nbd.append(checks[i])
        return cls(m, n, checks+bits)

    def shortstr(self):
        m, n = self.m, self.n
        rows = []
        for i in range(m): # checks
            row = ['.']*n
            p = self.points[i]
            for p1 in p.nbd:
                row[p1.idx-m] = '1'
            row = ''.join(row)
            rows.append(row)
        return '\n'.join(rows)


def from_ham(H):
    n = len(H)
    points = []
    for i in range(n):
        p = Point('(%s)'%H[i, i], i)
        points.append(p)
    for i in range(n):
      for j in range(n):
        if i==j:
            continue
        if H[i, j]:
            points[i].nbd.append(points[j])
    bag = Bag(n, n, points)
    return bag


def get_perm(m, n, fn):

    U = numpy.zeros((m, m), dtype=int)
    for i in range(m):
        j = fn[i]
        U[i, j] = 1

    V = numpy.zeros((n, n), dtype=int)
    for i in range(n):
        j = fn[i+m]-m
        V[j, i] = 1

    return U, V


def search_recursive(bag0, bag1, fn=None, depth=2):

    assert depth>0

    if fn is None:
        fn = {}
        assert len(bag0)==len(bag1)
        assert bag0 is not bag1

    orbits0 = bag0.get_orbits(depth)
    orbits1 = bag1.get_orbits(depth)

    if len(orbits0) != len(orbits1):
        return

    keys0 = orbits0.keys()
    keys1 = orbits1.keys()
    keys0.sort()
    keys1.sort()
    if keys0 != keys1:
        return

    idx = len(fn)

    # choose any uncoloured bag0 point
    p = bag0.points[idx]
    assert p.colour == ''

    key = p.get_desc(depth)
    orbit = orbits1[key]

    p.colour = str(idx)

    # go through each candidate in bag1
    for p1 in orbit:
        assert p1.colour == ''
    
        p1.colour = str(idx)
    
        assert fn.get(idx) is None
        fn[idx] = p1.idx
    
        if len(fn) == len(bag0):
            yield dict(fn)

        else:

            for _fn in search_recursive(bag0, bag1, fn, depth):
                yield _fn

        del fn[idx]
        assert len(fn) == idx

        p1.colour = ''

    p.colour = ''



"""
class Search(object):

    def __init__(self, bag0, bag1, depth=2):
        assert depth>0
        assert len(bag0)==len(bag1)
        assert bag0 is not bag1
        self.bags = bag0, bag1
        self.fn = {} # current state
        self.stack = []

    def search(self):
        self.done = False
        while not done:

    def push(self, state):
        self.stack.append(state)

    def pop(self):
        if len

    def attempt(self):

        fn = self.fn
        bag0, bag1 = self.bags
        assert len(fn) == len(stack)
"""


class Backtrack(Exception):
    pass


class State(object):
    def __init__(self, bag0, bag1, idx0, depth):
        orbits0 = bag0.get_orbits(depth) # map: desc -> list of points
        orbits1 = bag1.get_orbits(depth) # map: desc -> list of points
    
        if len(orbits0) != len(orbits1):
            raise Backtrack()
    
        keys0 = orbits0.keys()
        keys1 = orbits1.keys()
        keys0.sort()
        keys1.sort()
        if keys0 != keys1:
            raise Backtrack()

        # choose any uncoloured bag0 point
        p0 = bag0.points[idx0]
        assert p0.colour == ''
    
        key0 = p0.get_desc(depth)
        self.orbit1 = orbits1[key0]
        assert self.orbit1 # otherwise: wtf?
        self.idx0 = idx0 # source index
        self.idx1 = 0 # search target index
        self.p0 = p0
        self.p1 = None
    
    def do(self, fn):
        # make assignment: idx0 -> idx1
        p0 = self.p0
        assert p0.colour == ''
        p0.colour = str(self.idx0)

        p1 = self.orbit1[self.idx1]
        assert p1.colour == ''
        p1.colour = str(self.idx0)
    
        assert fn.get(self.idx0) is None
        fn[self.idx0] = p1.idx
        assert self.p1 is None
        self.p1 = p1

    def undo(self, fn):
        # undo assignment
        del fn[self.idx0]
        assert self.p1 is not None
        p0 = self.p0
        p1 = self.p1
        assert p1.colour==str(self.idx0)
        assert p0.colour==str(self.idx0)
        p0.colour = ''
        p1.colour = ''
        self.p1 = None

    def next(self):
        assert self.p1 is None
        self.idx1 += 1
        if self.idx1 >= len(self.orbit1):
            raise Backtrack()


def search(bag0, bag1, depth):

    fn = {}
    idx = 0
    state = State(bag0, bag1, idx, depth)
    stack = [state]

    while stack:

        #print "stack:", len(stack)

        state = stack[-1]
        state.do(fn)

        #print fn

        if len(fn) == len(bag0):
            #print "FOUND"
            yield fn

        else:
            # try to add another state
            try:
                idx = len(fn)
                _state = State(bag0, bag1, idx, depth)
                stack.append(_state)
                continue
    
            except Backtrack:
                # the above do() doesn't work
                pass

        # next
        while stack:
            state = stack[-1]
            #print "UNDO"
            state.undo(fn)
            try:
                #print "NEXT"
                state.next()
                break # ok, finished backtracking
            except Backtrack:
                #print "POP"
                stack.pop() # discard this guy



def all_autos(Gx):
    #Gx = parse(gcolor_gauge)
    m, n = Gx.shape

    bag0 = Tanner.build(Gx)
    bag1 = Tanner.build(Gx)

    for fn in search(bag0, bag1):
        U, V = get_perm(m, n, fn)
        yield U, V


def test():

    # Find rotation symmetry of the code. It's S_4 with order 24.

    Gx = parse(gcolor_gauge)
    m, n = Gx.shape

    bag0 = Tanner.build(Gx)
    bag1 = Tanner.build(Gx)

    #global search
    #search = search_recursive

    count = 0
    for fn in search(bag0, bag1, depth=2):
        #print "iso", fn
        bag = bag0.map(fn)
        #print bag.shortstr()
        U, V = get_perm(m, n, fn)
        Gx1 = numpy.dot(U, numpy.dot(Gx, V))
        assert numpy.abs(Gx-Gx1).sum()==0
        count += 1
    #print "count:", count
    assert count == 24

    # S_3 symmetry of cubical hamiltonian
    depth = 1
    H = parse(cube_ham)
    bag0 = from_ham(H)
    bag1 = from_ham(H)
    count = 0 
    for fn in search(bag0, bag1, depth=depth):
#        #write('.')
#        print [fn[i] for i in range(len(fn))]
#        for i in range(len(fn)):
#            assert bag0[i].get_desc(depth) == bag1[fn[i]].get_desc(depth)
#            print '\t', bag0[i].get_desc(depth),
#            print len(bag0[i].nbd)
#            print '\t', bag1[fn[i]].get_desc(depth)
#        print
        count += 1
    assert count == 6



from argv import Argv 
argv = Argv()


if __name__ == "__main__":

    test()

    print "OK"



