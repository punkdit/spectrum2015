#!/usr/bin/env python

import os, sys

import numpy
import networkx as nx
from matplotlib import pyplot

from solve import zeros2, enum2, row_reduce, span, shortstr, shortstrx, solve, rank, find_kernel
import isomorph
from isomorph import Bag, Point, write


from argv import Argv
argv = Argv()


def all_subsets(n):

    if n==0:
        yield []
        return

    if n==1:
        yield []
        yield [0]
        return

    for subset in all_subsets(n-1):
        yield subset
        yield subset + [n-1] # sorted !!

assert len(list(all_subsets(5))) == 2**5

def factorial(n):
    r = 1
    for i in range(1, n+1):
        r *= i
    return r

assert factorial(0) == 1
assert factorial(1) == 1
assert factorial(2) == 2
assert factorial(3) == 2*3
assert factorial(4) == 2*3*4


def choose(items, n):
    if n > len(items):
        return
    if n == 0:
        yield ()
        return
    if n == 1:
        for item in items:
            yield (item,)
        return
    for i, item in enumerate(items):
        for rest in choose(items[i+1:], n-1):
            yield (item,)+rest

assert len(list(choose(range(4), 1))) == 4
assert len(list(choose(range(4), 2))) == 6
assert len(list(choose(range(4), 3))) == 4


#
# http://mathworld.wolfram.com/q-BinomialCoefficient.html
#
def qbinomial(n, m, q=2):
    "n choose_q m"
    assert n>=m>=0, (n, m)

    if n==m or m==0:
        return 1

    if m==1:
        top = 1-q**n
        bot = 1-q
        assert top%bot == 0
        return top//bot

    return q**m * qbinomial(n-1, m) + qbinomial(n-1, m-1)

assert qbinomial(4, 2) == 35


class Geometry(object):

    def __init__(self, incidence, tpmap):
        """
            incidence: list of item pairs (i, j)
            tpmap: dict mapping each item to it's type
        """
        items = set() # all items
        nbd = {} # map item -> list of incident items
        tplookup = {} # map type -> list of items of that type
        for i, j in incidence:
            items.add(i)
            items.add(j)
        for i in items:
            nbd[i] = []
        for i, j in incidence:
            if j not in nbd[i]:
                nbd[i].append(j)
            if i not in nbd[j]:
                nbd[j].append(i)
        for jtems in nbd.values():
            assert len(set(jtems))==len(jtems), "nbd:%s"%nbd # uniq

        for i in items:
          for j in items:
            try:
                i==j
            except ValueError:
                print "cant compare %r and %r" % (i, j)
                raise
        for i in items:
            if i not in nbd[i]:
                nbd[i].append(i)
            for j in nbd[i]:
                if i==j:
                    continue

                assert tpmap[i] != tpmap[j]
                if i not in nbd[j]:
                    nbd[j].append(i)

        for t in tpmap.values():
            tplookup[t] = []
        for item, t in tpmap.items():
            tplookup[t].append(item)
        for jtems in nbd.values():
            assert len(set(jtems))==len(jtems), "nbd:%s"%nbd # uniq

        incidence = [] # rebuild this
        for item, jtems in nbd.items():
            for jtem in jtems:
                incidence.append((item, jtem))
        incidence.sort()

        self.incidence = set(incidence) # incidence relation: list of pairs (i, j)
        self.tpmap = dict(tpmap) # map item -> type
        self.tplookup = tplookup # map type -> list of items
        self.types = list(set(tpmap.values()))
        self.rank = len(self.types)
        self.items = items
        self.nbd = nbd
        self.check_geometry()

    def __str__(self):
        incidence = list(self.incidence)
        incidence.sort()
        ss = []
        for i, j in incidence:
            if i<j:
                ss.append("%s--%s" % (i, j))
        return "Geometry([%s], %s)"%(', '.join(ss), self.tpmap)
    __repr__ = __str__

    def __contains__(self, pair):
        return pair in self.incidence

    @classmethod
    def polygon(cls, n):
        verts = ['p%d'%i for i in range(n)]
        edges = ['l%d'%i for i in range(n)]
        incidence = []
        tpmap = {}
        for i in range(n):
            incidence.append((verts[i], edges[i]))
            incidence.append((verts[i], edges[(i+1)%n]))
            tpmap[verts[i]] = 'p'
            tpmap[edges[i]] = 'l'
        #print incidence
        #print tpmap
        return cls(incidence, tpmap)

    @classmethod
    def simplex(cls, dim):
        tps = range(dim+1)
        items = [tuple(idxs) for idxs in all_subsets(dim+1) if idxs] # skip empty set
        incidence = []
        tpmap = {}
        for i in items:
            tpmap[i] = len(i)-1
            for j in items:
                if set(i).intersection(set(j)) == set(i):
                    incidence.append((i, j))
        #print "simplex: %d"%dim
        #print "incidence:", incidence
        #print "tpmap:", tpmap
        return cls(incidence, tpmap)

    @classmethod
    def cube(cls):
        "the verts, edges and faces of a cube"
        verts = []
        edges = []
        faces = []
        incidence = []
        for i in [0, 1]:
         for j in [0, 1]:
          for k in [0, 1]:
            verts.append((i, j, k))
        for edge in choose(verts, 2):
            i, j = edge
            if sum([abs(i[ii]-j[ii]) for ii in range(3)])==1:
                e = (i, j)
                edges.append(e)
                incidence.append((i, e))
                incidence.append((j, e))
        assert len(edges)==12
        assert len(incidence)==24
        for idx in range(3):
          for face in choose(verts, 4):
            r = face[0][idx]
            for v in face:
                if v[idx] != r:
                    break
            else:
                faces.append(face)
                for v in face:
                    incidence.append((v, face))
                for edge in choose(face, 2):
                    if edge in edges:
                        incidence.append((edge, face))
        assert len(faces)==6
        #for (i,j) in incidence[24:]:
        #    print i, '---', j
        assert len(incidence)==24 + 6*4 + 6*4
        tpmap = {}
        for v in verts:
            tpmap[v] = 'v'
        for e in edges:
            tpmap[e] = 'e'
        for f in faces:
            tpmap[f] = 'f'
        g = Geometry(incidence, tpmap)
        return g

    def is_flag(self, flag):
        if not len(set(flag))==len(flag):
            print "X"
            return False
        for i in flag:
            if not i in self.items:
                print "I"
                return False
            for j in flag:
                if not (i, j) in self.incidence:
                    print "C"
                    return False
        #print "."
        return True

    def ordered_flags(self, flag=[]):
        #yield flag
        nbd = set(self.items)
        for i in flag:
            nbd = nbd.intersection(self.nbd[i])
        #print "all_flags:", flag, nbd
        for i in nbd:
            if i in flag:
                continue
            yield flag+[i]
            for _flag in self.ordered_flags(flag+[i]):
                yield _flag

    def all_flags(self):
        flags = [flag for flag in self.ordered_flags()]
        for flag in flags:
            flag.sort()
        flags = [tuple(flag) for flag in flags]
        flags = list(set(flags))
        flags.sort()
        return flags

    def maximal_flags(self):
        "all maximal flags"
        flags = self.all_flags()
        maximal = []
        for f in flags:
            #print "is maximal %s ?"%(f,)
            for g in flags:
                if f is g:
                    continue
                if set(g).intersection(f) == set(f):
                    #print "%s contains %s" % (g, f)
                    break
            else:
                #print "maximal:", f
                maximal.append(f)
        return maximal

    def flags_type(self, tps):
        "all flags of type tps"
        tps = list(tps)
        if len(tps)==0:
            return [()]
        if len(tps)==1:
            return [(i,) for i in self.tplookup[tps[0]]]
        flags = []
        incidence = self.incidence
        for flag in self.flags_type(tps[:-1]):
            #print "flags_type", flag
            tp = tps[-1]
            for i in self.tplookup[tp]:
                for j in flag:
                    #print "?", (i, j)
                    if (i, j) not in incidence:
                        break
                else:
                    flags.append(flag + (i,))
        return flags

    def check_geometry(self):
        # every maximal flag is a chamber
        for flag in self.maximal_flags():
            if len(flag)!=self.rank:
                print "flag %s is not a chamber (rank=%d)"%(flag, self.rank)
                assert 0

    def is_digon(self):
        if self.rank != 2:
            return False
        tpi = self.types[0]
        tpj = self.types[1]
        for i in self.tplookup[tpi]:
          for j in self.tplookup[tpj]:
            if (i, j) not in self.incidence:
                return False
        return True

    def residue(self, flag):
        items = []
        for i in flag:
            assert i in self.items
        for i in self.items:
            if i in flag:
                continue
            for j in flag:
                if (i, j) not in self.incidence:
                    break
            else:
                items.append(i)
        items = list(set(items))
        incidence = [(i, j) for i in items for j in items if (i, j) in self.incidence]
        assert len(set(incidence))==len(incidence), str(incidence)
        #print "residue:", incidence
        tpmap = dict((i, self.tpmap[i]) for i in items)
        g = Geometry(incidence, tpmap)
        return g

    def get_diagram(self):
        diagram = []
        for (i, j) in choose(self.types, 2):
            tps = [k for k in self.types if k!=i and k!=j]
            assert len(tps)<20
            #print "flag_type:", tps
            for flag in self.flags_type(tps):
                h = self.residue(flag)
                assert h.types == [i, j] or h.types == [j, i]
                if h.is_digon():
                    break
            else:
                diagram.append((i, j))
        diagram.sort()
        return diagram

    def is_projective_plane(self):
        if len(self.types) != 2:
            return False
        ptp, ltp = self.types
        points = self.tplookup[ptp]
        lines = self.tplookup[ltp]
        nbd = self.nbd
        for line in lines:
            for other in lines:
                if other is line:
                    continue
                ps = set(nbd[line]).intersection(set(nbd[other]))
                if len(ps)!=1:
                    #print line
                    #print other
                    #print ps
                    return False
                p = list(ps)[0]
                assert self.tpmap[p]==ptp
        for i in points:
          for j in points:
            if i==j:
                continue
            ilines = set(nbd[i]).intersection(set(nbd[j]))
            if len(ilines)!=1:
                print "L"
                return False
            line = list(ilines)[0]
            assert self.tpmap[line]==ltp
        return True

    def get_system(self):
        flags = list(self.maximal_flags())

        graphs = {} # map tp -> graph
        for tp in self.types:
            graph = nx.Graph()
            for flag in flags:
                graph.add_node(flag)
                #print "node:", tp, flag
            graphs[tp] = graph

        for i in flags:
            ii = set(i)
            for j in flags:
                jj = set(j)
                kk = ii.intersection(jj)
                if len(kk)==1:
                    tp = self.tpmap[list(kk)[0]]
                    #print "edge:", tp, (i, j)
                    graphs[tp].add_edge(i, j)
        tpmap = {}
        for tp, graph in graphs.items():
            equ = nx.connected_components(graph)
            tpmap[tp] = equ
            #print "equ:", tp
            #for items in equ:
            #    print items
        system = System(flags, tpmap)
        return system

    def get_graph(self):
        graph = nx.Graph()
        for item in self.items:
            graph.add_node(item)
        for i, j in self.incidence:
            if i!=j and i<j:
                graph.add_edge(i, j)
        return graph

    def rel(self, a, b):
        incidence = self.incidence
        A = zeros2(len(a), len(b))
        for i, aa in enumerate(a):
          for j, bb in enumerate(b):
            A[i, j] = (aa, bb) in incidence
        return A

    def get_bag(self):
        items = list(self.items)
        items.sort()
        tpmap = self.tpmap
        points = []
        lookup = {}
        for idx, item in enumerate(items):
            tp = tpmap[item]
            p = Point(str(tp), idx)
            #p = Point("XX", idx)
            lookup[item] = p
            points.append(p)
        #print points
        for a, b in self.incidence:
            lookup[a].nbd.append(lookup[b])
        bag = Bag(points)
        for p in points:
            print p, p.nbd
        return bag
        
    def get_symmetry(self):
        bag0 = self.get_bag()
        bag1 = self.get_bag()
        count = 0
        for f in isomorph.search(bag0, bag1):
            print f
            count += 1
        return count


class System(object):
    """
        Chamber system
    """
    def __init__(self, flags, tpmap):
        # for every type, we have an equivelance relation on the Chambers
        "tpmap : map type -> list of list of flags "
        self.flags = flags
        self.tpmap = tpmap

    def __str__(self):
        return "System(%s, %s)"%(self.flags, self.tpmap)

    def get_graph(self):
        graph = nx.Graph()
        flags = self.flags
        for i, c in enumerate(flags):
            graph.add_node(i)
        for tp, flagss in self.tpmap.items():
            for flags in flagss:
                for flag in flags:
                  i = self.flags.index(flag)
                  for _flag in flags:
                    _i = self.flags.index(_flag)
                    graph.add_edge(i, _i)
        return graph

    def get_bag(self):
        flags = self.flags
        tpmap = self.tpmap
        lookup = {}
        points = []
        for flag in flags:
            p = Point("f", len(lookup))
            lookup[flag] = p
            points.append(p)
        for tp, flagss in tpmap.items():
            p = Point("t", len(lookup))
            lookup[tp] = p
            points.append(p)
            for i, flags in enumerate(flagss):
                e = Point("e", len(lookup))
                lookup[tp, i] = e
                points.append(e)
                p.nbd.append(e)
                e.nbd.append(p)
                for flag in flags:
                    e.nbd.append(lookup[flag])
                    lookup[flag].nbd.append(e)
        bag = Bag(points)
        return bag

    def get_symmetry(self):
        bag0 = self.get_bag()
        bag1 = self.get_bag()
        count = 0
        for f in isomorph.search(bag0, bag1):
            #print f
            count += 1
        return count


def genidx(*shape):
    if len(shape)==0:
        yield ()
    else:
        for idx in range(shape[0]):
            for _idx in genidx(*shape[1:]):
                yield (idx,)+_idx


def tesellate_3d(n=3):
    "cubic tesselation of 3-torus by n*n*n cubes and verts"

    assert n<10
    cubes = {}
    verts = {}
    incidence = []
    tpmap = {}
    for i,j,k in genidx(n,n,n):
        cubes[i,j,k] = 'c%d%d%d'%(i,j,k)
        verts[i,j,k] = 'v%d%d%d'%(i,j,k)
    for i,j,k in genidx(n,n,n):
        c = cubes[i,j,k]
        for di,dj,dk in genidx(2,2,2):
            v = verts[(i+di)%n, (j+dj)%n, (k+dk)%n]
            incidence.append((c, v))
    for c in cubes.values():
        tpmap[c] = 'c'
    for v in verts.values():
        tpmap[v] = 'v'
    g = Geometry(incidence, tpmap)
    return g



def fano():
    # https://www.maa.org/sites/default/files/pdf/upload_library/22/Ford/Lam305-318.pdf
    points = range(1, 8)
    lines = [(1,2,4), (2,3,5), (3,4,6), (4,5,7), (5,6,1), (6,7,2), (7,1,3)]
    assert len(lines)==7
    incidence = []
    tpmap = {}
    for point in points:
        tpmap[point] = 'p'
    for line in lines:
        tpmap[line] = 'l'
        for point in line:
            incidence.append((point, line))
    g = Geometry(incidence, tpmap)
    assert g.is_projective_plane()
    return g


def projective(n, dim=2):
    # Take n-dim F_2-vector space
    # points are subspaces of dimension 1
    # lines are subspaces of dimension 2
    # etc.

    def get_key(L):
        vs = [str(v) for v in span(L) if v.sum()]
        vs.sort()
        key = ''.join(vs)
        return key

    assert n>1

    points = []
    for P in enum2(n):
        if P.sum()==0:
            continue
        points.append(P)
    #print "points:", len(points)

    lines = []
    lookup = {}
    for L in enum2(2*n):
        L.shape = (2, n)
        L = row_reduce(L)
        if len(L)!=2:
            continue
        key = get_key(L)
        if key in lookup:
            continue
        lines.append(L)
        lookup[key] = L
    #print "lines:", len(lines)

    spaces = []
    if n>3 and dim>2:
        m = 3
        lookup = {}
        for A in enum2(m*n):
            A.shape = (m, n)
            A = row_reduce(A)
            if len(A)!=m:
                continue
            key = get_key(A)
            if key in lookup:
                continue
            spaces.append(A)
            lookup[key] = A
        print "spaces:", len(spaces)

    incidence = []
    tpmap = {} 
    for point in points:
        point = str(point)
        tpmap[point] = 0
        #print point

    for L in lines:
        line = str(tuple(tuple(row) for row in L))
        tpmap[line] = 1
        for P in span(L):
            if P.sum():
                incidence.append((str(P), line))

    for A in spaces:
        space = freeze(A)
        tpmap[space] = 2
        for P in span(A):
            if P.sum()==0:
                continue
            incidence.append((str(P), space))
        for L in lines:
            B = solve(A.transpose(), L.transpose())
            if B is not None:
                line = str(tuple(tuple(row) for row in L))
                incidence.append((space, line))

    g = Geometry(incidence, tpmap)
    if dim==2:
        assert g.get_diagram() == [(0, 1)]
    elif dim==3:
        assert g.get_diagram() == [(0, 1), (1, 2)]
    return g


def test():

    # triangle:
    g = Geometry(
        "aB aC bA bC cB cA".split(), 
        {'a':'p', 'b':'p', 'c':'p', 'A':'l', 'B':'l', 'C':'l'})

    assert len(g.items) == 6
    #print g.flags_type('l') 
    #print g.flags_type('p') 

    #g = g.residue(["a"])

    for n in range(2, 7):
        g = Geometry.polygon(n)
        assert len(g.items)==2*n
        assert len(g.maximal_flags())==2*n

        assert len(g.flags_type(['p'])) == n
        assert len(g.flags_type(['l'])) == n
        assert len(g.flags_type(['p', 'l'])) == 2*n
        if n>3:
            assert not g.is_projective_plane()

    g = Geometry.simplex(2)
    h = g.residue([(0,), (0,1)])
    assert len(h.items)==1
    assert len(g.flags_type([0, 1])) == 6

    for dim in range(4):
        g = Geometry.simplex(dim)
        assert len(g.items)==2**(dim+1)-1
        assert len(g.maximal_flags())==factorial(dim+1)

        h = g.residue([(0,)])

        if dim>0:
            g.residue([(0,), (0,1)])

    g = Geometry.cube()
    d = g.get_diagram()
    assert d == [('e', 'v'), ('f', 'e')] # woo !

    assert Geometry.simplex(4).get_diagram() == [(0, 1), (1, 2), (2, 3)] # A_4 Dynkin diagram !

    #g = tesellate_3d()

    g = fano()
    assert g.get_diagram() == [('p', 'l')]

    n = argv.get("n", 3)
    assert n>=3
    g = projective(n)

    assert len(g.tplookup[0]) == qbinomial(n, 1)
    assert len(g.tplookup[1]) == qbinomial(n, 2)

    if argv.cube:
        g = Geometry.cube()
    elif argv.simplex:
        g = Geometry.simplex(2)
    elif argv.triangle:
        g = Geometry.polygon(3)
    elif argv.projective:
        dim = argv.get("dim", 2)
        g = projective(n, dim)

        if dim>2:
            assert len(g.tplookup[2]) == qbinomial(n, 3)
    else:
        return

    if argv.system:
        s = g.get_system()
        graph = s.get_graph()
        #nx.draw(graph)
        #pyplot.draw()
        #pyplot.show()

    if argv.symmetry:
        #print g.get_symmetry()
        s = g.get_system()
        s.get_symmetry()
        return

    flags = list(g.maximal_flags())
    N = len(flags)
    if N>26:
        names = dict((flags[i], chr(ord('A')+(i//26))+chr(ord('A')+(i%26))) for i in range(N))
    else:
        names = dict((flags[i], chr(ord('A')+i)) for i in range(N))

    def le(f, g):
        return numpy.alltrue(f<=g)

    print "flags:", N
    a = flags[0]
    rank = len(a)
    table = {} # calculate all geometrical relationships (triangles)
    unit = freeze(g.rel(a, a))
    for b in flags:
        AB = g.rel(a, b)
        Bs = []
        for c in flags:
            BC = g.rel(b, c)
            AC = g.rel(a, c)
            key = freeze(AB), freeze(BC)
            #table.setdefault(key, []).append(AC)
            table.setdefault(key, set([])).add(freeze(AC))
            #print shortstrx(AB, BC, AC)
            #print

    items = list(table.items())
    items.sort()

    mul = {} # Murphy product
    bruhat = {} # Bruhat order bruhat[a, b] : a <= b
    star = {} # anti-involution (..?)
    elements = set()
    for key, values in items:
        elements.add(key[0])
        elements.add(key[1])
        AB, BC = thaw(key[0]), thaw(key[1])
        bruhat[key] = le(BC, AB)
        values = [thaw(f) for f in values]
        x = values[0]
        for y in values:
            #assert le(AC, y)
            if le(y, x):
                x = y
        #mul.append((AB, BC, x))
        #assert numpy.allclose(numpy.dot(AB, BC), x)
        mul[key] = freeze(x) # Murphy product
    assert freeze(g.rel(a, a)) in elements
    print "relations:", len(elements)

    order = {}
    for a in elements:
      for b in elements:
        order[a, b] = False
    for a in elements:
      for b in elements:
        order[a, mul[a, b]] = True
        order[b, mul[a, b]] = True

    #assert order == bruhat # nope..
    bruhat = order    

    for a in elements:
      assert bruhat[a, a]
      for b in elements:
        assert bruhat[a, mul[a, b]] 
        assert bruhat[b, mul[a, b]] 
        if a != b:
            assert not bruhat[a, b] or not bruhat[b, a]

    elements = list(elements)
    def cmp(g, h):
        if g==h:
            return 0
        elif bruhat[h,g]:
            return +1
        return -1
    elements.sort(cmp)

    lhom = {}
    rhom = {}
    for c in elements:
      for b in elements:
        for a in elements:
          if bruhat[c, mul[a, b]] and (lhom.get((c, b)) is None or bruhat[a, lhom[c, b]]):
            lhom[c, b] = a
          if bruhat[b, mul[c, a]] and (rhom.get((c, b)) is None or bruhat[a, rhom[c, b]]):
            rhom[c, b] = a

#    for a in elements:
#     for b in elements:
#      for c in elements:
#        if bruhat[c, mul[a, b]]:
#            assert bruhat[lhom[c, b], a]
#        if bruhat[b, mul[c, a]]:
#            assert bruhat[rhom[c, b], a]

    for u in elements:
      for v in elements:
        for w in elements:

            lhs = bruhat[w, mul[u, v]]
            rhs = bruhat[lhom[w, v], u]
            assert lhs==rhs

            lhs = int(bruhat[u, mul[v, w]])
            rhs = int(bruhat[rhom[v, u], w])
            #assert lhs==rhs

    #poset_homology(elements, bruhat)

    I = elements[0]
    L = elements[1]
    P = elements[2]
    assert mul[I, L] == L
    assert mul[L, P] in elements
    LP = mul[L, P]
    PL = mul[P, L]
    LPL = elements[-1]
    assert mul[LP, LPL] == LPL
    names = {I:"I", L:"L", P:"P", LP:"LP", PL:"PL", LPL:"LPL"}

    assert len(elements)==6
    for a in elements:
        assert a in names

    if argv.magnitude:
        magnitude_homology(g, flags, elements, names, mul, bruhat)


    return

    if argv.yang_baxter:
        for f in elements:
          for g in elements:
            try:
                if mul[f, mul[g, f]] == mul[g, mul[f, g]]:
                    print "+",
                else:
                    print "-",
            except KeyError:
                pass
        print

    for f in elements:
        assert mul[f, unit] == f
        assert mul[unit, f] == f

    if argv.assoc or 1:
        # Check assoc
        for f in elements:
          for g in elements:
            #mul[f, g]
            #meet[f, g]
            for h in elements:
                try:
                    assert mul[mul[f,g], h] == mul[f, mul[g, h]]
                except KeyError:
                    pass
    
    if argv.quantale or 1:
        # Check quantale
        for f in elements:
          for g in elements:
            for h in elements:
                try:
                    assert mul[f, meet[g, h]] == meet[mul[f, g], mul[f, h]]
                    assert mul[meet[g, h], f] == meet[mul[g, f], mul[h, f]]
                except KeyError:
                    pass

    print "\nOK"


def uniqtuples(itemss):
    if len(itemss)==0:
        yield ()
    elif len(itemss)==1:
        items = itemss[0]
        for head in items:
            yield (head,)
        return
    for head in itemss[0]:
        for tail in uniqtuples(itemss[1:]):
            if head != tail[0]:
                yield (head,)+tail


def alltuples(itemss):
    if len(itemss)==0:
        yield ()
    elif len(itemss)==1:
        items = itemss[0]
        for head in items:
            yield (head,)
        return
    for head in itemss[0]:
        for tail in alltuples(itemss[1:]):
            yield (head,)+tail


def magnitude_homology(g, flags, monoid, names, mul, bruhat):

    print "magnitude_homology:"

    N = argv.get("N", 2)

    I = monoid[0]
    for b in monoid:
        assert mul[I, b] == b

    for f in flags:
        assert g.is_flag(f)

    def length(chain):
        for f in chain:
            assert g.is_flag(f)
        assert len(chain)
        if len(chain)==1:
            return I
        A = g.rel(chain[0], chain[1])
        A = freeze(A)
        assert A in monoid
        for i in range(1, len(chain)-1):
            B = freeze(g.rel(chain[i], chain[i+1]))
            A = mul[A, B]
        return A

    chains = {} # map l,n -> []
    lengths = list(monoid)
    if argv.get("l"):
        name = argv.l
        for l in monoid:
            if names[l] == name:
                lengths = [l]
                break
    for l in lengths:
        assert l in names
        for n in range(N):

            if n==0 and l==I:
                nchains = [(flag,) for flag in flags] # objects of X
    
            else:
                nchains = []
                for chain in alltuples((flags,)*(n+1)):
                    if length(chain) == l:
                        nchains.append(chain)
            chains[l,n] = nchains
            print "l = %s, |%d-chains|=%d" % (names[l], n, len(nchains))

        bdys = {} # map 
        for n in range(N-1):
            # bdy maps n+1 chains to n chains
            print "bdy: chains[%d] -> chains[%d]" % (n+1, n),
            bdy = {}
            source = chains[l, n+1]
            target = chains[l, n]
            #print "|chains[%d]|=%d |chains[%d]|=%d" % (n+1, len(source), n, len(target))
            #for chain in target:
            #    print "\t", chain
            for col, chain in enumerate(source):
                assert len(chain)==n+2
                for i in range(n+2):
                    chain1 = chain[:i] + chain[i+1:]
                    assert len(chain1)==n+1
                    if length(chain1) == l:
                        #print chain1
                        #assert len(set(chain1))==n+1 # uniq NOT !
                        row = target.index(chain1)
                        bdy[row, col] = bdy.get((row, col), 0) + (-1)**i
            print "nnz:", len(bdy), "range:", set(bdy.values())
            bdys[n+1] = bdy

        # bdys[n]: map n-chains to (n-1)-chains

        for i in range(1, N-1):
            b1, b2 = bdys[i], bdys[i+1]
            b12 = compose(b1, b2)
            #print "len(bdy*bdy):", len(b12.values())
            for value in b12.values():
                assert value == 0, value

            find_homology(b1, b2, len(chains[l, i-1]), len(chains[l, i]), len(chains[l, i+1]))

            #if i==2 and names[l]=="L":
            #    return


def find_homology(g, f, *dims):

    print "find_homology"
    print dims[2], "-f->", dims[1], "-g->", dims[0]

    F = numpy.zeros((dims[1], dims[2]))
    for row, col in f.keys():
        v = f[row, col]
        F[row, col] = v
    #print shortstr(F)

    G = numpy.zeros((dims[0], dims[1]))
    for row, col in g.keys():
        v = g[row, col]
        G[row, col] = v
    #print shortstr(G)

    GF = numpy.dot(G, F)
    #print shortstr(GF)
    assert numpy.abs(GF).sum() == 0

    import homology
    if f and g:
        d = homology.bettiNumber(G, F)
        print "MH:", d

# See also:
# http://stackoverflow.com/questions/15638650/is-there-a-standard-solution-for-gauss-elimination-in-python


def find_homology_2(g, f, *dims):

    print "find_homology"
    print dims[2], "-f->", dims[1], "-g->", dims[0]

    F = numpy.zeros((dims[1], dims[2]))
    for row, col in f.keys():
        v = f[row, col]
        F[row, col] = v % 2
    #print shortstr(F)

    G = numpy.zeros((dims[0], dims[1]))
    for row, col in g.keys():
        v = g[row, col]
        G[row, col] = v % 2
    #print shortstr(G)

    GF = numpy.dot(G, F) % 2
    #print shortstr(GF)
    assert numpy.abs(GF).sum() == 0

    F = F.astype(numpy.int32)
    G = G.astype(numpy.int32)

    print "rank(f)=%d, rank(g)=%d" % (rank(F), rank(G))
    if g:
        kerg = find_kernel(G)
        print "ker(g):", len(kerg)
    if f:
        kerf = find_kernel(F)
        print "ker(f):", len(kerf)




def compose(g, f):
    # first f then g
    h = {}
    #print g.keys()
    for row, c in g.keys():
        for r, col in f.keys():
            if c==r:
                #print row, col
                h[row, col] = h.get((row, col), 0) + g[row, c] * f[r, col]

    return h



def test_projective_plane():

    names = {}
    I, L, P = elements[:3]
    LP = mul[L, P]
    PL = mul[P, L]
    LPL = PLP = mul[LP, L]
    assert LPL == elements[-1]
    BOT = elements[-1]
    names[I] = "I  "
    names[L] = "L  "
    names[P] = "P  "
    names[mul[L, P]] = "LP "
    names[mul[P, L]] = "PL "
    names[BOT] = "_|_"
    assert len(names)==6

    star = {I:I, L:L, P:P, LP:PL, PL:LP, LPL:LPL}

#    for a in elements:
#      for b in elements:
#        assert rhom[a, b] == star[lhom[star[b], star[a]]]
#    for a in elements:
#        assert star[star[a]] == a

    show_table(elements, names, mul, 'mul')
    print

    for i in elements:
        assert bruhat[i, i]
    assert bruhat[I, L]
    assert bruhat[L, LP]
    assert bruhat[LP, LPL]

    show_table(elements, names, lhom, "lhom")
    print 

    show_table(elements, names, rhom, "rhom")


def show_table(elements, names, mul, desc):

    print desc[:3],
    for i in elements:
        print names[i],
    print
    for ii, i in enumerate(elements):
      print names[i],
      for jj, j in enumerate(elements):
        k = mul[i, j]
        #print "%2d"%elements.index(k),
        print "%3s"%names.get(k, "?"),
      print

def freeze(A):
    items = [A.shape]
    for idx in genidx(*A.shape):
        items.append((idx, A[idx]))
    items = tuple(items)
    return items

def thaw(items):
    #idx, value = items[0]
    shape = items[0]
    A = numpy.zeros(shape, dtype=numpy.int32)
    for (idx, value) in items[1:]:
        A[idx] = value
    return A


def all_chains(n, above, chain=()):
    elements = above.keys()

    if n==1:
        for a in elements:
            yield (a,)
    elif n==2:
        for a in elements:
          for b in above[a]:
            yield (a, b)
        return

    # build a nested for-loop
    lines = [
        "def get_all(elements, above):",
        "  for a0 in elements:",
    ]
    indent = '   '
    for i in range(n-1):
        lines.append(indent + "for a%d in above[a%d]:" % (i+1, i))
        indent += ' '
    tpl = "(" + ','.join('a%d'%i for i in range(n)) + ")"
    lines.append(indent + "yield %s"%tpl)
    #print '\n'.join(lines)
    exec ''.join(l+'\n' for l in lines)

    for chain in get_all(elements, above):
        yield chain


def poset_homology(elements, order):

    above = {}
    for a in elements:
        above[a] = [b for b in elements if order[a, b] and b!=a] # a <= b

    r = len(elements) # 0-th homology
    homology = [r]
    n = 1
    while 1:

        chains = set(all_chains(n+1, above))
        for chain in chains:
            assert len(set(chain))==len(chain)==n+1

        size = len(chains)
        if size==0:
            break
        homology.append(size)
        r += (-1)**n * size

        n += 1

    print "poset_homology:", homology, "=", r

    return r


"""
TWF links:

http://math.ucr.edu/home/baez/week162.html # jordan algebras

Lie groups & geometry, Borel groups, parabolic subgroups:
http://math.ucr.edu/home/baez/week178.html
http://math.ucr.edu/home/baez/week180.html

the finite case & q-polynomials:
http://math.ucr.edu/home/baez/week186.html

"""

"""
buildings as enriched categories:
https://golem.ph.utexas.edu/category/2010/02/3000_and_one_thing_to_think_ab.html#c032006
https://ncatlab.org/nlab/show/building#buildings_as_metric_spaces

magnitude of enriched category:
https://golem.ph.utexas.edu/category/2011/06/the_magnitude_of_an_enriched_c.html

https://golem.ph.utexas.edu/category/2016/08/monoidal_categories_with_proje.html#c051238

magnitude homology of a graph:
https://golem.ph.utexas.edu/category/2015/05/categorifying_the_magnitude_of.html

magnitude homology in general:
https://golem.ph.utexas.edu/category/2016/09/magnitude_homology.html#more
"""



if __name__ == "__main__":

    test()



