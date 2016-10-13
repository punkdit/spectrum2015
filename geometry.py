#!/usr/bin/env python

import os, sys

import networkx as nx
from matplotlib import pyplot

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
                print "%s not a chamber (rank=%d)"%(flag, self.rank)
                assert 0, self

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


def projective(n):
    # Take n-dim F_2-vector space
    # points are subspaces of dimension 1
    # lines are subspaces of dimension 2
    # etc.
    from solve import zeros2, enum2, row_reduce, span, shortstr

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

    # TOO SLOW:
#    for m in range(3, n+1):
#        spaces = []
#        lookup = {}
#        for A in enum2(m*n):
#            A.shape = (m, n)
#            A = row_reduce(A)
#            if len(A)!=m:
#                continue
#            key = get_key(A)
#            if key in lookup:
#                continue
#            spaces.append(A)
#            lookup[key] = A
#        print "spaces:", len(spaces)

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

    g = Geometry(incidence, tpmap)
    assert g.get_diagram() == [(0, 1)]
    return g


def test():

    desc = r"""

    a
    |\ B
    | \
   C|  c
    | /
    |/ A
    b

    """

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
    g = projective(n)

    assert len(g.tplookup[0]) == qbinomial(n, 1)
    assert len(g.tplookup[1]) == qbinomial(n, 2)

    s = g.get_system()
    #print s

    graph = s.get_graph()
    #nx.draw(graph)
    #pyplot.draw()
    #pyplot.show()


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


"""
https://ncatlab.org/nlab/show/building#buildings_as_metric_spaces
https://golem.ph.utexas.edu/category/2010/02/3000_and_one_thing_to_think_ab.html#c032006

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



