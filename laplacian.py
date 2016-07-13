#!/usr/bin/env python

import sys
from math import *

import numpy

from solve import enum2, zeros2, shortstr

from argv import Argv
argv = Argv()


EPSILON = 1e-8



def lstr2(xs, decimals=3):
    m, n = xs.shape
    lines = [lstr(xs[i, :], decimals) for i in range(m)]
    return '[%s]'%',\n '.join(line for line in lines)


def lstr(xs, decimals=4):
    if len(xs.shape)>1:
        return lstr2(xs, decimals)
    rxs = []
    for x in xs: 
        if abs(x.imag)<EPSILON:
            x = x.real
        if abs(x-round(x))<EPSILON:
            x = int(round(x))
        x = "%.*f"%(decimals, x)
        #if x.replace('0', '')=='.':
        #    x = '0.'
        rxs.append(x)
    s = '[%s]'%', '.join(x.rjust(decimals+2) for x in rxs)
    s = s.replace(" 0,", "  ,")
    s = s.replace(" 0]", "  ]")
    return s


class Graph(object):

    def __init__(self, verts, edges, potential={}, **deco):
        self.verts = verts
        self.n = len(verts)
        self.edges = edges
        self.lookup = dict((vert, i) for (i, vert) in enumerate(verts))
        #self.potential = dict((i,potential.get(v, 0)) for i,v in enumerate(verts))
        self.potential = dict(potential)
        bdy = dict((v, set()) for v in verts)
        for (v1, v2) in edges:
            bdy[v1].add(v2)
        self.bdy = bdy
        self.__dict__.update(deco)

    def __str__(self):
        verts, edges = self.verts, self.edges
        #return "Graph(%d:%s, %d:%s)"%(len(verts), verts, len(edges), edges)
        return "Graph(%d, %d)"%(len(verts), len(edges))

    def get_adj(self):
        verts, edges = self.verts, self.edges
        potential = self.potential
        lookup = self.lookup
        n = self.n
        A = zeros2(n, n)
        for v0, v1 in edges:
            i = lookup[v0]
            j = lookup[v1]
            A[i, j] = 1
            assert i!=j
        for i in range(n):
            #A[i, i] = potential.get(i, 0)
            A[i, i] = potential.get(verts[i], 0)
        return A

    def get_laplacian(self):
        verts, edges = self.verts, self.edges
        potential = self.potential
        lookup = self.lookup
        n = self.n
        A = self.get_adj().astype(numpy.float64)
        D = []
        I = numpy.zeros((n, n))
        for i in range(n):
            d = A[i].sum()
            D.append(d**-0.5)
            I[i, i] = 1.
        D = numpy.array(D)
        #print "D:", D
        D.shape = (n, 1)
        AD = A*D
        #print "A:", A
        D.shape = (1, n)
        DAD = D*AD
        #print "A:", A
        L = I - DAD
        return L

    def get_laplacian2(self):
        # appears to be the same as get_laplacian()
        verts, edges = self.verts, self.edges
        potential = self.potential
        lookup = self.lookup
        n = self.n
        A = self.get_adj().astype(numpy.float64)
        D = [] # degrees 
        for i in range(n):
            d = 1.*A[i].sum()
            D.append(d)
        D = numpy.array(D)
        L = numpy.zeros((n, n))
        for i in range(n):
            for j in range(n):
                if i==j:
                    L[i, j] = 1 - A[i, i] / D[i]
                elif A[i, j]:
                    L[i, j] = -A[i, j] / sqrt(D[i]*D[j])
        return L

#    def get_combinatorial_laplacian(self):
#        A = self.get_adj().astype(numpy.float64)
#        L = -A
#        for i in range(n):
#            d = 1.*A[i].sum()
#            D.append(d)

    @classmethod
    def cube(cls, d, w=0):
        verts = [shortstr(v, zero='0') for v in enum2(d)]
        verts.sort(key = lambda v : v.count('1'))
        print verts
        potential = {}
        for v in verts:
            r = 2*w*(d - v.count('1'))
            assert r>=0
            potential[v] = r
        edges = []
        for v1 in verts:
          for v2 in verts:
            if v1==v2:
              continue
            bv1 = eval('0b'+v1)
            bv2 = eval('0b'+v2)
            u = bv1 ^ bv2
            #print v1, v2, bin(u)
            u = bin(u)
            if u.count('1')==1:
                edges.append((v1, v2))
        return cls(verts, edges, potential)

    def dump(self):
        #print self
        #print len(self.verts)
        #print len(self.edges)
        print "Graph(%d, %d)"%(len(self.verts), len(self.edges))

        A = self.get_adj()

        #print lstr(A, 0)
        assert numpy.allclose(A, A.transpose())

        vals, vecs = numpy.linalg.eigh(A)
        print "eigvals:", lstr(vals[-20:])
        vec = vecs[:, -1]
        if vec[0] < 0:
            vec = -1*vec
        print "eigvec:", lstr(vec)
        vec = numpy.log(vec)
        print "log eigvec:", lstr(vec)
        gap = vals[-1] - vals[-2]
        print "eigvec_2:", lstr(vecs[:, -2])
        print "gap:", gap





"""

The spectral gap of the normalized laplacian shrinks like 1/N where
N is the degree of the graph. The spectral gap of the
adjacency matrix is constant. 

The laplacian eigenvalue bunches arithmetically with
the potential, whereas the adjacency matrix eigenvalue
bunches geometrically with the potential.

No idea how to leverage a lower bound on the laplacian
gap into a lower bound on the adjacency matrix gap.

"""


def main():

    d = argv.get("d", 3)
    w = argv.get("w", 0)

    graph = Graph.cube(d, w)

    print graph

    A = graph.get_adj()

    print "A:"
    print A
    vals, vecs = numpy.linalg.eigh(A)
    print list(vals)
    print vecs[:, -1]
    print

    #graph.dump()

    L = graph.get_laplacian()
    #print lstr2(L)

    L2 = graph.get_laplacian2()
    assert numpy.allclose(L, L2)

    vals, vecs = numpy.linalg.eigh(L)
    print list(vals)
    print vecs[:, 0]


if __name__ == "__main__":

    main()


