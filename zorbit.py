#!/usr/bin/env python

import sys, os

import numpy

from solve import shortstr, parse, dot2, zeros2, array2
from solve import row_reduce, RowReduction, span

from lanczos import write, show_eigs


def genidx(shape):
    if len(shape)==0:
        yield ()
    else:
        for idx in range(shape[0]):
            for _idx in genidx(shape[1:]):
                yield (idx,)+_idx



def build_gcolor(size):

    from qupy.ldpc import gcolor

    lattice = gcolor.Lattice(size)

    n = len(lattice.qubits)
    print lattice

    code = lattice.build_code(check=False)
    #Ex = lattice.Ex
    Gx, Gz = code.Gx, code.Gz
    Hx, Hz = code.Hx, code.Hz

    return Gx, Gz, Hx


def build_compass(l):

    n = l**2

    keys = [(i, j) for i in range(l) for j in range(l)]
    coords = {}  
    for i, j in keys:
        for di in range(-l, l+1):
          for dj in range(-l, l+1):
            coords[i+di, j+dj] = keys.index(((i+di)%l, (j+dj)%l))

    m = n 
    Gx = zeros2(m, n)
    Gz = zeros2(m, n)

    idx = 0 
    for i in range(l):
      for j in range(l):
        Gx[idx, coords[i, j]] = 1 
        Gx[idx, coords[i, j+1]] = 1 

        Gz[idx, coords[i, j]] = 1 
        Gz[idx, coords[i+1, j]] = 1 
        idx += 1

    assert idx == m

    mx = l-1
    Hx = zeros2(mx, n)
    for idx in range(l-1):
      for j in range(l):
        Hx[idx, coords[j, idx]] = 1
        Hx[idx, coords[j, idx+1]] = 1

    return Gx, Gz, Hx


def build_isomorph(Gx):
    from isomorph import Tanner, search
    bag0 = Tanner.build(Gx)
    bag1 = Tanner.build(Gx)
    m, n = Gx.shape
    count = 0
    perms = []
    keys = range(m, len(bag0))
    #print "keys:", keys
    for fn in search(bag0, bag1):
        #print fn
        perm = tuple(fn[i]-m for i in keys)
        perms.append(perm)
        count += 1
    print "isomorphisms:", count
    #print len(set(perms))
    return perms


def orbiham(H):
    from isomorph import from_ham, search
    import networkx as nx

    n = len(H)
    graph = nx.Graph()
    for i in range(n):
        graph.add_node(i)

    bag0 = from_ham(H)
    bag1 = from_ham(H)

    count = 0
    for fn in search(bag0, bag1):
        #print fn
        write('.')
        for i, j in fn.items():
            graph.add_edge(i, j)
        count += 1
    print

    equs = nx.connected_components(graph)
    m = len(equs)

    print "isomorphisms:", count
    print "components:", m

    P = numpy.zeros((n, m))
    Q = numpy.zeros((m, n))
    for i, equ in enumerate(equs):
        for j in equ:
            P[j, i] = 1
        Q[i, j] = 1

    #print shortstr(P)
    #print shortstr(Q)

    H = numpy.dot(Q, numpy.dot(H, P))
    return H


def sparse_orbiham(n, H):
    from isomorph import from_sparse_ham, search
    import networkx as nx

    graph = nx.Graph()
    for i in range(n):
        graph.add_node(i)

    bag0 = from_sparse_ham(n, H)
    bag1 = from_sparse_ham(n, H)

    print "search..."

    count = 0
    for fn in search(bag0, bag1):
        #print fn
        write('.')
        for i, j in fn.items():
            graph.add_edge(i, j)
        count += 1
    print

    equs = nx.connected_components(graph)
    m = len(equs)

    print "isomorphisms:", count
    print "components:", m

    P = numpy.zeros((n, m))
    Q = numpy.zeros((m, n))
    for i, equ in enumerate(equs):
        for j in equ:
            P[j, i] = 1
        Q[i, j] = 1

    #print shortstr(P)
    #print shortstr(Q)

    H = numpy.dot(Q, numpy.dot(H, P))
    return H


def sparse_orbiham_2(n, H):
    from isomorph import from_sparse_ham, search
    import networkx as nx

    graph = nx.Graph()
    for i in range(n):
        graph.add_node(i)

    bag0 = from_sparse_ham(n, H)
    bag1 = from_sparse_ham(n, H)

    print "search..."

    count = 0
    for fn in search(bag0, bag1):
        #print fn
        write('.')
        for i, j in fn.items():
            graph.add_edge(i, j)
        count += 1
    print

    equs = nx.connected_components(graph)
    m = len(equs)

    print "isomorphisms:", count
    print "components:", m

    P = numpy.zeros((n, m))
    Q = numpy.zeros((m, n))
    for i, equ in enumerate(equs):
        for j in equ:
            P[j, i] = 1
        Q[i, j] = 1

    #print shortstr(P)
    #print shortstr(Q)

    H = numpy.dot(Q, numpy.dot(H, P))
    return H


def main():

    if argv.gcolor:
        size = argv.get("size", 1)
        Gx, Gz, Hx = build_gcolor(size)

    elif argv.compass:
        l = argv.get('l', 3)
        Gx, Gz, Hx = build_compass(l)

    elif argv.gcolor2:

        from gcolor import build as build1
        Gx, Gz, Hx = build1()

    else:

        return

    #print shortstr(Gx)

    n = Gx.shape[1]

#    perms = [tuple(range(n))]
#    if argv.isomorph:
#        perms = build_isomorph(Gx)
#    #print perms[:10]

    print "Hx:"
    print shortstr(Hx)

    Hxr = row_reduce(Hx)

    print "Hxr:"
    print shortstr(Hxr)

    rr = RowReduction(Hx)

    #print "Gx:"
    #print shortstr(Gx)

    Gxr = []
    for g in Gx:
        g = rr.reduce(g)
        Gxr.append(g)
    Gxr = array2(Gxr)

    Gxr = row_reduce(Gxr, truncate=True)

    print "Gxr:", len(Gxr)
    print shortstr(Gxr)

    import networkx as nx
    graph = nx.Graph()

    verts = []
    lookup = {}
    for i, v in enumerate(span(Gxr)):
        lookup[v.tostring()] = i
        verts.append(v)
    print "span:", len(verts)

    mz = len(Gz)
    n = len(verts)

    if n <= 1024 and not argv.sparse:
        H = numpy.zeros((n, n))
        for i, v in enumerate(verts):
            count = dot2(Gz, v).sum()
            H[i, i] = mz - 2*count
            for g in Gx:
                v1 = (g+v)%2
                v1 = rr.reduce(v1)
                j = lookup[v1.tostring()]
                H[i, j] += 1
    
        print H
    
        vals, vecs = numpy.linalg.eig(H)
        show_eigs(vals)

        if argv.orbiham:
            H1 = orbiham(H)
            print "orbiham:"
            print H1
            vals, vecs = numpy.linalg.eig(H1)
            show_eigs(vals)

    else:
        print "building H..."
        H = {}

        for i, v in enumerate(verts):
            if i%1000==0:
                write('.')
            count = dot2(Gz, v).sum()
            H[i, i] = mz - 2*count
            for g in Gx:
                v1 = (g+v)%2
                v1 = rr.reduce(v1)
                j = lookup[v1.tostring()]
                H[i, j] = H.get((i, j), 0) + 1
    
        print "nnz:", len(H)

        if argv.orbiham:
            H1 = sparse_orbiham(len(verts), H)
            print "orbiham:"
            print H1
            vals, vecs = numpy.linalg.eig(H1)
            show_eigs(vals)



from argv import Argv
argv = Argv()

if __name__ == "__main__":

    main()





