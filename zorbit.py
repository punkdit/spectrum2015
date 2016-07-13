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


def sparse_orbiham_nx(n, H):
    import networkx as nx
    from networkx.algorithms.isomorphism import GraphMatcher

    bag = nx.Graph()
    for i in range(n):
        bag.add_node(i, syndrome=H[i,i])

    for i in range(n):
      for j in range(n):
        if i==j:
            continue
        if H.get((i, j)):
            bag.add_edge(i, j)

    def node_match(n0, n1):
        return n0['syndrome'] == n1['syndrome']

    matcher = GraphMatcher(bag, bag, node_match=node_match)

    print "search..."

    graph = nx.Graph()
    for i in range(n):
        graph.add_node(i)

    count = 0
    for iso in matcher.isomorphisms_iter(): # too slow :P
        #print iso
        write('.')
        for i, j in iso.items():
            graph.add_edge(i, j)
        count += 1
    print

    equs = nx.connected_components(graph)
    m = len(equs)

    print "isomorphisms:", count
    print "components:", m


def sparse_orbiham(n, H):
    from isomorph import from_sparse_ham, search
    import networkx as nx

    graph = nx.Graph()
    for i in range(n):
        graph.add_node(i)

    print "bag0"
    bag0 = from_sparse_ham(n, H)

    print "bag1"
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


def sparse_orbiham_nauty(n, degree, A, U):

    idxs = A.keys()
    idxs.sort()
    edges = n*degree

    rows = dict((i, []) for i in range(n))
    for (i, j) in idxs:
        rows[i].append(j)
        
    import pnauty
    print "pnauty.init_graph", n, edges, degree
    assert n*degree == edges
    pnauty.init_graph(n, edges, degree)

    i0 = None
    #for idx, weight in A.items():
    for idx in idxs: # sorted
        i, j = idx
        assert A[idx] > 0
        for count in range(A[idx]):
            if i != i0:
                edge = 0
                i0 = i
            else:
                edge += 1
            assert edge < degree
            #print "add_edge", (i, j, edge)
            pnauty.add_edge(i, j, edge)

    #pnauty.search()
    #return

    verts = range(n)
    verts.sort(key = lambda i : U[i])

    label = None
    for idx in range(n-1):
        i0 = verts[idx]
        i1 = verts[idx+1]
        if U[i0] == U[i1]:
            pnauty.set_partition(idx, i0, 1)
        else:
            pnauty.set_partition(idx, i0, 0)
    pnauty.set_partition(n-1, verts[-1], 0)

    print "pnauty.search"
    orbits = pnauty.search()
    canonical = list(set(orbits))
    canonical.sort()
    N = len(canonical)
    print "orbits:", N

    lookup = {}
    for i, idx in enumerate(orbits):
        lookup[i] = canonical.index(idx)

    H = numpy.zeros((N, N))
    for idx in range(N):
        i = canonical[idx]
        # i indexes a row of A
        for j in rows[i]:
            value = A[i, j]
            jdx = lookup[j]
            H[idx, jdx] += value

        H[idx, idx] += U[i]

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

    mr = len(Gxr)
    print "Gxr:", mr
    print shortstr(Gxr)

    v0 = None
    excite = argv.excite
    if excite is not None:
        v0 = zeros2(n)
        v0[excite] = 1

    verts = []
    lookup = {}
    for i, v in enumerate(span(Gxr)): # XXX does not scale well
        if v0 is not None:
            v = (v+v0)%2
        lookup[v.tostring()] = i
        verts.append(v)
    print "span:", len(verts)
    assert len(lookup) == len(verts)

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
        print "building H",
        A = {} # adjacency
        U = [] # potential

        for i, v in enumerate(verts):
            if i%1000==0:
                write('.')
            count = dot2(Gz, v).sum()
            #H[i, i] = mz - 2*count
            U.append(mz - 2*count)
            for g in Gx:
                v1 = (g+v)%2
                v1 = rr.reduce(v1)
                j = lookup[v1.tostring()]
                A[i, j] = A.get((i, j), 0) + 1
    
        print "\nnnz:", len(A)
        #print A
        degrees = {} # out-degree
        for key, value in A.items():
            i, j = key
            degrees[i] = degrees.get(i, 0) + value
        degrees = list(set(degrees.values()))
        print "degrees:", degrees
        assert len(degrees) == 1
        degree = degrees[0]

        #for value in A.values():
        #    assert value==1

        if not argv.orbiham:
            return

        H1 = sparse_orbiham_nauty(len(verts), degree, A, U)
        if H1 is None:
            return

        print "orbiham:"
        print H1
        if len(H1)<=1024 and 0:
            vals, vecs = numpy.linalg.eig(H1)
        else:
            from scipy.sparse.linalg import eigs
            vals, vecs = eigs(H1, k=min(len(H1)-5, 40), which="LM")

        show_eigs(vals)


from argv import Argv
argv = Argv()

if __name__ == "__main__":

    main()





