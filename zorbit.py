#!/usr/bin/env python

import sys, os
import math

import numpy
from numpy import log2
from scipy import sparse
from scipy.sparse.linalg import eigs, eigsh

from solve import shortstr, shortstrx, parse, eq2, dot2, zeros2, array2, identity2
from solve import row_reduce, RowReduction, span, get_reductor
from solve import u_inverse, find_logops, solve, find_kernel, linear_independent
from solve import System, Unknown, pseudo_inverse, enum2
from solve import find_errors, find_stabilizers, check_commute
from solve import minweightall

from isomorph import Tanner, search, from_sparse_ham, search_recursive, Backtrack, from_ham
from lanczos import write, show_eigs
from code import lstr2


def genidx(shape):
    if len(shape)==0:
        yield ()
    else:
        for idx in range(shape[0]):
            for _idx in genidx(shape[1:]):
                yield (idx,)+_idx



def build_isomorph(Gx):
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


def get_cycles(perm):
    f = dict((i, j) for i, j in enumerate(perm) if i!=j)
    cycles = []
    remain = set(f.keys())
    while remain:
        i = iter(remain).next()
        cycle = [i]
        while 1:
            remain.remove(i)
            j = f[i]
            if j == cycle[0]:
                break
            assert j not in cycle
            cycle.append(j)
            i = j
        cycles.append(cycle)
    return cycles


def write_gap(perms):
    name = argv.get("name", "zorbit.gap")
    print "open(%r, 'w')"%(name,)
    f = open(name, 'w')
    f.write("G:=Group(")
    for i, perm in enumerate(perms):
        #print len(perm)
        #print str(perm)
        cycles = get_cycles(perm)
        if not cycles:
            continue
        s = ''.join(str(tuple(c)) for c in cycles)
        s = s.replace(' ', '') 
        f.write(s)
        if i+1<len(perms):
            f.write(",\n  ")
    f.write(");\n")
    f.write("G1:= SmallerDegreePermutationRepresentation(G);;\n");
    f.write("Image(G1);\n");
    f.close()


def build_orbigraph(H, syndromes=None):
    import networkx as nx

    if syndromes is not None:
        return build_orbigraph_syndromes(H, syndromes)

    n = len(H)
    graph = nx.Graph()
    for i in range(n):
        graph.add_node(i)

    bag0 = from_ham(H)
    bag1 = from_ham(H)

    count = 0
    fs = set()
    for fn in search(bag0, bag1):
        f = [None]*n
        #print fn
        for i, j in fn.items():
            assert i<n and j<n
            f[i] = j
            graph.add_edge(i, j)
        #for i, j in enumerate(f):
        #    if i!=j:
        #        print (i,j),
        #print
        f = tuple(f)
        if f in fs:
            #write('/')
            continue # <---- continue
        fs.add(f)
        write('.')
        count += 1
    print
    print "isomorphisms:", count
    if argv.gap:
        write_gap(fs)

    equs = nx.connected_components(graph)
    equs = list(equs)
    m = len(equs)

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


def build_orbigraph_syndromes(H, syndromes):
    import networkx as nx

    print "build_orbigraph_syndromes"

    n = len(H)
    graph = nx.Graph()
    for i in range(n):
        graph.add_node(i)

    bag00 = from_ham(H)
    bag11 = from_ham(H)

    count = 0
    for f0 in search(bag00, bag11):

        bag0 = from_ham(H, syndromes)
        bag1 = from_ham(H, syndromes)

        try:
            for f1 in search(bag0, bag1, fn=dict(f0)):
                write('.')
                break
        except Backtrack:
            continue

        for i, j in f0.items():
            graph.add_edge(i, j)
        count += 1
    print
    print "isomorphisms:", count

    equs = nx.connected_components(graph)
    m = len(equs)

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


def sparse_orbigraph_nx(n, H):
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


def sparse_orbigraph(n, H):
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


def sparse_orbigraph_nauty(degree, A, U):

    n = len(U)

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


def do_lanczos(A, U):
    print "do_lanczos: building Hs"
    H = dict(A)
    for i in range(len(U)):
        H[i, i] = H.get((i, i), 0) + U[i]
    #H1 = sparse.lil_matrix((len(U), len(U)))
    keys = H.keys()
    keys.sort()
    data = []
    rows = []
    cols = []
    for idx in keys:
        #H1[idx] = H[idx]
        data.append(H[idx])
        rows.append(idx[0])
        cols.append(idx[1])
    H1 = sparse.coo_matrix((data, (rows, cols)), (len(U), len(U)))
    H1 = sparse.csr_matrix(H1, dtype=numpy.float64)

    print "do_lanczos: eigsh"
    vals, vecs = eigsh(H1, k=min(len(U)-5, 40), which="LM")

    return vals, vecs


def do_orbigraph(A, U):
    print "do_orbigraph"
    #print A
    degrees = {} # out-degree
    for key, value in A.items():
        i, j = key
        degrees[i] = degrees.get(i, 0) + value
    degrees = list(set(degrees.values()))
    #print "degrees:", degrees
    assert len(degrees) == 1
    degree = degrees[0]

    #for value in A.values():
    #    assert value==1

    H1 = sparse_orbigraph_nauty(degree, A, U)
    if H1 is None:
        return

    N = len(H1)
    rows = range(N)
    rows.sort(key = lambda i : -H1[i, i])

    if argv.sort:
        H1 = H1[rows, :]
        H1 = H1[:, rows]

    if argv.show and N < 40:
        print lstr2(H1, 0)
    elif argv.show:
        print "orbigraph:"
        for i in rows:
            print "%d:"%i,
            for j in rows:
                if H1[i, j]:
                    print "%d:%d"%(j, H1[i, j]),
            print

    if argv.sympy:
        from sympy import Matrix
        A = H1.astype(numpy.int32)
        A = Matrix(A)
        print A.eigenvals()

    if len(H1)<=1024:
        #print "numpy.linalg.eig"
        vals, vecs = numpy.linalg.eig(H1)
    else:
        vals, vecs = eigs(H1, k=min(len(H1)-5, 40), which="LM")

    return vals, vecs


class Code(object):
    def __init__(self, filename=None):
        self.filename = filename
        self.lines = []
        self.dent = 0
    def indent(self):
        self.dent += 1
    def dedent(self):
        self.dent -= 1
    def begin(self):
        self.append("{")
        self.indent()
    def end(self):
        self.dedent()
        self.append("}")
    def append(self, line):
        line = '    '*self.dent+line
        self.lines.append(line)
    def output(self):
        if self.filename is None:
            self.lines.append("")
            return '\n'.join(self.lines)
        f = open(self.filename, 'w')
        for line in self.lines:
            print >>f, line
        f.close()


def dense(Gx, Gz, Hx, Hz, Rx, Rz, Pxt, Qx, Pz, Tx, **kw):

    r, n = Rx.shape

    N = 2**r

    gz = len(Gz)
#    print "Hz:"
#    print shortstr(Hz)
    print "Hx|Tx:"
    print shortstrx(Hx, Tx)
    print "Hx:"
    for i, h in enumerate(Hx):
        print i, shortstr(h), h.sum()
    print "GzTx"
    GzTx = dot2(Gz, Tx.transpose())
    for i, h in enumerate(GzTx.transpose()):
        print i, shortstr(h), h.sum()

#    print "Rx:"
#    print shortstr(Rx)
    print "Tx:", len(Tx)
    #print shortstr(Tx)

    RR = dot2(Gz, Rx.transpose())
    PxtQx = dot2(Pxt, Qx)

    excite = argv.excite
    if excite is None:
        excite = kw.get("excite")

    if excite:
        excites = [excite]
    else:
        excites = genidx((2,)*len(Tx))

    vec0 = None

    for excite in excites:

        print "excite:", excite
        assert len(excite)==len(Tx)

        t = zeros2(n)
        for i, ex in enumerate(excite):
            if ex:
                t = (t + Tx[i])%2
        #print "t:", shortstr(t)
        Gzt = dot2(Gz, t)
        #print "Gzt:", shortstr(Gzt)

        if N<=1024:
            H = numpy.zeros((N, N))
        else:
            H = None
        A = {}
        U = []

        #for i in range(N):
        pos = neg = 0
        for i, v in enumerate(genidx((2,)*r)):
            v = array2(v)
            syndrome = (dot2(Gz, Rx.transpose(), v) + Gzt)%2
            value = gz - 2*syndrome.sum()
            #print shortstr(dot2(Rx.transpose(), v)), value
            if H is not None:
                H[i, i] = value
            U.append(value)

        for i, v in enumerate(genidx((2,)*r)):
            v = array2(v)
            #print shortstr(v),
            for g in Gx:
                u = (v + dot2(g, PxtQx))%2
                j = eval('0b'+shortstr(u, zero='0'))
                if H is not None:
                    H[i, j] += 1
                A[i, j] = A.get((i, j), 0) + 1

        #print H

        if argv.orbigraph:
            vals, vecs = do_orbigraph(A, U)
            show_eigs(vals)

        #if vec0 is not None:
        #    Hv = numpy.dot(H, vec0)
        #    print numpy.dot(vec0, Hv),
        #    Hv /= numpy.linalg.norm(Hv)
        #    print numpy.dot(vec0, Hv)

        if argv.solve:
            assert N <= 1024
            assert numpy.allclose(H, H.transpose())
            vals, vecs = numpy.linalg.eigh(H)
            if argv.show_eigs:
                show_eigs(vals)
            vals = list(vals)
            vals.sort()
            val0 = vals[-1] # top one is last
            assert vals[-2] < val0 - 1e-4
            print "excite:", excite,
            print "eigval:", val0
            vec0 = vecs[:,-1]
            if vec0[0] < 0:
                vec0 = -vec0

        #break
        if argv.plot:
            from pyx import graph

            xdata = U
            ydata = list(vec0)

            g = graph.graphxy(
                width=16,
                x=graph.axis.linear(reverse=True),
                y=graph.axis.log(min=0.8*vec0.min(), max=1.0))
            # either provide lists of the individual coordinates
            g.plot(graph.data.values(x=xdata, y=ydata))
            # or provide one list containing the whole points
            #g.plot(graph.data.points(list(zip(range(10), range(10))), x=1, y=2))
            g.writePDFfile("pic-groundstate.pdf")

            return



def show_delta(Gx, Gz, Hx, Hz, Rx, Rz, Pxt, Qx, Pz, Tx, **kw):

    r, n = Rx.shape

    N = 2**r

    gz = len(Gz)
    GzTx = dot2(Gz, Tx.transpose())
    RR = dot2(Gz, Rx.transpose())
    PxtQx = dot2(Pxt, Qx)

    if argv.excite:
        excites = [argv.excite]
    else:
        excites = genidx((2,)*len(Tx))

    for excite in excites:

        print "excite:", excite
        assert len(excite)==len(Tx)

        t = zeros2(n)
        for i, ex in enumerate(excite):
            if ex:
                t = (t + Tx[i])%2
        print "t:", shortstr(t)
        Gzt = dot2(Gz, t)
        #print "Gzt:", shortstr(Gzt)

        #for i in range(N):
        pos = neg = 0
        for i, v in enumerate(genidx((2,)*r)):
            v = array2(v)
            syndrome0 = dot2(Gz, Rx.transpose(), v)
            syndrome = (dot2(Gz, Rx.transpose(), v) + Gzt)%2
            delta = syndrome.sum() - syndrome0.sum()
            if delta > 0:
                pos += 1
            elif delta < 0:
                neg += 1
            if delta:
                print "%3d %3d" % (gz-2*syndrome0.sum(), gz-2*syndrome.sum())
        print pos, neg, "total:", i+1


def dense_full(Gx, Gz, Hx, Hz, Rx, Pxt, Qx, Pz, Tx, **kw):
    " find orbigraph for hamiltonian component Gamma "

    gz, n = Gz.shape

    if argv.excite:
        excites = [argv.excite]

    else:
        excites = genidx((2,)*len(Tx))

    for excite in excites:

        print "excite:", excite
        assert len(excite)==len(Tx)

        t = zeros2(n)
        for i, ex in enumerate(excite):
            if ex:
                t = (t + Tx[i])%2
        #print "t:", shortstr(t)
        Gzt = dot2(Gz, t)
        #print "Gzt:", shortstr(Gzt)

        # This is our basis
        Bx = array2([v+t for v in Rx] + [v+t for v in Hx])
        Bx %= 2
        r = len(Bx)
        N = 2**r
        Bx = row_reduce(Bx, truncate=True)
        assert len(Bx)==r # linearly independent rows
        Cx = u_inverse(Bx)

        if N<=1024:
            H = numpy.zeros((N, N))
        else:
            H = None
        A = {}
        U = []

        #for i in range(N):
        pos = neg = 0
        for i, v in enumerate(genidx((2,)*r)):
            v = array2(v)
            syndrome = dot2(Gz, Bx.transpose(), v)
            value = gz - 2*syndrome.sum()
            #print shortstr(dot2(Rx.transpose(), v)), value
            if H is not None:
                H[i, i] = value
            U.append(value)

        for i, v in enumerate(genidx((2,)*r)):
            v = array2(v)
            u0 = dot2(Bx.transpose(), v)
            #print shortstr(v),
            for g in Gx:
                u1 = (u0 + g) % 2
                v1 = dot2(Cx.transpose(), u1)
                assert v1.shape == v.shape
                j = eval('0b'+shortstr(v1, zero='0'))
                if H is not None:
                    H[i, j] += 1
                A[i, j] = A.get((i, j), 0) + 1

        #print H

        if argv.orbistab:
            hx = Hx[0]
            do_orbistab(Bx, Cx, H, hx)

        if argv.orbigraph:
            vals, vecs = do_orbigraph(A, U)
            show_eigs(vals)

        if argv.solve:
            assert N <= 1024
            assert numpy.allclose(H, H.transpose())
            vals, vecs = numpy.linalg.eigh(H)
            if argv.show_eigs:
                show_eigs(vals)
            #print vals
            vals = list(vals)
            vals.sort()
            val0 = vals[-1] # top one is last
            vec0 = vecs[:,-1]
            if vec0[0] < 0:
                vec0 = -vec0
            assert vals[-2] < val0 - 1e-4
            print "excite:", excite,
            print "eigval:", val0



def do_orbistab(Bx, Cx, H, sx):

    r, n = Bx.shape
    N = 2**r

    S = numpy.zeros((N, N))
    for i, v in enumerate(genidx((2,)*r)):
        v = array2(v)
        u = dot2(Bx.transpose(), v)

        u1 = (u+sx)%2
        v1 = dot2(Cx.transpose(), u1)
        j = eval('0b'+shortstr(v1, zero='0'))
        S[i, j] = 1.0
    assert numpy.allclose(S, S.transpose())

    P = 0.5*(numpy.eye(N) - S)
    H1 = numpy.dot(P, numpy.dot(H, P))

    vals, vecs = numpy.linalg.eigh(H1)
    show_eigs(vals)
    idx = numpy.argmax(vals)
    evec = vecs[:, idx]
    print idx
    print evec
    #return

    pairs = []
    for i, v in enumerate(genidx((2,)*r)):
        v = array2(v)
        u = dot2(Bx.transpose(), v)

        u1 = (u+sx)%2
        v1 = dot2(Cx.transpose(), u1)
        j = eval('0b'+shortstr(v1, zero='0'))
        assert i!=j
        if evec[i]>evec[j] or (evec[i]==evec[j] and i<j):
            pairs.append((i, j)) # uniq

    print "pairs:", len(pairs)
    print "N:", N
    print pairs

    P = numpy.zeros((N, N/2))
    Q = numpy.zeros((N/2, N))
    for idx, pair in enumerate(pairs):
        i, j = pair
        P[i, idx] = 1
        P[j, idx] = -1
        Q[idx, i] = 1

    H = numpy.dot(Q, numpy.dot(H, P))
    #print lstr2(H, 0)

    evec = numpy.dot(Q, evec)
    print evec
    print "val:", (numpy.dot(evec, numpy.dot(H, evec))) / (numpy.dot(evec, evec))

    vals, vecs = numpy.linalg.eig(H)
    show_eigs(vals.real)
    #idx = numpy.argmax(vals)
    print vals.real
    #print "idx:", idx
    #v = vecs[:, idx]
    #print "v:", v 
    


def getnum(v):
    x = 0
    for i in v:
        if i:
            x += 1
        x *= 2
    return x/2

def getnum(v):
    s = ''.join(str(i) for i in v)
    return '0b'+s


def slepc(Gx, Gz, Hx, Hz, Rx, Rz, Pxt, Qx, Pz, Tx, **kw):

    name = argv.get("name", "ex3.tmp.c")
    print "slepc: name=%s"%name

    r = len(Rx)
    n = 2**r
    assert (r<40), "ugh"

    #code = Code("body.h")
    code = Code()

    code.append("#define DIMS (%d)"%n)

    code.append("static void matmult(PetscScalar *py, const PetscScalar *px, long nx)")
    code.begin()
    code.append("assert(DIMS == %d);"%n)
    code.append("assert(nx == %d);"%n)
    code.append("memset(py, 0, sizeof(PetscScalar)*nx);")

    offset = argv.get("offset", None)

    mz = len(Gz)
    t = None
    #excite = argv.excite
    #if excite is None:
    excite = kw.get("excite") or argv.excite

    if excite is not None:
        print "excite:", excite
        if type(excite) is tuple:
            t = Tx[excite[0]]
            for i in range(1, len(excite)):
                t = (t + Tx[excite[i]])%2
        else:
            assert type(excite) in (int, long)
            t = Tx[excite]
        print "t:", shortstr(t)
        Gzt = dot2(Gz, t)
        print "Gzt:", shortstr(Gzt)

    weights = kw.get("weights")
    if weights is not None:
        assert len(weights)==len(Gx)

    RR = dot2(Gz, Rx.transpose())

    PxtQx = dot2(Pxt, Qx)
    gxs = [getnum(dot2(gx, PxtQx)) for gx in Gx]
    gxs.sort()
    uniq_gxs = list(set(gxs))
    uniq_gxs.sort()

    code.append("long v;")
    code.append("int k;")
    code.append("struct timeval t0, t1;")
    code.append("gettimeofday(&t0, NULL);")
    code.append("for(v=0; v<%d; v++)"%n)
    code.begin()
    code.append("double pxv = px[v];")
    if n >= 128:
        code.append(r'if((v+1) %% %d == 0)' % (n//128))
        code.begin()
        code.append("gettimeofday(&t1, NULL);")
        code.append("long delta = t1.tv_sec-t0.tv_sec;")
        code.append("if(delta>1)")
        code.append('{printf("[%lds]", delta);fflush(stdout);}')
        code.append('t0 = t1;')
        code.end()
    code.append("k = 0;")
    for i, row in enumerate(RR):
        if t is not None and Gzt[i]==1:
            code.append("k += (countbits_fast(v&%s)+1) %% 2;" % getnum(row))
        else:
            code.append("k += countbits_fast(v&%s) %% 2;" % getnum(row))
    cutoff = argv.cutoff
    if cutoff is not None:
        code.append("if(k>%d) continue; // <-------- continue" % cutoff)
    else:
        code.append("if(k>cutoff) continue; // <-------- continue")
    code.append("py[v] += pxv * (%d - 2*k);" % mz)

    if weights is None:
        for gx in uniq_gxs:
            s = '+'.join(['pxv']*gxs.count(gx))
            code.append("py[v^%s] += %s;" % (gx, s))
    else:
        gxs = [getnum(dot2(gx, PxtQx)) for gx in Gx]
        for i, gx in enumerate(gxs):
            code.append("py[v^%s] += %s*pxv;" % (gx, weights[i]))

    code.end()
    code.end()

    if name is None:
        return

    s = code.output()

    src = open("ex3.c").read()
    match = '\n#include "body.h"\n'
    assert match in src
    src = src.replace(match, s)
    assert name and name.endswith(".c")
    f = open(name, 'w')
    tag = hash(src)
    print("hash(src):", tag)
    f.write(src)
    f.close()

    import socket
    host = socket.gethostname()
    if host == "bucket":
        cmd = "gcc MATCH.c -O3 -o MATCH -I/home/simon/local/petsc/arch-linux2-c-debug/include -I/home/simon/local/petsc/include/petsc/mpiuni -I/home/simon/local/petsc/include -I/home/simon/local/slepc-3.7.1/include -I/home/simon/local/slepc-3.7.1/arch-linux2-c-debug/include/ -L/home/simon/local/petsc/arch-linux2-c-debug/lib -L/home/simon/local/slepc-3.7.1/arch-linux2-c-debug/lib -lpetsc -lslepc"
    elif host == "hero":
        cmd = "gcc MATCH.c -O3 -o MATCH -I/usr/include/openmpi -I/usr/include/petsc -I/usr/include/slepc -lpetsc -lslepc -lmpi"
    else:
        cmd = "gcc -O3 MATCH.c -I/suphys/sburton/include/ -o MATCH -lpetsc -L$PETSC_DIR/$PETSC_ARCH/lib -L$SLEPC_DIR/$PETSC_ARCH/lib -lslepc"

    cmd = cmd.replace("MATCH.c", name)
    stem = name[:-2]
    cmd = cmd.replace("MATCH", stem)

    rval = os.system(cmd)
    assert rval == 0
    #print("exec:", hash(open(stem).read()))

    nev = argv.get("nev", 1)
    cmd = "./%s -eps_nev %d -eps_ncv %d -eps_largest_real" 

    if argv.plot:
        cmd += " -eps_view_vectors binary:evec.bin "

    cmd = cmd%(stem, nev, max(2*nev, 4))
    eps_tol = argv.get("eps_tol", 1e-4)
    if eps_tol is not None:
        cmd += " -eps_tol %s "%eps_tol

    #cmd += " -eps_type arnoldi -info -eps_monitor -eps_tol 1e-3"
    print cmd
    #rval = os.system(cmd)
    #assert rval == 0
    f = os.popen(cmd)
    s = f.read()
    #print(s)
    lines = s.split('\n')
    vals = []
    for line in lines:
        line = line.strip()
        flds = line.split()
        #print("parse", flds)
        try:
            a, b = flds
            a = float(a)
            b = float(b)
            vals.append(a)
        except:
            pass

    if not argv.plot:
        print("vals:", vals)
        return vals

    assert argv.plot.endswith(".pdf")

    s = open("evec.bin").read()
    sz = 8*2**r

    if len(s)==sz+8:
        s = s[8:]
    elif len(s)==sz+16:
        s = s[16:]
    #elif len(s)==2*(sz+16): # got two vectors
    #    s = s[16:16+sz] # pick the first vector
    elif len(s)%(sz+16) == 0:
        count = len(s)/(sz+16)
#        s = s[16:16+sz] # pick the first vector
        ev_idx = argv.get("ev_idx", 0)
        s = s[16+ev_idx*(16+sz):(ev_idx+1)*(16+sz)]
    else:
        assert 0, "sz=%d but s=%s"%(sz, len(s))

    vec0 = numpy.fromstring(s, dtype=">d")
    assert len(vec0)==2**r

    assert excite is None

    print "graphing..."
    gz, n = Gz.shape
    xdata = []
    lookup = {}
    GzRxt = dot2(Gz, Rx.transpose())
    for i, v in enumerate(genidx((2,)*r)):
        v = array2(v)
        lookup[v.tostring()] = i
        syndrome = dot2(GzRxt, v)
        value = gz - 2*syndrome.sum()
        xdata.append(value)

    pdata = {}
    ndata = {}
    my = 20. # mul y
    EPSILON = argv.get("EPSILON", 1e-6)

    def yfunc(y):
        y = log2(abs(y))
        y = int(round(my*y))
        return y

    for i in range(len(vec0)):
        x = xdata[i] # integer
        y = vec0[i]
        if abs(y) < EPSILON:
            continue
        if y > 0.:
            y = yfunc(y)
            pdata[x, y] = pdata.get((x, y), 0) + 1
        else:
            y = yfunc(y)
            ndata[x, y] = ndata.get((x, y), 0) + 1

    from pyx import graph, canvas, path, trafo, color, deco, text
    
    north = [text.halign.boxcenter, text.valign.top]
    northeast = [text.halign.boxright, text.valign.top]
    northwest = [text.halign.boxleft, text.valign.top]
    south = [text.halign.boxcenter, text.valign.bottom]
    southeast = [text.halign.boxright, text.valign.bottom]
    southwest = [text.halign.boxleft, text.valign.bottom]
    east = [text.halign.boxright, text.valign.middle]
    west = [text.halign.boxleft, text.valign.middle]
    center = [text.halign.boxcenter, text.valign.middle]
    
    c = canvas.canvas()
    sx = 0.4
    sy = 1.4
    tr = trafo.scale(sx, sy)

    green = color.rgb(0.2,0.6,0.2)
    brown = color.rgb(0.8,0.2,0.2)
    grey = color.rgb(0.4,0.4,0.4)
    lgrey = color.rgb(0.8,0.8,0.8)

    W = 2*gz
    H = log2(EPSILON)
    dy = 0.5 * 1.2/my

    X0 = -gz
    Y0 = 0.

    def showp(gx, radius):
        v = dot2(gx, PxtQx)
        syndrome = dot2(GzRxt, v)
        x = gz - 2*syndrome.sum()
        i = lookup[v.tostring()]
        #print syndrome, syndrome.sum(), vec0[i]
        y = (1./my)*yfunc(vec0[i]) + 0.5*dy
        #y = 0.5*dy + log2(abs(vec0[i]))
        c.fill(path.circle(-x*sx, y*sy, radius), [lgrey])

    showp(zeros2(n), 0.8)
    for gx in Gx:
        showp(gx, 0.4)

    for gx0 in Gx:
      for gx1 in Gx:
        gx = (gx0+gx1)%2
        if gx.sum()==0:
            continue
        showp(gx, 0.2)

    for i in range(0, gz+1):
        x, y = X0+2*i, Y0
        c.stroke(path.line(x, y, x, y+H), [tr, grey])
        if i%2 == 0:
            c.text(x*sx, y*sy + 0.2, "%d"%i, south)

    c.stroke(path.line(X0, Y0, X0+1.0*W+3.5, Y0), [tr, deco.earrow(size=0.5)])
    c.stroke(path.line(X0, Y0, X0, Y0+1.0*H-0.5), [tr, deco.earrow(size=0.5)])

    y = 1.0
    i = 0
    while y > EPSILON:

        x = X0*sx
        y1 = sy*(1./my)*yfunc(y)
        c.stroke(path.line(x, y1, x-0.1, y1))

        c.text(x-0.3, y1, "%d"%i, east)

        y /= 2.
        i -= 1

    R = 0.15
    for key, value in pdata.items():
        x, y = key
        y = y/my
        x = -x
        value = 1 + math.log(value)
        r = R*value
        #c.stroke(path.circle(x, y, r))
        #c.stroke(path.line(x, y, x+r, y), [brown, tr])
        c.fill(path.rect(x, y, r, dy), [brown, tr])

    for key, value in ndata.items():
        x, y = key
        y = y/my
        x = -x
        value = 1 + math.log(value)
        r = R*value
        #c.stroke(path.circle(x, y, r))
        #c.stroke(path.line(x-r, y, x, y), [green, tr])
        c.fill(path.rect(x-r, y, r, dy), [green, tr])

    c.writePDFfile(argv.plot)


    if 0:
        print "graph.."
    
        g = graph.graphxy(
            width=16,
            x=graph.axis.linear(reverse=True),
            y=graph.axis.linear())
            #y=graph.axis.log(min=0.8*vec0.min(), max=1.2*vec0.max()))
    
        g.plot(graph.data.values(x=xdata, y=ydata))
        g.writePDFfile(argv.plot)



def get_perm(fn):
    n = len(fn)
    perm = {}
    for i in range(n):
        if i != fn[i]:
            perm[i] = fn[i]
    return perm


def do_symmetry(Gx, Gz, Hx, Hz):
    "find permutation symmetries of the code"

    bag0 = Tanner.build(Gx, Gz)
    bag1 = Tanner.build(Gx, Gz)

    #n = Gx.shape[1]
    #for i in range(n):
    #    print Gx[:, i].sum(),
    #print

    rows = [shortstr(h) for h in Hx]
    bits = [p for p in bag0 if p.desc=='b']
    #print bits

    count = 0
    perms = []
    #keys = range(m, len(bag0))
    #print "keys:", keys
    for fn in search(bag0, bag1):
        write('.')
        perm = [bag1[fn[p.idx]].row for p in bits]
        #print perm

        # this is the action of graph automorphism on the stabilizers
#        rows1 = [rows.index(shortstr(h[perm])) for h in Hx]
#        print rows1 
#        perms.append(rows1)

        count += 1
    print
    print "isomorphisms:", count

    print "stabilizer orbits:",
    marked = {}
    for i in range(len(Hx)):
        if i in marked:
            continue
        print i,
        marked[i] = True
        for perm in perms:
            marked[perm[i]] = True
    print


def do_chainmap(Gx, Gz):

    print "do_chainmap"

    gx, n = Gx.shape
    gz, n = Gz.shape

    Qx = Unknown(gx, gx)
    Qz = Unknown(gz, gz)
    P = Unknown(n, n)
    
    system = System(Qx, P, Qz)

    system.append(dot2(Qz, Gz), dot2(Gz, P))
    system.append(dot2(Qx, Gx), dot2(Gx, P))

    kern = system.kernel()
    print "kernel:", len(kern)

    T = system.solve()
    Qx, P, Qz = T

    # zero !
    #print shortstrx(Qx, P, Qz)


def find_ideals():
    import networkx as nx
    from models import build_model

    model = build_model()
    if argv.verbose:
        print model
        print
    #print(shortstrx(Hx, Hz))

    Gx, Gz = model.Gx, model.Gz
    Hx, Hz = model.Hx, model.Hz
    print shortstrx(Gx, Gz)

    Rx, Rz = model.Rx, model.Rz
    Rxt = Rx.transpose()
    Rzt = Rz.transpose()

    Px, Pz = model.Px, model.Pz
    Pxt = Px.transpose()
    Pzt = Pz.transpose()

    r, n = Rx.shape

    assert eq2(dot2(Rx, Pxt), Rx)

#    print shortstrx(Gx, Gz)
#    print

    PGx = dot2(Gx, Pxt)
    PGz = dot2(Gz, Pzt)
    
#    print shortstrx(PGx, PGz)
#    print

    if 0:
        ng = len(PGx)

        graph = nx.Graph()
        for i in range(ng):
            graph.add_node(i)
    
        for i, gi in enumerate(PGx):
            for j, gj in enumerate(PGz):
                if (gi*gj).sum()%2:
                    graph.add_edge(i, j)

    mx = len(Gx)
    mz = len(Gz)
    graph = nx.Graph()
    for i in range(mx+mz):
        graph.add_node(i)

    for ix in range(mx):
        for iz in range(mz):
            if (Gx[ix]*Gz[iz]).sum()%2:
                graph.add_edge(ix, mx+iz)

#    excite = argv.excite
#    if type(excite) is int: 
#        _excite = [0]*len(model.Tx)
#        _excite[excite] = 1
#        excite = tuple(_excite)
#    #elif excite is None:
#    print("excite:", excite)

    equs = nx.connected_components(graph)
    equs = list(equs)

    print(model)
    print("find_ideals:", len(equs))

    models = []
    for equ in equs:
        equ = list(equ)
        equ.sort()
        #print(equ)
        gxs = []
        gzs = []
        for i in equ:
            if i<mx:
                gxs.append(Gx[i])
            else:
                gzs.append(Gz[i-mx])
        _model = build_model(array2(gxs), array2(gzs))
        print(_model)
        models.append(_model)

    if not argv.solve:
        return

    if argv.excite:
        if argv.minweightall:
            excites = minweightall(Hz)
            #print("excite", excite)
        elif argv.exciteall:
            excites = list(enum2(len(Hz)))[1:]
        else:
            assert len(Hz), Hz
            excites = []
            for i in range(len(Hz)):
                excite = array2([0]*len(Hz))
                excite[i] = 1
                excites.append(excite)
        print("excites", (excites))
    else:
        excites = [None]

    if eq2(model.Hx, model.Hz):
        print "self dual stabilizers"

    _excite = None
    top = None
    for excite in excites:
        #if excite is not None:
        #    print "excite:", (excite)
        #    g = dot2(model.Hx.transpose(), excite)
        #    model.show_stabx(g)
        if excite is not None:
            tx = dot2(excite, model.Tx)

        total = 0.
        gaps = []

        for _model in models:
            #print(_model)
            if excite is not None:
                _tx = dot2(tx, _model.Hz.transpose(), _model.Tx)
                _excite = dot2(_tx, _model.Hz.transpose())
                #print(shortstrx(_excite))
                #print _excite
    
            r = len(_model.Rx)
            if r <= 12 and not argv.slepc and not argv.sparse:

                H = _model.build_ham(_excite) #weights=weights1)
                if argv.orbigraph:
                    H1 = build_orbigraph(H)
                
                vals, vecs = numpy.linalg.eigh(H)

            elif argv.sparse and r <= 19:

                vals = _model.sparse_ham_eigs(_excite) #weights=weights1)

            elif r <= 27:
                #for i in range(len(_excite)):
                #    _excite[i] = 1
                #print("excite:", tuple(_excite))
                #vals = slepc(excite=tuple(_excite), **_model.__dict__)
                if _excite is not None:
                    _excite = tuple(_excite)
                #vals = slepc(excite=_excite, **_model.__dict__)
                vals = _model.do_slepc(excite=_excite)
                #vals = [0, 0]

            else:
                assert 0, "r=%d too big"%r
            vals = list(vals)
            vals.sort(reverse=True)

            total += vals[0]
            gaps.append(vals[0]-vals[1])
    
        if top is None or total > top:
            top = total
        print("eval_1:", total)
        print("eval_2:", total-min(gaps))
    print("top:", top)


def main():

    import models

    assert not argv.orbiham, "it's called orbigraph now"

    if argv.find_ideals:
        find_ideals()
        return

    Gx, Gz, Hx, Hz = models.build()

    if argv.chainmap:
        do_chainmap(Gx, Gz)

    if argv.symmetry:
        do_symmetry(Gx, Gz, Hx, Hz)
        return
    
    #print shortstrx(Gx, Gz)
    if argv.report:
        print "Hz:"
        for i, h in enumerate(Hz):
            print i, shortstr(h), h.sum()
    #print shortstr(find_stabilizers(Gx, Gz))

    Lz = find_logops(Gx, Hz)
    Lx = find_logops(Gz, Hx)
    #print "Lz:", shortstr(Lz)

    if Lz.shape[0]*Lz.shape[1]:
        print Lz.shape, Gx.shape
        check_commute(Lz, Gx)
        check_commute(Lz, Hx)

    Px = get_reductor(Hx) # projector onto complement of rowspan of Hx
    Pz = get_reductor(Hz) 

    Rz = [dot2(Pz, g) for g in Gz]
    Rz = array2(Rz)
    Rz = row_reduce(Rz, truncate=True)
    rz = len(Rz)

    n = Gx.shape[1]
    print "n =", n
    if len(Lx):
        print "Lx Lz:"
        print shortstrx(Lx, Lz)
    print "Hx:", len(Hx), "Hz:", len(Hz)
    print "Gx:", len(Gx), "Gz:", len(Gz)

    Rx = [dot2(Px, g) for g in Gx]
    Rx = array2(Rx)

    Rx = row_reduce(Rx, truncate=True)
    rx = len(Rx)
    print "Rx:", rx, "Rz:", rz
    if argv.show:
        print shortstrx(Rx, Rz)

    Qx = u_inverse(Rx)
    Pxt = Px.transpose()
    assert eq2(dot2(Rx, Qx), identity2(rx))
    assert eq2(dot2(Rx, Pxt), Rx)

    #print shortstr(dot2(Pxt, Qx))
    PxtQx = dot2(Pxt, Qx)
    lines = [shortstr(dot2(g, PxtQx)) for g in Gx]
    lines.sort()
    #print "PxtQx:"
    #for s in lines:
    #    print s
    #print "RzRxt"
    #print shortstr(dot2(Rz, Rx.transpose()))

    offset = argv.offset

    if len(Hz):
        Tx = find_errors(Hz, Lz, Rz)
    else:
        Tx = zeros2(0, n)

    if argv.dense:
        dense(**locals())
        return

    if argv.dense_full:
        dense_full(**locals())
        return

    if argv.show_delta:
        show_delta(**locals())
        return

    if argv.slepc:
        slepc(**locals())
        return

#    if argv.orbigraph:
#        from linear import orbigraph
#        orbigraph(**locals())
#        return

    v0 = None

#    excite = argv.excite
#    if excite is not None:
#        v0 = zeros2(n)
#        v0[excite] = 1

    verts = []
    lookup = {}
    for i, v in enumerate(span(Rx)): # XXX does not scale well
        if v0 is not None:
            v = (v+v0)%2
            v = dot2(Px, v)
        lookup[v.tostring()] = i
        verts.append(v)
    print "span:", len(verts)
    assert len(lookup) == len(verts)

    mz = len(Gz)
    n = len(verts)

    if argv.lie:
        U = []
        for i, v in enumerate(verts):
            count = dot2(Gz, v).sum()
            Pxv = dot2(Px, v)
            assert count == dot2(Gz, Pxv).sum()
            U.append(mz - 2*count)
        uniq = list(set(U))
        uniq.sort(reverse=True)
        s = ', '.join("%d(%d)"%(val, U.count(val)) for val in uniq)
        print s
        print "sum:", sum(U)
        return
        

    if n <= 1024 and argv.solve:
        H = numpy.zeros((n, n))
        syndromes = []
        for i, v in enumerate(verts):
            syndromes.append(dot2(Gz, v))
            count = dot2(Gz, v).sum()
            Pxv = dot2(Px, v)
            assert count == dot2(Gz, Pxv).sum()
            H[i, i] = mz - 2*count
            for g in Gx:
                v1 = (g+v)%2
                v1 = dot2(Px, v1)
                j = lookup[v1.tostring()]
                H[i, j] += 1
    
        if argv.showham:
            s = lstr2(H, 0).replace(',  ', ' ')
            s = s.replace(' 0', ' .')
            s = s.replace(', -', '-')
            print s
    
        vals, vecs = numpy.linalg.eigh(H)
        show_eigs(vals)

        if argv.orbigraph:
            if argv.symplectic:
                H1 = build_orbigraph(H, syndromes)
            else:
                H1 = build_orbigraph(H)
            print "orbigraph:"
            print H1
            vals, vecs = numpy.linalg.eig(H1)
            show_eigs(vals)

    elif argv.sparse:
        print "building H",
        A = {} # adjacency
        U = [] # potential

        if offset is None:
            offset = mz + 1 # make H positive definite

        for i, v in enumerate(verts):
            if i%1000==0:
                write('.')
            count = dot2(Gz, v).sum()
            #H[i, i] = mz - 2*count
            U.append(offset + mz - 2*count)
            for g in Gx:
                v1 = (g+v)%2
                v1 = dot2(Px, v1)
                j = lookup[v1.tostring()]
                A[i, j] = A.get((i, j), 0) + 1
    
        print "\nnnz:", len(A)

        if argv.lanczos:
            vals, vecs = do_lanczos(A, U)

        elif argv.orbigraph:
            vals, vecs = do_orbigraph(A, U)

        else:
            return

        vals -= offset # offset doesn't change vecs

        show_eigs(vals)

    elif argv.orbigraph:

        assert n<=1024

        H = numpy.zeros((n, n))
        syndromes = []
        for i, v in enumerate(verts):
            syndromes.append(dot2(Gz, v))
            count = dot2(Gz, v).sum()
            Pxv = dot2(Px, v)
            assert count == dot2(Gz, Pxv).sum()
            H[i, i] = mz - 2*count
            for g in Gx:
                v1 = (g+v)%2
                v1 = dot2(Px, v1)
                j = lookup[v1.tostring()]
                H[i, j] += 1
    
        if argv.showham:
            s = lstr2(H, 0).replace(',  ', ' ')
            s = s.replace(' 0', ' .')
            s = s.replace(', -', '-')
            print s
    
        if argv.symplectic:
            H1 = build_orbigraph(H, syndromes)
        else:
            H1 = build_orbigraph(H)
        #print "orbigraph:"
        #print H1


from argv import Argv
argv = Argv()

if __name__ == "__main__":

    main()





