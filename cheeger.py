#!/usr/bin/env python

import sys, os

import numpy
from scipy import sparse

from solve import shortstr, shortstrx, parse, eq2, dot2, zeros2, array2, identity2
from solve import row_reduce, RowReduction, span, get_reductor
from solve import u_inverse, find_logops, solve, find_kernel, linear_independent

import models
from lanczos import write
from code import lstr2
from zorbit import get_reductor, find_errors, find_logops, check_commute, check_conjugate
from zorbit import genidx

from argv import Argv
argv = Argv()


def mksparse(H, N):
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
    H1 = sparse.coo_matrix((data, (rows, cols)), (N, N))
    H1 = sparse.csr_matrix(H1, dtype=numpy.float64)

    return H1


def xprojector(Bx, Cx, opx, sign=+1):
    r = len(Bx)
    N = 2**r
    P = {}
    for i in range(N):
        P[i, i] = 0.5
    for i, v in enumerate(genidx((2,)*r)):
        v = array2(v)
        u = dot2(Bx.transpose(), v)
        u1 = (u+opx)%2
        v1 = dot2(Cx.transpose(), u1)
        j = eval('0b'+shortstr(v1, zero='0'))
        assert i!=j, opx
        P[i, j] = +0.5*sign
    P = mksparse(P, N)
    return P


def mkxop(Bx, Cx, opx):
    r = len(Bx)
    N = 2**r
    S = {}
    for i, v in enumerate(genidx((2,)*r)):
        v = array2(v)
        u = dot2(Bx.transpose(), v)
        u1 = (u+opx)%2
        v1 = dot2(Cx.transpose(), u1)
        j = eval('0b'+shortstr(v1, zero='0'))
        S[i, j] = 1.0
    S = mksparse(S, N)
    return S


def mkzop(Bx, Cx, opz):
    r = len(Bx)
    N = 2**r
    S = {}
    for i, v in enumerate(genidx((2,)*r)):
        v = array2(v)
        u = dot2(Bx.transpose(), v)
        check = dot2(u, opz)
        assert check in (0, 1)
        S[i, i] = 1 - 2*check
    S = mksparse(S, N)
    return S


class Graph(object):
    def __init__(self, verts, edges):
        self.verts = set(verts)
        self.N = len(self.verts)
        self.edges = set(edges) # directed edges!
        nbd = dict((i, []) for i in self.verts) # map vert index -> list of neighbours
        for (i, j) in edges:
            nbd[i].append(j)
        self.nbd = nbd

    def __str__(self):
        return "Graph(%s, %s)"%(len(self.verts), len(self.edges))

    def get_dist(self, i, N):
        nbd = self.nbd
        dist = [None]*N
        bdy = [i]
        dist[i] = 0
        while bdy:
            _bdy = []
            for i in bdy:
                d = dist[i]
                assert d is not None
                for j in nbd[i]:
                    if dist[j] == None:
                        dist[j] = d+1
                        _bdy.append(j)
            bdy = _bdy
    
        assert None not in dist
        return dist


def mkop(s, lookup, Bx, r):
    op = {}
    for i, v in enumerate(genidx((2,)*r)):
        u = dot2(Bx.transpose(), v)
        u = (s+u)%2
        op[i] = lookup[u.tostring()]
    return op


def do_cut(Sx, Gx, Bx, Cx, A, r, N, lookup, **kw):
    print "do_cut"

    n = Bx.shape[1]

    assert len(A.keys()) == (2**r) * len(Gx)

    graph = Graph(range(N), A.keys())
    print graph

    u = zeros2(n)

    #s = Sx[0]

    paths = []
    edges = []
    for s in Sx:
    
        word = [mkop(g, lookup, Bx, r) for g in Gx if (s*g).sum()==2]
        print "word:", len(word)
    
        for i in graph.verts:
            path = []
            for op in word:
                j = op[i]
                edge = (i, j)
                path.append(edge)
                edges.append(edge)
                i = j
            #print path
            paths.append(path)
    
    # paths do not overlap
    assert len(edges)==len(set(edges)) # uniq

    print "path length:", len(paths[0])
    print "paths:", len(paths)
    print "paths/verts:", 1.0*len(paths)/N

    # path length == li
    # num paths == (lj-1)*N
    # paths per basis: (lj-1)

    return

    #######################################
    # Try to partition the basis vectors:

    i0 = lookup[u.tostring()]
    i1 = lookup[s.tostring()]
    dist0 = graph.get_dist(i0, N)
    dist1 = graph.get_dist(i1, N)

    count = 0
    for i, j in graph.edges:
        if dist0[i] == dist0[j] and dist1[i]==dist1[j]:
            count += 1
    print "indistinguishable edges:", count

    #print dist0
    #print dist1
    A = set([i for i in range(N) if dist0[i] <= dist1[i]])
    B = set([i for i in range(N) if dist0[i] >= dist1[i]])
    print len(A), len(B), N

    cut = []
    for i, j in graph.edges:
        if i in A and j in B:
            cut.append((i, j))
    print "cut:", len(cut)
    cut = set(cut)
    for i, j in cut:
        assert j, i in cut

    cut0 = cut
    for s in Sx[1:]:
        op = mkop(s, lookup, Bx, r)
        cut1 = set((op[i], op[j]) for (i, j) in cut)
        cut0 = cut0.intersection(cut1)
        print "intersection:", len(cut0)




def main():

    Gx, Gz, Sx, Sz = models.build()

    mx = len(Sx)

    Lz = find_logops(Gx, Sz)
    check_commute(Lz, Gx)
    check_commute(Lz, Sx)

    Lx = find_logops(Gz, Sx)
    check_commute(Lx, Gz)
    check_commute(Lx, Sz)

    check_conjugate(Lx, Lz)

    Px = get_reductor(Sx) # projector onto complement of rowspan of Sx
    Pz = get_reductor(Sz) 

    Rz = [dot2(Pz, g) for g in Gz]
    Rz = array2(Rz)
    Rz = row_reduce(Rz, truncate=True)
    rz = len(Rz)
    print "Rz:", rz

    Rx = [dot2(Px, g) for g in Gx]
    Rx = array2(Rx)
    Rx = row_reduce(Rx, truncate=True)
    rx = len(Rx)
    print "Rx:", rx

    Tx = find_errors(Sz, Lz, Rz)
    check_commute(Tx, Lz)
    check_commute(Tx, Rz)
    Tz = find_errors(Sx, Lx, Rx)
    check_commute(Tz, Lx)
    check_commute(Tz, Rx)

    #print Gx.shape
    #print shortstr(Gx)

    gz, n = Gz.shape

    excite = argv.excite

    print "excite:", excite

    t = zeros2(n)
    if excite is not None:
        assert len(excite)==len(Tx)
        for i, ex in enumerate(excite):
            if ex:
                t = (t + Tx[i])%2
        #print "t:", shortstr(t)

    Gzt = dot2(Gz, t)
    #print "Gzt:", shortstr(Gzt)

    # This is our basis
    Bx = array2([v+t for v in Rx] + [v+t for v in Sx])
    Bx %= 2
    r = len(Bx)
    N = 2**r
    Bx = row_reduce(Bx, truncate=True)
    assert len(Bx)==r # linearly independent rows
    Cx = u_inverse(Bx)

    H = {}
    A = {}
    U = []

    print "build U"
    pos = neg = 0
    lookup = {}
    for i, v in enumerate(genidx((2,)*r)):
        v = array2(v)
        u = dot2(Bx.transpose(), v)
        lookup[u.tostring()] = i
        syndrome = dot2(Gz, u)
        value = gz - 2*syndrome.sum()
        #print shortstr(dot2(Rx.transpose(), v)), value
        H[i, i] = value
        U.append(value)

    print "build A"
    for i, v in enumerate(genidx((2,)*r)):
        v = array2(v)
        u0 = dot2(Bx.transpose(), v)
        #print shortstr(v),
        for g in Gx:
            u1 = (u0 + g) % 2
            #v1 = dot2(Cx.transpose(), u1)
            #assert v1.shape == v.shape
            #j = eval('0b'+shortstr(v1, zero='0'))
            j = lookup[u1.tostring()]
            A[i, j] = A.get((i, j), 0) + 1
            H[i, j] = H.get((i, j), 0) + 1

    #print A

    if argv.cut:

        do_cut(**locals())
        return

    H1 = mksparse(H, N)
    vals, vecs = sparse.linalg.eigsh(H1, k=40, which="LA")
    idx = numpy.argmax(vals)
    print "eigval:", vals[idx]
    vec = vecs[:, idx]
    if vec.min() < -1e-4:
        vec *= -1
    vec *= 1./numpy.linalg.norm(vec)
    print vec

    def show(H, v, msg=""):
        Hv = H.dot(v)
        print msg,
        print v.dot(Hv),
        print v.dot(Hv)/numpy.linalg.norm(Hv)

    print shortstr(Sx)

    def project(H, syndrome):
        for i in range(mx):
            #sign = syndrome[i]
            sign = -1 if i in syndrome else 1
            P = xprojector(Bx, Cx, Sx[i], sign)
            H = P.dot(H.dot(P))
        return H

    PHP0 = project(H1, [0])
    PHP01 = project(H1, [0, 2])

    vals, vecs = sparse.linalg.eigsh(PHP0, k=40, which="LA")
    idx = numpy.argmax(vals)
    print "eigval:", vals[idx]
    vec0 = vecs[:, idx]

    vals, vecs = sparse.linalg.eigsh(PHP01, k=40, which="LA")
    idx = numpy.argmax(vals)
    print "eigval:", vals[idx]
    vec01 = vecs[:, idx]

    cut0 = show_cut(A, vec, vec0)
    cut01 = show_cut(A, vec, vec01)

    print "intersection:", len(cut0.intersection(cut01))



def show_cut(A, vec0, vec):
    N = len(vec)
    weight = 0.
    cut = set()
    for i, j in A:
        assert i!=j
        assert (j,i) in A

        vivj = vec[i] * vec[j]
        if vivj==0.:
            fail
        if vivj<0: #-1e-10:
            cut.add((i, j))
            w = vec0[i] * vec0[j] * A[i, j]
            #assert w>1e-8
            if w<-1e-10:
                print w,
            weight += abs(w)
    print
    print "cut size:", len(cut)
    print "cut / edges:", 1.*len(cut)/len(A)
    print "cut / basis:", 1.*len(cut) / N
    #print 1.*len(A)/len(cut)
    print "weight:", weight
    for i, j in cut:
        assert (j, i) in cut
    return cut







if __name__=="__main__":

    if argv.profile:
        import cProfile as profile
        profile.run("main()")
    else:

        main()

