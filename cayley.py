#!/usr/bin/env python

import sys, os

import numpy
from scipy import sparse

from solve import shortstr, shortstrx, parse, eq2, dot2, zeros2, array2, identity2
from solve import row_reduce, RowReduction, span, get_reductor
from solve import u_inverse, find_logops, solve, find_kernel, linear_independant

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


def mkprojector(Bx, Cx, opx, sign=+1):
    r = len(Bx)
    N = 2**r
    S = {}
    P = {}
    for i in range(N):
        P[i, i] = 0.5
    for i, v in enumerate(genidx((2,)*r)):
        v = array2(v)
        u = dot2(Bx.transpose(), v)

        u1 = (u+opx)%2
        v1 = dot2(Cx.transpose(), u1)
        j = eval('0b'+shortstr(v1, zero='0'))
        S[i, j] = 1.0
        assert i!=j, opx
        P[i, j] = +0.5*sign
    S = mksparse(S, N)
    P = mksparse(P, N)
    return P


def main():

    Gx, Gz, Hx, Hz = models.build()

    mx = len(Hx)

    Lz = find_logops(Gx, Hz)
    check_commute(Lz, Gx)
    check_commute(Lz, Hx)

    Lx = find_logops(Gz, Hx)
    check_commute(Lx, Gz)
    check_commute(Lx, Hz)

    check_conjugate(Lx, Lz)

    Px = get_reductor(Hx) # projector onto complement of rowspan of Hx
    Pz = get_reductor(Hz) 

    Rz = [dot2(Pz, g) for g in Gz]
    Rz = array2(Rz)
    Rz = row_reduce(Rz, truncate=True)
    rz = len(Rz)
    print "Rz:", rz

    Tx = find_errors(Hz, Lz, Rz)

    Rx = [dot2(Px, g) for g in Gx]
    Rx = array2(Rx)

    Rx = row_reduce(Rx, truncate=True)
    rx = len(Rx)
    print "Rx:", rx

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
    Bx = array2([v+t for v in Rx] + [v+t for v in Hx])
    Bx %= 2
    r = len(Bx)
    N = 2**r
    Bx = row_reduce(Bx, truncate=True)
    assert len(Bx)==r # linearly independant rows
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

    H1 = mksparse(H, N)
    vals, vecs = sparse.linalg.eigsh(H1, k=40, which="LA")
    idx = numpy.argmax(vals)
    vec0 = vecs[:, idx]
    if vec0.min() < -1e-4:
        vec0 *= -1
    vec0 *= 1./numpy.linalg.norm(vec0)
    print vec0

    #print "Hx:"
    #print shortstr(Hx)
    for i in range(mx):
        sign = -1 if i==0 else 1
        P = mkprojector(Bx, Cx, Hx[i], sign)
        H1 = P.dot(H1.dot(P))

    #P = mkprojector(Bx, Cx, Lx[0], +1) # no no !!! we already killed Lx
    #H1 = P.dot(H1.dot(P))

    print "eigsh"
    vals, vecs = sparse.linalg.eigsh(H1, k=min(len(U)-5, 40), which="LA")

    #print list(vals)
    idx = numpy.argmax(vals)
    print "eval:", vals[idx]

    vec = vecs[:, idx]

    print vec

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
            assert w>1e-8
            weight += w
    print "cut:", len(cut), len(A),
    print 1.*len(A)/len(cut)
    print "weight:", weight
    for i, j in cut:
        assert (j, i) in cut







if __name__=="__main__":

    if argv.profile:
        import cProfile as profile
        profile.run("main()")
    else:

        main()

