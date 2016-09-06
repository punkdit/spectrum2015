#!/usr/bin/env python

import sys, os
from fractions import gcd

import numpy
from numpy import concatenate

import networkx as nx

from solve import get_reductor, array2, row_reduce, dot2, shortstr, zeros2, shortstrx, eq2
from solve import u_inverse, find_kernel, find_logops, identity2, solve
from isomorph import write
from zorbit import check_conjugate, check_commute
from zorbit import find_errors

import models
from models import genidx

"""
Find the closure of the gauge operators under bracket.
"""

def bracket(a, b):
    a, b = a.copy(), b.copy()
    assert len(a)==len(b)
    assert len(a)%2==0
    n = len(a)//2
    a.shape = (n, 2)
    b.shape = (n, 2)
    #print a
    #print b
    b = b[:, [1, 0]]
    #print b
    c = a*b
    return c.sum() % 2


def is_zop(a):
    a = a.copy()
    a.shape = (len(a)//2, 2)
    assert a.sum()
    return a[:, 0].sum() == 0


def mkop(xop, zop):
    if xop is None:
        xop = zeros2(len(zop))
    if zop is None:
        zop = zeros2(len(xop))
    xop, zop = xop.view(), zop.view()
    n = len(xop)
    xop.shape = 1,n
    zop.shape = 1,n
    op = numpy.concatenate((xop, zop))
    op = op.transpose().copy()
    op.shape = (2*n,)
    return op


def closure(ops):
    "take closure under bracket operation (up to factor of 2)"
    #found = set(op.tostring() for op in ops)
    found = set(op.tostring() for op in ops)
    pairs = [(a, b) for a in ops for b in ops if a.tostring()!=b.tostring()]
    while 1:
        new = []
        for a, b in pairs:
            if bracket(a, b):
                c = (a+b)%2
                s = c.tostring()
                if s not in found:
                    found.add(s)
                    new.append(c)
        if not new:
            break
        pairs = [(a, b) for a in ops+new for b in new if a.tostring()!=b.tostring()]
        ops.extend(new)
        if argv.verbose:
            print len(ops)
    return ops


def check_sy(Lx, Hx, Tx, Rx, Lz, Hz, Tz, Rz, **kw):

    check_conjugate(Lx, Lz)
    check_commute  (Lx, Hz)
    check_commute  (Lx, Tz)
    check_commute  (Lx, Rz)

    check_commute  (Hx, Lz)
    check_conjugate(Hx, Tz)
    check_commute  (Hx, Hz)
    check_commute  (Hx, Rz)

    check_commute  (Tx, Lz)
    check_commute  (Tx, Tz)
    check_conjugate(Tx, Hz)
    check_commute  (Tx, Rz)

    check_commute  (Rx, Lz)
    check_commute  (Rx, Hz)
    check_commute  (Rx, Tz)
    check_conjugate(Rx, Rz)


def find_errors(Hx, Lx):
    "find inverse of Hx commuting with Lx"
    
    # find Tz
    n = Hx.shape[1]

    #print "find_errors:"
    #print shortstrx(Hx, Lx)
    #print
    Lx = row_reduce(Lx)
    k = len(Lx)
    mx = len(Hx)

    HL = row_reduce(concatenate((Lx, Hx)))
    #print shortstrx(Lx, Hx, HL)
    assert len(HL) == mx+k
    assert k+mx <= n, (k, mx, n)

    U = zeros2(mx+k, n)
    U[:mx] = Hx 
    U[mx:mx+k] = Lx 

    B = zeros2(mx+k, mx)
    B[:mx] = identity2(mx)

    #print shortstrx(U, B)
    #print

    Tz_t = solve(U, B)
    assert Tz_t is not None, "no solution"
    Tz = Tz_t.transpose()
    assert len(Tz) == mx

    check_conjugate(Hx, Tz)
    check_commute(Lx, Tz)

    return Tz


def test_model():

    Gx, Gz, Hx, Hz = models.build()
    n = Hx.shape[1]

    check_commute(Hx, Hz)
    check_commute(Gx, Hz)
    check_commute(Hx, Gz)

    #Px = get_reductor(concatenate((Lx, Hx)))
    #Pz = get_reductor(concatenate((Lz, Hz)))
    Px = get_reductor(Hx)
    Pz = get_reductor(Hz)

    # Lz = find_logops( Hx            , Hz            )
    #      find_logops( ............. , ............. )
    #                 ( commutes with , orthogonal to )
    #                 ( ............. , ............. )

    Lz = find_logops(Gx, Hz)

    if 0:
        PGz = get_reductor(Gz)
        Lz = dot2(Lz, PGz.transpose())
        Lz = row_reduce(Lz)
    
        print shortstrx(Lz, Gz, Hz)

    assert len(row_reduce(concatenate((Lz, Hz))))==len(Lz)+len(Hz)
    assert len(row_reduce(concatenate((Lz, Gz))))==len(Lz)+len(row_reduce(Gz))

    # Tz = find_errors( Hx            , Lx            )
    #      find_errors( ............. , ............. )
    #                 ( conjugate to  , commutes with )
    #                 ( ............. , ............. )

    Lx = find_errors(Lz, Gz) # invert Lz, commuting with Gz

    check_commute  (Lx, Gz)
    check_commute  (Lx, Hz)
    check_conjugate(Lx, Lz)
    check_commute  (Lz, Gx)
    check_commute  (Lz, Hx)


    # Lx | Lz
    # Hx | ?
    # ?  | Hz
    # ?  | ?
    #Rz = find_logops(concatenate((Lx, Hx)), Hz)
    Rz = dot2(Gz, Pz.transpose())
    Rz = row_reduce(Rz)

    check_commute  (Rz, Lx)
    check_commute  (Rz, Hx)

    Rx = dot2(Gx, Px.transpose())
    Rx = row_reduce(Rx)

    check_commute  (Rx, Lz)
    check_commute  (Rx, Hz)

    # Lx | Lz
    # Hx | ?
    # ?  | Hz
    # Rx'| Rz'

    Tz = find_errors(Hx, concatenate((Lx, Rx)))
    Tx = find_errors(Hz, concatenate((Lz, Rz, Tz)))

    assert len((concatenate((Lx, Hx, Tx, Rx)))) == n
    assert len((concatenate((Lz, Hz, Tz, Rz)))) == n
    assert len(row_reduce(concatenate((Lx, Hx, Tx, Rx)))) == n
    assert len(row_reduce(concatenate((Lz, Hz, Tz, Rz)))) == n

    check_commute  (Rz, Tx)

    Rx = find_errors(Rz, concatenate((Lz, Hz, Tz)))

    check_conjugate(Rx, Rz)
    check_commute  (Rx, Hz)
    check_commute  (Rx, Tz)
    check_commute  (Rx, Lz)

    Rxt = Rx.transpose()
    Rzt = Rz.transpose()

    Pxt = Px.transpose()
    Pzt = Pz.transpose()

    check_sy(Lx, Hx, Tx, Rx, Lz, Hz, Tz, Rz)

#    for gx in Gx:
#        gx1 = dot2(gx, Pxt)
#        print dot2(gx, Rzt), dot2(gx1, Rzt), dot2((gx+gx1)%2, Rzt)

    assert eq2(dot2(Gz, Rxt), dot2(Gz, Pzt, Rxt))
    assert eq2(dot2(Gx, Rzt), dot2(Gx, Pxt, Rzt))

#    print shortstrx(dot2(Rx, Pz), Rx)

    assert eq2(dot2(Rx, Pz), Rx)
    assert eq2(dot2(Rz, Px), Rz)

    assert len(find_kernel(dot2(Gz, Rx.transpose())))==0

    r, n = Rx.shape

    #print shortstrx(Rx, Rz)
    print "r =", len(Rx)

    ops = []
    found = set()
    for gx in Gx:
        #rx = dot2(Px, gx)
        #rx = dot2(gx, PxtRzt)
        rx = dot2(gx, Rzt)
        #print rx, dot2(gx, PxtRzt)
        assert rx.sum()
        rx = mkop(rx, None)
        s = rx.tostring()
        if s not in found:
            found.add(s)
            ops.append(rx)

    print
    cartan = []
    for gz in Gz:
        #rz = dot2(Pz, gz)
        #rz = dot2(gz, PztRxt)
        rz = dot2(gz, Rxt)
        assert rz.sum()
        rz = mkop(None, rz)
        s = rz.tostring()
        if s not in found:
            found.add(s)
            ops.append(rz)
            cartan.append(rz)

    #for op in ops:
    #    print op

    ops = closure(ops)
    print "algebra dimension:", len(ops)
    cartan = [op for op in ops if is_zop(op)]
    print "cartan dimension:", len(cartan)

    return

    lookup = {}
    for i, op in enumerate(ops):
        lookup[op.tostring()] = i

    N = len(ops)

    graph = nx.Graph()
    for i in range(N):
        graph.add_node(i)
    for i, A in enumerate(ops):
      for j, B in enumerate(ops):
        if bracket(A, B)==0:
            continue
        C = (A+B)%2
        k = lookup[C.tostring()]
        graph.add_edge(i, k)
        graph.add_edge(j, k)
    equs = nx.connected_components(graph)
    print "ideals:", len(equs), [len(equ) for equ in equs]

    graph = nx.Graph()
    for i in range(N):
        graph.add_node(i)

    H = []
    for z in cartan:
        A = zeros2(N, N)
        for i, op in enumerate(ops):
            if is_zop(op):
                continue
            c = bracket(z, op)
            if c==0:
                continue
            op1 = (z+op)%2
            j = lookup[op1.tostring()]
            #print "%s->%s" % (i, j),
            A[j, i] = 1
            graph.add_edge(j, i)
        #print
        #print shortstr(A)
        #print
        H.append(A)
    
    for A in H:
      for B in H:
        assert numpy.allclose(numpy.dot(A, B), numpy.dot(B, A))

    equs = nx.connected_components(graph)
    print "orbits:", len(equs), [len(equ) for equ in equs]
    trivial = len([equ for equ in equs if len(equ)==1])
    print "trivial:", trivial, "non-trivial", len(equs)-trivial

    #for irrep in genidx((2,)*r):
    #    print "irrep:", irrep

    return

    r, n = Rx.shape
    N = 2**r
    gz = len(Gz)

    RR = dot2(Gz, Rx.transpose())
    PxtRzt = dot2(Px.transpose(), Rzt)

    if N<=1024:
        H = numpy.zeros((N, N))
    else:
        H = None
    A = {}
    U = []

    basis = []
    lookup = {}
    for i, v in enumerate(genidx((2,)*r)):
        v = array2(v)
        lookup[v.tostring()] = i
        basis.append(v)

    ops = []

    # zops
    for gz in Gz:
        elems = {}
        for i, v in enumerate(basis):
            bit = dot2(gz, Rx.transpose(), v)
            elems[(i, i)] = 1 - 2*bit # send 0,1 -> +1,-1

    # xops
    for gx in Gx:
        elems = {}
        for i, v in enumerate(basis):
            u = (v+dot2(gx, PxtRzt))%2
            j = lookup[u.tostring()]
            elems[j, i] = 1
        for (i, j) in elems:
            assert elems[i, j] == elems[j, i]





from argv import Argv

argv = Argv()

if __name__ == "__main__":

    fn = argv.next()
    if argv.profile:
        import cProfile as profile
        profile.run("%s()"%fn)

    else:
        fn = eval(fn)
        fn()


