#!/usr/bin/env python

import sys, os
from fractions import gcd

import numpy

import networkx as nx

from solve import get_reductor, array2, row_reduce, dot2, shortstr, zeros2, shortstrx, eq2
from solve import u_inverse, find_kernel, find_logops
from isomorph import write

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
    pairs = [(a, b) for a in ops for b in ops if not numpy.allclose(a, b)]
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
        pairs = [(a, b) for a in ops+new for b in new if not numpy.allclose(a, b)]
        ops.extend(new)
        if argv.verbose:
            print len(ops)
    return ops


def test_model():

    Gx, Gz, Hx, Hz = models.build()

    Lz = find_logops(Gx, Hz)
    Lx = find_logops(Gz, Hx)

    #Px = get_reductor(numpy.concatenate((Lx, Hx)))
    #Pz = get_reductor(numpy.concatenate((Lz, Hz)))
    Px = get_reductor(Hx)
    Pz = get_reductor(Hz)
    Pxt = Px.transpose()
    Pzt = Pz.transpose()

    Rx = dot2(Gx, Px.transpose())
    Rx = row_reduce(Rx, truncate=True)

    Rz = dot2(Gz, Pz.transpose())
    Rz = row_reduce(Rz, truncate=True)

    Rxt = Rx.transpose()
    Rzt = Rz.transpose()

    #print shortstrx(dot2(Gz, Rxt), dot2(Gz, Pz, Rxt))

    assert dot2(Hx, Rzt).sum() == 0
    assert dot2(Hz, Rxt).sum() == 0

#    for gx in Gx:
#        gx1 = dot2(gx, Pxt)
#        print dot2(gx, Rzt), dot2(gx1, Rzt), dot2((gx+gx1)%2, Rzt)

    assert eq2(dot2(Gz, Rxt), dot2(Gz, Pzt, Rxt))
    assert eq2(dot2(Gx, Rzt), dot2(Gx, Pxt, Rzt))

#    print shortstrx(dot2(Rx, Pz), Rx)

    assert eq2(dot2(Rx, Pz), Rx)
    assert eq2(dot2(Rz, Px), Rz)

    #return

    Qx = u_inverse(Rx) # a.k.a Rz.transpose()
    PxtQx = dot2(Pxt, Qx)

    Qz = u_inverse(Rz) # a.k.a Rx.transpose()
    PztQz = dot2(Pzt, Qz)

    assert eq2(dot2(Gz, Qz), dot2(Gz, Pzt, Qz))
    assert eq2(dot2(Gx, Qx), dot2(Gx, Pxt, Qx))

    #print shortstrx(Rz, Qz)

    r, n = Rx.shape

    #print shortstrx(Rx, Rz)
    print "r =", len(Rx)

    ops = []
    found = set()
    for gx in Gx:
        #rx = dot2(Px, gx)
        #rx = dot2(gx, PxtQx)
        rx = dot2(gx, Qx)
        print rx, dot2(gx, PxtQx)
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
        #rz = dot2(gz, PztQz)
        rz = dot2(gz, Qz) # fail
        assert rz.sum()
        rz = mkop(None, rz)
        s = rz.tostring()
        if s not in found:
            found.add(s)
            ops.append(rz)
            cartan.append(rz)

    ops = closure(ops)
    print "size:", len(ops)
    cartan = [op for op in ops if is_zop(op)]
    print "cartan:", len(cartan)

    lookup = {}
    for i, op in enumerate(ops):
        lookup[op.tostring()] = i

    N = len(ops)

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
    print "orbits:", len(equs)
    for equ in equs:
        print len(equ),
    print

    #for irrep in genidx((2,)*r):
    #    print "irrep:", irrep

    return

    r, n = Rx.shape
    N = 2**r
    gz = len(Gz)

    Qx = u_inverse(Rx)
    RR = dot2(Gz, Rx.transpose())
    PxtQx = dot2(Px.transpose(), Qx)

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
            u = (v+dot2(gx, PxtQx))%2
            j = lookup[u.tostring()]
            elems[j, i] = 1
        for (i, j) in elems:
            assert elems[i, j] == elems[j, i]





from argv import Argv

argv = Argv()

if __name__ == "__main__":

    if argv.profile:
        import cProfile as profile
        profile.run("main()")

    else:
        fn = argv.next()
        fn = eval(fn)
        fn()


