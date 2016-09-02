#!/usr/bin/env python

import sys, os
from fractions import gcd

import numpy

import networkx as nx

from solve import get_reductor, array2, row_reduce, dot2, shortstr, zeros2, shortstrx
from solve import u_inverse
from isomorph import write

import models
from models import genidx


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


def test_model():

    Gx, Gz, Hx, Hz = models.build()

    Px = get_reductor(Hx) # projector onto complement of rowspan of Hx
    Pz = get_reductor(Hz)
    Rz = [dot2(Pz, g) for g in Gz] 
    Rz = array2(Rz)
    Rz = row_reduce(Rz, truncate=True)

    Rx = [dot2(Px, g) for g in Gx] 
    Rx = array2(Rx)
    Rx = row_reduce(Rx, truncate=True)

    n = Rx.shape[1]

    print shortstrx(Rx, Rz)

    ops = []
    for rx in Gx:
        ops.append(mkop(rx, None))
    for rz in Gz:
        ops.append(mkop(None, rz))

    bracket(mkop(rx,None), mkop(None, rz))

    found = set(op.tostring() for op in ops)
    while 1:
        new = []
        for a in ops:
          for b in ops:
            if bracket(a, b):
                c = (a+b)%2
                s = c.tostring()
                if s not in found:
                    found.add(s)
                    new.append(c)
        if not new:
            break
        ops.extend(new)
        print len(ops)
    print len(ops)

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


