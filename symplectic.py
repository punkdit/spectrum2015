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

zero = 0


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

    test()

    if argv.profile:
        import cProfile as profile
        profile.run("main()")

    else:
        fn = argv.next()
        fn = eval(fn)
        fn()


