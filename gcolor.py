#!/usr/bin/env python

import sys


import numpy

from qupy.ldpc.solve import zeros2, shortstr, write, dot2, parse, array2, shortstrx
from qupy.ldpc.solve import enum2, span, find_logops
from qupy.ldpc.css import CSSCode, check_commute



from argv import Argv
argv = Argv()

def mkop(n, ops):
    A = zeros2(len(ops), n)
    for i, op in enumerate(ops):
        for j in op:
            A[i, j] = 1
    return A


def build():

    n = 39
    m = 10

    delta = 19

    top = n-1 # top qubit

    # bottom faces
    bfaces = [
        [0, 1, 2, 3],
        [1, 4, 5, 6, 7, 2],
        [3, 2, 7, 8],
        [4, 9, 10, 5],
        [8, 7, 6, 11, 12, 13],
        [9, 14, 15, 10],
        [5, 10, 15, 16, 11, 6],
        [11, 16, 17, 12],
        [13, 12, 17, 18]]

    faces = list(bfaces) + [[i+delta for i in face] for face in bfaces]

    # bottom edges
    bedges = []
    for face in bfaces:
        f = len(face)
        for i in range(f):
            bedges.append([face[i], face[(i+1)%f]])

    # extrude bottom edges to make a face
    for edge in bedges:
        a, b = edge
        face = edge + [a+delta, b+delta]
        faces.append(face)

    stabs = []
    for face in bfaces:
        stabs.append(face + [i+delta for i in face])

    # top faces
    for face in [
        [0, 1, 4, 9, 14],
        [0, 3, 8, 13, 18],
        [14, 15, 16, 17, 18]]:
        face = [i+delta for i in face] + [top]
        faces.append(face)

    stabs.append([i+delta for i in range(19)] + [top])

    g = len(faces)
    #print "faces:", g

    for stab in stabs:
        assert len(stab)%2 == 0, stab

    for face in faces:
        assert len(face)%2 == 0, face

    Gx = mkop(n, faces)
    Gz = Gx.copy()

    Hz = mkop(n, stabs)
    Hx = Hz.copy()

    # bottom face
    Lx = mkop(n, [range(19)])
    Lz = Lx.copy()

    check_commute(Hx, Hz)
    check_commute(Hx, Gz)
    check_commute(Hz, Gx)
    check_commute(Gx, Lz)
    check_commute(Gz, Lx)
    check_commute(Hx, Lz)
    check_commute(Hz, Lx)

    #code = CSSCode(Hx=Hx, Gx=Gx, Hz=Hz, Gz=Gz, build=False)

    Lx = find_logops(Gz, Hx)

    #print "Lx:", shortstr(Lx)

    return Gx, Gz, Hx



if __name__=="__main__":
    build()


