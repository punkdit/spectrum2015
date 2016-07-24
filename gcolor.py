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

    # bottom faces: points must be adjacent in each face
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

    def check_faces():
        items = [list(face) for face in faces]
        for face in items:
            assert len(face)%2 == 0, face
            face.sort()
        assert len(set([tuple(face) for face in items]))==len(items)
    check_faces()

    # bottom edges
    bedges = []
    for face in bfaces:
        f = len(face)
        for i in range(f):
            bedges.append([face[i], face[(i+1)%f]])

    # edges are not yet unique..
    for edge in bedges:
        edge.sort()
    bedges = list(set([tuple(e) for e in bedges]))

    # extrude bottom edges to make a face
    for edge in bedges:
        edge = list(edge)
        a, b = edge
        face = edge + [a+delta, b+delta]
        faces.append(face)
    check_faces()

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
    check_faces()

    stabs.append([i+delta for i in range(19)] + [top])

    g = len(faces)
    #print "faces:", g

    for stab in stabs:
        assert len(stab)%2 == 0, stab

    #faces.sort()
    #for face in faces:
    #    print face

    Gx = mkop(n, faces)
    Gz = Gx.copy()

    rows = [shortstr(g) for g in Gx]
    #rows.sort()
    #for i, row in enumerate(rows):
    #    print row, faces[i]
    assert len(set(rows))==len(rows)

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


