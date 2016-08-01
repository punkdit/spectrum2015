#!/usr/bin/env python

import sys, os

import numpy

from solve import shortstr, shortstrx, parse, eq2, dot2, zeros2, array2, identity2
from solve import row_reduce, RowReduction, span, get_reductor
from solve import u_inverse, find_logops, solve, find_kernel, linear_independant

from lanczos import write
from code import lstr2

from argv import Argv
argv = Argv()



def genidx(shape):
    if len(shape)==0:
        yield ()
    else:
        for idx in range(shape[0]):
            for _idx in genidx(shape[1:]):
                yield (idx,)+_idx


def check_conjugate(A, B):
    if A is None or B is None:
        return
    assert A.shape == B.shape
    I = numpy.identity(A.shape[0], dtype=numpy.int32)
    assert eq2(dot2(A, B.transpose()), I)


def check_commute(A, B):
    if A is None or B is None:
        return
    C = dot2(A, B.transpose())
    assert C.sum() == 0, "\n%s"%shortstr(C)



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


def build_compass(li, lj=None):

    if lj is None:
        lj = li

    n = li*lj

    keys = [(i, j) for i in range(li) for j in range(lj)]
    coords = {}  
    for i, j in keys:
        for di in range(-li, li+1):
          for dj in range(-lj, lj+1):
            coords[i+di, j+dj] = keys.index(((i+di)%li, (j+dj)%lj))

    m = n 
    Gx = zeros2(m, n)
    Gz = zeros2(m, n)

    idx = 0 
    for i in range(li):
      for j in range(lj):
        Gx[idx, coords[i, j]] = 1 
        Gx[idx, coords[i, j+1]] = 1 

        Gz[idx, coords[i, j]] = 1 
        Gz[idx, coords[i+1, j]] = 1 
        idx += 1

    assert idx == m

    mx = lj-1
    Hx = zeros2(mx, n)
    for idx in range(mx):
      for i in range(li):
        Hx[idx, coords[i, idx]] = 1
        Hx[idx, coords[i, idx+1]] = 1

    mz = li-1
    Hz = zeros2(mz, n)
    for idx in range(mz):
      for j in range(lj):
        Hz[idx, coords[idx, j]] = 1
        Hz[idx, coords[idx+1, j]] = 1

    assert dot2(Hx, Hz.transpose()).sum() == 0

    return Gx, Gz, Hx, Hz


def build_xy(n):

    m = n
    Gx = zeros2(m, n)
    Gz = zeros2(m, n)
    for i in range(m):
        Gx[i, i] = 1
        Gx[i, (i+1)%n] = 1
        
        Gz[i, i] = 1
        Gz[i, (i+1)%n] = 1

    Hx = zeros2(1, n)
    Hz = zeros2(1, n)

    Hx[:] = 1
    Hz[:] = 1

    return Gx, Gz, Hx, Hz


def build_ising(n):

    m = n
    Gx = zeros2(m, n)
    Gz = zeros2(m, n)
    for i in range(m):
        Gz[i, i] = 1
        Gz[i, (i+1)%n] = 1

        Gx[i, i] = 1 # transverse field

    Hx = zeros2(1, n)
    Hz = zeros2(0, n)

    Hx[:] = 1

    return Gx, Gz, Hx, Hz



def mkop(n, ops):
    A = zeros2(len(ops), n)
    for i, op in enumerate(ops):
        for j in op:
            A[i, j] = 1
    return A


def build_gcolor2():

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





def build():

    if argv.gcolor:
        size = argv.get("size", 1)
        Gx, Gz, Hx = build_gcolor(size)
        Hz = Hx.copy()

    elif argv.compass:
        l = argv.get('l', 3)
        li = argv.get('li', l)
        lj = argv.get('lj', l)
        Gx, Gz, Hx, Hz = build_compass(li, lj)

    elif argv.gcolor2:
        Gx, Gz, Hx = build_gcolor2()
        Hz = Hx.copy()

    elif argv.xy:
        n = argv.get('n', 4)
        Gx, Gz, Hx, Hz = build_xy(n)

    elif argv.ising:
        n = argv.get('n', 4)
        Gx, Gz, Hx, Hz = build_ising(n)

    else:

        return

    return Gx, Gz, Hx, Hz



