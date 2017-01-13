#!/usr/bin/env python

import sys, os
from random import seed

import numpy
from numpy import concatenate

from solve import shortstr, shortstrx, parse, eq2, dot2, zeros2, array2, identity2
from solve import row_reduce, RowReduction, span, get_reductor
from solve import u_inverse, find_logops, solve, find_kernel, linear_independant
from solve import rand2, find_stabilizers, find_errors

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


def build_compass3(li, lj=None, lk=None):

    if lj is None:
        lj = li

    if lk is None:
        lk = li

    n = li*lj*lk

    keys = [(i, j, k) for i in range(li) for j in range(lj) for k in range(lk)]
    coords = {}  
    for i, j, k in keys:
        for di in range(-li, li+1):
          for dj in range(-lj, lj+1):
            for dk in range(-lk, lk+1):
              coords[i+di, j+dj, k+dk] = keys.index(((i+di)%li, (j+dj)%lj, (k+dk)%lk))

    m = 2*n 
    Gx = zeros2(m, n)
    Gz = zeros2(m, n)

    idx = 0 
    for i in range(li):
      for j in range(lj):
       for k in range(lk):
        Gx[idx, coords[i, j, k]] = 1 
        Gx[idx, coords[i+1, j, k]] = 1 

        Gz[idx, coords[i, j, k]] = 1 
        Gz[idx, coords[i, j+1, k]] = 1 
        idx += 1

        Gx[idx, coords[i, j, k]] = 1 
        Gx[idx, coords[i, j+1, k]] = 1 

        Gz[idx, coords[i, j, k]] = 1 
        Gz[idx, coords[i, j, k+1]] = 1 
        idx += 1

    assert idx == m

#    mx = lj-1
#    Hx = zeros2(mx, n)
#    for idx in range(mx):
#      for i in range(li):
#        Hx[idx, coords[i, idx]] = 1
#        Hx[idx, coords[i, idx+1]] = 1
#
#    mz = li-1
#    Hz = zeros2(mz, n)
#    for idx in range(mz):
#      for j in range(lj):
#        Hz[idx, coords[idx, j]] = 1
#        Hz[idx, coords[idx+1, j]] = 1
#
#    assert dot2(Hx, Hz.transpose()).sum() == 0

    Hx = Hz = None

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

    if n%2 == 0:
        Hx = zeros2(1, n)
        Hz = zeros2(1, n)
    
        Hx[:] = 1
        Hz[:] = 1

    else:
        Hx = Hz = None


    return Gx, Gz, Hx, Hz


def build_random(n):

    weight = argv.get("weight", 3)
    coweight = argv.get("coweight")

    p = argv.get("p", 0.3)
    m = argv.get("m", n)
    mx = argv.get("mx", m)
    mz = argv.get("mz", m)

    if coweight is not None:
        Gx = rand2(n, mx, weight=coweight).transpose()
        Gz = rand2(n, mz, weight=coweight).transpose()

    else:
        Gx = rand2(mx, n, p=p, weight=weight)
        Gz = rand2(mz, n, p=p, weight=weight)

    Hx = Hz = None

    Gx = Gx[[i for i in range(m) if Gx[i].sum()], :]
    Gz = Gz[[i for i in range(m) if Gz[i].sum()], :]

    li = argv.get("li", True)

    if li:
        Gx = linear_independant(Gx)
        Gz = linear_independant(Gz)

    return Gx, Gz, Hx, Hz


def build_random_selfdual(n):

    weight = argv.get("weight", 3)
    m = argv.get("m", n)
    h = argv.get("h", 0)

    while 1:
        Gx = rand2(m, n, weight=weight)
        Gz = Gx.copy()
    
        Hx = Hz = None
    
        Gx = Gx[[i for i in range(m) if Gx[i].sum()], :]
        Gx = linear_independant(Gx)
    
        if len(Gx)<m:
            write("m")
            continue

        Gz = Gx.copy()

        Hx = find_stabilizers(Gz, Gx)
        Hz = find_stabilizers(Gx, Gz)

        if len(Hx)==h and len(Hz)==h:
            break

        write("H(%d,%d)"%(len(Hx), len(Hz)))

    print
    return Gx, Gz, Hx, Hz


def build_random_nostabs(n):

    m = argv.get("m", n)
    mx = argv.get("mx", m)
    mz = argv.get("mz", m)
    h = argv.get("h", 0)
    hx = argv.get("hx", h)
    hz = argv.get("hz", h)
    while 1:

        Gx, Gz, Hx, Hz = build_random(n)

        if len(Gx)<mx or len(Gz)<mz:
            write("m")
            continue

        Hx = find_stabilizers(Gz, Gx)
        Hz = find_stabilizers(Gx, Gz)

        if len(Hx)==hx and len(Hz)==hz:
            break

        write("H(%d,%d)"%(len(Hx), len(Hz)))

    print

    return Gx, Gz, Hx, Hz


def build_pauli(n):

    m = 2**n
    Gx = zeros2(m, n)
    Gz = zeros2(m, n)
    for i, idx in enumerate(genidx((2,)*n)):
        for j in idx:
            Gx[i,j] = 1
            Gz[i,j] = 1

    Hx = zeros2(0, n)
    Hz = zeros2(0, n)

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



#def build_ising2():
#
#    l = argv.get("l", 6)
#    assert l%2 == 0
#
#    li = lj = l
#    n = l**2
#
#    keys = [(i, j) for i in range(li) for j in range(lj)]
#    coords = {}  
#    for i, j in keys:
#        for di in range(-li, li+1):
#          for dj in range(-lj, lj+1):
#            coords[i+di, j+dj] = keys.index(((i+di)%li, (j+dj)%lj))
#
#    assert n%4==0
#    m = n/4
#
#    Gx = zeros2(m, n)
#    Gz = zeros2(m, n)
#
#    idx = 0
#    for i1 in range(l//2):
#      for j1 in range(l//2):
#        i = 2*i1
#        j = 2*j1
#        Gx[idx, coords[i, j]] = 1
#        Gx[idx, coords[i+1, j]] = 1
#        Gx[idx, coords[i, j+1]] = 1
#        Gx[idx, coords[i+1, j+1]] = 1
#
#        Gz[idx, coords[i+1, j+1]] = 1
#        Gz[idx, coords[i+2, j+1]] = 1
#        Gz[idx, coords[i+1, j+2]] = 1
#        Gz[idx, coords[i+2, j+2]] = 1
#
#        idx += 1
#
#    return Gx, Gz, None, None


def build_hex(li, lj=None):
    if lj is None:
        lj = li
    n = li*lj

    keys = [(i, j) for i in range(li) for j in range(lj)]
    coords = {}  
    for i, j in keys:
        for di in range(-li, li+1):
          for dj in range(-lj, lj+1):
            coords[i+di, j+dj] = keys.index(((i+di)%li, (j+dj)%lj))

    Gx = []
    if argv.open:
        idxs = range(li-1)
        jdxs = range(lj-1)
    else:
        idxs = range(li)
        jdxs = range(lj)

    for i in idxs:
      for j in jdxs:
        g = zeros2(n)
        g[coords[i,   j]] = 1 
        g[coords[i,   j+1]] = 1 
        g[coords[i+1, j+1]] = 1 
        Gx.append(g)

        g = zeros2(n)
        g[coords[i,   j]] = 1 
        g[coords[i+1, j]] = 1 
        g[coords[i+1, j+1]] = 1 
        Gx.append(g)
    Gx = array2(Gx)

    Gz = Gx.copy()

    return Gx, Gz, None, None



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


def build_projective(n, dim=2):

    import geometry
    g = geometry.projective(n, dim)

    P = g.types[0]
    L = g.types[1]
    if dim==3:
        L = g.types[2]

    points = g.tplookup[P]
    lines = g.tplookup[L]

    #lines = lines[:-4] # throw one out
    #points = points[:-1] # throw one out

    n = len(points)
    m = len(lines)
    Gx = zeros2(m, n)
    for i, line in enumerate(lines):
        for j, point in enumerate(points):
            if (line, point) in g.incidence:
                Gx[i, j] = 1

    #print shortstr(Gx)

    Gz = Gx.copy()

    Hx = None
    Hz = None

    return Gx, Gz, Hx, Hz




def build(name=""):

    if name:
        setattr(argv, name, True) # hack this

    _seed = argv.get("seed")
    if _seed is not None:
        numpy.random.seed(_seed)
        seed(_seed)

    size = argv.get("size", 1)

    if argv.gcolor2 or (argv.gcolor and size==1.5):
        Gx, Gz, Hx = build_gcolor2()
        Hz = Hx.copy()

    elif argv.gcolor:
        Gx, Gz, Hx = build_gcolor(size)
        Hz = Hx.copy()

    elif argv.compass:
        l = argv.get('l', 3)
        li = argv.get('li', l)
        lj = argv.get('lj', l)
        Gx, Gz, Hx, Hz = build_compass(li, lj)

    elif argv.compass3:
        l = argv.get('l', 3)
        li = argv.get('li', l)
        lj = argv.get('lj', l)
        lk = argv.get('lk', l)
        Gx, Gz, Hx, Hz = build_compass3(li, lj, lk)

    elif argv.hex:
        l = argv.get('l', 3)
        li = argv.get('li', l)
        lj = argv.get('lj', l)
        Gx, Gz, Hx, Hz = build_hex(li, lj)

    elif argv.xy:
        n = argv.get('n', 4)
        Gx, Gz, Hx, Hz = build_xy(n)

    elif argv.ising:
        n = argv.get('n', 4)
        Gx, Gz, Hx, Hz = build_ising(n)

    elif argv.random:
        n = argv.get('n', 4)
        Gx, Gz, Hx, Hz = build_random(n)

    elif argv.random_nostabs:
        n = argv.get('n', 4)
        Gx, Gz, Hx, Hz = build_random_nostabs(n)

    elif argv.random_selfdual:
        n = argv.get('n', 4)
        Gx, Gz, Hx, Hz = build_random_selfdual(n)

    elif argv.pauli:
        n = argv.get('n', 2)
        Gx, Gz, Hx, Hz = build_pauli(n)

    elif argv.projective:
        n = argv.get('n', 3)
        dim = argv.get('dim', 2)
        Gx, Gz, Hx, Hz = build_projective(n, dim)

    elif argv.test:
        Gx, Gz, Hx, Hz = build_test()

    else:

        name = argv.next()
        try:
            fn = eval("build_%s"%name)
        except NameError:
            print "no model found"
            raise

        Gx, Gz, Hx, Hz = fn()

    if Hx is None:
        Hx = find_stabilizers(Gz, Gx)
    if Hz is None:
        Hz = find_stabilizers(Gx, Gz)

    if argv.flip:
        Gz, Gx = Gx, Gz
        Hz, Hx = Hx, Hz

    if argv.show:
        print "Gx Gz:"
        print shortstrx(Gx, Gz)
        if len(Hx):
            print "Hx Hz:"
            print shortstrx(Hx, Hz)

    return Gx, Gz, Hx, Hz


def build_reduced():

    Gx, Gz, Hx, Hz = build()

    Px = get_reductor(Hx) # projector onto complement of rowspan of Hx
    Pz = get_reductor(Hz)
    Rz = [dot2(Pz, g) for g in Gz]
    Rz = array2(Rz)
    Rz = row_reduce(Rz, truncate=True)

    Rx = [dot2(Px, g) for g in Gx]
    Rx = array2(Rx)
    Rx = row_reduce(Rx, truncate=True)

    return Rx, Rz



class Model(object):
    def __init__(self, attrs):
        self.__dict__.update(attrs)
        self.Qx = self.Rz.transpose() # backwards compat

    def __str__(self):
        return "Model(Lx/z: %d, Gx: %d, Gz: %d, Hx: %d, Hz: %d, Rx/z: %d)" % (
            len(self.Lx), len(self.Gx), len(self.Gz), 
            len(self.Hx), len(self.Hz), len(self.Rx))

    def build_ham(self, excite=None, weights=None, Jx=1., Jz=1.):
        Gx, Gz = self.Gx, self.Gz        
        Rx, Rz = self.Rx, self.Rz        
        Hx, Hz = self.Hx, self.Hz        
        Tx, Tz = self.Tx, self.Tz        
        gz = len(Gz)
        r = len(Rx)
        n = self.n

        if type(excite) is int:
            _excite = [0]*len(Tx)
            _excite[excite] = 1
            excite = tuple(_excite)

        if excite is not None:
            assert len(excite)==len(Tx)
    
            t = zeros2(n)
            for i, ex in enumerate(excite):
                if ex:
                    t = (t + Tx[i])%2
            #print "t:", shortstr(t)
            Gzt = dot2(Gz, t)

        else:
            Gzt = 0

        if weights is None:
            weights = [1.]*len(Gx)
        assert len(weights) == len(Gx), len(weights)

        H = numpy.zeros((2**r, 2**r))
        for i, v in enumerate(genidx((2,)*r)):
            v = array2(v)
            syndrome = (dot2(Gz, Rx.transpose(), v) + Gzt)%2
            value = gz - 2*syndrome.sum()
            #print shortstr(dot2(Rx.transpose(), v)), value
            H[i, i] = Jz*value
            #U.append(value)

        Pxt = self.Px.transpose()
        Qx = Rz.transpose()
        #print dot2(Rx, Qx)
        PxtQx = dot2(Pxt, Qx)
        for i, v in enumerate(genidx((2,)*r)):
            v = array2(v)
            #print shortstr(v),
            #for g in Gx:
            for j, g in enumerate(Gx):
                u = (v + dot2(g, PxtQx))%2
                k = eval('0b'+shortstr(u, zero='0'))
                H[i, k] += Jx*weights[j]
                #A[i, k] = A.get((i, k), 0) + 1

        return H

#        H = numpy.zeros((n, n))
#        syndromes = [] 
#        for i, v in enumerate(verts):
#            syndromes.append(dot2(Gz, v))
#            count = dot2(Gz, v).sum()
#            Pxv = dot2(Px, v)
#            assert count == dot2(Gz, Pxv).sum()
#            H[i, i] = mz - 2*count
#            for g in Gx:
#                v1 = (g+v)%2
#                v1 = dot2(Px, v1)
#                j = lookup[v1.tostring()]
#                H[i, j] += 1 




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



def build_model(Gx=None, Gz=None, Hx=None, Hz=None):

    if Gx is None:
        Gx, Gz, Hx, Hz = build()

    n = Gx.shape[1]

    if Hx is None:
        Hx = find_stabilizers(Gz, Gx)
    if Hz is None:
        Hz = find_stabilizers(Gx, Gz)

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
    assert Lz.shape[1] == n

    if 0:
        PGz = get_reductor(Gz)
        Lz = dot2(Lz, PGz.transpose())
        Lz = row_reduce(Lz)
    
        print shortstrx(Lz, Gz, Hz)

    if len(Lz):
        #print Lz.shape, Hz.shape
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

    assert eq2(dot2(Gz, Rxt), dot2(Gz, Pzt, Rxt))
    assert eq2(dot2(Gx, Rzt), dot2(Gx, Pxt, Rzt))

#    print shortstrx(dot2(Rx, Pz), Rx)

    assert eq2(dot2(Rx, Pz), Rx) 
    assert eq2(dot2(Rz, Px), Rz) 

    assert len(find_kernel(dot2(Gz, Rx.transpose())))==0

    model = Model(locals())
    return model



if __name__ == "__main__":

    Gx, Gz, Hx, Hz = build()

    model = build_model(Gx, Gz, Hx, Hz)

    print model

    print "Hx/Hz:"
    print shortstrx(model.Hx, model.Hz)
    print
    print "Gx/Gz:"
    print shortstrx(Gx, Gz)
    print
    print "Lx/Lz:"
    print shortstrx(model.Lx, model.Lz)

    if len(model.Lx):
        w = min([v.sum() for v in span(model.Lx) if v.sum()])
        print "distance:", w
        


