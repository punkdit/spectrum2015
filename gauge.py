#!/usr/bin/env python

import sys, os

import numpy

from pyx import canvas, path, deco, trafo, style, text, color, unit, epsfile, deformer
rgb = color.rgb
rgbfromhexstring = color.rgbfromhexstring
red, green, blue, yellow = (rgb.red,
    rgbfromhexstring("#008000"),
    rgb.blue, rgb(0.75, 0.75, 0))
blue = rgb(0., 0., 0.8)
lred = rgb(1., 0.4, 0.4)
lgreen = rgb(0.4, 0.8, 0.2)
    
from solve import shortstr, shortstrx, parse, eq2, dot2, zeros2, array2, identity2
from solve import row_reduce, RowReduction, span, get_reductor
from solve import u_inverse, find_logops, solve, find_kernel, linear_independant

from lanczos import write
from code import lstr2

from argv import Argv
argv = Argv()

from qupy.ldpc import gcolor


def distance(graph, src, tgt):
    if type(src) is int:
        src=[src]
    if type(tgt) is int:
        tgt=[tgt]
    src = set(src)
    tgt = set(tgt)
    d = 0
    bdy = set(src)
    while not bdy.intersection(tgt):
        assert bdy, "src tgt in different components"
        _bdy = set()
        for i in bdy:
            for j in graph[i]:
                if j not in src:
                    _bdy.add(j)

        src.update(_bdy)
        bdy = _bdy

        d += 1

    items = list(bdy.intersection(tgt))
    assert len(items)==1
    return 1.*d, items[0]


def conv(r, c0, c1):
    assert 0.<=r<=1.
    c = tuple( (1.-r)*c0[i] + r*c1[i] for i in range(3) )
    return c

def conv3(r0, r1, r2, c0, c1, c2):
    assert 0.<=r0<=1.
    assert 0.<=r1<=1.
    assert 0.<=r2<=1.
    c = tuple( r0*c0[i] + r1*c1[i] + r2*c2[i] for i in range(3) )
    return c

def sub2(p0, p1):
    return (p0[0]-p1[0], p0[1]-p1[1])

def cross2(p0, p1):
    return (p0[0]*p1[1] - p0[1]*p1[0])
    


def orient(ps): # FAIL
    ps = list(ps)
    p0 = ps[0]
    i = 1
    while i+1<len(ps):
        p1 = sub2(ps[i], p0)
        p2 = sub2(ps[i+1], p0)
        if cross2(p1, p2) < 0:
            ps[i], ps[i+1] = ps[i+1], ps[i]
        i += 1
    return ps


def orient(graph, idxs):
    N = len(idxs)
    idxs = list(idxs)
    i0 = idxs.pop(0)
    result = [i0]
    while idxs:
        for i1 in graph[i0]:
            if i1 in idxs:
                result.append(i1)
                idxs.remove(i1)
                i0 = i1
                break
    assert len(result)==N
    return result


def main():

    size = argv.get("size", 3)


    lattice = gcolor.Lattice(size)

    n = len(lattice.qubits)
    print lattice

    code = lattice.build_code(check=False)
    #Ex = lattice.Ex
    Gx, Gz = code.Gx, code.Gz
    Hx, Hz = code.Hx, code.Hz

    CORNERS = 4
    EDGES = 6
    FACES = 4

    corner_idxs = []
    edge_idxs = []
    face_idxs = []
    bulk_idxs = []
    items = [corner_idxs, edge_idxs, face_idxs, bulk_idxs]

    for i in range(n):
        w = Hx[:, i].sum()
        assert 1<=w<=4
        items[w-1].append(i)

    assert len(corner_idxs)==CORNERS
    
    graph = dict((i, []) for i in range(n)) # map qubit -> nbd
    for g1 in Gx:
      for g2 in Gx:
        if str(g1)==str(g2):
            continue
        g = g1*g2
        for i in numpy.where(g)[0]:
          for j in numpy.where(g)[0]:
            if i==j:
              continue
            if i not in graph[j]:
              graph[j].append(i)
            if j not in graph[i]:
              graph[i].append(j)

    #print graph
    for i in range(n):
        assert len(graph[i])==CORNERS or i in corner_idxs and len(graph[i])==3
    

    coords = {} # map qubit -> (x, y, z)

    # corners of a simplex:
    R = 5.
    cs = [
        (R, R, R),
        (R, -R, -R),
        (-R, -R, R),
        (-R, R, -R)]
    for ii, i in enumerate(corner_idxs):
        coords[i] = cs[ii]
        
    # map pair of corners to list of edge qubits
    cedges = dict(((i, j), []) for i in corner_idxs for j in corner_idxs if i<j)
    assert len(cedges)==EDGES

    for i in edge_idxs:
        # find closest two corners
        djs = [distance(graph, i, j)[0] for j in corner_idxs]
        idxs = range(CORNERS)
        idxs.sort(key = lambda idx : djs[idx])
        #print idxs, djs
        i0, i1 = idxs[:2]
        dj0, dj1 = djs[i0], djs[i1]
        c0, c1 = coords[corner_idxs[i0]], coords[corner_idxs[i1]]
        r = dj0 / (dj0+dj1)
        coords[i] = conv(r, c0, c1)
        ci0, ci1 = corner_idxs[i0], corner_idxs[i1]
        if ci0>ci1:
            ci0,ci1=ci1,ci0
        cedges[ci0, ci1].append(i)

    print coords
    print cedges.values()
    match = [corner_idxs[i] for i in [0,1,3]]
    match.sort()
    for i in face_idxs:

        # find three closest corners
        djs = [distance(graph, i, j) for j in corner_idxs]
        djs.sort()
        assert djs[3][0] > djs[2][0]
        djs = djs[:3]
        cidxs = [dj[1] for dj in djs] # corner idxs
        cidxs.sort()
        if cidxs!=match:
            continue
        c0, c1, c2 = cidxs
        eidxss = [cedges[c0, c1], cedges[c0,c2], cedges[c1,c2]] # edges
        djs = [distance(graph, i, eidxs) for eidxs in eidxss]

        cs = [coords[dj[1]] for dj in djs]
        #print cs

        ds = [dj[0] for dj in djs]
        d = 2*(ds[0]+ds[1]+ds[2])
        coord = conv3(
            (ds[1]+ds[2])/d, 
            (ds[0]+ds[2])/d, 
            (ds[0]+ds[1])/d, 
            *cs)
        coords[i] = coord

    #for i in bulk_idxs:
    #    djs = [distance(graph, i, face) 


    def tran(x, y, z):
        return x+0.4*z, y+0.1*z
        return (x + 0.6*y + 0.6*z), (y + 0.2*z)
    
    c = canvas.canvas()

    tr = trafo.scale(sx=1, sy=1)

    for g in Gx:
        idxs = [i for i in numpy.where(g)[0]]
        if None in [coords.get(i) for i in idxs]:
            continue
        idxs = orient(graph, idxs)
        ps = [tran(*coords[i]) for i in idxs]
        ps = [path.moveto(*ps[0])]+[path.lineto(*p) for p in ps[1:]]+[path.closepath()]
        c.fill(path.path(*ps), [red, tr])
    
    for i, (x0, y0, z0) in coords.items():
        x0, y0 = tran(x0, y0, z0)
        for j in graph[i]:
            if j not in coords:
                continue
            x1, y1 = tran(*coords[j])
            c.stroke(path.line(x0, y0, x1, y1), [tr])

    for i, (x, y, z) in coords.items():
        x, y = tran(x, y, z)

        radius = 0.20 + 0.00*z
        c.fill(path.circle(x, y, radius), [tr])

    c.writePDFfile("pic-qubits.pdf")






if __name__=="__main__":
    main()


