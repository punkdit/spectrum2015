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
colors = [red, green, blue, yellow]

    
from lanczos import write

from argv import Argv
argv = Argv()

from qupy.ldpc import gcolor
import models


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

def conv4(r0, r1, r2, r3, c0, c1, c2, c3):
    assert 0.<=r0<=1.
    assert 0.<=r1<=1.
    assert 0.<=r2<=1.
    assert 0.<=r3<=1.
    c = tuple( r0*c0[i] + r1*c1[i] + r2*c2[i] + r3*c3[i] for i in range(3) )
    return c

def sub2(p0, p1):
    return (p0[0]-p1[0], p0[1]-p1[1])

def cross2(p0, p1):
    return (p0[0]*p1[1] - p0[1]*p1[0])
    


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


class Cell(object):
    dim = None
    mark = None # the "colour" of this cell: 0,1,2 or 3
    cellmap = {}
    def __init__(self, idxs, coords, deco=[]):
        idxs = list(idxs)
        self.idxs = idxs
        key = list(idxs)
        key.sort()
        key = tuple(key)
        self.cellmap[key] = self
        coords = [coords.get(i) for i in idxs]
        assert len(idxs)==len(coords)
        self.deco = deco
        self.ready = None not in coords
        if self.ready:
#            self.z = 1.*sum(coords[i][2] for i in range(len(coords))) / len(coords)
            self.z = min(coords[i][2] for i in range(len(coords)))
        else:
            self.z = 0.
        self.coords = coords
        self.bdy = [] # boundary
        self.cobdy = [] # co-boundary

    @classmethod
    def get(cls, *idxs):
        key = list(idxs)
        key.sort()
        key = tuple(key)
        cell = cls.cellmap.get(key)
        return cell

    def add(self, child):
        assert self.dim == child.dim+1
        for idx in child.idxs:
            assert idx in self.idxs
        if child not in self.bdy:
            self.bdy.append(child)
        if self not in child.cobdy:
            child.cobdy.append(self)

    def tran(self, x, y, z):
        return x+0.4*z, y+0.1*z
        return (x + 0.6*y + 0.6*z), (y + 0.2*z)
    
    def render(self, c):
        pass


class Vertex(Cell):
    dim = 0
    radius = 0.20
    def render(self, c, deco=[]):
        ps = [self.tran(*coord) for coord in self.coords]
        assert len(ps)==1
        x, y = ps[0]
        c.fill(path.circle(x, y, self.radius), deco+self.deco)


class Edge(Cell):
    dim = 1
    def render(self, c, deco=[]):
        ps = [self.tran(*coord) for coord in self.coords]
        assert len(ps)==2
        x0, y0 = ps[0]
        x1, y1 = ps[1]
        c.stroke(path.line(x0, y0, x1, y1), deco+self.deco)


class Face(Cell):
    dim = 2

    def render(self, c, deco=[]):
        cs = [self.tran(*coord) for coord in self.coords]
        assert len(cs)>2
        ps = [path.moveto(*cs[0])]+[path.lineto(*p) for p in cs[1:]]+[path.closepath()]
        c.fill(path.path(*ps), deco+self.deco)
        for x, y in cs:
            c.fill(path.circle(x, y, Vertex.radius), deco)
        c.stroke(path.path(*ps), deco)


class Body(Cell):
    dim = 3



def main():

    size = argv.get("size", 3)
    if size==1.5:
        Gx, Gz, Hx, Hz = models.build("gcolor2")
    else:
        Gx, Gz, Hx, Hz = models.build("gcolor")

    gx, n = Gx.shape

    CORNERS = 4
    EDGES = 6
    FACES = 4

    # -----------------------------------------------------
    # Where does each qubit live in relation to boundary ?
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
    
    # ---------------------------------------
    # Coordinatize qubits

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
    #match = [corner_idxs[i] for i in [0,1,3]]
    #match.sort()
    cfaces = {}
    for key in [(0,1,2), (0,1,3), (0,2,3), (1,2,3)]:
        key = tuple(corner_idxs[k] for k in key)
        cfaces[key] = []
    assert len(cfaces)==FACES
    for i in face_idxs:

        # find three closest corners
        djs = [distance(graph, i, j) for j in corner_idxs]
        djs.sort()
        assert djs[3][0] > djs[2][0]
        djs = djs[:3]
        cidxs = [dj[1] for dj in djs] # corner idxs
        cidxs.sort()
        #if cidxs!=match:
        #    continue
        c0, c1, c2 = cidxs
        cfaces[tuple(cidxs)].append(i)
        eidxss = [cedges[c0, c1], cedges[c0,c2], cedges[c1,c2]] # edges
        djs = [distance(graph, i, eidxs) for eidxs in eidxss]

        cs = [coords[dj[1]] for dj in djs]
        ds = [dj[0] for dj in djs]
        d = 2*(ds[0]+ds[1]+ds[2])
        coord = conv3((ds[1]+ds[2])/d, (ds[0]+ds[2])/d, (ds[0]+ds[1])/d, *cs)
        coords[i] = coord

    print "cfaces:", cfaces
    for i in bulk_idxs:
        djs = [distance(graph, i, face) for face in cfaces.values()]
        cs = [coords[dj[1]] for dj in djs]
        ds = [dj[0] for dj in djs]
        d = 3*sum(ds)
        coord = conv4((ds[1]+ds[2]+ds[3])/d, (ds[0]+ds[2]+ds[3])/d, (ds[0]+ds[1]+ds[3])/d, (ds[0]+ds[1]+ds[2])/d, *cs)
        coords[i] = coord


    # ------------------------
    # Build Cells & connect

    cells = []

    for i in range(n):
        vertex = Vertex([i], coords)
        cells.append(vertex)
        for j in graph[i]:
            if i>j:
                continue
            edge = Edge([i, j], coords)
            cells.append(edge)

    for i in range(n):
        vi = Cell.get(i)
        for j in graph[i]:
            if i>j:
                continue
            vj = Cell.get(j)
            e = Cell.get(i, j)
            e.add(vi)
            e.add(vj)

    bodys = [] # sic
    for h in Hx:
        idxs = [i for i in numpy.where(h)[0]]
        body = Body(idxs, coords)
        bodys.append(body)

    for g in Gx:
        idxs = [i for i in numpy.where(g)[0]]
        #if None in [coords.get(i) for i in idxs]:
        #    continue
        idxs = orient(graph, idxs)
        face = Face(idxs, coords, [red])
        cells.append(face)
        for ii in range(len(idxs)):
            i = idxs[ii]
            j = idxs[(ii+1)%len(idxs)]
            e = Cell.get(i, j)
            assert e
            face.add(e)
        for body in bodys:
            if len(set(body.idxs).intersection(face.idxs)) == len(face.idxs):
                body.add(face)

    cells.extend(bodys)

    # --------- mark ---------

    for i, idx in enumerate(corner_idxs):
        vertex = Cell.get(idx)
        face = vertex.cobdy[0].cobdy[0]
        assert len(face.cobdy)==1
        body = face.cobdy[0]
        body.mark = i

    # -------- render ---------
    
    c = canvas.canvas()
    tr = trafo.scale(sx=1, sy=1)

    cells.sort(key = lambda cell : -cell.z)
    for cell in cells:
        if not cell.ready:
            continue
        if cell.dim != 2:
            continue
            
        deco = [color.transparency(0.2)]
        cell.render(c, deco)


    c.writePDFfile("pic-qubits.pdf")






if __name__=="__main__":
    main()


