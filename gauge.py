#!/usr/bin/env python

import sys, os

import numpy

from pyx import canvas, path, deco, trafo, style, text, color, unit, epsfile, deformer
rgb = color.rgb
rgbhex = color.rgbfromhexstring
red, green, blue, yellow = (rgb.red,
    rgbhex("#008000"),
    rgb.blue, rgb(0.75, 0.75, 0))
blue = rgb(0., 0., 0.8)
lred = rgb(1., 0.4, 0.4)
lgreen = rgb(0.4, 0.8, 0.2)
colors = [red, green, blue, yellow]


north = [text.halign.boxcenter, text.valign.top]
northeast = [text.halign.boxright, text.valign.top]
northwest = [text.halign.boxleft, text.valign.top]
south = [text.halign.boxcenter, text.valign.bottom]
southeast = [text.halign.boxright, text.valign.bottom]
southwest = [text.halign.boxleft, text.valign.bottom]
east = [text.halign.boxright, text.valign.middle]
west = [text.halign.boxleft, text.valign.middle]
center = [text.halign.boxcenter, text.valign.middle]



    
from lanczos import write
from solve import shortstrx

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
    #assert len(items)==1
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
        self.idxs = idxs # qubit support
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

    def __str__(self):
        return "%s%s" % (self.__class__.__name__, tuple(self.idxs))

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

    R = argv.get("R", 2.5)

    corner_coords = [(R, R, R), (R, -R, -R), (-R, -R, R), (-R, R, -R)]
    # z points away from viewpoint
    BACK, FRONT = [0, 2], [1, 3]
    def tran(self, x, y, z):
        return x+0.4*z, y+0.1*z

    corner_coords = [(0,-R,-R), (R,-R,R), (-R,-R,R), (0,0.8*R,0)]
    BACK, FRONT = [1, 2], [3, 0]
    def tran(self, x, y, z):
        return x+0.3*z, y+0.4*z

        #return (x + 0.6*y + 0.6*z), (y + 0.2*z)
    
    def render(self, c):
        pass


class Vertex(Cell):
    dim = 0
    radius = 0.08
    def render(self, c, deco=[]):
        ps = [self.tran(*coord) for coord in self.coords]
        assert len(ps)==1
        x, y = ps[0]
        c.fill(path.circle(x, y, self.radius), deco+self.deco)

        if argv.debug:
            c.text(x-0.2, y-0.4, "%s"%self.idxs[0], northeast)


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

    def render(self, c, deco=[], edges=True, verts=True, face=True):
        cs = [self.tran(*coord) for coord in self.coords]
        assert len(cs)>2
        ps = [path.moveto(*cs[0])]+[path.lineto(*p) for p in cs[1:]]+[path.closepath()]

        if face:
            c.fill(path.path(*ps), deco+self.deco)
        if verts:
            for x, y in cs:
                c.fill(path.circle(x, y, Vertex.radius)) #, deco)
        if edges:
            c.stroke(path.path(*ps)) #, deco)

#    def birender(self, c, deco1=[], deco2=[], edges=True, verts=True):
#        cs = [self.tran(*coord) for coord in self.coords]
#        assert len(cs) in [4, 6]
#        if len(cs)==4:
#            cs1 = [cs[i] for i in [0, 2, 1, 3]]
#            cs2 = [cs[i] for i in [0, 2, 3, 1]]
#        if len(cs)==6:
#            cs1 = [cs[i] for i in [0, 3, 4, 1, 2, 5]]
#            cs2 = [cs[i] for i in [0, 3, 2, 5, 4, 1]]
#        
#        ps = [path.moveto(*cs[0])]+[path.lineto(*p) for p in cs[1:]]+[path.closepath()]
#        if verts:
#            for x, y in cs:
#                c.fill(path.circle(x, y, Vertex.radius)) #, deco)
#        if edges:
#            c.stroke(path.path(*ps), [style.linejoin.round])
#
#        for (cs, deco) in [(cs1, deco1), (cs2, deco2)]:
#            ps = [path.moveto(*cs[0])]+[path.lineto(*p) for p in cs[1:]]+[path.closepath()]
#            c.fill(path.path(*ps), deco+self.deco)

    def birender(self, c, deco1=[], deco2=[], edges=True, verts=True):
        cs = [self.tran(*coord) for coord in self.coords]
        n = len(cs)
        c0 = (sum(c[0] for c in cs)/n, sum(c[1] for c in cs)/n)

        deco = deco1
        for (i, deco) in [(0, deco1), (1, deco2)]:
            while i < len(cs):
                cs1 = [c0, cs[i], cs[(i+1)%n]]
                ps = [path.moveto(*cs1[0])]+[path.lineto(*p) for p in cs1[1:]]+[path.closepath()]
                c.fill(path.path(*ps), deco+self.deco)
                i += 2

        ps = [path.moveto(*cs[0])]+[path.lineto(*p) for p in cs[1:]]+[path.closepath()]
        if verts:
            for x, y in cs:
                c.fill(path.circle(x, y, 0.7*Vertex.radius)) #, deco)
        if edges:
            c.stroke(path.path(*ps), [style.linejoin.round])



class Body(Cell):
    dim = 3



def mark_stabs(Hx, marks={}):

    marks = dict(marks)
    COLORS = 4

    mx, n = Hx.shape

    #print shortstrx(Hx)

    from pyfinder.expr import Sort, Function, Theory, Constant
    from pyfinder.solver3 import Solver

    color = Sort("color", 4)

    sorts = [color]
    funcs = [Function("c%d"%i, [], color) for i in range(mx)]
    exprs = []
    for i in range(mx):
        if i in marks:
            c = marks[i]
            c = Constant(c, color)
            exprs.append(funcs[i]() == c)
        else:
            for j in range(i+1, mx):
                if (Hx[i]*Hx[j]).sum():
                    e = funcs[i]() != funcs[j]()
                    exprs.append(e)

    theory = Theory(exprs)
    solver = Solver(theory, sorts, funcs)
    print "marks:", marks
    for position in solver.solve(max_count=1, verbose=False):
        for i, idx in enumerate(position):
            assert marks.get(i) in (None, idx)
            marks[i] = idx
    print "marks:", marks
    return marks




def make(Gx, Gz, Hx, Hz, **kw):

    gx, n = Gx.shape

    NAME = "gcolor2" if n==39 else "gcolor"

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

#    if NAME == "gcolor2" and 0:
#        corner_idxs.extend([0, 14, 18, 38])
#        edge_idxs.extend([1,4,9,3,8,13,15,16,17,19,33,37])
#        face_idxs.extend([2,5,6,7,10,11,12])
#        face_idxs.extend([i+19 for i in [1,4,9,3,8,13,15,16,17]])
#        bulk_idxs.extend([i+19 for i in [2,5,6,7,10,11,12]])
#        assert len(corner_idxs+edge_idxs+face_idxs+bulk_idxs)==n
#

    print "graph..."

    for i in range(n):
        w = Hx[:, i].sum()
        assert 1<=w<=4
        items[w-1].append(i)

    assert len(corner_idxs)==CORNERS

    # Find edges between qubits
    # This is defined by two gauge operators intersecting.
    graph = dict((i, []) for i in range(n)) # map qubit -> nbd

    GG = numpy.dot(Gx, Gx.transpose())
    idxs, jdxs = numpy.where(GG)
    for i, j in zip(idxs, jdxs):
        if i==j:
            continue
        g = Gx[i]*Gx[j]
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

    print "coordinates..."

    coords = {} # map qubit -> (x, y, z)

    # corners of a simplex:
    for ii, i in enumerate(corner_idxs):
        coords[i] = Cell.corner_coords[ii]
        
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

    #print coords
    #print cedges.values()

    #match = [corner_idxs[i] for i in [0,1,3]]
    #match.sort()
    cfaces = {}
    for key in [(0,1,2), (0,1,3), (0,2,3), (1,2,3)]:
        key = tuple(corner_idxs[k] for k in key)
        cfaces[key] = []
    assert len(cfaces)==FACES

    if NAME == "gcolor2":
      for i in face_idxs:
        j = i-19
        c3 = corner_idxs[-1]
        if 0<i<18:
            cidxs = (0, 14, 18)
        elif j in (1, 4, 9):
            cidxs = (0, 14, c3)
        elif j in (3, 8, 13):
            cidxs = (0, 18, c3)
        elif j in (15, 16, 17):
            cidxs = (14, 18, c3)
        else:
            assert 0, i

        c0, c1, c2 = cidxs
        cfaces[tuple(cidxs)].append(i)
        eidxss = [cedges[c0, c1], cedges[c0,c2], cedges[c1,c2]] # edges
        djs = [distance(graph, i, eidxs) for eidxs in eidxss]

        cs = [coords[dj[1]] for dj in djs]
        ds = [dj[0] for dj in djs]
        d = 2*(ds[0]+ds[1]+ds[2])
        coord = conv3((ds[1]+ds[2])/d, (ds[0]+ds[2])/d, (ds[0]+ds[1])/d, *cs)
        coords[i] = coord


    else:
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

    #print "cfaces:", cfaces
    for i in bulk_idxs:
        djs = [distance(graph, i, face) for face in cfaces.values()]
        cs = [coords[dj[1]] for dj in djs]
        ds = [dj[0] for dj in djs]
        d = 3*sum(ds)
        coord = conv4((ds[1]+ds[2]+ds[3])/d, (ds[0]+ds[2]+ds[3])/d, (ds[0]+ds[1]+ds[3])/d, (ds[0]+ds[1]+ds[2])/d, *cs)
        coords[i] = coord


    if NAME == "gcolor2":
        _, y0, _ = coords[19]
        x1, y1, z1 = coords[38]
        for i in range(19,38):
            x0, y0, z0 = coords[i-19]
            #y = -0.5*Cell.R
            coords[i] = conv(0.3, (x0, y0, z0), (x1, y1, z1))

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
        face = Face(idxs, coords) #, [red])
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

    print "colors..."

    stab_marks = {}
    ext_mark = {} # map exterior face qubit -> opposite corner mark
    for i, idx in enumerate(corner_idxs):

        h = Hx[:, idx]
        assert h.sum() == 1 # corner touches one stabilizer
        j = numpy.where(Hx[:, idx])[0][0]
        stab_marks[j] = i

        vertex = Cell.get(idx)
        face = vertex.cobdy[0].cobdy[0]
        assert len(face.cobdy)==1
        body = face.cobdy[0]
        body.mark = i # mark corner stabilizer

        for key, idxs in cfaces.items():
          if idx not in key:
            for idx1 in idxs:
                assert ext_mark.get(idx1) is None
                ext_mark[idx1] = i
    #print "ext_mark:", ext_mark

    stab_marks = mark_stabs(Hx, stab_marks)
    #print stab_marks
    for i, mark in stab_marks.items():
        idxs = numpy.where(Hx[i])[0]
        body = Cell.get(*idxs)
        assert body.mark is None or body.mark == mark
        body.mark = mark

    # Exterior faces get a mark
    exterior_faces = []
    interior_faces = []
    front_faces = []
    back_faces = []
    for cell in cells:
        if cell.dim != 2:
            continue
        if len(cell.cobdy)!=1:
            interior_faces.append(cell)
            mark = [c.mark for c in cell.cobdy]
            mark.sort()
            cell.mark = tuple(mark)
            continue
        exterior_faces.append(cell)
        mark = None
        for idx in cell.idxs:
            m = ext_mark.get(idx)
            if m is None:
                continue
            assert mark is None or m==mark, (mark, m)
            mark = m

        assert mark is not None
        cell.ext_mark = mark

        if cell.ext_mark in Cell.BACK:
            front_faces.append(cell)
        elif cell.ext_mark in Cell.FRONT:
            back_faces.append(cell)
        else:
            assert 0, mark

        if mark == cell.cobdy[0].mark:
            print "BUG "*5
        mark = [mark, cell.cobdy[0].mark]
        mark.sort()
        cell.mark = tuple(mark)


    #front_faces = [face for face in exterior_faces if face.mark in Cell.BACK]
    #back_faces = [face for face in exterior_faces if face.mark in Cell.FRONT]

    print "done."

    return Model(locals())


    
class Model(object):
    def __init__(self, attrs):
        self.__dict__.update(attrs)
        self.c = canvas.canvas()
        self.X = 0.
        self.Y = 0.

    def render(self, **kw):

        c = canvas.canvas()
        tr = trafo.scale(sx=1, sy=1)
    
        cells = self.cells
        cells.sort(key = lambda cell : -cell.z)
    
        #deco = [color.transparency(kw.get("transparency", 0.4))]
    
        verts = kw.get("verts", False)
        edges = kw.get("edges", False)
    
        # much brighter without these guys!
        deco = [color.transparency(0.4)]
        for cell in self.back_faces:
            cell.render(c, deco, verts=verts, edges=edges, face=False)
    
        deco = [color.transparency(0.4)]
        for cell in self.interior_faces:
            cell.render(c, deco, verts=verts, edges=edges)
    
        colors = [rgbhex("ff1e1e"), rgbhex("3ee53e"), rgbhex("011fff"), rgbhex("fdff57"), ]
        deco = [color.transparency(0.3)]
        for cell in self.front_faces:
            #cl = colors[cell.mark]
            cl = colors[cell.cobdy[0].mark]
            cell.render(c, deco + [cl], verts=verts, edges=edges)
    
        show_idxs = kw.get("show_idxs", [])
        for idx in show_idxs:
            cell = Cell.get(idx)
            cell.render(c)
    
        if kw.get("label", False):
            c.text(0.1*Cell.R, -1.5*Cell.R, "$n=%d$"%self.n, north)

        self.c.insert(c, [trafo.translate(self.X, self.Y)])
        self.X += 2.3 * Cell.R

    def render_ideal(self, **kw):

        c = canvas.canvas()
        tr = trafo.scale(sx=1, sy=1)
    
        cells = self.cells
        cells.sort(key = lambda cell : -cell.z)
    
        verts = kw.get("verts", False)
        edges = kw.get("edges", False)

        faces = [cell for cell in cells if cell.dim==2]
            
        deco = [rgbhex("e1d9e0")]
            #color.transparency(0.8)]
        for cell in self.back_faces:
            cell.render(c, deco, verts=verts, edges=edges)

        deco = []
        #colors = [rgbhex("bf3f3f"), rgbhex("7fbf7f"), rgbhex("3f3fa5"), rgbhex("dfdf7f")]
        colors = [rgbhex("ff1e1e"), rgbhex("3ee53e"), rgbhex("011fff"), rgbhex("fdff57"), ]
        for cell in faces:
            #print cell.mark
            if cell.mark != (0,3) and cell.mark != (1,2):
                continue
            #cell.render(c, deco, verts=verts, edges=edges)
            deco1 = [colors[cell.mark[0]]]+deco
            deco2 = [colors[cell.mark[1]]]+deco
            cell.birender(c, deco1, deco2, verts=True, edges=edges)
    
        if kw.get("label", False):
            c.text(0.1*Cell.R, -1.5*Cell.R, label, north)

        self.c.insert(c, [trafo.translate(self.X, self.Y)])
        self.X += 2.3 * Cell.R

    def newline(self):
        self.X = 0.
        self.Y -= 2.3 * Cell.R

    def save(self, filename):
        print "saving", filename
        self.c.writePDFfile(filename)
        self.c = canvas.canvas()
        self.X = 0.
        self.Y = 0.



def main():

    size = argv.get("size", 2)
    transparency = argv.get("transparency", 0.4)

    Gx, Gz, Hx, Hz = models.build("gcolor")

    filename = argv.get("filename", "pic-gcolor.pdf")

    model = make(Gx, Gz, Hx, Hz)

    if argv.ideal:
        model.render_ideal(transparency=transparency, label="", verts=False, edges=True,)
    else:
        model.render( transparency=transparency, label=True, verts=True, edges=True,)

    model.save(filename)


    




if __name__=="__main__":
    main()


