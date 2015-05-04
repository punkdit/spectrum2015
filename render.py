#!/usr/bin/env python

import sys
import math
from random import *
from math import *

from pyx import canvas, path, deco, trafo, style, text, color, unit, epsfile, deformer



rgb = color.rgb
rgbfromhexstring = color.rgbfromhexstring

red, green, blue, yellow = (rgb.red,
    rgbfromhexstring("#008000"),
    rgb.blue, rgb(0.75, 0.75, 0))

blue = rgb(0., 0., 0.8)
lred = rgb(1., 0.4, 0.4)
lgreen = rgb(0.4, 0.8, 0.2)


text.set(mode="latex") 
#text.set(docopt="10pt")
text.preamble(r'\usepackage{amsfonts}')
#text.preamble(r"\def\I{\mathbb{I}}")


c = canvas.canvas()


l = 3
my = 1.2
mx = 0.9

coords = {}
for i in range(l):
  for j in range(l):
   for k in range(2):
    x = 2*j + i
    y = 1.7*i + k
    y *= my
    x *= mx
    coords[i, j, k] = (x, y)
    #for di in (-l, 0, l):
    #  for dj in (-l, 0, l):
    #    coords[i+di, j+dj, k] = (x, y)



radius = 0.45
roff = 0.25
ta = [text.halign.boxcenter, text.valign.middle]
xcoords = {}
ycoords = {}
zcoords = {}
for i in range(l):
  for j in range(l):

    k = 0
    x, y = coords[i, j, k]
    xcoords[i, j, k] = xc, yc = (x, y+roff)
    ycoords[i, j, k] = xc, yc = (x+0.75**0.5*roff, y-0.5*roff)
    zcoords[i, j, k] = xc, yc = (x-0.75**0.5*roff, y-0.5*roff)

    k = 1
    x, y = coords[i, j, k]
    xcoords[i, j, k] = xc, yc = (x, y-roff)
    ycoords[i, j, k] = xc, yc = (x-0.75**0.5*roff, y+0.5*roff)
    zcoords[i, j, k] = xc, yc = (x+0.75**0.5*roff, y+0.5*roff)


la = [style.linewidth(0.3), style.linecap.round]
for i in range(l):
  for j in range(l):

    x0, y0 = xcoords[i, j, 0]
    x1, y1 = xcoords[i, j, 1]
    c.stroke(path.line(x0, y0, x1, y1), la+[lgreen])

    x0, y0 = ycoords[i, j, 0]
    if (i-1, j+1, 0) in ycoords:
        x1, y1 = ycoords[i-1, j+1, 1]
        c.stroke(path.line(x0, y0, x1, y1), la+[lgreen])

    x0, y0 = zcoords[i, j, 0]
    if (i-1, j, 0) in ycoords:
        x1, y1 = zcoords[i-1, j, 1]
        c.stroke(path.line(x0, y0, x1, y1), la+[lgreen])


for i in range(l):
  for j in range(l):

    k = 0
    x, y = coords[i, j, k]
    c.stroke(path.circle(x, y, radius))

    xc, yc = xcoords[i, j, k]
    c.text(xc, yc, "$X$", ta)
    xc, yc = ycoords[i, j, k]
    c.text(xc, yc, "$Y$", ta)
    xc, yc = zcoords[i, j, k]
    c.text(xc, yc, "$Z$", ta)

    k = 1
    x, y = coords[i, j, k]
    c.stroke(path.circle(x, y, radius))

    xc, yc = xcoords[i, j, k]
    c.text(xc, yc, "$X$", ta)
    xc, yc = ycoords[i, j, k]
    c.text(xc, yc, "$Y$", ta)
    xc, yc = zcoords[i, j, k]
    c.text(xc, yc, "$Z$", ta)



c.writeEPSfile("fig_00")
c.writePDFfile("fig_00")

sys.exit(0)



y = h-2*dr

c.text(w/2, y,
    r"{\huge Error Correction in a Fibonacci Fusion Code}",
    [text.halign.boxcenter, text.valign.top, text.size(1)])

y -= 2*dr
c.text(w/2, y,
    r"{\large Simon Burton$^1$, Courtney G. Brell$^{1,2}$, Steven T. Flammia$^1$}",
    [text.halign.boxcenter, text.valign.top])

#c.text(w/2, y,
#    r"{\large Simon Burton{\tiny(1)}, Courtney G. Brell\(1,2\), Steven T. Flammia{(1)}}",
#    [text.halign.boxcenter, text.valign.top])

y -= dr
c.text(w/2, y,
    r"{\it $^1$ Centre for Engineered Quantum Systems, School of Physics, The University of Sydney, Sydney, Australia}",
    [text.halign.boxcenter, text.valign.top])

y -= dr
c.text(w/2, y,
    r"{\it $^2$ Institut f\"ur Theoretische Physik, Leibniz Universit\"at Hannover, Appelstra\ss e 2, 30167 Hannover, Germany}",
    [text.halign.boxcenter, text.valign.top])


#c.text(0, h-10*dr,
#"""We consider an encoded qubit stored in the fusion space
#of four Fibonacci anyons on a lattice.
#We use a phenomenological model of short-time 
#thermalization processes including anyon pair-creation, 
#hopping, exchange, and decoherence to simulate the interaction 
#of this system with a thermal bath. We 
#conduct monte-carlo simulations to investigate the ability of 
#an active error-correction scheme to restore the system 
#to its ground state and preserve the encoded 
#information. 
#The efficient simulation of these processes is 
#not in general expected to be possible on 
#a classical computer, but we demonstrate that in 
#some regimes relevant to error-correction (particularly well below 
#the error-correction threshold), such simulation is indeed possible.""",
#    [text.parbox(dw), ])


pip = 0.08

def rect(x, y, dx, dy, intxt=None, bottxt=None):
    "x, y: lower left"
    c.stroke(path.rect(x, y, dx, dy), [style.linestyle.dashed, style.linewidth.thick, green])
    if intxt:
        c.text(x+dx/2., y+dy/2., intxt,
            [text.valign.middle, text.halign.boxcenter])
    if bottxt:
        c.text(x+dx-pip, y-pip, bottxt,
            [text.valign.top, text.halign.boxright])

bw = 0.8
def triple(x, y, left, right, total, bw=bw, pip=pip):
    rect(x, y, bw, bw, left)
    rect(x+bw+pip, y, bw, bw, right)
    rect(x-pip, y-pip, 2*bw+3*pip, bw+2*pip, bottxt=total)

def arrow(x0, y0, x1, y1, extra=[]):
    c.stroke(path.line(x0, y0, x1, y1),
        extra+[deco.earrow(size=0.2)])


# -------------------------------------------------------------------

# Column 1 ----------------------------------------------------------

x = m
y = h-9*dr

#for i in range(30):
#    y0 = y-i*dr
#    if y0 < 0.:
#        break
#    cross(0, y0)

box = c.text(x, y, r"""
Here we consider a two dimensional system with Fibonacci
anyon excitations. %What does this mean?

The total charge enclosed in a region may be measured with possible
outcomes $\tau$ which indicates a Fibonacci anyon charge, or ${\mathbb{I}}$
which indicates vacuum:
""".strip(), [text.parbox(dw-2*m), ])
#print repr(box.height)
#print unit.tocm(box.height)

y -= 5*dr
rect(x+dw/2-1., y-1., 1., 1., r"$\I$ or $\tau$")

y -= 3.*dr
box = c.text(x, y, r"""
Two such measurements commute when the two regions are either disjoint
or one is wholly contained within the other.

The following rules constrain measurment outcomes:
""".strip(), [text.parbox(dw-2*m), ])

y -= 6*dr
triple(x+m, y, r"$\I$", r"$\I$", r"$\I$")
triple(x+m+3.*bw, y, r"$\I$", r"$\tau$", r"$\tau$")
#triple(x+5*bw, y, r"$\tau$", r"$\I$", r"$\tau$")

y -= 3*dr
triple(x+m, y, r"$\tau$", r"$\tau$", r"$\I$ or $\tau$")

y -= 2*dr
box = c.text(x, y, r"""
The last rule indicates the non-local degree of freedom
associated with two Fibonacci anyons.

Continuing in this way, any region may be decomposed into a 
nested sequence of such pairwise fusions (disjoint regions):
""".strip(), [text.parbox(dw-2*m), ])

y -= 9*dr

N = 5
x0 = x + N*pip + m
y0 = y + N*pip
bw0 = bw - N*pip
for i in range(N):

    rw, rh = bw + i*(bw+2*pip), 1.4*bw-2*(N-i-1)*pip
    txt = r"$\tau$" if i==0 else None
    rect(x0, y0, rw, rh, txt)

    if i < N-1:
        txt = r"$\tau$" if i in [0, 2] else r"$\I$"
        rect(x0+rw+pip, y0, bw-pip, rh, txt)

    if i:
        tpl = (x0+2*rw/3, y-3*pip, x0+2*rw/3, y0)
        arrow(*tpl)
        c.text(tpl[0], tpl[1]-4*pip, ("$a$" if i<=2 else "$b$"),
            [text.valign.bottom, text.halign.boxcenter])

    x0 -= pip
    y0 -= pip


y -= 2.*dr
c.text(x, y,
#$a$ and $b$ denote $\tau$ or $\I$.
r"""$(a,b)$ is either of
$(\tau,\tau), (\tau,\I)\hbox{\ or\ } (\I,\I).$
We use a tree notation:""",
[text.parbox(dw-2*m), ])

#def tree(x, y, dx, dy, shape):
#    "x, y: bottom-left corner"
#    assert type(shape) is tuple
#    c, l, r = shape
#    if type(l) is not tuple and type(r) is not tuple:
#        coords = 
#    elif type(l) is not tuple:
#    elif type(r) is not tuple:


def tree(x, y, dx, dy, branches, leaves, braid=False, braid2=False, clr=None):

    for i, leaf in enumerate(leaves):
        if leaf:
            c.text(x + i*dx, y+pip, leaf, [text.halign.boxcenter])

    ldeco = [style.linewidth.Thick, clr or blue]
    if braid:

        a = 0.4
        c.stroke(path.line(x, y, x+a*dx, y-a*2*dy), ldeco)
        b = 0.6
        c.stroke(path.line(x+b*dx, y-b*2*dy, x+dx, y-2*dy), ldeco)

        c.stroke(path.line(x+dx, y, x, y-2*dy), ldeco)
        y -= 2*dy

    if braid2:

        col = lred
        c.stroke(path.line(x+3*dx, y, x+2*dx, y-2*dy), ldeco)
        c.fill(path.circle(x+2.5*dx, y-1*dy, pip), [color.rgb.white])
        c.stroke(path.line(x+4*dx, y, x+3*dx, y-2*dy), ldeco)
        c.fill(path.circle(x+3.25*dx, y-1.45*dy, pip), [color.rgb.white])
        c.stroke(path.line(x+1*dx, y, x+4*dx, y-2*dy), ldeco)
        c.fill(path.circle(x+1.75*dx, y-0.45*dy, pip), [color.rgb.white])
        c.stroke(path.line(x+2*dx, y, x+1*dx, y-2*dy), ldeco)
        c.stroke(path.line(x, y, x, y-2*dy), ldeco)
        c.stroke(path.line(x+5*dx, y, x+5*dx, y-2*dy), ldeco)

        y -= 2*dy

    for branch in branches:
        row, col, lr = branch[:3]
        x0 = x + col*dx + 0.5*row*dx + lr*dx
        x1 = x + col*dx + 0.5*row*dx + 0.5*dx
        y0 = y -row*dy
        y1 = y -(row+1)*dy
        c.stroke(path.line(x0, y0, x1, y1), ldeco)
        if len(branch)==4:
            extra = [text.valign.top]
            extra += [text.halign.boxleft] if lr else [text.halign.boxright] 
            c.text(x1, y1-pip, branch[3], extra)
                


y -= 2*dr
tau, vac = r"$\tau$ $\I$".split()
tdx, tdy = 0.6, 0.4
tree(x + dw/3, y, tdx, tdy,
    [
        (0,0,0,"$a$"), (0,0,1), (0,1,1), (0,2,1), (0,3,1),
        (1,0,0,"$a$"), (1,0,1), (1,1,1), (1,2,1),
        (2,0,0,"$b$"), (2,0,1), (2,1,1),
        (3,0,0,"$b$"), (3,0,1),
        (4,0,0),
    ], 
    [tau, tau, vac, tau, vac])


#y -= 2.5*dr
#c.text(x, y,
#r""" """,
#[text.parbox(dw-2*m), ])


y -= 5*dr

c.text(x, y,
r"""The so-called F-moves tell us how
to re-associate trees, for example:""",
[text.parbox(dw-2*m),])

y -= 2*dr
tree(x, y, tdx, tdy,
    [
        (0,0,0,tau), (0,0,1), (0,1,1),
        (1,0,0,tau), (1,0,1),
        (2,0,0,tau),
    ], 
    [tau, tau, tau])

c.text(x+m+2.2*tdx, y-1.5*dr, r"$=\ \ c_1$")

tree(x + 4*tdx, y, tdx, tdy,
    [
        (0,1,0), (0,2,0), (0,2,1, tau),
        (1,1,0), (1,1,1, tau),
        (2,0,1, tau),
    ], 
    ['', tau, tau, tau])

c.text(x+m+6.6*tdx, y-1.5*dr, r"$+\ \ c_2$")

tree(x + 8*tdx, y, tdx, tdy,
    [
        (0,1,0), (0,2,0), (0,2,1, tau),
        (1,1,0), (1,1,1, vac),
        (2,0,1, tau),
    ], 
    ['', tau, tau, tau])


y -= 4*dr
c.text(x, y,
r"""Particle exchange (braiding) picks up
a phase, for example:""",
[text.parbox(dw-2*m),])

y -= 2*dr
tree(x, y, tdx, tdy,
    [
        (0,0,0,), (0,0,1),
        (1,0,0,tau),
    ], 
    [tau, tau], braid=True)

c.text(x+1.5*tdx, y-2.*dr, r"$=\ \ c_3$")

tree(x+3*tdx, y-1*dr, tdx, tdy,
    [
        (0,0,0,), (0,0,1),
        (1,0,0,tau),
    ], 
    [tau, tau])

c.text(x+4.5*tdx, y-2.*dr, r",")

tree(x+6*tdx, y, tdx, tdy,
    [
        (0,0,0,), (0,0,1),
        (1,0,0,vac),
    ], 
    [tau, tau], braid=True)

c.text(x+7.5*tdx, y-2.*dr, r"$=\ \ c_4$")

tree(x+9*tdx, y-1*dr, tdx, tdy,
    [
        (0,0,0,), (0,0,1),
        (1,0,0,vac),
    ], 
    [tau, tau])


# -------------------------------------------------------------------

# Column 2

x = m + dw
y = h-9*dr


#y -= 5*dr
c.text(x, y,
r"""We store the qubit of our code in the
fusion space of four anyons:""",
[text.parbox(dw-2*m),])


y -= 3*dr
c.text(x, y-2.*dr, r"$|\psi\rangle \ = \ \alpha$")

tree(x + 3*tdx, y, tdx, tdy,
    [
        (0,0,0,tau), (0,0,1), (0,2,0), (0,2,1,tau), 
        (1,0,0), (1,1,1), 
        (2,0,0), (2,0,1),
        (3,0,0,vac),
    ], 
    [tau, tau, tau, tau], clr=yellow)

c.text(x+6.5*tdx, y-2.*dr, r"$+\ \beta$")

tree(x + 8*tdx, y, tdx, tdy,
    [
        (0,0,0,vac), (0,0,1), (0,2,0), (0,2,1,vac), 
        (1,0,0), (1,1,1), 
        (2,0,0), (2,0,1),
        (3,0,0,vac),
    ], 
    [tau, tau, tau, tau], clr=yellow)


# -------------------------------------------------------------------

# Column 3

x = m + 2*dw
y = h-9*dr

#y -= 5*dr

c.text(x, y,
r"""The fusion tree calculus is
inherently 1-dimensional, so we choose
some ordering for the anyons in our 
2-dimensional system:""",
[text.parbox(dw-2*m),])

y -= 3*dr

from pyx import metapost
from pyx.metapost.path import beginknot, endknot, smoothknot, tensioncurve

x0, y0 = tx0, ty0 = x, y

dx, dy = 1.5, 1.5

# Grid
for i in range(3):
    c.stroke(path.line(x0+0.5*dx, y0-(i+1.0)*dy, x0+3.5*dx, y0-(i+1.0)*dy))
    c.stroke(path.line(x0+(i+1.0)*dx, y0-0.5*dy, x0+(i+1.0)*dx, y0-3.5*dy))

c.stroke(path.rect(1.*dx + pip, -2.*dy + pip, dx-2*pip, dy-2*pip),
    [green, style.linewidth.Thick, style.linestyle.dashed,
    trafo.translate(tx0, ty0)])

def pair((x0, y0), (x1, y1), clr=None):
    extra = [clr or blue, style.linewidth.Thick, style.linecap.round, trafo.translate(tx0, ty0)]
    c.stroke(path.line(x0, y0, x1, y1), extra)
    c.stroke(path.circle(x0, y0, 1.3*pip), extra)
    c.stroke(path.circle(x1, y1, 1.3*pip), extra)


points0 = [
    (1.3*dx, -1.7*dy), (1.3*dx, -2.5*dy),
    (2.6*dx, -1.4*dy), (3.4*dx, -1.4*dy),
    (1.7*dx, -1.7*dy), (1.7*dx, -2.5*dy)]

points = list(points0)

idx = 0 # start
points.insert(idx, (1.3*dx, -0.7*dy)); idx+=3 # next pair
points.insert(idx, (1.3*dx, -2.7*dy)); idx+=1
points.insert(idx, (1.5*dx, -2.7*dy)); idx+=1
points.insert(idx, (1.5*dx, -2.*dy)); idx+=1
points.insert(idx, (1.5*dx, -1.4*dy)); idx+=3 # next pair
points.insert(idx, (3.7*dx, -1.4*dy)); idx+=1
points.insert(idx, (3.7*dx, -1.6*dy)); idx+=1
points.insert(idx, (1.7*dx, -1.6*dy)); idx+=3 # next pair
points.insert(idx, (1.7*dx, -3.3*dy)); # end

items = []
for i, p in enumerate(points):
    if i==0:
        items.append(path.moveto(*p))
    else:
        items.append(path.lineto(*p))

c.stroke(path.path(*items),
    [trafo.translate(tx0, ty0), lred, deformer.smoothed(2.0),
    style.linewidth.THICk, deco.earrow(size=0.3)])

for (s, t) in [(0, 1), (2, 3), (4, 5)]:
    pair(points0[s], points0[t])


y -= 4.*dy
c.text(x, y,
r"""When we yank the red line straight
we find out what braid moves to do:
""",
[text.parbox(dw-2*m),])

# Linear order:
ty0 -= 4.5*dy + 2*dr
dx = 0.9

c.stroke(path.line(0, 0, 7*dx, 0),
    [trafo.translate(tx0, ty0), lred,
    style.linewidth.THICk, deco.earrow(size=0.3)])

pair((1.*dx, 0), (2*dx, 0))
pair((3.*dx, 0), (4*dx, 0))
pair((5.*dx, 0), (6*dx, 0))

points = [
    (0.8, 0.),
    (1.5, 0.4),
    (2.6, 0),
    (3.5, -0.2),
    (4.8, 0.),
    (5.0, 0.2),
    (5.2, 0.),
    (3.5, -0.4),
    (2.3, 0.),
    (1.5, 0.2),
    (1.2, 0.),
    (1., -0.2),
]

items = []
for i, p in enumerate(points):
    px, py = p
    items.append(smoothknot(px*dx, py*dy))
    items.append(tensioncurve())

c.stroke(metapost.path.path(items),
    [trafo.translate(tx0, ty0), green, style.linewidth.Thick,
        style.linestyle.dashed])

y -= 6*dr
c.text(x, y,
r"""
In tree notation:
""",
[text.parbox(dw-2*m),])


y -= 2*dr
tdx = 0.8
tree(x+2*m, y, tdx, tdy,
    [
        (0,0,0,vac), (0,0,1), (0,2,0), (0,2,1,vac), (0,4,0), (0,4,1,vac),
        (1,0,0), (1,1,1), (1,3,1),
        (2,0,0), (2,0,1), (2,2,1),
        (3,0,0,), (3,1,1),
        (4,0,0), (4,0,1),
        (5,0,0),
    ], 
    [tau, tau, tau, tau, tau, tau], braid2=True)



# -------------------------------------------------------------------

# Error correction
# Column 2

x = 2*m + dw
y = h-9*dr

y -= 9*dr
c.text(x, y,
r"""
The algorithm for error correction
works by clustering nearby charges together.
In this case, error correction succeeds:
""",
[text.parbox(dw-2*m),])

y -= 4*dr
grid = {}

class Pair(object):
    def __init__(self, (i0, j0), (i1, j1)):
        self.i0 = i0
        self.j0 = j0
        self.i1 = i1
        self.j1 = j1


#seed(0)

dx, dy = 0.4, -0.5
l = 16

def pair((i0, j0), (i1, j1), clr=None):
    extra = [clr or blue, style.linewidth.Thick, style.linecap.round, trafo.translate(x, y)]
    extra.pop(1) # not Thick ?
    c.stroke(path.line(i0*dx, j0*dy, i1*dx, j1*dy), extra)
    c.stroke(path.circle(i0*dx, j0*dy, 1.3*pip), extra)
    c.stroke(path.circle(i1*dx, j1*dy, 1.3*pip), extra)

def rect(i0, j0, di=1, dj=1, txt=None):
    extra = [green, style.linewidth.Thick, style.linestyle.dashed, trafo.translate(x, y)]
    c.stroke(path.rect((i0-0.5)*dx, (j0-0.5)*dy, di*dx, dj*dy), extra)
    if txt:
        c.text((i0+di-0.5)*dx, (j0+dj-0.5)*dy-pip, txt,
            [trafo.translate(x, y), text.halign.boxright, text.valign.top])

pair((0, l/2), (l/2, l-1), clr=yellow)
pair((l-l/2-1, 0), (l-1, l-l/2-1), clr=yellow)

if 0:
    for _ in range(30):
        i0, j0 = randint(0, l-1), randint(0, l-1)
        if random() < 0.5:
            i1, j1 = i0+1, j0
        else:
            i1, j1 = i0, j0+1
        print "pair((%d, %d), (%d, %d))"%(i0, j0, i1, j1)
        pair((i0, j0), (i1, j1))

else:
    pair((10, 0), (11, 0))
    pair((11, 14), (11, 15))
    pair((1, 12), (1, 13))
    pair((11, 4), (11, 5))
    pair((11, 6), (11, 7))
    pair((12, 14), (13, 14))
    pair((12, 4), (13, 4))
    pair((12, 8), (13, 8))
    pair((13, 12), (14, 12))
    pair((13, 14), (13, 15))
    pair((13, 3), (13, 4))
    pair((13, 3), (14, 3))
    pair((13, 8), (13, 9))
    pair((14, 13), (14, 14))
    pair((14, 8), (15, 8))
    #pair((3, 15), (3, 16))
    pair((4, 12), (4, 13))
    pair((4, 14), (4, 15))
    pair((4, 3), (5, 3))
    pair((4, 8), (5, 8))
    pair((7, 13), (8, 13))
    pair((7, 1), (8, 1))
    pair((8, 11), (8, 12))
    pair((8, 14), (9, 14))
    pair((9, 10), (10, 10))
    pair((9, 14), (9, 15))
    pair((9, 6), (9, 7))
    pair((9, 9), (10, 9))

rect(0, 8, 1, 1, tau)

#rect(8, 15)
rect(7, 11, 3, 5, tau)

rect(7, 0, 2, 2, tau)
rect(10, 0, 2, 1, vac)
rect(4, 3, 2, 1, vac)
rect(4, 8, 2, 1, vac)

#rect(15, 7)
rect(11, 3, 5, 7, tau)
rect(9, 6, 1, 2, vac)
rect(9, 9, 2, 2, vac)

rect(11, 12, 4, 4, vac)
rect(1, 12, 1, 2, vac)
rect(4, 12, 1, 4, vac)


y -= 18*dr
c.text(x, y,
r"""
In this case error correction
fails:
""",
[text.parbox(dw-2*m),])

y -= 2*dr

pair((0, l/2), (l/2, l-1), clr=yellow)
pair((l-l/2-1, 0), (l-1, l-l/2-1), clr=yellow)

#seed(24)

if 0:
    for _ in range(40):
        i0, j0 = randint(0, l-2), randint(0, l-2)
        if random() < 0.5:
            i1, j1 = i0+1, j0
        else:
            i1, j1 = i0, j0+1
        print "pair((%d, %d), (%d, %d))"%(i0, j0, i1, j1)
        pair((i0, j0), (i1, j1))

pair((0, 6), (0, 7))
pair((0, 7), (1, 7))
pair((1, 11), (1, 12))
pair((1, 12), (2, 12))
pair((1, 9), (2, 9))
pair((2, 11), (3, 11))
pair((2, 6), (3, 6))
pair((2, 7), (3, 7))
pair((2, 9), (3, 9))
pair((3, 10), (3, 11))
pair((3, 12), (4, 12))
pair((4, 0), (4, 1))
pair((4, 7), (5, 7))
pair((5, 0), (5, 1))
pair((5, 7), (5, 8))
pair((5, 9), (6, 9))
pair((6, 9), (7, 9))
pair((7, 13), (8, 13))
pair((7, 2), (8, 2))
pair((8, 1), (9, 1))
pair((8, 7), (8, 8))
pair((8, 8), (8, 9))
pair((8, 8), (9, 8))
pair((9, 10), (9, 11))
pair((9, 12), (9, 13))
pair((9, 7), (9, 8))
pair((10, 0), (10, 1))
pair((10, 10), (10, 11))
pair((10, 13), (11, 13))
pair((10, 7), (11, 7))
pair((11, 7), (12, 7))
pair((12, 3), (13, 3))
pair((13, 0), (14, 0))
pair((13, 10), (13, 11))
pair((13, 1), (13, 2))
pair((13, 5), (14, 5))
pair((13, 6), (14, 6))
pair((13, 7), (14, 7))
pair((14, 1), (15, 1))

rect(4, 0, 2, 2, vac)
rect(7, 0, 4, 3, tau)
rect(12, 0, 4, 4, vac)
rect(0, 5, 16, 9, "%s or %s"%(vac, tau))
rect(8, 15, 1, 1, tau)


#x -= 2*m

# -------------------------------------------------------------------

# Column 3 ----------------------------------------------------------

x = m + 2*dw

y = 20*dr

c.text(x, y,
r"""Braiding with Fibonacci anyons
can realize universal quantum computing, so
we don't expect to be able to simulate
every instance in a reasonable amount of time.
Instead we bound the compute time, and any
instance that exceeds this bound is considered
``missing'' data:""",
[text.parbox(dw-2*m),])

y -= 18*dr
try:
    c.insert(epsfile.epsfile(x, y, "qip/fibonacci.eps", width=(dw-2*m)))
except IOError, e:
    print e


# -------------------------------------------------------------------


the_can = canvas.canvas()
the_can.insert(c, [trafo.scale(SCALE)])

the_can.writeEPSfile("poster")
the_can.writePDFfile("poster")

