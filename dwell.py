#!/usr/bin/env python

import sys
from math import *
from random import *

from pyx import canvas, path, deco, trafo, style, text, color, deformer
from pyx.color import rgb

import numpy
from code import lstr2


text.set(mode="latex") 
text.set(docopt="12pt")
text.preamble(r"\usepackage{amsmath,amsfonts,amssymb}")

text.preamble(r"\def\ket #1{|#1\rangle}")


def dopath(points, extra=[], fill=[], closepath=False, smooth=0.0):
    ps = [path.moveto(*points[0])]+[path.lineto(*p) for p in points[1:]]
    if closepath:
        ps.append(path.closepath())
    p = path.path(*ps)
    if fill:
        c.fill(p, [deformer.smoothed(smooth)]+extra+fill)
    c.stroke(p, [deformer.smoothed(smooth)]+extra)


class Turtle(object):
    def __init__(self, x, y, theta):
        self.x = x
        self.y = y
        self.theta = theta
        self.ps = [(x, y)]

    def fwd(self, d):
        self.x += d*sin(self.theta)
        self.y += d*cos(self.theta)
        self.ps.append((self.x, self.y))
        return self

    def reverse(self, d):
        self.fwd(-d)
        return self

    def right(self, dtheta, r=0.):
        theta = self.theta
        self.theta += dtheta
        if r==0.:
            return self
        N = 20
        x, y = self.x, self.y
        x0 = x - r*sin(theta-pi/2)
        y0 = y - r*cos(theta-pi/2)
        for i in range(N):
            theta += (1./(N))*dtheta
            x = x0 - r*sin(theta+pi/2)
            y = y0 - r*cos(theta+pi/2)
            self.ps.append((x, y))
        self.x = x
        self.y = y
        return self

    def left(self, dtheta, r=0.):
        self.right(-dtheta, -r)
        return self

    def stroke(self, deco=[], closepath=False):
        dopath(self.ps, deco, closepath)
        return self


black = rgb(0., 0., 0.) 
blue = rgb(0., 0., 0.8)
lred = rgb(1., 0.4, 0.4)
white = rgb(1., 1., 1.) 

#shade = rgb(0.75, 0.55, 0)

grey = rgb(0.75, 0.75, 0.75)
shade = grey
shade0 = rgb(0.75, 0.75, 0.75)
shade1 = rgb(0.80, 0.80, 0.80)
shade2 = rgb(0.85, 0.85, 0.85)

light_shade = rgb(0.85, 0.65, 0.1)
light_shade = rgb(0.9, 0.75, 0.4)


north = [text.halign.boxcenter, text.valign.top]
northeast = [text.halign.boxright, text.valign.top]
northwest = [text.halign.boxleft, text.valign.top]
south = [text.halign.boxcenter, text.valign.bottom]
southeast = [text.halign.boxright, text.valign.bottom]
southwest = [text.halign.boxleft, text.valign.bottom]
east = [text.halign.boxright, text.valign.middle]
west = [text.halign.boxleft, text.valign.middle]
center = [text.halign.boxcenter, text.valign.middle]


st_dashed = [style.linestyle.dashed]
st_dotted = [style.linestyle.dotted]
st_round = [style.linecap.round]

st_thick = [style.linewidth.thick]
st_Thick = [style.linewidth.Thick]
st_THICK = [style.linewidth.THICK]
print dir(style.linewidth)

def arc(*args):
    return path.path(path.arc(*args))

def arcn(*args):
    return path.path(path.arcn(*args))

###############################################################################


#n = int(sys.argv[1])
n = 12


H = numpy.zeros((n, n))
E = numpy.zeros(n)

H[0, 0] = 2.
if "double" in sys.argv or 1:
    H[n-1, n-1] = 2.

for i in range(n):
    if i>0:
        H[i, i-1] = 1
    if i+1<n:
        H[i, i+1] = 1

#print H

vals, vecs = numpy.linalg.eigh(H)

#print vals
a = a0 = vals[n-1]
a1 = vals[n-2]
print "a0:", a0
print "a1:", a1

v = v0 = vecs[:, n-1]
if v0[0] < 0.:
    v0 *= -1
v1 = vecs[:, n-2]
if v1[0] < 0.:
    v1 *= -1

U = H.copy()
for i in range(n):
    U[i, i] -= a

for i in range(1, n):
    c = U[i, i-1]
    U[i] -= (c/U[i-1, i-1])*U[i-1]

#print lstr2(U)
for i in range(3):
    print U[i, i], U[i, i]+a

w = numpy.dot(U, v)
assert numpy.abs(w).sum() < 1e-8

#print -a+2
#print -a - 1./(-a+2)
#print -a - 1./(-a-1./(-a+2))

E[0] = 2
for i in range(1, n):
    E[i] = 1./(a-E[i-1])
#print -a+E[0]
#print -a+E[1]
#print -a+E[2]
print E


for i in range(1, n):
    print i, v[i]/v[i-1], v1[i]/v1[i-1]




c = canvas.canvas()

dx = 1.0
dy = 3.0
x, y = 0., 0.

mx = 0.5*dx
c.stroke(path.line(x-1.*dx, y, x+(n+1)*dx, y), 
    [deco.earrow(size=0.3)])
c.stroke(path.line(x-mx, y-0.8*dy, x-mx, y+0.8*dy), 
    [deco.earrow(size=0.3)])

x1 = x + 0.5*(n-1)*dx
c.stroke(path.line(x1, y-0.5*dy, x1, y+0.5*dy),
    st_dotted+st_thick)
c.text(x1+0.2, y+0.5*dy, r"\emph{cut}", west)

for i in range(n):

    x1 = x+i*dx
    c.fill(path.circle(x1, y, 0.08))
    c.text(x1, y-0.4, r"$\ket{%d}$"%(i+1), north)

thick = [style.linewidth.THIck]
pts = [(x+i*dx, v0[i]*dy) for i in range(n)]
dopath(pts, extra=thick)

c.text(x+i*dx-0.4, v0[i]*dy, r"$\ket{v_1}$", east)

green = rgb(0.2, 0.8, 0.2)
red = rgb(1.0, 0.2, 0.2)
pts = [(x+i*dx, v1[i]*dy) for i in range(n)]
dopath(pts, extra=st_dashed+[red]+thick)

c.text(x+i*dx-0.4, v1[i]*dy, r"$\ket{v_2}$", east)

c.writePDFfile("pic-dwell.pdf")



