#!/usr/bin/env python2.7

import sys
from math import *
from random import *

from pyx import canvas, path, deco, trafo, style, text, color, deformer
from pyx.color import rgb


text.set(mode="latex") 
text.set(docopt="12pt")
text.preamble(r"\usepackage{amsmath,amsfonts,amssymb}")



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

st_thick = [style.linewidth.thick]
st_Thick = [style.linewidth.Thick]





c = canvas.canvas()


w = 2.0
h = 4.0

hk = 0.2 * h
hm = 0.6 * h
hr = 0.4 * h

dm = 0.2

def axis(x, y0, y1, label):

    y = (y0+y1)*0.5
    c.stroke(path.line(x, y0, x, y1))
    #c.stroke(path.line(x, y0-0.5*dm, x, y1+0.5*dm))
    c.stroke(path.line(x-0.1, y0-0.5*dm, x+0.1, y0-0.5*dm))
    c.stroke(path.line(x-0.1, y1+0.5*dm, x+0.1, y1+0.5*dm))

    c.fill(path.circle(x, y, 0.20), [white])
    c.text(x, y, label, center)


x, y = 0., 0.

c.stroke(path.rect(x, y, 2*w, hr))
c.stroke(path.line(x+w, y, x+w, y+hr), st_dotted)

c.text(x+0.5*w, y+0.5*hr, "$R$", center)

axis(x-0.5, y, y+hr, "$r$")

y += hr + dm

c.stroke(path.rect(x, y, w-0.5*dm, hm))
c.text(x+0.5*w, y+0.5*hm, "$S$", center)

c.stroke(path.rect(x+w+0.5*dm, y, w-0.5*dm, hm))
c.text(x+1.5*w, y+0.5*hm, "$T$", center)

axis(x-0.5, y, y+hm, "$m$")

y += hm + dm

c.stroke(path.rect(x, y, 2*w, hk))
c.stroke(path.line(x+w, y, x+w, y+hk), st_dotted)
c.text(x+0.5*w, y+0.5*hk, "$L$", center)

axis(x-0.5, y, y+hk, "$k$")

axis(x+2*w+0.5, 0., hr+hm+hk+2*dm, "$n$")

x, y = 0., 0.

dm2 = 0.5*dm
t = Turtle(x-dm2, y, 0.)
t.fwd(hr+dm+hm).right(pi/2, dm2).fwd(w-dm2).right(pi/2, dm2)
t.fwd(hm).left(pi/2, dm2).fwd(w-dm2).right(pi/2, dm2).fwd(hr)
t.right(pi/2, dm2).fwd(2*w).right(pi/2, dm2)
t.stroke()

c.text(1.6*w, -0.38, "$G$", center)

c.writePDFfile("pic-canonical.pdf")



###############################################################################



c = canvas.canvas()


w = 2.0
h = 1.0

hk = 1.0 * h
hm = 2.0 * h
hr = 1.0 * h

dm = 0.2


def mkop(row, col, op):

    y = row * 1.13*h + 0.5*h
    x = col * w + 0.25*w

    for i, opp in enumerate(op):
        dx = i * 0.16*w
        c.text(x+dx, y, "$%s$"%opp, center)


x, y = 0., 0.

c.stroke(path.rect(x, y, 2*w, hr))
c.stroke(path.line(x+w, y, x+w, y+hr), st_dotted)

#c.text(x+0.5*w, y+0.5*hr, "$R$", center)

ops = """
XXII ZIZI 
ZZZZ IIIX
XXXX ZZZI
XIXI ZZII
""".strip().split()

for col in [0,1]:
 for row in range(4):
    mkop(row, col, ops[col+2*row])

#axis(x-0.5, y, y+hr, "$r$")

y += hr + dm

c.stroke(path.rect(x, y, w-0.5*dm, hm))
#c.text(x+0.5*w, y+0.5*hm, "$S$", center)

c.stroke(path.rect(x+w+0.5*dm, y, w-0.5*dm, hm))
#c.text(x+1.5*w, y+0.5*hm, "$T$", center)

#axis(x-0.5, y, y+hm, "$m$")

y += hm + dm

c.stroke(path.rect(x, y, 2*w, hk))
c.stroke(path.line(x+w, y, x+w, y+hk), st_dotted)
#c.text(x+0.5*w, y+0.5*hk, "$L$", center)

#axis(x-0.5, y, y+hk, "$k$")

#axis(x+2*w+0.5, 0., hr+hm+hk+2*dm, "$n$")

x, y = 0., 0.

dm2 = 0.5*dm
t = Turtle(x-dm2, y, 0.)
t.fwd(hr+dm+hm).right(pi/2, dm2).fwd(w-dm2).right(pi/2, dm2)
t.fwd(hm).left(pi/2, dm2).fwd(w-dm2).right(pi/2, dm2).fwd(hr)
t.right(pi/2, dm2).fwd(2*w).right(pi/2, dm2)
t.stroke()

c.text(-0.5, 0.5*hr, "$R$", center)
c.text(-0.5, hr+dm+0.5*hm, "$S$", center)
c.text(2*w+0.5, hr+dm+0.5*hm, "$T$", center)
c.text(-0.5, hr+dm+hm+dm+0.5*hk, "$L$", center)

c.text(1.6*w, -0.38, "$G$", center)

c.writePDFfile("pic-gauge4.pdf")



