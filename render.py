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


l = 4
my = 1.2
mx = 0.9

coords = {}
for i in range(-1, l):
  for j in range(-1, l):
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


for i in range(l-1):

    x0, y0 = coords[i, 0, 0]
    x1, _ = coords[i, l-2, 0]

    c.stroke(path.line(x0-3*radius, y0, x1+0*radius, y0), [style.linestyle.dashed])
    c.text(x0-3*radius, y0+0.1, "$i=%s$"%(i+1))

    c.text(x1+2*radius, y0+0.1, "$k=1$", [text.valign.top])

    _, y0 = coords[i, 0, 1]
    #c.stroke(path.line(x0, y0, x1+3*radius, y0), [style.linestyle.dashed])
    c.text(x1+2*radius, y0+0.1, "$k=2$", [text.valign.top])



for j in range(l-1):
    x0, y0 = coords[0, j, 0]
    x1, y1 = coords[l-1, j, 0]

    c.stroke(path.line(x0, y0, x1, y1), [style.linestyle.dotted])
    c.text(x1-1.0, y1-0.1, "$j=%s$"%(j+1))


for i in range(l-1):
  for j in range(l-1):

    k = 0
    x, y = coords[i, j, k]
    c.fill(path.circle(x, y, 1.2*radius), [color.rgb.white])

    k = 1
    x, y = coords[i, j, k]
    c.fill(path.circle(x, y, 1.2*radius), [color.rgb.white])


xcoords = {}
ycoords = {}
zcoords = {}
for i in range(l):
  for j in range(l-1):

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
for i in range(l-1):
  for j in range(l-1):

    x0, y0 = xcoords[i, j, 0]
    x1, y1 = xcoords[i, j, 1]
    c.stroke(path.line(x0, y0, x1, y1), la+[lgreen])

    x0, y0 = ycoords[i, j, 0]
    if (i-1, j+1, 1) in ycoords:
        x1, y1 = ycoords[i-1, j+1, 1]
        c.stroke(path.line(x0, y0, x1, y1), la+[lgreen])

    x0, y0 = zcoords[i, j, 0]
    if (i-1, j, 1) in zcoords:
        x1, y1 = zcoords[i-1, j, 1]
        c.stroke(path.line(x0, y0, x1, y1), la+[lgreen])


for i in range(l-1):
  for j in range(l-1):

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



c.writePDFfile("fig_00")


c = canvas.canvas()


X, Y, Z = "$X$ $Y$ $Z$".split()

xoff = 0.25
yoff = 0.05

op0 = [
    (0, 0, 0, X),
    (0, 0, 1, Y),
    (1, 0, 0, Y),
    (1, 0, 1, Y),
    (2, 0, 0, Y),
    (2, 0, 1, Y),
    (3, 0, 0, X),
    (2, 1, 1, X),
    (3, 1, 0, X),
    (2, 2, 1, X),
    (3, 2, 0, X),
    (2, 3, 1, Y),
]

for idx in range(len(op0)-1):
    x0, y0 = coords[op0[idx][:3]]
    x1, y1 = coords[op0[idx+1][:3]]
    c.stroke(path.line(x0-xoff, y0+yoff, x1-xoff, y1+yoff), la+[lgreen])

for idx in range(len(op0)):
    i, j, k, opc = op0[idx]
    x, y = coords[i, j, k]
    c.text(x-xoff, y+yoff, opc, ta)


op1 = [
    (0, 0, 0, X),
    (0, 0, 1, Y),
    (1, 0, 0, Y),
    (1, 0, 1, Y),
    (2, 0, 0, X),
    (1, 1, 1, X),
    (2, 1, 0, X),
    (1, 2, 1, X),
    (2, 2, 0, X),
    (1, 3, 1, X),
    (2, 3, 0, Z),
]

for idx in range(len(op1)-1):
    x0, y0 = coords[op1[idx][:3]]
    x1, y1 = coords[op1[idx+1][:3]]
    c.stroke(path.line(x0+xoff, y0-yoff, x1+xoff, y1-yoff), la+[lgreen])

for idx in range(len(op1)):
    i, j, k, opc = op1[idx]
    x, y = coords[i, j, k]
    c.text(x+xoff, y-yoff, opc, ta)


for i in range(l):
  for j in range(l):

    if (i,j)==(l-1,l-1):
        continue

    k = 0
    x, y = coords[i, j, k]
    c.stroke(path.circle(x, y, radius))

    if i==l-1:
        continue

    k = 1
    x, y = coords[i, j, k]
    c.stroke(path.circle(x, y, radius))

x, y = coords[l-2,l-1,0]
c.text(x+1.2*radius, y-2.*radius, "$g_{341}$", ta)
x, y = coords[l-2,l-1,1]
c.text(x+1.0*radius, y+2.*radius, "$g_{342}$", ta)

c.writePDFfile("fig_01")




