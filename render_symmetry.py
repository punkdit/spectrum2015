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



from qupy.ldpc.solve import *

if 0:
    n = 9
    Sx = """
    11.......
    ...11....
    ......11.
    .11......
    ....11...
    .......11
    1.1......
    ...1.1...
    ......1.1
    """
    Sx = parse(Sx)
    print Sx
    
    Sz = """
    1..1.....
    .1..1....
    ..1..1...
    ...1..1..
    ....1..1.
    .....1..1
    1.....1..
    .1.....1.
    ..1.....1
    """
    Sz = parse(Sz)
    print Sz


#Sz = []
#for gx in Sx:
#    gz = gx.copy()
#    gz.shape = (3,3)
#    gz = gz.transpose().copy()
#    gz.shape = (n,)
#    Sz.append(gz)
#Sz = array2(Sz)

if 0:
    n = 6
    Sx = """
    11....
    ..11..
    ....11
    """
    Sx = parse(Sx)
    print Sx
    
    Sz = """
    1.1...
    .1.1..
    ..1.1.
    ...1.1
    1...1.
    .1...1
    """
    Sz = parse(Sz)
    print Sz



if 1:
    n = 6
    Sz = """
    11....
    ..11..
    ....11
    """
    Sz = parse(Sz)
    print Sz
    
    Sx = """
    1.1...
    .1.1..
    ..1.1.
    ...1.1
    1...1.
    .1...1
    """
    Sx = parse(Sx)
    print Sx


print shortstr(Sz)


#states = set(['.'*n])
states = ['.'*n]
for _ in range(3):
  for g in Sx:
    for s in list(states):
        op = parse(s)
        op.shape = (n,)
        op = (g+op)%2
        #states.add(shortstr(op))
        op = shortstr(op)
        if op not in states:
            states.append(op)

#print states
#assert len(states)==64

states = list(states)
states = [parse(v) for v in states]
def weight(v):
    v.shape = (n,)
    w = dot2(Sz, v).sum()
    w = len(Sz) - 2*w
    return w

#states.sort(key = lambda v : -weight(v))

states.sort(key = lambda v : str(v))

#i = 2
#while i+1 < len(states):
#    states[i], states[i+1] = states[i+1], states[i]
#    i += 1


c = canvas.canvas()

dx, dy = 2., 6.
radius = 0.6


def render(v, x0, y0, dx=0.25, dy=0.25):
    #c.fill(path.circle(x, y, radius), [color.rgb.white])
    #c.stroke(path.circle(x, y, radius))
    v = v.view()
    m, n = 3, 2
    dm = 0.1
    c.fill(path.rect(x0-dm-dx*m/2., y0-dm-dy*n/2., m*dx+2*dm, n*dy+2*dm), [color.rgb.white])
    v.shape = (m, n)
    x0 -= dx*m/2.
    y0 -= dy*n/2.
    #c.fill(path.rect(x0-dm, y0+dm, m*dx+2*dm, n*dy+2*dm), [color.rgb.white])
    #c.stroke(path.rect(x0, y0, m*dx, n*dy))
    for i in range(m):
     for j in range(n):
      p = path.rect(x0+i*dx, y0+j*dy, dx, dy)
      if v[i,j]==1:
       c.stroke(p)
       c.fill(p)
      else:
       #c.fill(p, [color.rgb.white])
       c.stroke(p)


nodes = {}
layout = {}
N = len(states)
R0 = 4.0
R1 = 4.7
for i, state in enumerate(states):
    print "state:", state
    theta = 2*pi*i/N
    x = R0*sin(theta)
    y = R0*cos(theta)
    nodes[str(state)] = x, y
    x = R1*sin(theta)
    y = R1*cos(theta)
    layout[str(state)] = x, y


for v0 in states:
  for u in Sx:
    v1 = (v0+u) % 2
    a = nodes[str(v0)]
    b = nodes[str(v1)]
    c.stroke(path.line(a[0], a[1], b[0], b[1]))


for v in states:
    x, y = nodes[str(v)]
    c.fill(path.circle(x, y, 0.08), [color.rgb.white])
    c.fill(path.circle(x, y, 0.05))

    x, y = layout[str(v)]
    render(v, x, y)


c.writePDFfile("fig_symmetry.pdf")

