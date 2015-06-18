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


states = set(['.'*n])
for _ in range(3):
  for g in Sx:
    for s in list(states):
        op = parse(s)
        op.shape = (n,)
        op = (g+op)%2
        states.add(shortstr(op))

#print states
#assert len(states)==64

states = list(states)
states = [parse(v) for v in states]
def weight(v):
    v.shape = (n,)
    w = dot2(Sz, v).sum()
    w = len(Sz) - 2*w
    return w

states.sort(key = lambda v : -weight(v))

#print states
ws = set([weight(v) for v in states])
ws = list(ws)
ws.sort(reverse=True)

rows = []
for w in ws:
    row = [v for v in states if weight(v)==w]
    rows.append(row)

order = {} # map str(v) -> left-right order
for i, v in enumerate(rows[0]):
    for u in Sx:
        v1 = (u+v)%2
        key = str(v1)
        if key not in order:
            order[key] = i
rows[1].sort(key = lambda v : order[str(v)])

if len(rows)>2:
    order = dict((str(v), []) for v in states) # map str(v) -> left-right order
    for i, v in enumerate(rows[1]):
        for u in Sx:
            v1 = (u+v)%2
            key = str(v1)
            #i0 = order.get(key, [])
            order[key].append(1.*i)
    rows[2].sort(key = lambda v : sum(order[str(v)]))

c = canvas.canvas()

dx, dy = 2., 6.
radius = 0.6

W = max([dx*len(row) for row in rows])

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

layout = {} # map str(v) -> (x, y)
circ = {}
y = dy * len(rows)
Y0 = dy * (len(rows)+1)
X0 = W/2.
for i, row in enumerate(rows):
    #x = (W - dx*len(row))/2.
    x = 0.
    scalex = W / (dx*len(row))
    for j, v in enumerate(row):
        key = str(v)
        layout[key] = (x, y)
        r0 = (i+1)*dy
        theta = 0.6*pi + 0.8 * pi * 1.*j / (len(row)-1)
        circ[key] = (-r0*sin(theta), 0.8*r0*cos(theta)-0.6)
        if i==0:
            circ[key] = (-r0*sin(theta), r0*cos(pi/2))
        x += dx*scalex
    y -= dy

layout = circ

for v0 in states:
  for u in Sx:
    v1 = (v0+u) % 2
    a = layout[str(v0)]
    b = layout[str(v1)]
    c.stroke(path.line(a[0], a[1], b[0], b[1]))


for v in states:
    x, y = layout[str(v)]
    render(v, x, y)


c.writePDFfile("fig_compass")

