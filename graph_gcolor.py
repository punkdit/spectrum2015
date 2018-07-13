#!/usr/bin/env python

import sys, os

#import numpy
#from numpy import concatenate

from pyx import canvas, path, deco, trafo, style, text, color, unit, epsfile, deformer

from argv import argv
from solve import *
import models



verbose = argv.verbose

Gx, Gz, Hx, Hz = models.build("gcolor")
#model = models.build_model(Gx, Gz, Hx, Hz)
#print model
#print


gx, n = Gx.shape
print shortstr(Gx)


class Node(object):
    def __init__(self, v):
        self.v = v
        self.s = str(v.data)
        self.Gzv = dot2(Gz, v)
        self.count = sum(self.Gzv)
        self.nbd = []

    def __hash__(self):
        return hash(self.s)

    def __eq__(self, other):
        return self.s == other.s

    def __ne__(self, other):
        return self.s != other.s

    def __str__(self):
        return "Node(%s, %d)"%(self.v, self.count)
    

v = zeros2(n)
node = Node(v)


nodes = {node : node}
ordered = list(nodes.keys())
bdy = set(nodes.keys())


while bdy:
    _bdy = set()
    for node in bdy:
        for g in Gx:
            v = (node.v + g)%2
            n1 = Node(v)
            if n1 not in nodes:
                nodes[n1] = n1
                ordered.append(n1)
                _bdy.add(n1)
            else:
                n1 = nodes[n1]
            node.nbd.append(n1)
    bdy = _bdy

#nodes = list(nodes)
#nodes.sort(key = lambda n : n.count)

rows = {}
nodes = ordered
for node in nodes:
    #print node
    assert len(node.nbd)==gx

    row = rows.get(node.count)

    if row is None:
        row = []
        rows[node.count] = row

    row.append(node)

counts = rows.keys()
counts.sort()
print counts


c = canvas.canvas()
r = 0.05
dx = 1.
dy = 1. / 3
W = 40.
H = 20.
dy = H / max(counts)

for count in counts:
    row = rows[count]

    dx = W / len(row)
    for i, node in enumerate(row):
        x = dx * (i + 0.5)
        y = dy * count
        node.x = x
        node.y = y

for node in nodes:
    x, y = node.x, node.y
    c.fill(path.circle(x, y, r))

    for n in node.nbd:
        c.stroke(path.line(
            x, y, n.x, n.y))

c.writePDFfile("pic-cayley.pdf")





