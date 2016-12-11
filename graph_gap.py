#!/usr/bin/env python

import sys, os
from math import *

from pyx import graph



g = graph.graphxy(
    width=14, height=6,
    x=graph.axis.linear(min=0., title="$n$"), 
    y=graph.axis.linear(min=0., title="gap"))

# Compass
xs = [16, 25, 36]
ys = [0.644, 0.452, 0.316]
g.plot(graph.data.values(x=xs, y=ys))
g.plot(graph.data.values(x=xs, y=ys), [graph.style.line()])

# 3D Compass
xs = [27]
ys = [0.538]
g.plot(graph.data.values(x=xs, y=ys))

# gauge color 
xs = [15, 65, 175]
ys = [3.241, 1.694, 1.049]
g.plot(graph.data.values(x=xs, y=ys))
g.plot(graph.data.values(x=xs, y=ys), [graph.style.line()])

g.writePDFfile("pic-gap.pdf")



