#!/usr/bin/env python

import sys, os
from math import *

from pyx import graph



g = graph.graphxy(
    width=14, height=6,
    x=graph.axis.linear(min=0., title=r"qubits: $n$"), 
    y=graph.axis.linear(min=0., title=r"gap: $\lambda_1-\lambda_2$"))


# XY
# note: n=25 has no stabilizers: 31.851942-31.349618
xs = [16, 24]
ys = [20.503324-20.109358, 30.645190-30.383016]
g.plot(graph.data.values(x=xs, y=ys))
g.plot(graph.data.values(x=xs, y=ys), [graph.style.line()])


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



