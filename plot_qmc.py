#!/usr/bin/env python3

import sys, os

import matplotlib.pyplot as pyplot

f = open(sys.argv[1])
ys = []
xs = []
i = 0
for line in f:
    line = line.strip()
    if not line:
        break

    ys.append(float(line))
    xs.append(i)
    i += 1
        

#pyplot.plot(xs, ys, 'bo')
pyplot.plot(xs, ys, 'bo')
pyplot.show()


