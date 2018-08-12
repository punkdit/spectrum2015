#!/usr/bin/env python3

import sys, os

import numpy
from scipy import stats

import matplotlib.pyplot as pyplot

from argv import argv

f = open(argv.next())
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

if argv.plot:
    pyplot.plot(xs, ys, 'bo')
    pyplot.show()

xs = numpy.array(xs)
ys = numpy.array(ys)

N = len(ys)
        
slope, intercept, r_value, p_value, std_err = stats.linregress(xs, ys)
print("samples:", N)
print("slope:", slope)
print("total slope:", slope*N)

ys -= ys.mean()

print("range:", ys.min(), ys.max())

rs = []
for i in range(1, 200):

    r = 0.
    for j in range(0, N-i):
        r += ys[j] * ys[j+i]
    r /= N
    rs.append(r)
#    print("%.1f"%(r/N), end=" ")
#print()

pyplot.plot(rs)
pyplot.show()


