#!/usr/bin/env python

import sys, os
from math import *

import numpy
from numpy import concatenate

#from solve import get_reductor, array2, row_reduce, dot2, shortstr, zeros2, shortstrx, eq2
#from solve import u_inverse, find_kernel, find_logops, identity2, solve
#from solve import check_conjugate, check_commute

from isomorph import write
import models
from models import genidx

from pyx import graph


def main():

    model = models.build_model()

    N = argv.get("N", 10)
    xs = []
    ys1 = []
    ys2 = []
    for i in range(N):
        theta = 0.5*pi*i/(N-1)

        Jz = cos(theta)
        Jx = sin(theta)
        H = model.build_ham(Jx=Jx, Jz=Jz)

        vals, vecs = numpy.linalg.eigh(H)

        vals = list(vals)
        vals.sort(reverse=True)
        print vals[:5]

        xs.append(theta)
        ys1.append(vals[0])
        ys2.append(vals[1])

    g = graph.graphxy(width=16, x=graph.axis.linear(), y=graph.axis.linear())
    g.plot(graph.data.values(x=xs, y=ys1))
    g.plot(graph.data.values(x=xs, y=ys2))
    #g.plot(graph.data.points(list(zip(range(10), range(10))), x=1, y=2))
    g.writePDFfile("pic-phase.pdf")



from argv import Argv

argv = Argv()

if __name__ == "__main__":

    fn = argv.next()
    if argv.profile:
        import cProfile as profile
        profile.run("%s()"%fn)

    else:
        fn = eval(fn)
        fn()


