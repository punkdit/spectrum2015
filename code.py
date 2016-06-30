#!/usr/bin/env python

import sys

import numpy

from qupy.ldpc.solve import *
from qupy.gauge import lstr2

ra.seed(0)

def reflect(H, k):

    assert k < len(H)
    H1 = H[:k, :k]

    for i in range(k):
        H1[i, i] += H[i, k:].sum()
    return H1


def build_compass(l):

    n = l**2

    keys = [(i, j) for i in range(l) for j in range(l)]
    coords = {}  
    for i, j in keys:
        for di in range(-l, l+1):
          for dj in range(-l, l+1):
            coords[i+di, j+dj] = keys.index(((i+di)%l, (j+dj)%l))

    m = n
    Gx = zeros2(m, n)
    Gz = zeros2(m, n)

    idx = 0
    for i in range(l):
      for j in range(l):
        Gx[idx, coords[i, j]] = 1
        Gx[idx, coords[i, j+1]] = 1

        Gz[idx, coords[i, j]] = 1
        Gz[idx, coords[i+1, j]] = 1
        idx += 1

    assert idx == m

    return Gx, Gz



def test():

    if argv.compass:

        l = argv.get('l', 3)
        Gx, Gz = build_compass(l)

    else:
        m = argv.get("m", 5)
        mx = argv.get("mx", m)
        mz = argv.get("mz", m)
        n = argv.get("n", 10)
    
        Gx = rand2(m, n)
        Gz = rand2(m, n)

    print "Gx:", len(Gx)
    #print shortstr(Gx)

    print "Gz:", len(Gz)
    #print shortstr(Gz)

    G = list(span(Gx))
    #print "<Gx>:", len(G)
    
    states = {} # map shortstr -> idx
    #v = zeros2(n, 1)
    #print v

    G.sort(key = lambda v : dot2(Gz, v).sum())
    k = len(G)

    H = zeros2(k, k)
    for idx, v in enumerate(G):
        sv = shortstr(v)
        states[shortstr(v)] = idx

    #print
    mz = len(Gz)
    _count = None
    for idx, v in enumerate(G):
        count = mz - 2*dot2(Gz, v).sum()
        if count != _count:
            print "idx:", idx, count
            _count = count
        #count = 2*(mx - dot2(Gz, v).sum())
        H[idx, idx] = count
        #print shortstr(v), count
        #print (idx, count),
        for g in Gx:
            v1 = (g+v)%2
            jdx = states[shortstr(v1)]
            H[idx, jdx] += 1
    #print
    print "G:", len(G)

    #print shortstr(H)

    i = argv.truncate
    if i:
        print "truncate:", i
        H = H[:i, :i]

    #for i in range(k):
    i = argv.reflect
    if i:

        print "reflect:", H[i-2, i-2], H[i-1, i-1]
        H = reflect(H, i)
        print "reflect:", len(H)

        #for i in range(len(H)):
        #    print H[i, i],
        #print
        #print shortstr(H1)

        idxs = range(len(H))
        idxs.sort(key = lambda idx : (-H[idx, idx], -H[idx].sum()))
        H = H[:, idxs]
        H = H[idxs, :]

        #for i in range(len(H)):
        #    print H[i, i],
        #print

    u, v = numpy.linalg.eigh(H)
    print "eigval:", u[-1]
    v = numpy.abs(v[:, -1])
    
    print "eigvec:", 
    #_x = 1.0
    for idx, x in enumerate(v):
        print "%.2f"%x,
        for jdx in range(idx):
            if H[idx, jdx]:
                if v[jdx] < v[idx]-1e-6:
                    print "*",
        #if x > _x+1e-4:
        #    print "*", H[idx-1,idx-1], H[idx, idx],
        #_x = x
    print

    if argv.sort:
        idxs = range(len(v))
        idxs.sort(key = lambda idx : -v[idx])
        #print v[idxs]
    
        H = H[:, idxs]
        H = H[idxs, :]
        print shortstr(H)
        
        u, v = numpy.linalg.eigh(H)
        print u
        v = numpy.abs(v[:, -1])
        print v




from qupy.tool.argv import Argv 
argv = Argv()

if __name__ == "__main__":

    test()






