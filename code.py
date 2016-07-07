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


def build_cube(n):

    m = n
    Gx = zeros2(m, n)
    Gz = zeros2(m, n)

    for i in range(m):
        Gx[i, i] = 1
        Gz[i, i] = 1

    return Gx, Gz



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


EPSILON = 1e-8

def lstr2(xs, decimals=4):
    m, n = xs.shape
    lines = [lstr(xs[i, :], decimals) for i in range(m)]
    return '[%s]'%',\n '.join(line for line in lines)


def lstr(xs, decimals=4, nozero=False):
    if len(xs.shape)>1:
        return lstr2(xs, decimals)
    rxs = [] 
    for x in xs:
        if abs(x.imag)<EPSILON:
            x = x.real
        if abs(x-round(x))<EPSILON:
            x = int(round(x))
        x = "%.*f"%(decimals, x)
        #if x.replace('0', '')=='.':
        #    x = '0.'
        rxs.append(x)
    s = '[%s]'%', '.join(x.rjust(decimals+2) for x in rxs) 
    if nozero:
        s = s.replace(" 0,", "  ,")
        s = s.replace(" 0]", "  ]")
    return s


def texstr(H, align='r'):
    s = lstr(H, 0)
    s = s.replace(', ', ' & ')
    s = s.replace(',\n', ' \\cr\n')
    s = s.replace('[', '')
    s = s.replace(']', '')
    s = r"""
\left(\begin{array}{%s}
%s
\end{array}\right)
    """ % (align*H.shape[1], s)
    return s



gcolor_gauge = """ 
1111...........
11..11.........
1.1.1.1........
..11..11.......
.1.1.1.1.......
....1111.......
11......11.....
1.1.....1.1....
........1111...
..11......11...
.1.1.....1.1...
1...1...1...1..
........11..11.
.1...1...1...1.
....11......11.
........1.1.1.1
..1...1...1...1
....1.1.....1.1
"""

gcolor_stab = """
11111111.......
1111....1111...
11..11..11..11.
1.1.1.1.1.1.1.1
"""




def test():

    if argv.compass:

        l = argv.get('l', 3)
        Gx, Gz = build_compass(l)

    elif argv.gcolor:

        from isomorph import parse
        Gx = parse(gcolor_gauge)
        Gz = Gx.copy()

    elif argv.cube:
        n = argv.get('n', 3)
        Gx, Gz = build_cube(n)

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
    print "G:", len(G)
    
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

    from isomorph import write
    if argv.isomorph:

        depth = argv.get('depth', 1)

        import isomorph
        bag0 = isomorph.from_ham(H)
        bag1 = isomorph.from_ham(H)
        count = 0
        for fn in isomorph.isos(bag0, bag1, depth=depth):
#            #write('.')
#            print [fn[i] for i in range(len(fn))]
#            for i in range(len(fn)):
#                assert bag0[i].get_desc(depth) == bag1[fn[i]].get_desc(depth)
#                print '\t', bag0[i].get_desc(depth),
#                print len(bag0[i].nbd)
#                print '\t', bag1[fn[i]].get_desc(depth)
#            print
            count += 1
        print
        print count

    if len(H)<10:
        #print texstr(H)
        print shortstr(H)

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


    if argv.eigs:

        u, v = numpy.linalg.eigh(H)
        print "eigs:", u[-20:]
        
#        for i in range(20):
#            print "eigval: %.6f" % u[i],
#            evec = v[:, i]
#            #print "eigvec:", 
#            #for idx, x in enumerate(evec):
#            #    print "%.2f"%x,
#                #for jdx in range(idx):
#                #    if H[idx, jdx]:
#                #        if evec[jdx] < evec[idx]-1e-6:
#                #            print "*",
#            print

    if 0:
        A = numpy.array([[3,3,0,0],[1,1,2,0],[0,2,-1,1],[0,0,3,-3]])
        print texstr(A)
        u, v = numpy.linalg.eig(A)
        for i in range(len(A)):
            print u[i], v[:, i]

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






