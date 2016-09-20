#!/usr/bin/env python

import sys, os


f = open(sys.argv[1])


rows = []

for line in f:
    line = line.strip()
    idxs, val = line.split(" eval: ")
    idxs = eval(idxs)
    val = float(val)

    rows.append((idxs, val))

EPSILON = 1e-8
n = len(rows)

def lt(idxs, jdxs):
    assert idxs != jdxs
    for idx, jdx in zip(idxs, jdxs):
        if idx > jdx:
            return False
    return True


for i in range(n):
  for j in range(i+1, n):

    if lt(rows[i][0], rows[j][0]):
        if rows[i][1] < rows[j][1]-EPSILON:
            print rows[i]
            print rows[j]
            assert 0, (i, j)


print len(rows)
print rows[-1]



    

