#!/usr/bin/env python

import sys, os


f = open(sys.argv[1])


rows = []

for line in f:
    line = line.strip()
    TOK = " -- "
    if TOK not in line:
        continue
    flds = line.split(TOK)
    assert len(flds) == 3
    weight, idxs, val = flds
    weight = int(weight)
    idxs = eval(idxs)
    val = float(val)

    rows.append((weight, idxs, val))

EPSILON = 1e-8
n = len(rows)
print len(rows), "rows"

def lt(idxs, jdxs):
    assert idxs != jdxs
    for idx, jdx in zip(idxs, jdxs):
        if idx > jdx:
            return False
    return True


for i in range(n):
  for j in range(i+1, n):

    if lt(rows[i][1], rows[j][1]):
        if rows[i][2] < rows[j][2]-EPSILON:
            print "%6d"%i, rows[i]
            print "%6d"%j, rows[j]
            print
            #assert 0, (i, j)
            assert sum(rows[i][1])>1





    

