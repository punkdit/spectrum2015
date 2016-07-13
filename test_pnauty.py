#!/usr/bin/env python

import sys


import pnauty


n = int(sys.argv[1])

verts = n
edges = n
degree = 2

print "init_graph"
pnauty.init_graph(verts, edges, degree)

for i in range(n):
    pnauty.add_edge(i, (i+1)%n, 0)
    pnauty.add_edge(i, (i-1)%n, 1)

for i in range(n):
    pnauty.set_partition(i, i, 1)


print "search:"

pnauty.search()




