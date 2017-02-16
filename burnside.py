#!/usr/bin/env python

r"""
Calculate the Burnside Ring corresponding to a particular G-set.

Klein geometry: 
(1) a set of "points" and a group G acting on the set.
(2) (conjugacy classes of) subgroups of G are the "geometric properties" 

For example, a triangle.
This has symmetry group of order 6, and lattice
of subgroups:

      6
     / \
    /   \
    3    2
    \   /
     \ /
      1

We name the corresponding properties as:

      Frames
     / \
    /   \
Points  Orientations
    \   /
     \ /
      Nothing



"""


from geometry import Geometry
import isomorph
from action import Perm, Action

from argv import argv


def main():

    n = argv.get("n", 3)
    poly = Geometry.polygon(n)

    bag0 = poly.get_bag()
    bag1 = poly.get_bag()

    items = range(len(bag0))

    perms = [] # group of perms
    for f in isomorph.search(bag0, bag1):
        print f
        perm = Perm(f, items)
        perms.append(perm)
    G = Action(perms, items)

    print "order(G):", len(G)

    orbits = G.orbits()

    print "orbits:", len(orbits)

    for orbit in orbits:

        orbit = list(orbit)
        perms = [g.restrict(orbit) for g in G]
        G1 = Action(perms, orbit)



if __name__ == "__main__":

    main()



