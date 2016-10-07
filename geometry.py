#!/usr/bin/env python

import os, sys



class Geometry(object):

    def __init__(self, incidence, tpmap):
        """
            incidence: list of item pairs (i, j)
            tpmap: dict mapping each item to it's type
        """
        items = set() # all items
        nbd = {} # map item -> list of incident items
        tplookup = {} # map type -> list of items of that type
        for i, j in incidence:
            items.add(i)
            items.add(j)
        for i in items:
            nbd[i] = []
        for i, j in incidence:
            nbd[i].append(j)
            nbd[j].append(i)

        for i in items:
            if i not in nbd[i]:
                nbd[i].append(i)
            for j in nbd[i]:
                if i==j:
                    continue

                assert tpmap[i] != tpmap[j]
                if i not in nbd[j]:
                    nbd[j].append(i)

        for t in tpmap.values():
            tplookup[t] = []
        for item, t in tpmap.items():
            tplookup[t].append(item)
        for jtems in nbd.values():
            assert len(set(jtems))==len(jtems) # uniq

        items = list(items)
        items.sort()

        incidence = [] # rebuild this
        for item, jtems in nbd.items():
            for jtem in jtems:
                incidence.append((item, jtem))
        incidence.sort()
        self.incidence = list(incidence) # incidence relation: list of pairs (i, j)
        self.tpmap = dict(tpmap) # map item -> type
        types = set(tpmap.values())
        self.rank = len(types)
        self.items = items
        self.nbd = nbd

        self.check()

    def __str__(self):
        return "Geometry(%s, %s)"%(self.incidence, self.tpmap)

    def ordered_flags(self, flag=[]):
        #yield flag
        nbd = set(self.items)
        for i in flag:
            nbd = nbd.intersection(self.nbd[i])
        #print "all_flags:", flag, nbd
        for i in nbd:
            if i in flag:
                continue
            yield flag+[i]
            for _flag in self.ordered_flags(flag+[i]):
                yield _flag

    def all_flags(self):
        flags = [flag for flag in self.ordered_flags()]
        for flag in flags:
            flag.sort()
        flags = [tuple(flag) for flag in flags]
        flags = list(set(flags))
        flags.sort()
        return flags

    def maximal_flags(self):
        flags = self.all_flags()
        maximal = []
        for f in flags:
            #print "is maximal %s ?"%(f,)
            for g in flags:
                if f is g:
                    continue
                if set(g).intersection(f) == set(f):
                    #print "%s contains %s" % (g, f)
                    break
            else:
                #print "maximal:", f
                maximal.append(f)
        return maximal

    def check(self):
        # every maximal flag is a chamber
        for flag in self.maximal_flags():
            assert len(flag)==self.rank, "%s not a chamber %d"%(flag, self.rank)

    @classmethod
    def simplex(cls, dim):
        tps = range(dim+1)
        


def main():

    desc = r"""

    a
    |\ B
    | \
   C|  c
    | /
    |/ A
    b

    """

    triangle = Geometry(
        "aB aC bA bC cB cA".split(), 
        {'a':'vertex', 'b':'vertex', 'c':'vertex', 'A':'edge', 'B':'edge', 'C':'edge'})

    assert len(triangle.items) == 6

    print desc
    print triangle
    #for flag in triangle.all_flags():
    #    print flag

    print "maximal:"
    for flag in triangle.maximal_flags():
        print flag



if __name__ == "__main__":

    main()



