#!/usr/bin/env python

import sys, os
import subprocess
from fractions import gcd

import numpy
from numpy import concatenate

import networkx as nx
from networkx.algorithms import bipartite

from solve import get_reductor, array2, row_reduce, dot2, shortstr, zeros2, shortstrx, eq2
from solve import u_inverse, find_kernel, find_logops, identity2, solve
from solve import check_conjugate, check_commute
from lanczos import write, show_eigs

import models
from models import genidx

import gauge


def test():

    Gx, Gz, Hx, Hz = models.build("gcolor")
    model = models.build_model(Gx, Gz, Hx, Hz)
    print model

    if 0:
        H = model.build_ham()
        print H
        vals, vecs = numpy.linalg.eigh(H)
        show_eigs(vals)

    Rx, Rz = model.Rx, model.Rz
    Rxt = Rx.transpose()
    Rzt = Rz.transpose()

    Px, Pz = model.Px, model.Pz
    Pxt = Px.transpose()
    Pzt = Pz.transpose()

    r, n = Rx.shape

    #print shortstrx(Gx, Gz)

    ng = len(Gx)
    graph = nx.Graph()
    for i in range(ng):
        graph.add_node(i)

    for i, gi in enumerate(Gx):
        for j, gj in enumerate(Gz):
            if (gi*gj).sum()%2:
                graph.add_edge(i, j)

    equs = nx.connected_components(graph)
    assert len(equs) == 3
    assert bipartite.is_bipartite(graph)

    color = bipartite.color(graph)

    total = 0.
    gap = None

    for equ in equs:
        GIx = [Gx[i] for i in equ if color[i] == 0]
        GIz = [Gz[i] for i in equ if color[i] == 1]
        GIx = array2(GIx)
        GIz = array2(GIz)

        print
        #print shortstrx(GIx, GIz)
        model = models.build_model(GIx, GIz)
        print model

        if len(model.Rx)<12 and 0:
            H = model.build_ham()
            print H.shape
            vals, vecs = numpy.linalg.eigh(H)
            #vals = list(vals)
            #show_eigs(vals)
            val = numpy.max(vals)
            total += 2*val

        else:
            vals = do_slepc(model)
            _gap = vals[0]-vals[1]
            print vals, "gap:", _gap
            gap = _gap if gap is None or _gap<gap else gap
            total += 2*vals[0]

    print "eval:", total
    print "gap:", gap
    print "eval_2:", total-gap


def do_slepc(model):
    from zorbit import slepc

    name = "ex3.tmp"
    os.unlink(name)

    slepc(**model.__dict__)

    #subprocess.check_output(args, *, stdin=None, stderr=None, shell=False, universal_newlines=False)

    args = ("./%s -eps_nev 2  -eps_hermitian   -eps_largest_real"%name).split()
    s = subprocess.check_output(args)

    vals = []
    lines = s.split("\n")
    flag = 0
    for line in lines:
        line = line.strip()
        flds = line.split()
        if line.startswith("-----"):
            flag += 1
        elif flag == 1:
            assert flds[0] == 'k'
        elif flag == 2:
            #print flds
            val = float(flds[1])
            vals.append(val)

    return vals



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


