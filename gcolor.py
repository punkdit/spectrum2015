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
from code import lstr2
from lanczos import write, show_eigs

import models
from models import genidx


def test():

    verbose = argv.verbose

    Gx, Gz, Hx, Hz = models.build("gcolor")
    model = models.build_model(Gx, Gz, Hx, Hz)
    if argv.verbose:
        print model
        print

    #print shortstrx(Hx)

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

    assert eq2(dot2(Rx, Pxt), Rx)

#    print shortstrx(Gx, Gz)
#    print

    PGx = dot2(Gx, Pxt)
    PGz = dot2(Gz, Pzt)
    
#    print shortstrx(PGx, PGz)
#    print

    ng = len(PGx)

    graph = nx.Graph()
    for i in range(ng):
        graph.add_node(i)

    for i, gi in enumerate(PGx):
        for j, gj in enumerate(PGz):
            if (gi*gj).sum()%2:
                graph.add_edge(i, j)

    equs = nx.connected_components(graph)
    assert len(equs) == 3, len(equs)
    assert bipartite.is_bipartite(graph)

    color = bipartite.color(graph)

    assert dot2(model.Hx, Pxt).sum() == 0

    mx = len(Hx)

    A = solve(concatenate((Hx, Rx)).transpose(), Gx.transpose())
    #print "A:"
    #print shortstrx(A)
    assert A.shape == (mx+r, ng)
    assert len(A[0]) == len(Gx)

    if argv.enum:

      #weightss = [weights]
      def weightss():
        for idxs in genidx((2,)*mx):
            weights = numpy.array([1.]*ng)
            print idxs,
            for i, idx in enumerate(idxs):
                if idx:
                    weights *= (1-2*A[i])
            yield weights
      weightss = weightss()

    else:
    
        weights = numpy.array([1.]*ng)
        k = argv.k
        if type(k) is tuple:
            ks = k
            hs = Hx[ks, :]
            print "stabilizer weights:", hs.sum(axis=1)
            h = Hx[ks, :].sum(axis=0) % 2
            print "stabilizer weight:", h.sum()
            for k in ks:
                weights *= (1-2*A[k])
        elif k is not None:
            weights = 1-2*A[k]
            print "stabilizer:", Hx[k].sum()
        if verbose:
            print "weights:", list(weights)
        weightss = [weights]

    for weights in weightss:
     total = 0.
     gap = None

     for equ in equs:
      for flag in [0, 1]:
        GIx = [PGx[i] for i in equ if color[i] == flag]
        if weights is not None:
            weights1 = [weights[i] for i in equ if color[i] == flag]
            #print "weights:", weights1
        GIz = [PGz[i] for i in equ if color[i] == 1-flag]
        GIx = array2(GIx)
        GIz = array2(GIz)

        #print shortstrx(GIx, GIz)
        model = models.build_model(GIx, GIz)
        if verbose:
            print model

        r1 = len(model.Rx)

        if r1<12 and not argv.slepc:
    
            H = model.build_ham(weights=weights1)

            vals, vecs = numpy.linalg.eigh(H)
            vals = list(vals)
            vals.sort(reverse=True)

        else:
            vals = do_slepc(model, weights=weights1)

        _gap = vals[0]-vals[1]
        if verbose:
            print "vals:", vals[:5], "gap:", _gap
        gap = _gap if gap is None or _gap<gap else gap
        total += vals[0]

        #return

     print "eval:", total,
     if verbose:
        print "gap:", gap,
        print "eval_2:", total-gap
     else:
        print
     sys.stdout.flush()



def do_slepc(model, weights=None):
    from zorbit import slepc

    name = "ex3.tmp"
    os.unlink(name)

    slepc(weights=weights, **model.__dict__)

    #subprocess.check_output(args, *, stdin=None, stderr=None, shell=False, universal_newlines=False)

    if weights is None:
        args = ("./%s -eps_nev 2  -eps_hermitian   -eps_largest_real"%name).split()
    else:
        args = ("./%s -eps_nev 2 -eps_ncv 40 -eps_type arnoldi -eps_largest_real"%name).split()

    #print args
    s = subprocess.check_output(args)

    #print s

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


