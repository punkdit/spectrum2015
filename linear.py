#!/usr/bin/env python

import sys, os

import numpy
from scipy import sparse
from scipy.sparse.linalg import eigs, eigsh

from solve import shortstr, shortstrx, parse, eq2, dot2, zeros2, array2, identity2
from solve import row_reduce, RowReduction, span, get_reductor
from solve import u_inverse, find_logops, solve, find_kernel, linear_independant
from solve import System, Unknown, pseudo_inverse

import isomorph
from lanczos import write, show_eigs
from code import lstr2


def genidx(shape):
    if len(shape)==0:
        yield ()
    else:
        for idx in range(shape[0]):
            for _idx in genidx(shape[1:]):
                yield (idx,)+_idx


class LinSet(object):
    """
        A linear independent set of rows of G.
    """
    def __init__(self, G, idxs=None, left=None):
        m, n = G.shape
        self.G = G
        if idxs is None:
            idxs = [] # choose these
            left = range(m) # further choices
        self.idxs = idxs
        self.left = left

    def __len__(self):
        return len(self.idxs)

    def __getitem__(self, i):
        return self.idxs[i]

    def __repr__(self):
        return "%s(%s)"%(self.__class__.__name__, self.idxs)
    __str__ = __repr__

    def get(self):
        U = self.G[self.idxs, :]
        assert U.shape == (len(self.idxs), self.G.shape[1])
        return U

    def iter(self):
        "iterate  LI sets that contain this one"
        G, idxs, left = self.G, self.idxs, self.left
        for ii, i in enumerate(left):
            idxs1 = idxs+[i]
            G1 = G[idxs1, :]
            if len(row_reduce(G1, truncate=True)) < len(G1):
                continue
            li = LinSet(G, idxs1, left[ii+1:])
            yield li

    def mul(self, v):
        G, idxs = self.G, self.idxs
        u = zeros2(G.shape[1])
        for ii, i in enumerate(v):
            if i:
                u += G[idxs[ii]]
        u %= 2
        return u


class LinList(LinSet):
    "ordered LinSet"
    def iter(self):
        "iterate ordered LI sets that contain this one"
        G, idxs, left = self.G, self.idxs, self.left
        for ii, i in enumerate(left):
            idxs1 = idxs+[i]
            G1 = G[idxs1, :]
            if len(row_reduce(G1, truncate=True)) < len(G1):
                continue
            left1 = left[:ii] + left[ii+1:]
            li = LinList(G, idxs1, left1)
            yield li


def allowed(G, lookup, src, tgt):
    # lookup: map vector string to row index in G
    r = len(tgt)
    assert r>1
    Q = {} # permutation
    #for v in genidx((2,)*r):
    for v in genidx((2,)*(r-1)):
        v = v+(1,) # always include the newest vector
        u = src.mul(v)
        idx = lookup.get(u.tostring())
#        print "src.mul:", v, "=", u, idx
        if idx is None: # not in G so we don't care
            continue
        u = tgt.mul(v)
        jdx = lookup.get(u.tostring())
#        print "->", u, jdx
        if jdx is None:
            # must send a G row vector to another G row vector
            assert v[-1] == 1
            return False
        kdx = Q.get(idx)
        if kdx is None:
            Q[idx] = jdx # OK found a perm
        elif kdx != jdx:
            assert v[-1] == 1
            assert 0
            #print "X"
            return False # Not a perm
#        print "->", u, jdx
#    print "allowed:", Q
    return True


def search(G, src0=None, tgt0=None):
    m, n = G.shape
    lookup = dict((g.tostring(), idx) for (idx, g) in enumerate(G))

    if src0 is None:
        src0 = LinSet(G)
        while len(src0) < n:
            for src0 in src0.iter():
                break
        tgt0 = LinList(G)

    src = src0
    U = src.get()
    U = U.transpose()
    U1 = pseudo_inverse(U)
    assert eq2(dot2(U, U1), identity2(n))

    for tgt in tgt0.iter():

        #assert len(src)==len(tgt)>0

        if len(tgt)==1:
            pass
        elif not allowed(G, lookup, src, tgt):
            continue

        if len(tgt)==n:
            assert allowed(G, lookup, src, tgt)
            V = tgt.get().transpose()
            P = dot2(V, U1)
            #print shortstr(dot2(V, U1))
            #yield src, tgt
            yield P

        else:
    
            for solution in search(G, src, tgt):
                yield solution


def find_isos(Gz, Rx):

    RR = dot2(Gz, Rx.transpose())
    rows = [shortstr(r) for r in RR]
    for row in rows:
        print rows.count(row),
    print

    uniq = list(set(rows))
    uniq.sort(reverse=True)
    print "uniq:", len(set(rows))

    idxs = [rows.index(row) for row in uniq]
    G = RR[idxs,:]

    print shortstrx(G)
    #print '\n'.join(uniq)

    normal = [shortstr(row) for row in G]
    normal.sort()
    normal = '\n'.join(normal)

    count = 0
    found = set()
    Ps = []
    for P in search(G):
        #print src, tgt
        #write('.')
        GP = dot2(G, P.transpose())
        #print
        #print shortstrx(P, G, GP)
        #print
        rows = [shortstr(row) for row in GP]
        rows.sort()
        rows = '\n'.join(rows)
        assert rows == normal

        s = P.transpose().tostring()
        if s not in found:
            found.add(s)
            Ps.append(P)
            write(count)
            count += 1
        else:
            write('/')
    print
    print "count:", count
    return Ps, found


def orbigraph(Gx, Gz, Hx, Hz, Rx, Rz, Pxt, Qx, Pz, Tx, **kw):

    r, n = Rx.shape
    N = 2**r

    gz = len(Gz)
    RR = dot2(Gz, Rx.transpose())
    PxtQx = dot2(Pxt, Qx)

    Ps, found = find_isos(Gz, Rx)

    if not argv.check:
        return

    H = {}
    A = {}
    U = []

    lookup = {}
    for i, v in enumerate(genidx((2,)*r)):
        v = array2(v)
        lookup[v.tostring()] = i
        syndrome = dot2(Gz, Rx.transpose(), v)
        value = gz - 2*syndrome.sum()
        #print shortstr(dot2(Rx.transpose(), v)), value
        H[i, i] = value
        U.append(value)

    vectors = []
    for i, v in enumerate(genidx((2,)*r)):
        v = array2(v)
        vectors.append(v)
        #print shortstr(v),
        for g in Gx:
            v1 = (v + dot2(g, PxtQx))%2
            #j = eval('0b'+shortstr(v1, zero='0'))
            j = lookup[v1.tostring()]
            H[i, j] = H.get((i, j), 0) + 1
            A[i, j] = A.get((i, j), 0) + 1

    print len(H)
    #return

    bag0 = isomorph.from_sparse_ham(N, H)
    bag1 = isomorph.from_sparse_ham(N, H)

    count = 0
    for fn in isomorph.search(bag0, bag1):

        P = zeros2(r, r)
        v0 = zeros2(r)
        for i in range(r):
            v0[i] = 1
            j = eval('0b'+shortstr(v0, zero='0'))
            idx = fn[j]
            v1 = vectors[idx]
            #print shortstr(v0), "->", shortstr(v1)
            P[:, i] = v1
            v0[i] = 0

        assert P.tostring() in found

        for i, v0 in enumerate(vectors):
            v1 = dot2(P, v0)
            j = fn[i]
            if not eq2(v1, vectors[j]):
                write('/')
                break
            assert shortstr(vectors[j])==shortstr(v1)
        else:
            count += 1
            write('.')

        print
        print shortstrx(P, RR, dot2(RR, P))
        print

    print
    print "count:", count

#    size = 1
#    for j in range(r):
#        size *= 2**r - 2**j
#    print "|GL_%d(F_2)| = %d"%(r, size)
#
#    size = 1
#    for j in range(1, r):
#        size *= 2**j
#    size = size*size
#    for j in range(1, r+1):
#        size *= j
#    print "|BWB| =", size



def main():

    from zorbit import find_logops, find_errors, check_commute, check_conjugate

    import models
    Gx, Gz, Hx, Hz = models.build()

    assert not argv.orbiham, "it's called orbigraph now"

    Lz = find_logops(Gx, Hz)
    #print "Lz:", shortstr(Lz)

    if Lz.shape[0]*Lz.shape[1]:
        print Lz.shape, Gx.shape
        check_commute(Lz, Gx)
        check_commute(Lz, Hx)

    Px = get_reductor(Hx) # projector onto complement of rowspan of Hx
    Pz = get_reductor(Hz) 

    Rz = [dot2(Pz, g) for g in Gz]
    Rz = array2(Rz)
    Rz = row_reduce(Rz, truncate=True)
    rz = len(Rz)
    print "Rz:", rz

    n = Gx.shape[1]
    print "Hx:", len(Hx)
    print "Gx:", len(Gx)
    print "Gz:", len(Gz)

    Rx = [dot2(Px, g) for g in Gx]
    Rx = array2(Rx)

    Rx = row_reduce(Rx, truncate=True)
    rx = len(Rx)
    print "Rx:", rx

    Qx = u_inverse(Rx)
    Pxt = Px.transpose()
    assert eq2(dot2(Rx, Qx), identity2(rx))
    assert eq2(dot2(Rx, Pxt), Rx)

    #print shortstr(dot2(Pxt, Qx))
    PxtQx = dot2(Pxt, Qx)

    offset = argv.offset

    if len(Hz):
        Tx = find_errors(Hz, Lz, Rz)
    else:
        Tx = zeros2(0, n)

    if argv.orbigraph:
        orbigraph(**locals())
        return




from argv import Argv
argv = Argv()

if __name__ == "__main__":

    main()





