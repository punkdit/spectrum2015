#!/usr/bin/env python

"""
Gaussian elimination over the rationals.
"""

import sys, os
from random import randint
from fractions import Fraction


import numpy
from numpy import dot

from argv import argv



def zeros(m, n):
    A = numpy.empty((m, n), dtype=object)
    A[:] = 0
    return A


def identity(m):
    I = zeros(m, m)
    for i in range(m):
        I[i, i] = 1
    return I


def eq(A, B):
    r = numpy.abs(A-B).sum()
    return r==0


def shortstr(A):
    return str(A)


def shortstrx(*args):
    return '\n'.join(str(A) for A in args)



def swap_row(A, j, k):
    row = A[j, :].copy()
    A[j, :] = A[k, :]
    A[k, :] = row

def swap_col(A, j, k):
    col = A[:, j].copy()
    A[:, j] = A[:, k]
    A[:, k] = col


def plu_reduce(A, truncate=False, check=False, verbose=False):
    """
    Solve PLU = A, st. P is permutation, L is lower tri, U is upper tri.
    Remove zero rows from U if truncate=True.
    """

    m, n = A.shape
    P = identity(m)
    L = identity(m)
    U = A.copy()

    assert m*n, (m, n)

    if verbose:
        print "plu_reduce:"
        print "%d rows, %d cols" % (m, n)

    i = 0
    j = 0
    while i < m and j < n:
        if verbose:
            print "i, j = %d, %d" % (i, j)
            print "P, L, U:"
            print shortstrx(P, L, U)

        assert i<=j
        if i and check:
            assert U[i:,:j].max() == 0 # XX rm

        # first find a nonzero entry in this col
        for i1 in range(i, m):
            if U[i1, j]:
                break
        else:
            j += 1 # move to the next col
            continue # <----------- continue ------------

        if i != i1:
            if verbose:
                print "swap", i, i1
            swap_row(U, i, i1)
            swap_col(P, i, i1)
            swap_col(L, i, i1)
            swap_row(L, i, i1)

        if check:
            A1 = dot(P, dot(L, U))
            assert eq(A1, A)

        r = U[i, j]
        assert r != 0
        for i1 in range(i+1, m):
            s = U[i1, j]
            if s != 0:
                t = Fraction(s, r)
                if verbose: 
                    print "add %s to %s" % (i, i1)
                L[i1, i] = t
                U[i1, :] -= t*U[i, :]
                assert U[i1, j] == 0

        if check:
            A1 = dot(P, dot(L, U))
            assert eq(A1, A)

        i += 1
        j += 1

    if check:
      for i in range(m):
        for j in range(i):
            assert U[i, j] == 0
      for i in range(m):
        for j in range(i+1, m):
            assert L[i, j] == 0

    if truncate:
        m = U.shape[0]-1
        #print "sum:", m, U[m, :], U[m, :].sum()
        while m>=0 and U[m, :].sum()==0:
            m -= 1
        U = U[:m+1, :]

    return P, L, U


def u_inverse(U, check=False, verbose=False):
    """invert a row reduced U
    """

    m, n = U.shape

    #items = []
    leading = []
    for row in range(m):
        cols = numpy.where(U[row, :])[0]
        if not len(cols):
            break
        col = cols[0]
        leading.append(col)

    U1 = zeros(n, m)

    #print shortstrx(U, U1)

    # Work backwards
    i = len(leading)-1 # <= m
    while i>=0:

        j = leading[i]
        #print "i=", i, "j=", j
        U1[j, i] = 1

        #print shortstrx(U, U1)

        k = i-1
        while k>=0:
            #print "dot"
            #print shortstr(U[k,:])
            #print shortstr(U1[:,i])
            if dot(U[k, :], U1[:, i]):
                j = leading[k]
                #print "set", j, i
                U1[j, i] = 1
            #print shortstr(U1[:,i])
            assert dot(U[k, :], U1[:, i]) == 0
            k -= 1
        i -= 1

    return U1


def l_inverse(L, check=False, verbose=False):
    """invert L (lower triangular, 1 on diagonal)
    """

    m, n = L.shape
    assert m==n
    L1 = identity(m)

    # Work forwards
    for i in range(m):
        #u = L1[:, i]
        for j in range(i+1, m):
            r = dot(L[j, :], L1[:, i])
            if r:
                L1[j, i] = 1
            assert dot(L[j, :], L1[:, i]) == 0

    assert eq(dot(L, L1), identity(m))
    return L1






def pseudo_inverse(A, check=False):
    m, n = A.shape
    if m*n == 0:
        A1 = zeros(n, m)
        return A1
    P, L, U = plu_reduce(A, verbose=False, check=check)
    L1 = l_inverse(L, check=check)
    U1 = u_inverse(U, check=check)
    A1 = dot(U1, dot(L1, P.transpose()))
    return A1


def solve(H, u, force=False, debug=False):
    "Solve Hv = u"
    A = pseudo_inverse(H)
    v = dot(A, u)
    if eq(dot(H, v), u) or force:
        return v



def test():

    m, n = 4, 5
    A = zeros(m, n)

    for i in range(m):
      for j in range(n):
        A[i, j] = Fraction(randint(-3, 3), randint(1, 3))

    P, L, U = plu_reduce(A, check=True, verbose=True)

    print P
    print L
    print U


if __name__ == "__main__":

    test()



