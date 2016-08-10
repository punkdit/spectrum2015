#!/usr/bin/env python

import sys

import numpy
from numpy import linalg as la

from code import texstr
from lanczos import show_eigs

from argv import Argv
argv = Argv()


I = numpy.array([[1, 0], [0, 1.]])
X = numpy.array([[0, 1], [1, 0.]])
Z = numpy.array([[1, 0], [0, -1.]])


def test():
    H = numpy.zeros((4, 4))
    H[0:2, 0:2] = X + Z + I
    H[2:, 2:] = X + Z - I
    
    H[0:2, 2:] = I
    H[2:, 0:2] = I
    
    #w, v = la.eigh(H)
    
    H2 = numpy.zeros((8, 8))
    H2[:4, :4] = H
    H2[4:, 4:] = H
    
    #print H2
    
    for i in range(4):
    
        H2[i,i+4] = 1
        H2[i+4,i] = 1
        H2[i+4,i+4] = H2[i,i]-4
    
    print H2
    w, v = la.eigh(H2)
    
    print w
    v = v[:,-1]
    print v
    
    print
    for i in range(4):
        print v[i] / v[i+4]
    

def reflect(A, n):
    assert n < len(A)
    B = A[:n, :n]
    B[n-1, n-1] += A[n-1, n]
    return B


def test_gcolor():

    # This computes energy of the groundstate of the n=15 gauge color code

    n = 7

    A = numpy.zeros((n, n))
    
    for k in range(n):
        A[k, k] = 3*(6-2*k)
        if k>0:
            # from k-1 to k
            A[k, k-1] = 3*k
        if k<6:
            # from k+1 to k
            A[k, k+1] = 3*(6-k)
    
    if argv.transpose():
        A = A.transpose()

    print A # not hermitian! not normal !
    #print texstr(A)

    u, v = la.eig(A)
    print "eigs:",
    for x in u:
        print x,
    print
    print "v[:, 0]:"
    print v[:, 0]

    return

    B = A[:n-1, :n-1]
    print
    print B

    u, v = la.eig(B)
    print u[0]
    print v[:, 0]


    B = reflect(A, n-1)
    #B[n-2, n-2] += 3

    for count in range(10):
        print
        print B

        u, v = la.eig(B)
        print "eigs:", u
        for i in range(len(u)):
            print v[:, i]
    
        B[2, 2] += 1


def test_lie():

    # This computes energy of the groundstate of the n=15 gauge color code

    n = argv.get("n", 7)

    A = numpy.zeros((n, n))

    m = 1
    
    for k in range(n):
        A[k, k] = m*(n-1-2*k)
        if k>0:
            # from k-1 to k
            A[k, k-1] = m*k
        if k<n-1:
            # from k+1 to k
            A[k, k+1] = m*(n-1-k)

    assert abs(A.trace()) < 1e-10, A.trace()

    print A # not hermitian! not normal !
    #print texstr(A)

    r2 = 2**0.5
    vals, vecs = la.eig(A)
    print "eigs/(2**0.5):"
    show_eigs(vals / r2)


def test_compass():
    A = numpy.zeros((3, 3))
    
    # self interactions
    A[1,1] += 4
    A[2,2] += 3
    
    # 1 --> 0
    A[0, 1] = 9
    
    # 0 --> 1
    A[1, 0] = 1
    
    # 2 --> 1
    A[1, 2] = 4
    
    # 1 --> 2
    A[2, 1] = 6
    
    print A
    
    # G_z
    A[0,0] += 9
    A[1,1] += 1
    A[2,2] += -3
    
    print A
    u, v = la.eig(A)
    print u
    print v#[:, 0]
    
    # eigvals: [ 11.21110255   6.          -3.21110255]

    print texstr(A)

    A = A.astype(numpy.int32)
    print A
    from sympy import Matrix
    A = Matrix(A)
    print A
    print A.eigenvals()


def test_orbits():

    A = [
    [3, 1, 1, 1, 0, 0, 0, 0],
    [1, 1, 0, 0, 1, 1, 0, 0],
    [1, 0, 1, 0, 1, 0, 1, 0],
    [1, 0, 0, 1, 0, 1, 1, 0],
    [0, 1, 1, 0, -1, 0, 0, 1],
    [0, 1, 0, 1, 0, -1, 0, 1],
    [0, 0, 1, 1, 0, 0, -1, 1],
    [0, 0, 0, 0, 1, 1, 1, -3]]

    A = numpy.array(A)
    #u, v = numpy.linalg.eigh(A)

    P = numpy.array([
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 1, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 1, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1]])
    #print texstr(P)

    AP = numpy.dot(A, P)
    print AP

    Q = numpy.array([
        [1, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 1, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 1]])

    #print texstr(Q)

    QAP = numpy.dot(Q, AP)
    print QAP

    #PAP = numpy.dot(P.transpose(), AP)
    #u, v = la.eigh(A)
    #print u
    #u1, v1 = la.eigh(PAP)
    #print (u1[-1] / u[-1])**0.5


def test_circulant():

    n = argv.get("n", 8)

    A = numpy.zeros((n, n))

    for i in range(n):
        A[i, (i+1)%n] = 1
        A[i, (i-1)%n] = 1

    vals, vecs = la.eigh(A)
    print list(vals)


if argv.gcolor:
    test_gcolor()
elif argv.lie:
    test_lie()
elif argv.compass:
    test_compass()
elif argv.orbits:
    test_orbits()
elif argv.circulant:
    test_circulant()


    



