#!/usr/bin/env python

import sys

import numpy
from numpy import linalg as la



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
    
    print A # not hermitian! not normal !

    u, v = la.eig(A)
    print u[0]
    print v[:, 0]

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



test_gcolor()


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



