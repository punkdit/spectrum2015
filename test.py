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
    


def gcolor_1():

    # This computes energy of the groundstate of the n=15 gauge color code

    A = numpy.zeros((7, 7))
    
    for k in range(7):
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

#gcolor_1()


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



