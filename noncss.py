#!/usr/bin/env python

import sys, os

import numpy

A = numpy.zeros((4, 4))

A[0, 1] = 1
A[1, 0] = 1
A[2, 3] = 1
A[3, 2] = 1

A[0, 2] = -1
A[1, 3] = -1
A[2, 0] = -1
A[3, 1] = -1

print A

vals, vecs = numpy.linalg.eigh(A)

print vals

print vecs


