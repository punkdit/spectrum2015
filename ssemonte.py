#!/usr/bin/env python

"""
SSE quantum monte-carlo
"""

from __future__ import print_function

import sys, os
from random import randint, random

import numpy
from numpy import kron, dot, allclose, log
from scipy.linalg import expm
from scipy import stats

from argv import argv
import models
from code import lstr2
from solve import eq2, dot2, zeros2, shortstrx


def dotx(items):
    idx = 0 
    A = items[idx]
    while idx+1 < len(items):
        B = items[idx+1]
        A = dot(A, B)
        idx += 1 
    return A

def powop(A, d):
    if d==0:
        n = len(A)
        return numpy.identity(n)
    if d==1:
        return A
    return dotx([A]*d)


def kronx(items):
    idx = 0 
    A = items[idx]
    while idx+1 < len(items):
        B = items[idx+1]
        A = kron(A, B)
        idx += 1 
    return A



H = (1./numpy.sqrt(2.))*numpy.array([[1., 1], [1, -1]])
I = numpy.array([[1., 0], [0, 1]])
X = numpy.array([[0., 1], [1, 0]])
Z = numpy.array([[1., 0], [0, -1]])

eZ = numpy.array([[numpy.e, 0], [0, 1./numpy.e]])
eX = dotx([H, eZ, H])


    
def get_dense(model):
    n = model.n
    N = 2**n
    H = numpy.zeros((N, N))
    for g in model.Gz:
        ops = [(Z if i else I) for i in g]
        op = kronx(ops)
        H += op

    for g in model.Gx:
        ops = [(X if i else I) for i in g]
        op = kronx(ops)
        H += op
    return H

def get_partition(model, beta):
    H = get_dense(model)
    H = beta * H
    A = expm(H)
    x = numpy.trace(A)
    return x


class Hamiltonian(object):
    def __init__(self, Gx, Gz, offset=0.):
        mx, n = Gx.shape
        self.n = n
        self.mx = mx
        assert Gz.shape[1] == n
        self.Gx = Gx
        self.Gz = Gz
        self.offset = offset # add multiple of the identity

    def get(self, u, v):
        assert u.shape == (self.n,)
        assert v.shape == (self.n,)
        if eq2(u, v):
            w = self.offset
            for g in self.Gz:
                w += 2.*dot2(g, u) - 1
        else:
            w = 0.
            for g in self.Gx:
                u1 = (g + u)%2
                if eq2(u1, v):
                    w += 1.
                    assert eq2((g+v)%2, u)
        return w

    def getx(self, config):
        assert len(config)
        if len(config)==1:
            u = config[0]
            return self.get(u, u)

        #config = list(config) + [config[0]]
        N = len(config)
        w = 1.
        for idx in range(N):
            w *= self.get(config[idx], config[(idx+1)%N])
            if w==0.:
                break
        return w

    def get_config(self, d):
        # path of length d
        assert d>0
        n = self.n
        if d==1:
            config = [rand2(n)]
        elif d==2:
            v = rand2(n)
            idx = randint(0, self.mx)
            if idx==self.mx:
                u = v.copy()
            else:
                u = (v + self.Gx[idx])%2
            config = [v, u]
        else:
            assert 0

        return config

    def get_size(self, d):
        assert d>0
        N = 2**self.n
        if d==1:
            size = N
        elif d==2:
            size = N * (self.mx + 1)
        else:
            assert 0

        return float(size)
    

def rand2(n):
    v = numpy.random.binomial(1, 0.5, (n,))
    return v

def factorial(n):
    assert n>=0
    if n < 4:
        return float([1, 1, 2, 6][n])
    return factorial(n-1) * n

assert factorial(3) == 6
assert factorial(4) == 24
assert factorial(5) == 120


def metropolis(model, offset=0.):
    n = model.n

    assert n <= 10
    N = 2**n

    H = Hamiltonian(model.Gx, model.Gz, offset)

    # d = 0
    pfunc = N

    counts = argv.get("counts", 100)
    trials = argv.get("trials", 1000)

    for d in range(1, 3):
#        print("d =", d)
        accept = 0
        samples = []

        for count in range(counts):
            #config = [rand2(n) for i in range(d)]
            config = H.get_config(d)
            assert len(config) == d
            W = H.getx(config)
            #print("config:", config)
#            print(W, end=" ")
            for trial in range(trials):
                config1 = H.get_config(d)
                W1 = H.getx(config1)
                #print(W1)
                if W==0 or random() <= W1/W:
                    config = config1
                    W = W1
                    accept += 1
#                    print(W, end=" ")
#            print()
#            print(shortstrx(*config))
            samples.append(W)
        if counts*trials:
            print("accept:", 1.*accept / (counts*trials))
        samples = numpy.array(samples)
        size = H.get_size(d)
#        scale = size * ((beta ** d) / factorial(d)) / (numpy.e ** offset)
        scale = size * ((beta ** d) / factorial(d))
        #scale = ((beta ** d) / factorial(d))
        samples *= scale
        avg, std = (numpy.mean(samples), numpy.std(samples))
        print("d = %d, val = %f, dev = %f" % (d, avg, std))
        pfunc += avg
        

    print("Z* = ", pfunc)

    # -----------------------------------------


def main():

    beta = argv.get("beta", 1.0)
    M = argv.get("M", 3)
    delta = beta / M

    verbose = argv.verbose

    model = models.build_model()
    print(model)

    n = model.n
    N = 2**n

    print("N =", N)

    #offset = len(model.Gz) # much harder!
    offset = 0.
    #print(H.min())
    #return

    H = get_dense(model)
    H += offset * numpy.identity(N)

    vals = numpy.linalg.eigvalsh(H)
    vals = list(vals)
    vals.sort()
    #print(vals[-3:])
    print("eigval:", vals[-1])

    print("-----------------")

    D = argv.get("D", 10)

    Z = 0.
    dZ = 0.
    for d in range(D):
        if d==0:
            val = N
            x = 0.
        else:
            tr = numpy.trace(powop(H, d))
            val = ((beta**d) / factorial(d)) * tr
            x = ((beta**(d-1)) / factorial(d-1)) * tr
        print("d = %d, val = %f, x = %f" %(d, val, x))
        #print()
        
        Z += val
        dZ += x

    print("Z =", Z)
    print("<val> = ", dZ/Z)

    print("-----------------")

    H = beta * H
    A = expm(H)
    x = numpy.trace(A)

    print("Z  =", x)

    if argv.slope or argv.plot:
        ys = []
        #xs = numpy.arange(1., 100., 1.)
        #xs = numpy.arange(0.01, 1., 0.01)
        #xs = numpy.arange(1., 10., 1.)
        #xs = numpy.arange(10., 20., 1.)
        xs = numpy.arange(0.9*beta, 1.1*beta, 0.01*beta)
        for beta in xs:
            y = get_partition(model, beta)
            y = log(y)
            ys.append(y)
    
        ys = numpy.array(ys)
    
        slope, intercept, r_value, p_value, std_err = stats.linregress(xs, ys)
        print("slope:", slope)
    
        if argv.plot:
            import matplotlib.pyplot as pyplot
            pyplot.plot(xs, ys, 'bo')
            pyplot.show()


    




if __name__ == "__main__":
    main()





