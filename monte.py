#!/usr/bin/env python

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


def dotx(items):
    idx = 0 
    A = items[idx]
    while idx+1 < len(items):
        B = items[idx+1]
        A = dot(A, B)
        idx += 1 
    return A


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


def mkzop(g, delta):
    ops = []
    for i in g:
        if i==1:
            A = Z
        elif i==0:
            A = I
        else:
            assert 0, str(g)
        #A = delta * A
        #A = expm(A)
        ops.append(A)
    A = kronx(ops)
    A = expm(delta * A) # SLOOOW !!
    return A
    

def mkxop(g, delta):
    A = mkzop(g, delta)
    n = len(g)
    HH = kronx([H]*n)
    A = dotx([HH, A, HH])
    return A


def test():

    delta = 0.1

    A = mkzop([0, 1], delta)
    B = expm(delta * kron(I, Z))
    assert allclose(A, B)

    A = mkzop([1, 1], delta)
    B = expm(delta * kron(Z, Z))
    assert allclose(A, B)

    A = mkxop([0, 1], delta)
    B = expm(delta * kron(I, X))
    assert allclose(A, B)

    A = mkxop([1, 1], delta)
    B = expm(delta * kron(X, X))
    assert allclose(A, B)

    #print(B)

    
def get_partition(model, beta):
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

    H = beta * H
    A = expm(H)
    x = numpy.trace(A)
    return x


test()
    

def main():

    beta = argv.get("beta", 1.0)
    M = argv.get("M", 3)
    delta = beta / M

    verbose = argv.verbose

    model = models.build_model()
    print(model)
    n = model.n

    assert n <= 10
    N = 2**n

    HA = numpy.identity(N)
    for g in model.Gz:
        A = mkzop(g, delta)
        HA = dot(HA, A)

    if verbose:
        s = lstr2(HA)
        s = s.replace("0.0000", "      ")
        print(s)

    HB = numpy.identity(N)
    for g in model.Gx:
        B = mkxop(g, delta)
        HB = dot(HB, B)

    if verbose:
        s = lstr2(HB)
        s = s.replace("0.0000", "      ")
        print(s)

    x = HB.min()
    if x < 0.:
        HB -= 2*x

    def get_weight(config):
        assert len(config) == 2*M
        W = 1.
        for idx in range(M):
            W *= HA[config[2*idx], config[2*idx+1]]
            W *= HB[config[2*idx+1], config[(2*idx+2)%(2*M)]]
        return W

    #config = [randint(0, N-1) for i in range(2*M)]
    def get_rand():
        config = [randint(0, N-1)]
        while len(config) < 2*M:
            while 1:
                idx = randint(0, N-1)
                if HA[config[-1], idx]>0.:
                    break
            config.append(idx)
            if len(config) == 2*M:
                break
            while 1:
                idx = randint(0, N-1)
                if HB[config[-1], idx]>0.:
                    break
            config.append(idx)
        return config

    trials = argv.get("trials", 10000)
    counts = argv.get("counts", 100)

    samples = []
    for count in range(counts):
      config = get_rand()
      config1 = list(config)
      W0 = get_weight(config)
      assert W0 > 0., (W0, trial)
      accept = 0
      for trial in range(trials):

        idx = 2*randint(0, M-1)
        val = randint(0, N-1)
        config1[idx] = val
        config1[idx+1] = val
        W1 = get_weight(config1)
        #if W1==0.:
        #    print(".", end="", flush=True)

        if random() <= W1 / W0:
            #config = config1
            config[idx] = val
            config[idx+1] = val
            #print(config1)
            W0 = W1
            accept += 1
            #print(W1)
        else:
            val0 = config[idx]
            config1[idx] = val0
            config1[idx+1] = val0

      #print("accept:", accept)

      #print("W0 =", W0)
      samples.append(W0)
    #W = sum(samples) / counts
    #print("W = ", W)
    #print("Z* =", W*N)
    samples = numpy.array(samples)
    samples *= N
    print("mean=%s, std=%s" % (numpy.mean(samples), numpy.std(samples)))
    #print(samples)
    
    # -----------------------------------------

    if not argv.exact:
        return

    A = numpy.identity(N)
    for i in range(M):
        A = dot(A, HA)
        A = dot(A, HB)

    print("Z  =", numpy.trace(A))
    print("Z  =", get_partition(model, beta))

    ys = []
#    xs = numpy.arange(1., 100., 1.)
    #xs = numpy.arange(0.01, 1., 0.01)
    #xs = numpy.arange(1., 10., 1.)
    xs = numpy.arange(10., 20., 1.)
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





