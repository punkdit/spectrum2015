#!/usr/bin/env python

"""
SSE quantum monte-carlo
"""

from __future__ import print_function

import sys, os
from random import randint, random, shuffle

import numpy
from numpy import kron, dot, allclose, log
from scipy.linalg import expm
from scipy import stats

from argv import argv
import models
from code import lstr2
from solve import eq2, dot2, zeros2, shortstrx

def ustr(u):
    s = ''.join([('1' if i else '.') for i in u])
    s = '('+s+')'
    return s


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

def rand2(n):
    v = numpy.random.binomial(1, 0.5, (n,))
    return v

def factorial(n):
    r = 1
    for i in range(1, n+1):
        r *= i
    return r

assert factorial(3) == 6
assert factorial(4) == 24
assert factorial(5) == 120



    
def get_dense(model, offset):
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

    if offset != 0.:
        H += offset * numpy.identity(N)

    return H

def get_partition(model, beta, offset):
    H = get_dense(model, offset)
    H = beta * H
    A = expm(H)
    x = numpy.trace(A)
    return x


    

class Metro(object):
    def __init__(self, Gx, Gz, beta, offset=0.):
        mx, n = Gx.shape
        self.n = n
        self.mx = mx
        self.mz = len(Gz)
        assert Gz.shape[1] == n
        self.Gx = Gx
        self.Gz = Gz
        self.beta = beta
        self.offset = offset # add multiple of the identity
        I = zeros2(n)
        letters = [I] + list(Gx)
        #print("letters:", [l for l in letters])
        self.letters = letters

    def get(self, u, v):
        assert u.shape == (self.n,)
        assert v.shape == (self.n,)
        #print("get", u, v)
        if eq2(u, v):
            w = self.offset + self.mz - 2*(dot2(self.Gz, u).sum())
        else:
            w = 0.
            for g in self.Gx:
                u1 = (g + u)%2
                if eq2(u1, v):
                    w += 1.
                    assert eq2((g+v)%2, u)
            assert w > 0.
        #print("   w =", w)
        return w

    def getx(self, spins):

        #print("getx", spins)
        assert len(spins) > 1
        assert eq2(spins[0], spins[-1])

        #spins = list(spins) + [spins[0]]
        N = len(spins)
        w = 1.
        for idx in range(N-1):
            val = self.get(spins[idx], spins[idx+1])
            #assert val > 0
            w *= val
            if w==0:
                break
        #print(" w =", w)
        return w

    def evaluate(self, state):
        letters = self.letters
        u, word = state
        #print("evaluate", word)
        assert len(word)
        spins = [u]
        for lidx in word:
            op = letters[lidx]
            u = (u + op)%2
            spins.append(u)
        r = self.getx(spins)
        nword = len(word)
        r *= (self.beta**nword) / factorial(nword) # CHECK THIS
        return r

    def get_init(self):
        u = rand2(self.n)
        nword = 1 + int(round(self.beta) * len(self.Gx))
        half = nword//2 or 1
        word = [randint(0, len(self.letters)-1) for i in range(half)]
        word = word + word # double each letter
        shuffle(word)
        state = u, word
        return state

    def move(self, state0):
        u, word = state0
        word = list(word)
        nlet = len(self.letters)
        nword = len(word)

        i = randint(0, 5)

        #print("move", i)

        if i==0:
            # add I
            idx = randint(0, nword)
            word.insert(idx, 0)

        elif i==1 and 0 in word:
            # remove I
            while 1:
                idx = randint(0, nword-1)
                if word[idx] == 0:
                    word.pop(idx)
                    break

            if len(word)==0:
                word.append(0)

        elif i==2:
            # add a pair
            idx = randint(0, nword)
            lidx = randint(0, nlet-1)
            word.insert(idx, lidx)
            word.insert(idx, lidx)

        elif i==3 and nword>1:
            # swap
            idx = randint(0, nword-2)
            word[idx], word[idx+1] = word[idx+1], word[idx]

        elif i==4 and nword>1:
            # remove pair
            for _ in range(nword): # tune this?
                idx = randint(0, nword-2)
                jdx = randint(idx+1, nword-1)
                if word[idx]==word[jdx]:
                    word.pop(jdx)
                    word.pop(idx)
                    break
            if len(word)==0:
                word.append(0)

        elif i==5:
            # translate
            u = u.copy()
            idx = randint(0, len(u)-1)
            u[idx] = 1 - u[idx] # flip

#        elif i==6:
#            # rotate
#        else:
#            assert 0

        return u, word


    def get_sample(self):
        n = self.n
    
        assert n <= 10
        N = 2**n
    
        state = self.get_init()

        W = self.evaluate(state)
        #print(statestr(state), W)

        trials = argv.get("trials", 100)
        for trial in range(trials):

            #print(state)
            u, word = state

            state1 = self.move(state)
#            print(statestr(state1), end=" ")
            
            #print(state1)
            W1 = self.evaluate(state1)

#            print("(%s)"%W1, end=" ")

            if W==0 or random() <= W1 / W:
                state = state1
                W = W1
                #print(statestr(state), W)
                #print("W =", W)

        n = len(state[1])
        if argv.verbose:
            print(statestr(state), W)
        return n


def statestr(state):
    u, word = state
    s = ustr(u) + "[%s]"%(''.join(str(i) for i in word))
    return s

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

    offset = len(model.Gz) # much harder!
    #offset = 0.
    #print(H.min())
    #return

    if argv.metro:
    
        print("-----------------")

        metro = Metro(model.Gx, model.Gz, beta, offset)
    
        counts = argv.get("counts", 100)
        ns = []
        for count in range(counts):
            n = metro.get_sample()
            ns.append(n)

        ns = numpy.array(ns)
        avg = numpy.mean(ns)
        print("<n> = %s, <lambda> = %s" % (avg, avg/beta))


    if not argv.exact:
        return

    print("-----------------")

    H = get_dense(model, offset)

    vals = numpy.linalg.eigvalsh(H)
    vals = list(vals)
    vals.sort()
    #print(vals[-3:])
    print("eigval:", vals[-1])
    print("gap:", vals[-1] - vals[-2])

    Z = 0.
    dZ = 0.
    for val in vals:
        x = numpy.exp(beta * val)
        Z += x
        dZ += val * x
    
    print("Z =", Z)
    print("slope =", dZ/Z)

    if argv.taylor:
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

        # 1/Z(dZ/dbeta) : should be close to eigval for large beta
        print("<val> = ", dZ/Z)

    if argv.slope or argv.plot:

        print("-----------------")
    
        H = beta * H
        A = expm(H)
        x = numpy.trace(A)
    
        print("Z  =", x)

        ys = []
        #xs = numpy.arange(1., 100., 1.)
        #xs = numpy.arange(0.01, 1., 0.01)
        #xs = numpy.arange(1., 10., 1.)
        xs = numpy.arange(10., 20., 1.)
        #xs = numpy.arange(0.9*beta, 1.1*beta, 0.01*beta)
        for beta in xs:
            y = get_partition(model, beta, offset)
            y = log(y)
            ys.append(y)
    
        ys = numpy.array(ys)
    
        slope, intercept, r_value, p_value, std_err = stats.linregress(xs, ys)
        print("slope:", slope) # not working very well...
    
        if argv.plot:
            import matplotlib.pyplot as pyplot
            pyplot.plot(xs, ys, 'bo')
            pyplot.show()


    




if __name__ == "__main__":
    main()





