#!/usr/bin/env python

"""
SSE quantum monte-carlo
"""

from __future__ import print_function

import sys, os
from random import randint, random, shuffle, seed

import numpy
from numpy import kron, dot, allclose, log
from scipy.linalg import expm
from scipy import stats

from argv import argv
import models
from code import lstr2
from solve import eq2, dot2, zeros2, shortstrx, array2

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



    
def get_dense(model, offset, Jx=1.0, Jz=1.0):
    n = model.n
    N = 2**n
    H = numpy.zeros((N, N))
    for g in model.Gz:
        ops = [(Z if i else I) for i in g]
        op = kronx(ops)
        H += Jz*op

    for g in model.Gx:
        ops = [(X if i else I) for i in g]
        op = kronx(ops)
        H += Jx*op

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
    def __init__(self, Gx, Gz, beta, offset=0., Jx=1.0, Jz=1.0):
        mx, n = Gx.shape
        self.n = n
        self.mx = mx
        self.mz = len(Gz)
        assert Gz.shape[1] == n
        self.Gx = Gx
        self.Gz = Gz
        self.Jx = Jx
        self.Jz = Jz
        assert Jx >= 0.
        assert Jz >= 0.
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
        Jx, Jz = self.Jx, self.Jz
        if eq2(u, v):
            w = self.offset + Jz*(self.mz - 2*(dot2(self.Gz, u).sum()))
            assert w>=0.
        elif 1:
            # faster to add u & v, then count occurances in Gx ?
            w = 0.
            for g in self.Gx:
                u1 = (g + u)%2
                if eq2(u1, v):
                    w += Jx
                    assert eq2((g+v)%2, u)
            assert w > 0.
            assert w==Jx
        else:
            w = Jx
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
        if len(word)==0:
            return 1.0
        spins = [u]
        for lidx in word:
            op = letters[lidx]
            u = (u + op)%2
            spins.append(u)
        r = self.getx(spins)
        nword = len(word)
        #print("evaluate:", r, nword)
        r *= 1. * (self.beta**nword) / factorial(nword) # CHECK THIS
        return r

    def get_init(self):
        u = zeros2(self.n)
        word = []
        return (u, word)

    def get_init_rand(self):
        u = rand2(self.n)
        nword = 1 + int(round(2*self.beta*len(self.Gx)))
        half = nword//2 or 1
        word = [randint(0, len(self.letters)-1) for i in range(half)]
        word = word + word # double each letter
        shuffle(word)
        state = u, word
        return state

    def move(self, state0):
        global debug
        debug = None

        u, word = state0
        word = list(word)
        nlet = len(self.letters)
        nword = len(word)

        C = 6 # choices
        i = randint(0, C-1)
        base = 1./C

        #print("move", i)

        ratio = 1.0

        if i==0:
            # add I
            idx = randint(0, nword)

#            count = 1
#            jdx = idx
#            while jdx<len(word) and word[jdx]==0:
#                jdx += 1
#                count += 1
#            jdx = idx-1
#            while jdx>=0 and word[jdx]==0:
#                jdx -= 1
#                count += 1
#
            total = word.count(0)

#            fwd = base * (count+1) / (nword+1)
#            rev = base * (count+1) / (total+1)

            ratio = 1. * (nword+1) / (total+1)

            word.insert(idx, 0)

        elif i==1 and 0 in word:
            # remove I
            while 1:
                idx = randint(0, nword-1)
                if word[idx] == 0:
                    break

            total = word.count(0)
            ratio = 1.*total / nword

            word.pop(idx)

        elif i==2:
            # add a pair
            idx = randint(0, nword)
            jdx = randint(0, nword+1)
            lidx = randint(0, nlet-1)
            word.insert(idx, lidx)
            word.insert(jdx, lidx)

            counts = [word.count(i) for i in range(nlet)]
            total = 0 # total number of pairs in word
            for count in counts:
                if count:
                    total += count * (count-1) / 2 # pairs for this letter
            p = 1./total

            fwd = base*2./(nword+1)/(nword+2)/nlet
            rev = base*p
            debug = fwd, rev
            ratio = rev/fwd

        elif i==3 and nword>1:
            # remove pair
            counts = [word.count(i) for i in range(nlet)]

            total = 0 # total number of pairs in word
            for count in counts:
                if count:
                    total += count * (count-1) / 2 # pairs for this letter
            p = 1./total

            check = 0
            while 1:
                idx = randint(0, nword-1)
                #jdx = randint(idx+1, nword-1)
                jdx = randint(0, nword-1)
                if idx!=jdx and word[idx]==word[jdx]:
                    if jdx<idx:
                        jdx, idx = idx, jdx
                    word.pop(jdx)
                    word.pop(idx)
                    break
                check += 1
                assert check < 1e6, "wup... nword=%d"%nword

            fwd = base*p
            rev = base*2./nword/(nword-1)/nlet
            debug = fwd, rev
            ratio = rev/fwd

        elif i==4 and nword>1:
            # swap
            idx = randint(0, nword-2)
            word[idx], word[idx+1] = word[idx+1], word[idx]

        elif i==5:
            # translate
            u = u.copy()
            idx = randint(0, len(u)-1)
            u[idx] = 1 - u[idx] # flip

        return ratio, (u, word)

    def get_sample(self, verbose):
        n = self.n
    
        assert n <= 10
        N = 2**n
    
        state = self.get_init()

        W = self.evaluate(state)
        if verbose>1:
            indent = "    "
            print(indent+statestr(state), W)

        trials = argv.get("trials", 100)
        for trial in range(trials):

            ratio, state1 = self.move(state)
            W1 = self.evaluate(state1)
#            if verbose>1:
#                print(2*indent+statestr(state1), W1)
#            try:
#                ratio = W1*rev/(W*fwd)
#            except ZeroDivisionError:
#                W = 0.
            if W==0 or random() <= (ratio*W1/W):
                state = state1
                W = W1
                if verbose>1:
                    print(indent+statestr(state), W)

        nword = len(state[1])
        if verbose:
            print(statestr(state), W)
        return nword

    def sample_fwd(self, state0, state1, trials=10000):
        # sample transition probability
        count = 0.
        for trial in range(trials):
            ratio, state = self.move(state0)
            if eq2(state[0], state1[0]) and state[1]==state1[1]:
                count += 1.
        return count / trials



def statestr(state):
    u, word = state
    s = ustr(u) + "[%s]"%(''.join(str(i) for i in word))
    return s

def eqstate(state0, state1):
    return eq2(state0[0], state1[0]) and state0[1]==state1[1]


def test():
    EPSILON = 1e-14

    n = 4
    offset = 4
    Gx, Gz, _, _ = models.build_ising(n)
    beta = 1.0
    metro = Metro(Gx, Gz, beta, offset)

    state = (zeros2(n), [0])
    r = metro.evaluate(state)
    assert r==8.0

    state = (zeros2(n), [0, 0])
    r = metro.evaluate(state)
    assert r==32.0

    state = array2([1,0,0,0]), [0, 0]
    r = metro.evaluate(state)
    assert r==8.0
    #print(statestr(state), r)

    state = array2([0,0,0,0]), [0,1,1,0]
    r = metro.evaluate(state)
    #print(statestr(state), r)
    assert abs(r - (64./factorial(4))) < EPSILON

test()


def test_balance():
    n = 4
    offset = 4
    Gx, Gz, _, _ = models.build_ising(n)
    beta = 1.0
    metro = Metro(Gx, Gz, beta, offset)

    state0 = (zeros2(n), [])
    state1 = (zeros2(n), [0])

    trials = argv.get("trials", 100000)

    state0 = metro.get_init_rand()
    state1 = state0
    while eqstate(state0, state1):
        ratio, state1 = metro.move(state0)
    print(statestr(state0))
    print(statestr(state1))
    if debug:
        fwd, rev = debug
        print("fwd = %s" % fwd)
        print("rev = %s" % rev)
    fwd = metro.sample_fwd(state0, state1, trials)
    rev = metro.sample_fwd(state1, state0, trials)
    print("sample fwd = %s" % fwd)
    print("sample rev = %s" % rev)
    print("ratio = %s, sampled ratio = %s" % (ratio, rev/fwd))



def main():


    beta = argv.get("beta", 1.0)

    verbose = argv.get("verbose", 0)

    Jx = argv.get("Jx", 1.0)
    Jz = argv.get("Jz", 1.0)

#    theta = argv.get("theta")

    model = models.build_model()
    print(model)
    print("Gx Gz:")
    print(shortstrx(model.Gx, model.Gz))
    print("Jx = %s, Jz = %s" % (Jx, Jz))

    simseed = argv.get("simseed", 0)
    for i in range(simseed):
        random()

    #offset = 0.
    #offset = len(model.Gz) # much harder!
    offset = argv.get("offset", len(model.Gz))

    n = model.n
    N = 2**n

    print("N =", N)
    print("offset =", offset)


    if argv.metro:
    
        print("-----------------")

        metro = Metro(model.Gx, model.Gz, beta, offset, Jx, Jz)
    
        counts = argv.get("counts", 100)
        ns = []
        for count in range(counts):
            n = metro.get_sample(verbose)
            ns.append(n)

        ns = numpy.array(ns)
        avg = numpy.mean(ns)
        dev = numpy.std(ns) / (len(ns)**0.5)
        print("<n> = %s, <lambda> = %s +/- %.2f" % (avg, avg/beta, dev))


    if not argv.exact and not argv.taylor:
        return

    print("-----------------")

    H = get_dense(model, offset, Jx, Jz)

    vals = numpy.linalg.eigvalsh(H)
    vals = list(vals)
    vals.sort()
    #print(vals[-3:])
    print("eigval:", vals[-1] - offset)
    print("gap:", vals[-1] - vals[-2])

    Z = 0.
    dZ = 0.
    for val in vals:
        x = numpy.exp(beta * val)
        Z += x
        dZ += val * x
    
    #print("Z =", Z) # FIX: correct for offset
    print("<H> =", dZ/Z - offset)

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
            print("n = %d, val = %f, x = %f" %(d, val, x))
            #print()
            
            Z += val
            dZ += x
    
        #print("Z =", Z) # FIX: correct for offset

        # 1/Z(dZ/dbeta) : should be close to eigval for large beta
        print("<H> = ", dZ/Z - offset)

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

    _seed = argv.get("seed")
    if _seed is not None:
        numpy.random.seed(_seed)
        seed(_seed)

    if argv.test:
        test_balance()
    else:
        main()





