#!/usr/bin/env python

"""
Find the closure of the gauge operators under bracket.
"""

import sys, os
from fractions import gcd, Fraction
from random import choice

import numpy
from numpy import concatenate

import networkx as nx

from solve import get_reductor, array2, row_reduce, dot2, shortstr, zeros2, shortstrx, eq2
from solve import u_inverse, find_kernel, find_logops, identity2, solve, rank
from solve import check_conjugate, check_commute
from isomorph import write

import models
from models import genidx
import cbracket
cbracket = cbracket.bracket


def cross(itemss):
    if len(itemss)==0:
        yield ()
    else:
        for head in itemss[0]:
            for tail in cross(itemss[1:]):
                yield (head,)+tail




def bracket(a, b):
    #return cbracket(a.tostring(), b.tostring())
    c1 = cbracket(a.tostring(), b.tostring())
    a, b = a.copy(), b.copy()
    assert len(a)==len(b)
    assert len(a)%2==0
    n = len(a)//2
    a.shape = (n, 2)
    b.shape = (n, 2)
    #print a
    #print b
    b = b[:, [1, 0]]
    #print b
    c = a*b
    c = c.sum() % 2
    assert c==c1
    return c


def halfbracket(a, b):
    a, b = a.copy(), b.copy()
    assert len(a)==len(b)
    assert len(a)%2==0
    n = len(a)//2
    a.shape = (n, 2)
    a = a[:, 0]
    b.shape = (n, 2)
    b = b[:, 1]
    c = a*b
    c = c.sum() % 2
    return c


def is_zop(a):
    a = a.copy()
    a.shape = (len(a)//2, 2)
    assert a.sum()
    return a[:, 0].sum() == 0


def get_xop(a):
    a = a.copy()
    a.shape = (len(a)//2, 2)
    assert a[:, 1].sum()==0
    a = a[:, 0]
    a.shape = (len(a),)
    return a

def get_zop(a):
    a = a.copy()
    a.shape = (len(a)//2, 2)
    assert a[:, 0].sum()==0
    a = a[:, 1]
    a.shape = (len(a),)
    return a


def mkop(xop, zop):
    if xop is None:
        xop = zeros2(len(zop))
    if zop is None:
        zop = zeros2(len(xop))
    xop, zop = xop.view(), zop.view()
    n = len(xop)
    xop.shape = 1,n
    zop.shape = 1,n
    op = numpy.concatenate((xop, zop))
    op = op.transpose().copy()
    op.shape = (2*n,)
    return op


def find_ideal(gen, ops=None):
    "find ideal generated by gen, inside lie algebra generated by ops"
    if ops is None:
        ops = list(gen)
    gen = list(gen)
    found = set(op.tostring() for op in gen)
    while 1:
        new = []
        for a in gen:
          _a = a.tostring()
          for b in ops:
            #if bracket(a, b):
            if cbracket(_a, b.tostring()):
                c = (a+b)%2
                s = c.tostring()
                if s not in found:
                    found.add(s)
                    new.append(c)
        if not new:
            break
        gen.extend(new)
        if argv.verbose:
            print len(gen),
            sys.stdout.flush()
    if argv.verbose:
        print
    return gen

closure = find_ideal


class Operator(object):
    """ Sparse vector (with integer components) representing a sum of Pauli operators.
        Define addition and multiplication, so these things form an algebra.
    """
    def __init__(self, v={}):
        if type(v) is dict:
            v = dict(v)
        elif type(v) is list:
            v = {array2(v).tostring() : 1}
        else:
            v = {v.tostring() : 1}
        self.v = v
        n = 0
        if v:
            k = iter(v.keys()).next()
            k = numpy.fromstring(k, dtype=numpy.int32)
            n = len(k)
            assert n%2==0
        self.n = n//2 # number of qubits

    def __str__(self):
        items = []
        for key, value in self.v.items():
            if value != 0:
                #if isinstance(value, Fraction):
                items.append("%s*%s"%(str(value), numpy.fromstring(key, dtype=numpy.int32)))
        s = "Operator(%s)"%('+'.join(items))
        s = s.replace("+-", "-")
        return s
    __repr__ = __str__

    def __eq__(self, other):
        assert self.n==0 or other.n==0 or self.n==other.n
        op = self - other
        v = op.v
        for value in v.values():
            if value != 0:
                return False
        return True

    def __ne__(self, other):
        return not (self==other)

    def is_zero(self):
        for value in self.v.values():
            if value != 0:
                return False
        return True

    def __hash__(self):
        v = self.v
        keys = list(v.keys())
        keys.sort()
        data = tuple((key, v[key]) for key in keys)
        return hash(data)

    def __add__(self, other):
        assert self.n==0 or other.n==0 or self.n==other.n
        v = dict(self.v)
        for key, value in other.v.items():
            v[key] = v.get(key, 0) + value
        op = Operator(v)
        return op
        
    def __sub__(self, other):
        assert self.n==0 or other.n==0 or self.n==other.n
        v = dict(self.v)
        for key, value in other.v.items():
            v[key] = v.get(key, 0) - value
        op = Operator(v)
        return op

    def __rmul__(self, r):
        v = {}
        assert Fraction(r) == r
        r = Fraction(r)
        for key, value in self.v.items():
            v[key] = r*value
        op = Operator(v)
        return op

    def __div__(self, r):
        assert Fraction(r) == r
        r = Fraction(r)
        return (1/r)*self

    def __neg__(self):
        return -1*self

    def __mul__(self, other):
        assert self.n==0 or other.n==0 or self.n==other.n
        v = {}
        for k1, v1 in self.v.items():
          if v1==0:
            continue
          k1 = numpy.fromstring(k1, dtype=numpy.int32)
          for k2, v2 in other.v.items():
            if v2==0:
                continue
            k2 = numpy.fromstring(k2, dtype=numpy.int32)
            sign = -1 if halfbracket(k1, k2) else +1
            k = (k1+k2)%2 # product of Pauli operators (up to sign)
            #print "__mul__", v1, k1, v2, k2, sign*v1*v2, k
            k = k.tostring()
            v[k] = v.get(k, 0) + sign*v1*v2
        op = Operator(v)
        return op

    def reduce(self):
        values = [abs(v) for v in self.v.values()]
        values = list(set(values))
        if not values:
            return self
        factor = reduce(gcd, values)
        v = {}
        for key, value in self.v.items():
            assert value % factor == 0
            v[key] = value//factor
        op = Operator(v)
        return op

    def bracket(self, other):
        return self*other - other*self

    def commutes(self, other):
        return self*other == other*self

    def anticommutes(self, other):
        return self*other != other*self

    def dot(self, other):
        assert self.n==0 or other.n==0 or self.n==other.n
        r = 0
        for k1, v1 in self.v.items():
            v2 = other.v.get(k1, 0)
            r += v1*v2
        return r

    def scale(self, other):
        "find r such that r*self == other"
        if other.is_zero():
            return 0
        assert self.n==0 or other.n==0 or self.n==other.n
        r = None
        for k1, v1 in self.v.items():
            if v1 == 0:
                continue
            v2 = other.v.get(k1, 0)
            assert v2
            r = Fraction(v2, v1)
            break
        assert r*self == other
        if r.denominator==1:
            r = r.numerator
        return r
        


def test_algebra():

    II = Operator([0,0,0,0])
    XI = Operator([1,0,0,0])
    XX = Operator([1,0,1,0])
    IX = Operator([0,0,1,0])
    ZI = Operator([0,1,0,0])
    ZZ = Operator([0,1,0,1])
    IZ = Operator([0,0,0,1])
    ops = [II, XI, XX, IX, ZI, ZZ, IZ]

    assert II==II
    assert II!=XI
    assert (II+XI) == (XI+II)
    assert (XI+XI) == 2*XI
    assert (XI+XI)/2 == XI

    for op in ops:
        assert op*op == II
    
    assert (XI * IX) == XX
    assert (XI * ZI) == -(ZI*XI)
    assert (XI * IZ) == (IZ * XI)

    assert (II - XI) * ZZ == ZZ - (XI*ZZ)

    assert hash(XI * IZ) == hash(IZ * XI)
test_algebra()


class Algebra(object):
    """ List of Operator objects that form the basis for an algebra.
    """
    def __init__(self, ops):
        self.ops = ops # basis

    def __str__(self):
        return "%s([%s])"%(
            self.__class__.__name__,
            ', '.join(str(op) for op in self.ops))

    def __getitem__(self, idx):
        return self.ops[idx]

    def __len__(self):
        return len(self.ops)

    def components(self, g):
        v = []
        for h in self.ops:
            r = h.dot(g)
            v.append(r)
        return tuple(v)


class Cartan(Algebra):
    "Cartan subalgebra, represented by a basis of operators."

    def get_eig(self, g):
        r""" find \lambda, a root, such that:
            [h, g] = \lambda(h) g for h in Cartan
        """
        root = []
        for h in self.ops:
            g1 = h.bracket(g)
            if g1.is_zero():
                root.append(0)
            elif g1 == 2*g:
                root.append(2)
            elif g1 == -2*g:
                root.append(-2)
            else:
                assert 0, "not an eigenvalue"
        root = tuple(root)
        return root

    def act(self, root, h):
        #print "act:", root, h
        v = self.components(h)
        #print "v:", v
        r = 0
        for i in range(len(self)):
            r += v[i] * root[i]
        return r


def build_roots(ops, hamiltonian):
    """ build all roots using orbits of X-type stabilizers.
    """

    assert ops
    I = Operator([0]*len(ops[0]))
    zero = Operator()

    cartan = Cartan([Operator(op) for op in ops if is_zop(op)])
    gops = [Operator(op) for op in ops if not is_zop(op)]

    print "build_roots", len(cartan), len(gops)

    eigs = set()
    roots = set()
    root_space = {} # map root -> operator
    for g in gops:
        hs = [h for h in cartan if h.anticommutes(g)]

        #print "hs:", len(hs)

        for signs in cross([(+1, -1)]*len(hs)):
            #print signs
            g1 = g
            for i, h in enumerate(hs):
                assert h != zero
                #assert h.anticommutes(g1)
                if signs[i] == 1:
                    g1 = (I+h)*g1
                else:
                    assert signs[i] == -1
                    g1 = (I-h)*g1
                #assert g1 != zero, signs # hmmm...
                #if g1 == zero:
                #    break

            if g1 == zero:
                #print "zero"
                #break # <-----------------
                continue

            if g1 in eigs or -g1 in eigs:
                continue

            eigs.add(g1)
            root = cartan.get_eig(g1)
    
            if root not in roots:
                roots.add(root)
                root_space[root] = g1
                r = sum(r**2 for r in root)
                print r, root
                if r==0:
                    print g
                    print g1
                    print g2


    print "eigs:", len(eigs)
    print "roots:", len(roots)

    roots = list(roots)
    lookup = dict((root, i) for i, root in enumerate(roots))

    algebra = Algebra(eigs)
    hc = algebra.components(hamiltonian)
    hc = [(-r if r<0 else r) for r in hc]
    print "Hamiltonian:", hc

    from action import Perm, Group

    perms = []
    for alpha in roots:
        nalpha = tuple(-a for a in alpha)

        X = root_space[alpha]
        Y = root_space[nalpha]
        H = X.bracket(Y)
        #print "X:", X
        #print "Y:", Y
        #print "H:", H
    
        r = X.scale(H.bracket(X)) # r*X == H.bracket(X)
        r = Fraction(2, r)
        H = r*H
        assert H.bracket(X) == 2*X
        assert H.bracket(Y) == -2*Y
        #print X.bracket(Y) == H # doesn't matter
    
        #print cartan
        #print cartan.components(H)
        assert cartan.act(alpha, H) == 2
    
        alpha = numpy.array(alpha)
        #print "alpha:", alpha
        perm = []
        for beta in roots:
            c = cartan.act(beta, H)
            assert int(c)==c
            c = int(c)
            #print beta, c,
            root = numpy.array(beta) - c*alpha
            root = tuple(root)
            #print root
            assert root in roots
            perm.append(lookup[root])
        #print [hc[i] for i in perm]
        #if hc == [hc[i] for i in perm]:
        #    print "*"
        assert len(set(perm))==len(roots)
        perms.append(perm)

    print "perms:", len(perms)
    # Note: Weyl group is generated by reflections of the simple roots.
    # So the nodes in the Dynkin diagram correspond to the simple roots.

    N = len(roots)
    if N <= 24:
        items = range(N)
        gen = [Perm(tuple(perm), items) for perm in perms]
        group = Group.generate(gen)
        print "Weyl group:", len(group)
        count = 0
        for g in group:
            if [hc[g[i]] for i in items] == hc:
                #print g
                count += 1
        print "Hamiltonian stabilizer subgroup:", count

    else:
        items = range(N)
        gen = [Perm(tuple(perm), items) for perm in perms]
        #group = Group.generate(gen)
        #print "Weyl group:", len(group)
        count = 0
        for g in gen:
            if [hc[g[i]] for i in items] == hc:
                #print g
                count += 1
        print "Hamiltonian stabilizer gen:", count
        
    lengths = {}
    for root in roots:
        r = sum(r**2 for r in root)
        lengths[r] = lengths.get(r, 0) + 1
    if len(lengths)>1:
        keys = lengths.keys()
        assert len(keys)==2, lengths
        keys.sort()
        short, long_ = keys
        print "short roots:", lengths[keys[0]]
        print "long  roots:", lengths[keys[1]]


def search_roots(ops, hamiltonian):
    """ Search for roots as in build_roots but use a backtracking
        algorithm. Does not find all roots.
    """

    assert ops
    I = Operator([0]*len(ops[0]))
    zero = Operator()

    cartan = Cartan([Operator(op) for op in ops if is_zop(op)])
    gops = [Operator(op) for op in ops if not is_zop(op)]

    print "search_roots", len(cartan), len(gops)

    eigs = set()
    roots = set()
    root_space = {} # map root -> operator
    for g in gops:
        hs = [h for h in cartan if h.anticommutes(g)]

        #print "hs:", len(hs)
        #for signs in cross([(+1, -1)]*len(hs)):

        signs = [+1]
        stack = [g]
        idx = 0
        while 1:

            #print "idx =", idx

            assert len(stack)==idx+1
            assert len(signs)==idx+1
            g = stack[idx]

            h = hs[idx]
            if signs[idx] == 1:
                g1 = (I+h)*g
            else:
                assert signs[idx] == -1
                g1 = (I-h)*g

            if g1.is_zero():

                if signs[idx] == 1:
                    signs[idx] = -1
                    continue

                # backtrack
                while idx>=0 and signs[idx] == -1:
                    idx -= 1
                    stack.pop()
                    signs.pop()

                if idx>=0:
                    signs[idx] = -1
                    continue

                break # finished search

            if idx+1 < len(hs):
                stack.append(g1)
                signs.append(+1)
                idx += 1
                continue

            #if g1 in eigs or -g1 in eigs:
            #    continue

            if g1 not in eigs and -g1 not in eigs:
                eigs.add(g1)
                root = cartan.get_eig(g1)

                assert root not in roots
                roots.add(root)
                root_space[root] = g1
                r = sum(r**2 for r in root)
                print "%d: %s"%(r, ' '.join("%2d"%ri for ri in root))
                if r==0:
                    print g
                    print g1
                    print g2

            if not argv.all_roots:
                break

            if 1:
                if signs[idx] == 1:
                    signs[idx] = -1
                    continue

                # backtrack
                while idx>=0 and signs[idx] == -1:
                    idx -= 1
                    stack.pop()
                    signs.pop()

                if idx>=0:
                    signs[idx] = -1
                    continue

                break # finished search

    print "eigs:", len(eigs)
    print "roots:", len(roots)



def report_lie(cartan, ideal):
    n = len(cartan) # rank
    m = len(ideal)
    print "ideal:", len(ideal), "cartan:", len(cartan)
    if m == n**2 + 2*n:
        print "A_%d = sl_%d" % (n, n+1)
    elif m == 2*n**2 + n:
        print "B_%d = so_%d, *OR*" % (n, 2*n+1)
        print "C_%d = sp_%d" % (n, 2*n)
    elif m == 2*n**2 - n:
        print "D_%d = so_%d" % (n, 2*n)
    else:
        print "Unkown lie algebra"
    if n==2 and m==14:
        print "G_2"
    elif n==4 and m==52:
        print "F_4"
    elif n==6 and m==78:
        print "E_6"
    elif n==7 and m==133:
        print "E_7"
    elif n==8 and m==248:
        print "E_8"



def test_model():

    Gx, Gz, Hx, Hz = models.build()
    model = models.build_model(Gx, Gz, Hx, Hz)

    Rx, Rz = model.Rx, model.Rz
    Rxt = Rx.transpose()
    Rzt = Rz.transpose()

    Px, Pz = model.Px, model.Pz
    Pxt = Px.transpose()
    Pzt = Pz.transpose()

    r, n = Rx.shape

    #print shortstrx(Rx, Rz)
    print shortstrx(Gx, Gz)
    print "Gx:", len(Gx), "Gz:", len(Gz), "Rx/z:", len(Rx)

    print shortstrx(Hx, Hz)
    print "Hx:", len(Hx), "Hz:", len(Hz)

    #RR = dot2(Rz, Rx.transpose()) == I
    #assert rank(RR)==r

    ops = []
    found = set()
    gops = []
    hamiltonian = []
    for gx in Gx:
        rx = dot2(gx, Rzt)
        rx = mkop(rx, None)
        hamiltonian.append(rx.tostring())
        if rx.sum()==0: # stabilizer..
            continue
        s = rx.tostring()
        if s not in found:
            found.add(s)
            ops.append(rx)
            gops.append(rx)

    cartan = []
    for gz in Gz:
        rz = dot2(gz, Rxt)
        rz = mkop(None, rz)
        hamiltonian.append(rz.tostring())
        if rz.sum()==0: # stabilizer..
            continue
        s = rz.tostring()
        if s not in found:
            found.add(s)
            ops.append(rz)
            cartan.append(rz)

    hamiltonian = Operator(dict((op, 1) for op in hamiltonian))

    if argv.closure:
        ops = closure(ops)
        print "algebra dimension:", len(ops)

    if argv.closure and argv.test:
        n = len(ops)
        lookup = dict((ops[i].tostring(), i) for i in range(n))
        As = []
        for h in cartan:
          A = numpy.zeros((n, n))
          for i, a in enumerate(ops):
            if cbracket(h, a):
                b = (a+h)%2
                j = lookup[b.tostring()]
                A[i, j] = 1.
          As.append(A)
          #print shortstr(A)
        for A in As:
          for B in As:
            assert numpy.allclose(numpy.dot(A, B), numpy.dot(B, A))

    if argv.render_ideals:
        import gauge
        model = gauge.make(Gx, Gz, Hx, Hz)

        found = set()
        i = 0
        for gx in Gx:
            #rx = dot2(gx, Rzt)
            gx = dot2(gx, Pxt)
            s = gx.tostring()
            if s in found:
                continue
            found.add(s)
            idxs = numpy.where(gx)[0]
            model.render(show_idxs=idxs)
            i += 1
            if i%6==0:
                model.newline()
        model.save("pic-gcolor-gauge.pdf")

        remain = list(gops)
        i = 0
        while remain:
            #print "remain:", len(remain)
            op = remain.pop(0)
            ideal = find_ideal([op], ops)
            print "ideal:", len(ideal)
            count = 0
            s_ideal = set(str(op) for op in ideal)

            for op in gops:
                if str(op) in s_ideal:
                    rx = get_xop(op)
                    gx = dot2(rx, Rx)
                    print "   ", shortstr(gx)
                    idxs = numpy.where(gx)[0]
                    model.render(show_idxs=idxs)

            #model.save("pic-ideal-%d.pdf"%i)
            model.newline()
            i += 1

            print "count:", len([1 for op in gops if str(op) in s_ideal])
            print
            ideal = set(op.tostring() for op in ideal)
            remain = [op for op in remain if not op.tostring() in ideal]
        model.save("pic-gcolor-ideals.pdf")

        return

    found = set()
    i = 0
    for gx in Gx:
        #rx = dot2(gx, Rzt)
        gx = dot2(gx, Pxt)
        s = gx.tostring()
        if s in found:
            continue
        found.add(s)


    # ----------- Find the ideals ------------------

    nh = 0
    remain = list(gops)
    i = 0
    while remain:
        #print "remain:", len(remain)
        op = remain.pop(0)
        ideal = find_ideal([op], ops)
        cartan = [op for op in ideal if is_zop(op)]
        report_lie(cartan, ideal)

        if argv.roots or argv.build_roots:
            build_roots(ideal, hamiltonian)

        if argv.search_roots:
            search_roots(ideal, hamiltonian)

        nh += len(cartan)
        count = 0
        s_ideal = set(str(op) for op in ideal)

        gxs = []
        for gx in Gx:
            rx = dot2(gx, Rzt)
            rx = mkop(rx, None)
            if str(rx) in s_ideal:
                #rx = get_xop(op)
                #gx = dot2(rx, Rx)
                gxs.append(rx)
        gzs = []
        for gz in Gz:
            rz = dot2(gz, Rxt)
            rz = mkop(rz, None)
            if str(rz) in s_ideal:
                #rz = get_zop(op)
                #gz = dot2(rz, Rz)
                gzs.append(rz)

        if argv.show:
            gxs = array2(gxs)
            gzs = array2(gzs)
            print "Gx:"
            print shortstrx(gxs)
            print "Gz:"
            print shortstrx(gzs)
            A = dot2(gxs, gzs.transpose())
            #print shortstrx(gxs, gzs, A)
            print "A.transpose():"
            print shortstr(A.transpose())

        ideal = set(op.tostring() for op in ideal)
        remain = [op for op in remain if not op.tostring() in ideal]

        if argv.single_ideal:
            break

    return

    lookup = {}
    for i, op in enumerate(ops):
        lookup[op.tostring()] = i

    N = len(ops)

#    graph = nx.Graph()
#    for i in range(N):
#        graph.add_node(i)
#    for i, A in enumerate(ops):
#      for j, B in enumerate(ops):
#        if bracket(A, B)==0:
#            continue
#        C = (A+B)%2
#        k = lookup[C.tostring()]
#        graph.add_edge(i, k)
#        graph.add_edge(j, k)
#    equs = nx.connected_components(graph)
#    print "ideals:", len(equs), [len(equ) for equ in equs]

    if argv.stop2:
        return

    graph = nx.Graph()
    for i in range(N):
        graph.add_node(i)

    H = []
    for z in cartan:
        A = zeros2(N, N)
        for i, op in enumerate(ops):
            if is_zop(op):
                continue
            c = bracket(z, op)
            if c==0:
                continue
            op1 = (z+op)%2
            j = lookup[op1.tostring()]
            #print "%s->%s" % (i, j),
            A[j, i] = 1
            graph.add_edge(j, i)
        #print
        #print shortstr(A)
        #print
        H.append(A)
    
    for A in H:
      for B in H:
        assert numpy.allclose(numpy.dot(A, B), numpy.dot(B, A))

    equs = nx.connected_components(graph)
    print "orbits:", len(equs), [len(equ) for equ in equs]
    trivial = len([equ for equ in equs if len(equ)==1])
    print "trivial:", trivial, "non-trivial", len(equs)-trivial

    #for irrep in genidx((2,)*r):
    #    print "irrep:", irrep

    return

    r, n = Rx.shape
    N = 2**r
    gz = len(Gz)

    RR = dot2(Gz, Rx.transpose())
    PxtRzt = dot2(Px.transpose(), Rzt)

    if N<=1024:
        H = numpy.zeros((N, N))
    else:
        H = None
    A = {}
    U = []

    basis = []
    lookup = {}
    for i, v in enumerate(genidx((2,)*r)):
        v = array2(v)
        lookup[v.tostring()] = i
        basis.append(v)

    ops = []

    # zops
    for gz in Gz:
        elems = {}
        for i, v in enumerate(basis):
            bit = dot2(gz, Rx.transpose(), v)
            elems[(i, i)] = 1 - 2*bit # send 0,1 -> +1,-1

    # xops
    for gx in Gx:
        elems = {}
        for i, v in enumerate(basis):
            u = (v+dot2(gx, PxtRzt))%2
            j = lookup[u.tostring()]
            elems[j, i] = 1
        for (i, j) in elems:
            assert elems[i, j] == elems[j, i]





from argv import Argv

argv = Argv()

if __name__ == "__main__":

    fn = "test_model"

    if fn is None:
        pass

    elif argv.profile:
        import cProfile as profile
        profile.run("%s()"%fn)

    else:
        fn = eval(fn)

        while 1:
            fn()
            if not argv.forever:
                break


