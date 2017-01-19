#!/usr/bin/env python

import sys
from heapq import heappush, heappop, heapify
from random import random

import numpy
import scipy.sparse.linalg as la


try:
    from pyx import canvas, path, deco, trafo, style, text, color, deformer
    from pyx.color import rgb 
    
    if 0:
        text.set(mode="latex") 
        text.set(docopt="12pt")
        text.preamble(r"\usepackage{amsmath,amsfonts,amssymb}")
        
        text.preamble(r"\def\ket #1{|#1\rangle}")

    black = rgb(0., 0., 0.)
    blue = rgb(0., 0., 0.8)
    lred = rgb(1., 0.4, 0.4)
    white = rgb(1., 1., 1.)
    grey = rgb(0.75, 0.75, 0.75)
    shade = grey
    shade0 = rgb(0.75, 0.75, 0.75)
    shade1 = rgb(0.80, 0.80, 0.80)
    shade2 = rgb(0.85, 0.85, 0.85)
    
    light_shade = rgb(0.85, 0.65, 0.1)
    light_shade = rgb(0.9, 0.75, 0.4)
    
    north = [text.halign.boxcenter, text.valign.top]
    northeast = [text.halign.boxright, text.valign.top]
    northwest = [text.halign.boxleft, text.valign.top]
    south = [text.halign.boxcenter, text.valign.bottom]
    southeast = [text.halign.boxright, text.valign.bottom]
    southwest = [text.halign.boxleft, text.valign.bottom]
    east = [text.halign.boxright, text.valign.middle]
    west = [text.halign.boxleft, text.valign.middle]
    center = [text.halign.boxcenter, text.valign.middle]
    
    st_dashed = [style.linestyle.dashed]
    st_dotted = [style.linestyle.dotted]
    st_round = [style.linecap.round]
    
    st_thick = [style.linewidth.thick]
    st_Thick = [style.linewidth.Thick]
    st_THICK = [style.linewidth.THICK]
    

except ImportError:
    print "no pyx module found"



from code import texstr

EPSILON = 1e-8
scalar = numpy.float64

def write(s):
    sys.stdout.write(str(s)+' ')
    sys.stdout.flush()


class Op(la.LinearOperator):

    verbose = False
    def __init__(self, n):
        self.n = n
        la.LinearOperator.__init__(self, 
            (n, n),
            dtype=scalar, 
            matvec=self.matvec)

    def __len__(self):
        return self.n

    def __mul__(self, other):
        return MulOp(self, other)

    def __add__(self, other):
        return SumOp(self.n, [self, other])

    def __sub__(self, other):
        return SumOp(self.n, [self, -1.0*other])

    def __rmul__(self, alpha):
        return RMulOp(self, alpha)

    def matvec(self, v, w=None, verbose=False):
        if self.verbose:
            write('.')
        assert len(v) == self.n
        if w is None:
            return v
        w += v
        return w

    def todense(self):
        n = self.n
        A = numpy.zeros((n, n))
        v = numpy.zeros(n)
        for i in range(n):
            v[i] = 1
            A[:, i] = self.matvec(v)
            v[i] = 0
        return A


class IdOp(Op):
    pass


class XOp(Op):
    def __init__(self, n, idxs, alpha=1.0):
        " X operator to specified qubits "

        if type(idxs) in (int, long):
            idxs = [idxs]
        for idx in idxs:
            assert 0<=idx<n
        assert len(set(idxs))==len(idxs), idxs
        perm = []
        assert n<=30
        ns = range(n)
        mask = [1 if i in idxs else 0 for i in ns]
        #print "XOp:", mask
        for i in range(2**n):
            #print i,
            bits = []
            for flip in mask:
                bit = (1 - i%2) if flip else i%2
                bits.append(bit)
                i >>= 1
            #print bits,
            j = 0
            for bit in reversed(bits):
                j += bit
                j <<= 1
            j >>= 1
            #print j
            perm.append(j)
        self.perm = perm
        self.idxs = list(idxs)
        self.alpha = alpha
        Op.__init__(self, 2**n)

    def __str__(self):
        return "XOp(%d, %s, %s)"%(self.n, self.idxs, self.alpha)

    def matvec(self, v, w=None, verbose=False):
        if self.verbose:
            write('.')
        assert len(v) == len(self.perm)
        if w is None:
            w = numpy.zeros(len(v))
        #for i, j in enumerate(self.perm):
        #    w[j] += v[i]
        v = v[self.perm]
        w += self.alpha * v
        return w


class ZOp(Op):
    def __init__(self, n, idxs, alpha=1.0):
        " Z operator to specified qubits "
        if type(idxs) in (int, long):
            idxs = [idxs]
        for idx in idxs:
            assert 0<=idx<n
        assert len(set(idxs))==len(idxs), idxs
        phases = []
        assert n<=30
        ns = range(n)
        mask = [1 if i in idxs else 0 for i in ns]
        #print "XOp:", mask
        for i in range(2**n):
            #print i,
            phase = 1
            for flip in mask:
                if flip and i%2:
                    phase *= -1
                i >>= 1
            phases.append(phase)
        #print "ZOp:", idxs, phases
        self.phases = numpy.array(phases)
        self.idxs = list(idxs)
        self.alpha = alpha
        Op.__init__(self, 2**n)

    def __str__(self):
        return "ZOp(%d, %s, %s)"%(self.n, self.idxs, self.alpha)

    def matvec(self, v, w=None, verbose=False):
        if self.verbose:
            write('.')
        assert len(v) == len(self.phases)
        if w is None:
            w = numpy.zeros(len(v))
        #for i, phase in enumerate(self.phases):
        #    w[i] += phase * v[i]
        w += self.alpha * self.phases * v
        return w


def mkop(tp, s):
    n = len(s)
    idxs = []
    if type(s) is str:
        for i, x in enumerate(s):
            if s[i]=='1':
                idxs.append(i)
    else:
        for i, x in enumerate(s):
            assert s[i] in [0,1]
            if s[i]:
                idxs.append(i)
    return tp(n, idxs)


class SumOp(Op):
    def __init__(self, n, ops):
        self.ops = ops
        self.n = n
        self.count = 0
        Op.__init__(self, self.n)

    def matvec(self, v, w=None, verbose=True):
        #print "SumOp.matvec"
        if self.verbose:
            write('.')
        assert len(v)==self.n
        if w is None:
            w = numpy.zeros(len(v))
        for op in self.ops:
            op.matvec(v, w)
        self.count += 1
        return w


class RMulOp(Op):
    def __init__(self, op, alpha=1.0):
        self.op = op
        self.alpha = alpha
        Op.__init__(self, op.n)

    def matvec(self, v, w=None, verbose=True):
        if self.verbose:
            write('.')
        assert len(v)==self.n
        if w is None:
            w = numpy.zeros(len(v))
        w += self.alpha * self.op.matvec(v, verbose=verbose)
        return w


class MulOp(Op):
    def __init__(self, a, b):
        # first do b then a !!
        self.ops = a, b
        assert a.n == b.n, (a.n, b.n)
        Op.__init__(self, a.n)

    def matvec(self, v, w=None, verbose=True):
        if self.verbose:
            write('.')
        assert len(v)==self.n
        if w is None:
            w = numpy.zeros(len(v))
        a, b = self.ops
        v = b.matvec(v, verbose=verbose)
        v = a.matvec(v, verbose=verbose)
        w += v
        return w


from smap import SMap
def strop(l, idx):
    keys = [(i, j) for i in range(l) for j in range(l)]
    s = SMap()
    bits = []
    for i in range(l):
      for j in range(l):
        s[i, j] = 'X' if idx%2 else '.'
        idx //= 2
    return str(s)
    

class Model(object):
    def get_syndrome(self, v):
        syndrome = []
        for op in self.xstabs+self.zstabs:
            v1 = op.matvec(v)
            r = numpy.dot(v, v1) # real values
            syndrome.append(r)
        return syndrome


from models import build_model
class CodeModel(Model):

    def __init__(self, code=None):
        if code is None:
            code = build_model()
    
        self.n = code.n
        self.xops = [mkop(XOp, g) for g in code.Gx]
        self.zops = [mkop(ZOp, g) for g in code.Gz]
        self.xstabs = [mkop(XOp, g) for g in code.Hx]
        self.zstabs = [mkop(ZOp, g) for g in code.Hz]
        self.A = SumOp(2**code.n, self.xops+self.zops)




def plot(model, v):

    import networkx as nx

    graph = nx.Graph()

    dim = len(v)
    syndrome = numpy.zeros(dim)
    for op in model.zops:
        syndrome += op.phases

    idxs = [i for i in range(dim) if abs(v[i])>EPSILON]
    print "nnz:", len(idxs)

    for i in idxs:
        graph.add_node(i, value=v[i])

    for i in idxs:
        for xop in model.xops:
            j = xop.perm[i]
            graph.add_edge(i, j)
            graph.add_edge(j, i)

    #idxs.sort(key = lambda i : -syndrome[i]) # stable sort
    idxs.sort(key = lambda i : -abs(v[i]))

    row = []
    rows = [row]
    val = v[idxs[0]]
    for idx in idxs:
        val1 = v[idx]
        #print idx, val1
        if abs(abs(val)-abs(val1))>1e-6:
            assert row
            val = val1
            row = []
            rows.append(row)
            #print idx, val1
        row.append((idx, val1))

    print "species:", len(rows)
    for row in rows:
        print len(row),
    print

    row = rows[0]
    n0 = len(row)
    row.sort(key = lambda (idx,val):-val)

    paths = {}
    for idx,_ in row:
        assert idx in graph, idx
        ps = nx.shortest_path(graph, target=idx) # paths to idx
        ps = dict((jdx, len(path)) for jdx,path in ps.items())
        #print idx, ps
        paths[idx] = ps

    for row in rows[1:]:
        n = len(row)
        closest = []
        targets = {}
        for i in range(n):
            idx,_ = row[i]
            jdxs = [jdx for jdx,_ in rows[0]]
            #jdxs.sort(key = lambda jdx,idx=idx : len(nx.shortest_path(graph, idx, jdx)))
            jdxs.sort(key = lambda jdx,idx=idx : paths[jdx].get(idx, 9999))
            targets[idx] = jdxs[0] # closest guy in top row
        #print targets
        row.sort(key = lambda (idx,val) : targets[idx])

    r = 0.05
    dy = 0.2
    dx = 1.2*r

    c = canvas.canvas()

    x, y = 0., 0.
    for row in rows:
        x = -len(row)*0.5*dx
        for idx, val in row:
            #p = path.circle(x, y, r)
            p = path.rect(x, y, r, r)
            if val > 0.:
                #print "F",
                c.fill(p)
            c.stroke(p)
            x += dx
        y -= dy

    c.writePDFfile("pic-wavefunction.pdf")


def save(v, fname):
    a = v.tostring()
    f = open(fname, 'w')
    f.write(a)
    f.close()


def load(fname):
    f = open(fname)
    s = f.read()
    f.close()
    a = numpy.fromstring(s)
    return a


def commutes(A, B, sign=-1):
    n = A.n
    v = numpy.random.normal(size=n)

    Av = A.matvec(v)
    BAv = B.matvec(Av)

    Bv = B.matvec(v)
    ABv = A.matvec(Bv)

    r = numpy.abs(BAv + sign * ABv).sum()

    return r < 1e-6


def anticommutes(A, B):
    return commutes(A, B, +1)


def show_eigs(vals):
    vals = list(vals)
    vals.sort(reverse=True)
    counts = {}
    val0 = vals[0]-10
    for val in vals:
        if abs(val-val0) < 1e-6:
            counts[val0] += 1
        else:
            val0 = val
            counts[val0] = 1
    vals = counts.keys()
    vals.sort(reverse=True)
    print "eigval, degeneracy:"
    for val in vals:
        print '\t', val, counts[val]
    


def test():

    norm = lambda v : (v**2).sum()**0.5

    model = CodeModel()

    if argv.test:
        stabs = model.xstabs+model.zstabs
        for A in stabs:
            assert commutes(A, model.xlogop)
            assert commutes(A, model.zlogop)
            for B in stabs:
                assert commutes(A, B)
        #assert anticommutes(model.xlogop, model.zlogop)
        print "OK"
        return

    if argv.threads:
        os.environ['OMP_NUM_THREADS'] = str(argv.threads)

    dim = model.A.n

    v0 = None

    if argv.stabdist:
        v0 = numpy.zeros(len(model.A))
        v0[0] = 1
    
    k = argv.get("k", 2)

    if k=="all":
        k = A.n

    A = model.A

    if argv.perturb:
#        # doesn't work very well...
#        alpha = 0.0001
#        perturb = [alpha * op for op in model.xstabs]
#        A = SumOp(2**model.n, A.ops+perturb)

        alpha = argv.get("alpha", 1e-4)

        def rnd(alpha):
            return alpha * (2*random() - 1.)

        ops = []
        for i in range(model.n):
            ops.append(rnd(alpha)*XOp(model.n, [i]))
            ops.append(rnd(alpha)*ZOp(model.n, [i]))
        A = SumOp(2**model.n, A.ops + ops)


    projs = []
    I = IdOp(A.n)
    flip = argv.get("flip", 0)
    zflip = argv.get("zflip", 0)
    if zflip:
        stabs = model.zstabs + model.xstabs
        flip = zflip
    else:
        stabs = model.xstabs + model.zstabs

    for op in stabs:
        if flip:
            op = 0.5 * (I - op)
            flip -= 1
        else:
            op = 0.5 * (I + op)
        projs.append(op)

    if projs:
        P = projs[0]
        for i in range(1, len(projs)):
            P = P * projs[i]

    if argv.stabilize:
        A = P*A*P

    if argv.logop:
        P = 0.5 * (I + model.xlogop)
        A = P*A*P

    if argv.exact:
        assert A.n < 2**14
        A = A.todense()
        vals, vecs = numpy.linalg.eigh(A)

    elif argv.power:

        v = numpy.random.normal(size=A.n)
        v /= norm(v)
    
        sigma = 1.0
    
        while 1:
    
            u = A.matvec(v)
            eigval = numpy.dot(v, u)
    
            u = u + sigma * v
    
            r = norm(u)
            if r>EPSILON:
                u /= r
    
            err = norm(u-v)
            print "delta:", err 
            print "eig:", eigval
    
            if err < 1e-4:
                break
    
            v = u
    

    else:

        A.verbose = True
        which = argv.get("which", 'LA')
        vals, vecs = la.eigsh(A, k=k, v0=v0, which=which, maxiter=None) #, tol=1e-8)
        #vals, vecs = la.eigs(A, k=k, v0=v0, which='LR', maxiter=None) #, tol=1e-8)
        print
    
    # vals go from smallest to highest
    show_eigs(vals)
    #print "iterations:", model.A.count


    if argv.verbose:
        for i in range(k):
            print i, "syndrome", model.get_syndrome(vecs[:, i])

    return

    idx = argv.get("idx", 0)
    k = vecs.shape[1]

    v0 = vecs[:, k-1-idx]

    pos, = numpy.where(v0>+EPSILON)
    neg, = numpy.where(v0<-EPSILON)
    print "pos: %d, neg: %d"%(len(pos), len(neg))
    print "v[0]", v0[0]
    print v0[pos[:10]], v0[neg[:10]]

#    v0 = numpy.abs(v0)
#
#    count = 0
#    idxs = []
#    values = []
#    for i in range(dim):
#        if abs(v0[i]) > EPSILON:
#            #write(i)
#            idxs.append(i)
#            values.append(v0[i])
#    print "nnz:", len(idxs)

    if argv.plot:
        plot(model, v0)

    return

    propagate = SumOp(2**model.n, model.xops)

    syndrome = numpy.zeros(dim)
    for op in model.zops:
        syndrome += op.phases

    idxs = [i for i in range(dim) if v0[i]>EPSILON]
    idxs.sort(key = lambda i : -syndrome[i])

    value = v0[idxs[0]]
    best = []
    for idx in idxs:
        if v0[idx] >= value-EPSILON:
            best.append(idx)
        else:
            break

    print "best:", best
    marked = set(best)
    for i in idxs:
        for xop in model.xops:
            j = xop.perm[i]
            if j in marked:
                continue
            marked.add(j)
            continue # <---------
            if syndrome[j] > syndrome[i]:
                #print "*", i, '->', j, 'and', syndrome[i], syndrome[j]
                print strop(l, i), syndrome[i]
                print "->"
                print strop(l, j), syndrome[j]
                print v0[i], "->", v0[j]
                print

    for i in idxs:
        for xop in model.xops:
            j = xop.perm[i]
            if syndrome[j] > syndrome[i] and v0[j]<v0[i]:
                #print "*", i, '->', j, 'and', syndrome[i], syndrome[j]
                print strop(l, i), syndrome[i]
                print "->"
                print strop(l, j), syndrome[j]
                print v0[i], "->", v0[j]
                print



from argv import Argv 
argv = Argv()


if __name__ == "__main__":

    if argv.search:
        search()
    else:
        test()




