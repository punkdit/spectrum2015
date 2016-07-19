#!/usr/bin/env python

import sys, os

import numpy
from scipy import sparse
from scipy.sparse.linalg import eigs, eigsh

from solve import shortstr, parse, eq2, dot2, zeros2, array2, identity2
from solve import row_reduce, RowReduction, span, get_reductor
from solve import u_inverse, find_logops, solve

from lanczos import write, show_eigs
from code import lstr2


def genidx(shape):
    if len(shape)==0:
        yield ()
    else:
        for idx in range(shape[0]):
            for _idx in genidx(shape[1:]):
                yield (idx,)+_idx


def check_conjugate(A, B):
    if A is None or B is None:
        return
    assert A.shape == B.shape
    I = numpy.identity(A.shape[0], dtype=numpy.int32)
    assert eq2(dot2(A, B.transpose()), I)


def check_commute(A, B):
    if A is None or B is None:
        return
    C = dot2(A, B.transpose())
    assert C.sum() == 0, "\n%s"%shortstr(C)





def build_gcolor(size):

    from qupy.ldpc import gcolor

    lattice = gcolor.Lattice(size)

    n = len(lattice.qubits)
    print lattice

    code = lattice.build_code(check=False)
    #Ex = lattice.Ex
    Gx, Gz = code.Gx, code.Gz
    Hx, Hz = code.Hx, code.Hz

    return Gx, Gz, Hx


def build_compass(l):

    n = l**2

    keys = [(i, j) for i in range(l) for j in range(l)]
    coords = {}  
    for i, j in keys:
        for di in range(-l, l+1):
          for dj in range(-l, l+1):
            coords[i+di, j+dj] = keys.index(((i+di)%l, (j+dj)%l))

    m = n 
    Gx = zeros2(m, n)
    Gz = zeros2(m, n)

    idx = 0 
    for i in range(l):
      for j in range(l):
        Gx[idx, coords[i, j]] = 1 
        Gx[idx, coords[i, j+1]] = 1 

        Gz[idx, coords[i, j]] = 1 
        Gz[idx, coords[i+1, j]] = 1 
        idx += 1

    assert idx == m

    mx = l-1
    Hx = zeros2(mx, n)
    for idx in range(l-1):
      for j in range(l):
        Hx[idx, coords[j, idx]] = 1
        Hx[idx, coords[j, idx+1]] = 1

    mz = l-1
    Hz = zeros2(mz, n)
    for idx in range(l-1):
      for j in range(l):
        Hz[idx, coords[idx, j]] = 1
        Hz[idx, coords[idx+1, j]] = 1

    assert dot2(Hx, Hz.transpose()).sum() == 0

    return Gx, Gz, Hx, Hz


def build_isomorph(Gx):
    from isomorph import Tanner, search
    bag0 = Tanner.build(Gx)
    bag1 = Tanner.build(Gx)
    m, n = Gx.shape
    count = 0
    perms = []
    keys = range(m, len(bag0))
    #print "keys:", keys
    for fn in search(bag0, bag1):
        #print fn
        perm = tuple(fn[i]-m for i in keys)
        perms.append(perm)
        count += 1
    print "isomorphisms:", count
    #print len(set(perms))
    return perms


def build_orbiham(H):
    from isomorph import from_ham, search
    import networkx as nx

    n = len(H)
    graph = nx.Graph()
    for i in range(n):
        graph.add_node(i)

    bag0 = from_ham(H)
    bag1 = from_ham(H)

    count = 0
    for fn in search(bag0, bag1):
        #print fn
        write('.')
        for i, j in fn.items():
            graph.add_edge(i, j)
        count += 1
    print

    equs = nx.connected_components(graph)
    m = len(equs)

    print "isomorphisms:", count
    print "components:", m

    P = numpy.zeros((n, m))
    Q = numpy.zeros((m, n))
    for i, equ in enumerate(equs):
        for j in equ:
            P[j, i] = 1
        Q[i, j] = 1

    #print shortstr(P)
    #print shortstr(Q)

    H = numpy.dot(Q, numpy.dot(H, P))
    return H


def sparse_orbiham_nx(n, H):
    import networkx as nx
    from networkx.algorithms.isomorphism import GraphMatcher

    bag = nx.Graph()
    for i in range(n):
        bag.add_node(i, syndrome=H[i,i])

    for i in range(n):
      for j in range(n):
        if i==j:
            continue
        if H.get((i, j)):
            bag.add_edge(i, j)

    def node_match(n0, n1):
        return n0['syndrome'] == n1['syndrome']

    matcher = GraphMatcher(bag, bag, node_match=node_match)

    print "search..."

    graph = nx.Graph()
    for i in range(n):
        graph.add_node(i)

    count = 0
    for iso in matcher.isomorphisms_iter(): # too slow :P
        #print iso
        write('.')
        for i, j in iso.items():
            graph.add_edge(i, j)
        count += 1
    print

    equs = nx.connected_components(graph)
    m = len(equs)

    print "isomorphisms:", count
    print "components:", m


def sparse_orbiham(n, H):
    from isomorph import from_sparse_ham, search
    import networkx as nx

    graph = nx.Graph()
    for i in range(n):
        graph.add_node(i)

    print "bag0"
    bag0 = from_sparse_ham(n, H)

    print "bag1"
    bag1 = from_sparse_ham(n, H)

    print "search..."

    count = 0
    for fn in search(bag0, bag1):
        #print fn
        write('.')
        for i, j in fn.items():
            graph.add_edge(i, j)
        count += 1
    print

    equs = nx.connected_components(graph)
    m = len(equs)

    print "isomorphisms:", count
    print "components:", m

    P = numpy.zeros((n, m))
    Q = numpy.zeros((m, n))
    for i, equ in enumerate(equs):
        for j in equ:
            P[j, i] = 1
        Q[i, j] = 1

    #print shortstr(P)
    #print shortstr(Q)

    H = numpy.dot(Q, numpy.dot(H, P))
    return H


def sparse_orbiham_nauty(degree, A, U):

    n = len(U)

    idxs = A.keys()
    idxs.sort()
    edges = n*degree

    rows = dict((i, []) for i in range(n))
    for (i, j) in idxs:
        rows[i].append(j)
        
    import pnauty
    print "pnauty.init_graph", n, edges, degree
    assert n*degree == edges
    pnauty.init_graph(n, edges, degree)

    i0 = None
    #for idx, weight in A.items():
    for idx in idxs: # sorted
        i, j = idx
        assert A[idx] > 0
        for count in range(A[idx]):
            if i != i0:
                edge = 0
                i0 = i
            else:
                edge += 1
            assert edge < degree
            #print "add_edge", (i, j, edge)
            pnauty.add_edge(i, j, edge)

    #pnauty.search()
    #return

    verts = range(n)
    verts.sort(key = lambda i : U[i])

    label = None
    for idx in range(n-1):
        i0 = verts[idx]
        i1 = verts[idx+1]
        if U[i0] == U[i1]:
            pnauty.set_partition(idx, i0, 1)
        else:
            pnauty.set_partition(idx, i0, 0)
    pnauty.set_partition(n-1, verts[-1], 0)

    print "pnauty.search"
    orbits = pnauty.search()
    canonical = list(set(orbits))
    canonical.sort()
    N = len(canonical)
    print "orbits:", N

    lookup = {}
    for i, idx in enumerate(orbits):
        lookup[i] = canonical.index(idx)

    H = numpy.zeros((N, N))
    for idx in range(N):
        i = canonical[idx]
        # i indexes a row of A
        for j in rows[i]:
            value = A[i, j]
            jdx = lookup[j]
            H[idx, jdx] += value

        H[idx, idx] += U[i]

    return H


def do_lanczos(A, U):
    print "do_lanczos: building Hs"
    H = dict(A)
    for i in range(len(U)):
        H[i, i] = H.get((i, i), 0) + U[i]
    #H1 = sparse.lil_matrix((len(U), len(U)))
    keys = H.keys()
    keys.sort()
    data = []
    rows = []
    cols = []
    for idx in keys:
        #H1[idx] = H[idx]
        data.append(H[idx])
        rows.append(idx[0])
        cols.append(idx[1])
    H1 = sparse.coo_matrix((data, (rows, cols)), (len(U), len(U)))
    H1 = sparse.csr_matrix(H1, dtype=numpy.float64)

    print "do_lanczos: eigsh"
    vals, vecs = eigsh(H1, k=min(len(U)-5, 40), which="LM")

    return vals, vecs


def do_orbiham(A, U):
    #print A
    degrees = {} # out-degree
    for key, value in A.items():
        i, j = key
        degrees[i] = degrees.get(i, 0) + value
    degrees = list(set(degrees.values()))
    print "degrees:", degrees
    assert len(degrees) == 1
    degree = degrees[0]

    #for value in A.values():
    #    assert value==1

    H1 = sparse_orbiham_nauty(degree, A, U)
    if H1 is None:
        return

    N = len(H1)
    rows = range(N)
    rows.sort(key = lambda i : -H1[i, i])

    if 0:
        print "orbiham:"
        for i in rows:
            print "%d:"%i,
            for j in rows:
                if H1[i, j]:
                    print "%d:%d"%(j, H1[i, j]),
            print

    if len(H1)<=1024:
        vals, vecs = numpy.linalg.eig(H1)
    else:
        vals, vecs = eigs(H1, k=min(len(H1)-5, 40), which="LM")

    return vals, vecs

class Code(object):
    def __init__(self, filename=None):
        self.filename = filename
        self.lines = []
        self.dent = 0
    def indent(self):
        self.dent += 1
    def dedent(self):
        self.dent -= 1
    def begin(self):
        self.append("{")
        self.indent()
    def end(self):
        self.dedent()
        self.append("}")
    def append(self, line):
        line = '    '*self.dent+line
        self.lines.append(line)
    def output(self):
        if self.filename is None:
            self.lines.append("")
            return '\n'.join(self.lines)
        f = open(self.filename, 'w')
        for line in self.lines:
            print >>f, line
        f.close()


def getnum(v):
    x = 0
    for i in v:
        if i:
            x += 1
        x *= 2
    return x/2

def getnum(v):
    s = ''.join(str(i) for i in v)
    return '0b'+s


def slepc(Gx, Gz, Hx, Hz, Rx, Rz, Pxt, Qx, Pz, Tx, **kw):

    print "slepc"

    r = len(Rx)
    n = 2**r
    assert (r<40), "ugh"

    #code = Code("body.h")
    code = Code()

    code.append("#define DIMS (%d)"%n)

    code.append("static void matmult(PetscScalar *py, const PetscScalar *px, long nx)")
    code.begin()
    code.append("assert(DIMS == %d);"%n)
    code.append("assert(nx == %d);"%n)
    code.append("memset(py, 0, sizeof(PetscScalar)*nx);")

    offset = argv.get("offset", 0)

    mz = len(Gz)
#    print "Hz:"
#    print shortstr(Hz)
#    print "Hx:"
#    print shortstr(Hx)
#    print "Rx:"
#    print shortstr(Rx)
    print "Tx:", len(Tx)
    #print shortstr(Tx)
    t = None
    excite = argv.excite
    if excite is not None:
        print "excite:", excite
        if type(excite) is tuple:
            t = Tx[excite[0]]
            for i in range(1, len(excite)):
                t = (t + Tx[excite[i]])%2
        else:
            t = Tx[excite]
        print "t:", shortstr(t)
        Gzt = dot2(Gz, t)
        print "Gzt:", shortstr(Gzt)

    RR = dot2(Gz, Rx.transpose())

    PxtQx = dot2(Pxt, Qx)
    gxs = [getnum(dot2(gx, PxtQx)) for gx in Gx]
    gxs.sort()
    uniq_gxs = list(set(gxs))
    uniq_gxs.sort()

    code.append("long v;")
    code.append("int k;")
    code.append("struct timeval t0, t1;")
    code.append("gettimeofday(&t0, NULL);")
    code.append("for(v=0; v<%d; v++)"%n)
    code.begin()
    code.append("double pxv = px[v];")
    if n >= 128:
        code.append(r'if((v+1) %% %d == 0)' % (n//128))
        code.begin()
        code.append("gettimeofday(&t1, NULL);")
        code.append("long delta = t1.tv_sec-t0.tv_sec;")
        code.append("if(delta>1)")
        code.append('{printf("[%lds]", delta);fflush(stdout);}')
        code.append('t0 = t1;')
        code.end()
    code.append("k = 0;")
    for i, row in enumerate(RR):
        if t is not None and Gzt[i]==1:
            code.append("k += (countbits_fast(v&%s)+1) %% 2;" % getnum(row))
        else:
            code.append("k += countbits_fast(v&%s) %% 2;" % getnum(row))
    cutoff = argv.cutoff
    if cutoff is not None:
        code.append("if(k>%d) continue; // <-------- continue" % cutoff)
    else:
        code.append("if(k>cutoff) continue; // <-------- continue")
    code.append("py[v] += pxv * (%d - 2*k);" % mz)
    for gx in uniq_gxs:
        s = '+'.join(['pxv']*gxs.count(gx))
        code.append("py[v^%s] += %s;" % (gx, s))
    code.end()
    code.end()

    s = code.output()

    src = open("ex3.c").read()
    match = '\n#include "body.h"\n'
    assert match in src
    src = src.replace(match, s)
    name = argv.get("name")
    assert name and name.endswith(".c")
    f = open(name, 'w')
    f.write(src)
    f.close()

    import socket
    host = socket.gethostname()
    if host == "bucket":
        cmd = "gcc MATCH.c -O3 -o MATCH -I/home/simon/local/petsc/arch-linux2-c-debug/include -I/home/simon/local/petsc/include/petsc/mpiuni -I/home/simon/local/petsc/include -I/home/simon/local/slepc-3.7.1/include -I/home/simon/local/slepc-3.7.1/arch-linux2-c-debug/include/ -L/home/simon/local/petsc/arch-linux2-c-debug/lib -L/home/simon/local/slepc-3.7.1/arch-linux2-c-debug/lib -lpetsc -lslepc"
    else:
        cmd = "gcc -O3 MATCH.c -I/suphys/sburton/include/ -o MATCH -lpetsc -L$PETSC_DIR/$PETSC_ARCH/lib -L$SLEPC_DIR/$PETSC_ARCH/lib -lslepc"
    cmd = cmd.replace("MATCH.c", name)
    cmd = cmd.replace("MATCH", name[:-2])

    print cmd
    rval = os.system(cmd)
    assert rval == 0


def find_errors(Hx, Lx, Rx):
    
    # find Tz
    n = Hx.shape[1]

    k = len(Lx)
    r = len(Rx)
    mx = len(Hx)
    assert k+r+mx <= n
    #mz = n - (k+r+mx)

    U = zeros2(mx+k+r, n)
    U[:mx] = Hx 
    U[mx:mx+k] = Lx 
    U[mx+k:mx+k+r] = Rx 
    B = zeros2(mx+k+r, mx)
    B[:mx] = identity2(mx)

    Tz_t = solve(U, B)
    Tz = Tz_t.transpose()
    assert len(Tz) == mx

    check_conjugate(Hx, Tz)
    check_commute(Lx, Tz)
    check_commute(Rx, Tz)

    return Tz


def main():

    if argv.gcolor:
        size = argv.get("size", 1)
        Gx, Gz, Hx = build_gcolor(size)
        Hz = Hx.copy()

    elif argv.compass:
        l = argv.get('l', 3)
        Gx, Gz, Hx, Hz = build_compass(l)

    elif argv.gcolor2:

        from gcolor import build as build1
        Gx, Gz, Hx = build1()
        Hz = Hx.copy()

    else:

        return

    Lz = find_logops(Gx, Hz)

    check_commute(Lz, Gx)
    check_commute(Lz, Hx)
    
    Px = get_reductor(Hx) # projector onto complement of rowspan of Hx
    Pz = get_reductor(Hz) 

    Rz = [dot2(Pz, g) for g in Gz]
    Rz = array2(Rz)
    Rz = row_reduce(Rz, truncate=True)
    rz = len(Rz)
    print "Rz:", rz

    Tx = find_errors(Hz, Lz, Rz)

    n = Gx.shape[1]
    print "Hx:", len(Hx)
    print "Gx:", len(Gx)
    print "Gz:", len(Gz)

    Rx = [dot2(Px, g) for g in Gx]
    Rx = array2(Rx)

    Rx = row_reduce(Rx, truncate=True)
    rx = len(Rx)
    print "Rx:", rx

    Qx = u_inverse(Rx)
    Pxt = Px.transpose()
    assert eq2(dot2(Rx, Qx), identity2(rx))
    assert eq2(dot2(Rx, Pxt), Rx)

    #print shortstr(dot2(Pxt, Qx))
    PxtQx = dot2(Pxt, Qx)
    lines = [shortstr(dot2(g, PxtQx)) for g in Gx]
    lines.sort()
    #print "PxtQx:"
    #for s in lines:
    #    print s
    #print "RzRxt"
    #print shortstr(dot2(Rz, Rx.transpose()))

    if argv.slepc:
        slepc(**locals())
        return

    v0 = None

#    excite = argv.excite
#    if excite is not None:
#        v0 = zeros2(n)
#        v0[excite] = 1

    verts = []
    lookup = {}
    for i, v in enumerate(span(Rx)): # XXX does not scale well
        if v0 is not None:
            v = (v+v0)%2
            v = dot2(Px, v)
        lookup[v.tostring()] = i
        verts.append(v)
    print "span:", len(verts)
    assert len(lookup) == len(verts)

    mz = len(Gz)
    n = len(verts)

    if n <= 1024 and argv.solve:
        H = numpy.zeros((n, n))
        for i, v in enumerate(verts):
            count = dot2(Gz, v).sum()
            Pxv = dot2(Px, v)
            assert count == dot2(Gz, Pxv).sum()
            H[i, i] = mz - 2*count
            for g in Gx:
                v1 = (g+v)%2
                v1 = dot2(Px, v1)
                j = lookup[v1.tostring()]
                H[i, j] += 1
    
        if argv.showham:
            s = lstr2(H, 0).replace(',  ', ' ')
            s = s.replace(' 0', ' .')
            print s
    
        vals, vecs = numpy.linalg.eigh(H)
        show_eigs(vals)

        if argv.orbiham:
            H1 = build_orbiham(H)
            print "orbiham:"
            print H1
            vals, vecs = numpy.linalg.eig(H1)
            show_eigs(vals)

    elif argv.sparse:
        print "building H",
        A = {} # adjacency
        U = [] # potential

        offset = mz + 1 # make H positive definite

        for i, v in enumerate(verts):
            if i%1000==0:
                write('.')
            count = dot2(Gz, v).sum()
            #H[i, i] = mz - 2*count
            U.append(offset + mz - 2*count)
            for g in Gx:
                v1 = (g+v)%2
                v1 = dot2(Px, v1)
                j = lookup[v1.tostring()]
                A[i, j] = A.get((i, j), 0) + 1
    
        print "\nnnz:", len(A)

        if argv.lanczos:
            vals, vecs = do_lanczos(A, U)

        elif argv.orbiham:
            vals, vecs = do_orbiham(A, U)

        else:
            return

        vals -= offset # offset doesn't change vecs

        show_eigs(vals)


from argv import Argv
argv = Argv()

if __name__ == "__main__":

    main()





