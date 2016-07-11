#!/usr/bin/env python

import sys, os

import numpy

from solve import shortstr, parse, dot2, zeros2
from lanczos import write, show_eigs


def genidx(shape):
    if len(shape)==0:
        yield ()
    else:
        for idx in range(shape[0]):
            for _idx in genidx(shape[1:]):
                yield (idx,)+_idx


def gen_code(n, Gx, Hx, perms):
    Gx = shortstr(Gx)
    Hx = shortstr(Hx)

    output = open("c_gorbits.c", "w")

    print >>output, r"""
#include "Python.h"
    """

    if 1:
        # char bitvectors
        print >>output, r"""

typedef struct _bv
{char bits[%(n)s+1];}
Bitvec;

static Bitvec s_uniq;

#define BV_ASSERT(bv) { int i; for(i=0;i<%(n)s;i++) assert((bv).bits[i]=='0' || (bv).bits[i]=='1'); assert((bv).bits[i]==0); }

//#define BV_ASSERT(bv) 

//#define BV_CLR(bv)  {(bv).bits[0] = 0;}
#define BV_CLR(bv)  {memset(&bv, 0, sizeof(Bitvec));}
#define BV_SET(tgt, src)  {memcpy((tgt).bits, (src).bits, sizeof(tgt));}
#define BV_CMP(tgt, src)  strncmp((tgt).bits, (src).bits, %(n)s+1)
#define BV_FLIP(tgt, i) {(tgt).bits[(i)] = '0'+'1'-(tgt).bits[(i)];}
#define BV_SETBIT(tgt, src, i) {(tgt).bits[(i)] = (src).bits[(i)];}
#define BV_SETMAP(tgt, src, i, j) {(tgt).bits[(i)] = (src).bits[(j)];}
#define BV_SETFLIP(tgt, src, i) {(tgt).bits[(i)] = '0'+'1'-(src).bits[(i)];}


void
BV_TOSTRING(char *str, Bitvec bv)
{
    BV_ASSERT(bv);
    strncpy(str, bv.bits, %(n)s+1);
}

void
BV_FROMSTRING(Bitvec *bv, const char *str)
{
    strncpy(bv->bits, str, %(n)s+1);
    BV_ASSERT(*bv);
}

"""%locals()

    else:
        # int bitvectors
        print >>output, r"""

typedef struct _bv
{int bits[%(n)s];}
Bitvec;

static Bitvec s_uniq;

#define BV_ASSERT(bv) { int i; for(i=0;i<%(n)s;i++) assert((bv).bits[i]==0 || (bv).bits[i]==1); }

#define BV_CLR(bv)  {memset(&bv, 0, sizeof(Bitvec));}
#define BV_SET(tgt, src)  {memcpy(&tgt, &src, sizeof(Bitvec));}
#define BV_CMP(tgt, src)  memcmp(&tgt, &src, sizeof(Bitvec))
#define BV_FLIP(tgt, i) {(tgt).bits[(i)] = 1-(tgt).bits[(i)];}
#define BV_SETBIT(tgt, src, i) {(tgt).bits[(i)] = (src).bits[(i)];}
#define BV_SETMAP(tgt, src, i, j) {(tgt).bits[(i)] = (src).bits[(j)];}
#define BV_SETFLIP(tgt, src, i) {(tgt).bits[(i)] = 1-(src).bits[(i)];}


void
BV_TOSTRING(char *str, Bitvec bv)
{
    BV_ASSERT(bv);
    int i;
    for(i=0;i<%(n)s;i++)
        str[i] = '0' + bv.bits[i];
    str[i] = 0;
}

void
BV_FROMSTRING(Bitvec *bv, const char *str)
{
    int i;
    for(i=0;i<%(n)s;i++)
        bv->bits[i] = str[i]-'0';
    BV_ASSERT(*bv);
}

"""%locals()


    print >>output, r"""
static void
list_append(PyObject *items, Bitvec bv)
{
    char str[%(n)s+1];
    BV_ASSERT(bv);
    BV_TOSTRING(str, bv);
    PyObject *string = PyString_FromString(str);
    PyList_Append(items, string);
    Py_DECREF(string);
}

"""%{'n':n}


    print >>output, r"""
static PyObject *
get_uniq(PyObject *self, PyObject *args)
{
    Bitvec src;
    const char *s_arg;
    if (!PyArg_ParseTuple(args, "s", &s_arg))
        return NULL;

    //printf("get_uniq %%s\n", s_arg);
    BV_FROMSTRING(&src, s_arg);
    BV_ASSERT(src);

    Bitvec target;
    Bitvec target1;

    BV_CLR(s_uniq);
    BV_CLR(target1);
    BV_SET(s_uniq, src);
    """%{'n':n}

    stabs = Hx.replace('.', '0').strip().split()
    m = len(stabs)
    ops = Gx.replace('.', '0').strip().split()
    ng = len(ops)

    for i in range(m):
        print >>output, "    int i_%d;"%i
    print >>output, "    BV_ASSERT(src);"
    print >>output, "    BV_SET(target, src);"
    print >>output, "    BV_ASSERT(src);"
    print >>output, "    BV_ASSERT(target);"
    for i in range(m):
        print >>output, "    for(i_%d=0; i_%d<2; i_%d++)" % (i, i, i)
        print >>output, "    {"

    print >>output, "    BV_ASSERT(target);"
    for perm in perms:
        for i, j in enumerate(perm):
            print >>output, "    BV_SETMAP(target1, target, %d, %d);"%(i, j)
        print >>output, "    BV_ASSERT(target1);"
        print >>output, "    BV_ASSERT(s_uniq);"
        print >>output, "    if(BV_CMP(s_uniq, target1)>0) BV_SET(s_uniq, target1);"
        print >>output, "    BV_ASSERT(s_uniq);"

    for i in range(m):
        stab = stabs[i]
        for j in range(n):
            if stab[j]=='1':
                print >>output, "    BV_FLIP(target, %d);"%j
        print >>output, "    }"

    print >>output, "    assert(BV_CMP(target, src)==0);"
    print >>output, "    char str[%s+1];"%n
    print >>output, "    BV_ASSERT(s_uniq);"
    print >>output, "    BV_TOSTRING(str, s_uniq);"
    print >>output, "    return PyString_FromString(str);"
    print >>output, "}"


    print >>output, """

static Bitvec all_targets[%(ng)s];

void
fill_gauge(Bitvec src)
{
    """%locals()

    for j, op in enumerate(ops):
        print >>output, "    // ", op
        #print >>output, "    strncpy(target, src, %d);"%(n+1)
        print >>output, "    BV_SET(all_targets[%d], src);"%j
        
        for i in range(n):
            flip = op[i]=='1'
            if flip:
                print >>output, "    BV_FLIP(all_targets[%d], %d)"%(j, i)
            #else:
            #    print >>output, "    BV_SETBIT(all_targets[%d], src, %d)"%(j, i)

#        print >>output, "    list_append(items, &target);"
#        
#    print >>output, "    return items;"
    print >>output, "}"

    print >>output, r"""

static PyObject *
get_gauge(PyObject *self, PyObject *args)
{
    Bitvec src;
    const char *s_arg;
    if (!PyArg_ParseTuple(args, "s", &s_arg))
        return NULL;

    //printf("get_gauge: s_arg = %%s\n", s_arg);
    assert(strlen(s_arg)==%(n)s);
    BV_FROMSTRING(&src, s_arg);

    fill_gauge(src);

    PyObject *items = PyList_New(0);
    int i;
    for(i=0; i<%(ng)s; i++)
        list_append(items, all_targets[i]);
    return items;
}

#define NGAUGE (%(ng)s)
#define NBITS (%(n)s)
#include "c_gorbits.h"

static PyMethodDef OrbitsMethods[] = 
{
    //{"get_orbit",  get_orbit, METH_VARARGS, "get the orbit"},
    {"get_gauge",  get_gauge, METH_VARARGS, "get the gauge"},
    {"get_uniq",  get_uniq, METH_VARARGS, "get the uniq"},
    {"all_orbits",  all_orbits, METH_VARARGS, "get all orbits"},
    //{"clear_uniq",  clear_uniq, METH_VARARGS, "clear the uniq"},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

PyMODINIT_FUNC
initc_gorbits(void)
{
    (void) Py_InitModule("c_gorbits", OrbitsMethods);
}

    """%locals()

    output.flush()
    output.close()


def build():
    cmd = "i686-linux-gnu-gcc -O1 -pthread -fno-strict-aliasing -g -fwrapv -Wall -Wstrict-prototypes -fPIC -I/usr/include/python2.7 -I/suphys/sburton/include/python2.7 -c c_gorbits.c -o c_gorbits.o"
    print cmd
    rval = os.system(cmd)
    assert rval==0
    cmd = "i686-linux-gnu-gcc -pthread -shared -Wl,-Bsymbolic-functions -Wl,-Bsymbolic-functions -Wl,-z,relro -fno-strict-aliasing -g -fwrapv -Wall -Wstrict-prototypes -D_FORTIFY_SOURCE=2 -g -fstack-protector --param=ssp-buffer-size=4 -Wformat -Werror=format-security c_gorbits.o -o c_gorbits.so"
    print cmd
    rval = os.system(cmd)
    assert rval==0



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

    return Gx, Gz, Hx


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


def orbiham(H):
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

    print shortstr(P)
    print shortstr(Q)

    H = numpy.dot(Q, numpy.dot(H, P))
    return H




def main():

    if argv.gcolor:
        size = argv.get("size", 1)
        Gx, Gz, Hx = build_gcolor(size)
    elif argv.compass:
        l = argv.get('l', 3)
        Gx, Gz, Hx = build_compass(l)

    else:

        return

    n = Gx.shape[1]

    perms = [tuple(range(n))]
    if argv.isomorph:
        perms = build_isomorph(Gx)
    print perms[:10]

    if argv.compile:
        gen_code(n, Gx, Hx, perms)
        try:
            os.unlink("c_gorbits.so")
        except:
            pass
        build()

        if not argv.run:
            return

    import c_gorbits as cg

    s0 = "0"*n
    s0 = cg.get_uniq(s0)
    assert s0.count('0')+s0.count('1') == len(s0), repr(s0)
    orbits = {s0:1}
    bdy = [s0]
    print "orbits:", len(orbits)

    while bdy:
        _bdy = []
        for s0 in bdy:
            assert s0.count('0')+s0.count('1') == len(s0), repr(s0)
            for s1 in cg.get_gauge(s0):
                assert s1.count('0')+s1.count('1') == len(s1), repr(s1)
                #print "s1:", s1
                s2 = cg.get_uniq(s1)
                assert s2.count('0')+s2.count('1') == len(s2), repr(s2)
                if s2 not in orbits:
                    orbits[s2] = 1
                    _bdy.append(s2)
        # new boundary
        bdy = _bdy
        print "orbits:", len(orbits), "bdy:", len(bdy)

        if len(orbits)>5000:
            break

    print "orbits:", len(orbits)
    orbits = list(orbits.keys())
    orbits.sort()

    mz = len(Gz)

    n = len(orbits)
    H = numpy.zeros((n, n))
    for i, s0 in enumerate(orbits):
        v = parse(s0).transpose()
        count = dot2(Gz, v).sum()
        H[i, i] = mz - 2*count
        row = {}
        for s1 in cg.get_gauge(s0):
            s1 = cg.get_uniq(s1)
            j = orbits.index(s1)
            H[i, j] += 1

    print H

    vals, vecs = numpy.linalg.eig(H)
    show_eigs(vals)

    if argv.orbiham:
        H1 = orbiham(H)
        print "orbiham:"
        print H1
        vals, vecs = numpy.linalg.eig(H1)
        show_eigs(vals)


from argv import Argv
argv = Argv()

if __name__ == "__main__":

    main()





