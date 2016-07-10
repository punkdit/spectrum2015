#!/usr/bin/env python

import sys, os

#import numpy



def genidx(shape):
    if len(shape)==0:
        yield ()
    else:
        for idx in range(shape[0]):
            for _idx in genidx(shape[1:]):
                yield (idx,)+_idx


def gen_code(n, gcolor_gauge, gcolor_stab):

    output = open("c_gorbits.c", "w")

    print >>output, r"""

#include "Python.h"

#define SWAP(a, b)  {tmp=a; a=b; b=tmp;}

typedef struct _bv
{char bits[%(n)s+1];}
Bitvec;

static Bitvec s_uniq;

#define BV_CLR(bv)  {(bv).bits[0] = 0;}
#define BV_SET(tgt, src)  {memcpy((tgt).bits, (src).bits, sizeof(tgt));}
#define BV_CMP(tgt, src)  strncmp((tgt).bits, (src).bits, %(n)s+1)
#define BV_FLIP(tgt, i) {(tgt).bits[(i)] = '0'+'1'-(tgt).bits[(i)];}
#define BV_SETBIT(tgt, src, i) {(tgt).bits[(i)] = (src).bits[(i)];}
#define BV_SETFLIP(tgt, src, i) {(tgt).bits[(i)] = '0'+'1'-(src).bits[(i)];}




static void
list_append(PyObject *items, Bitvec *s)
{
    PyObject *string = PyString_FromString(s->bits);
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
    strncpy(src.bits, s_arg, %(n)s+1);
    assert(strlen(src.bits)==%(n)s);

    Bitvec target;
    //target.bits[%(n)s] = 0;

    BV_CLR(s_uniq);
    """%{'n':n}

    stabs = gcolor_stab.replace('.', '0').strip().split()
    m = len(stabs)
    ops = gcolor_gauge.replace('.', '0').strip().split()
    ng = len(ops)

    for i in range(m):
        print >>output, "    int i_%d;"%i
    print >>output, "    BV_SET(target, src);"
    for i in range(m):
        print >>output, "    for(i_%d=0; i_%d<2; i_%d++)" % (i, i, i)
        print >>output, "    {"

    print >>output, "    if(BV_CMP(s_uniq, target)<0) BV_SET(s_uniq, target);"

    for i in range(m):
        stab = stabs[i]
        for j in range(n):
            if stab[j]=='1':
                print >>output, "    BV_FLIP(target, %d);"%j
        print >>output, "    }"

    print >>output, "    assert(strncmp(target.bits, src.bits, %d+1)==0);"%n
    print >>output, "    return PyString_FromString(s_uniq.bits);"
    print >>output, "}"


    print >>output, """

static Bitvec all_targets[%(ng)s];

void
fill_gauge(Bitvec src)
{
    Bitvec target;

    """%locals()

    for j, op in enumerate(ops):
        print >>output, "    // ", op
        #print >>output, "    strncpy(target, src, %d);"%(n+1)
        for i in range(n):
            flip = op[i]=='1'
            if flip:
                print >>output, "    BV_SETFLIP(all_targets[%d], src, %d)"%(j, i)
            else:
                print >>output, "    BV_SETBIT(all_targets[%d], src, %d)"%(j, i)

#        print >>output, "    list_append(items, &target);"
#        
#    print >>output, "    return items;"
    print >>output, "}"

    print >>output, """

static PyObject *
get_gauge(PyObject *self, PyObject *args)
{
    Bitvec src;
    const char *s_arg;
    if (!PyArg_ParseTuple(args, "s", &s_arg))
        return NULL;

    assert(strlen(s_arg)==%(n)s);
    strncpy(src.bits, s_arg, %(n)s+1);

    fill_gauge(src);

    PyObject *items = PyList_New(0);
    int i;
    for(i=0; i<%(ng)s; i++)
        list_append(items, &all_targets[i]);
    return items;
}


static PyMethodDef OrbitsMethods[] = 
{
    //{"get_orbit",  get_orbit, METH_VARARGS, "get the orbit"},
    {"get_gauge",  get_gauge, METH_VARARGS, "get the gauge"},
    {"get_uniq",  get_uniq, METH_VARARGS, "get the uniq"},
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
    cmd = "i686-linux-gnu-gcc -O3 -pthread -fno-strict-aliasing -DNDEBUG -g -fwrapv -Wall -Wstrict-prototypes -fPIC -I/usr/include/python2.7 -I/suphys/sburton/include/python2.7 -c c_gorbits.c -o c_gorbits.o"
    print cmd
    rval = os.system(cmd)
    assert rval==0
    cmd = "i686-linux-gnu-gcc -O3 -pthread -shared -Wl,-Bsymbolic-functions -Wl,-Bsymbolic-functions -Wl,-z,relro -fno-strict-aliasing -DNDEBUG -g -fwrapv -Wall -Wstrict-prototypes -D_FORTIFY_SOURCE=2 -g -fstack-protector --param=ssp-buffer-size=4 -Wformat -Werror=format-security c_gorbits.o -o c_gorbits.so"
    print cmd
    rval = os.system(cmd)
    assert rval==0




def main():

    from qupy.ldpc import gcolor
    size = argv.get("size", 1)
    lattice = gcolor.Lattice(size)

    n = len(lattice.qubits)
    print lattice

    code = lattice.build_code(check=False)
    #Ex = lattice.Ex
    Gx, Gz = code.Gx, code.Gz
    Hx, Hz = code.Hx, code.Hz

    Gx = gcolor.shortstr(Gx)
    Hx = gcolor.shortstr(Hx)

    if "compile" in sys.argv:
        gen_code(n, Gx, Hx)
        try:
            os.unlink("c_gorbits.so")
        except:
            pass
        build()

    import c_gorbits as cg

    s = "0"*n
    s = cg.get_uniq(s)
    orbits = {s:1}
    bdy = [s]
    print "orbits:", len(orbits)

    while bdy:
        _bdy = []
        for s0 in bdy:
            #print "bdy:", s0
            for s1 in cg.get_gauge(s0):
                #print "s1:", s1
                s2 = cg.get_uniq(s1)
                if s2 not in orbits:
                    orbits[s2] = 1
                    _bdy.append(s2)
        # new boundary
        bdy = _bdy
        print "orbits:", len(orbits), "bdy:", len(bdy)

        if len(orbits)>5000:
            break

    print "marked:", len(marked)


from argv import Argv
argv = Argv()

if __name__ == "__main__":

    main()





