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


gcolor_gauge = """ 
1111...........
11..11.........
1.1.1.1........
..11..11.......
.1.1.1.1.......
....1111.......
11......11.....
1.1.....1.1....
........1111...
..11......11...
.1.1.....1.1...
1...1...1...1..
........11..11.
.1...1...1...1.
....11......11.
........1.1.1.1
..1...1...1...1
....1.1.....1.1
"""

gcolor_stab = """ 
11111111.......
1111....1111...
11..11..11..11.
1.1.1.1.1.1.1.1
"""

gcolor_logop = "........1111111"


def gen_code(n, gcolor_gauge, gcolor_stab):

    output = open("c_gorbits.c", "w")

    print >>output, """

#include "Python.h"

#define SWAP(a, b)  {tmp=a; a=b; b=tmp;}

static char s_uniq[%(n)s+1];

static void
clear_uniq(void)
{
    s_uniq[0] = 0;
}

static void
set_uniq(char *s)
{
    if(strncmp(s_uniq, s, %(n)s) < 0)
        strncpy(s_uniq, s, %(n)s+1);
}

static void
list_append(PyObject *items, char *s)
{
    PyObject *string = PyString_FromString(s);
    PyList_Append(items, string);
    Py_DECREF(string);
}

"""%{'n':n}


    print >>output, """
static PyObject *
get_uniq(PyObject *self, PyObject *args)
{
    const char *src;
    if (!PyArg_ParseTuple(args, "s", &src))
        return NULL;

    assert(strlen(src)==%(n)s);
    //print "orbit", src

    char target[%(n)s+1];
    //char target1[%(n)s+1];
    //char tmp;

    target[%(n)s] = 0;
    //target1[%(n)s] = 0;

    clear_uniq();
    """%{'n':n}

    stabs = gcolor_stab.replace('.', '0').strip().split()
    m = len(stabs)
    for mask in genidx((2,)*m):
        print >>output, "    // ", mask
        #print >>output, "    strncpy(target, src, %d);"%(n+1)
        for i in range(n):
            flip = 0
            for j, k in enumerate(mask):
                if k and stabs[j][i]=='1':
                    flip = 1-flip
            if flip:
                print >>output, "    target[%d] = '0'+'1'-src[%d];"%(i, i)
            else:
                print >>output, "    target[%d] = src[%d];"%(i, i)
        print >>output, "    set_uniq(target);"
    print >>output, "    return PyString_FromString(s_uniq);"
    print >>output, "}"


    print >>output, """

static PyObject *
get_gauge(PyObject *self, PyObject *args)
{
    const char *src;
    if (!PyArg_ParseTuple(args, "s", &src))
        return NULL;

    assert(strlen(src)==%(n)s);
    //print "orbit", src
    PyObject *items = PyList_New(0);

    char target[%(n)s+1];
    target[%(n)s] = 0;
    """%{'n':n}

    ops = gcolor_gauge.replace('.', '0').strip().split()
    for op in ops:
        print >>output, "    // ", mask
        #print >>output, "    strncpy(target, src, %d);"%(n+1)
        for i in range(n):
            flip = op[i]=='1'
            if flip:
                print >>output, "    target[%d] = '0'+'1'-src[%d];"%(i, i)
            else:
                print >>output, "    target[%d] = src[%d];"%(i, i)

        print >>output, "    list_append(items, target);"
        
    print >>output, "    return items;"
    print >>output, "}"

    print >>output, """

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

    """
    output.flush()
    output.close()


def build():
    cmd = "i686-linux-gnu-gcc -pthread -fno-strict-aliasing -DNDEBUG -g -fwrapv -Wall -Wstrict-prototypes -fPIC -I/usr/include/python2.7 -I/suphys/sburton/include/python2.7 -c c_gorbits.c -o c_gorbits.o"

    print cmd
    rval = os.system(cmd)
    assert rval==0
    cmd = "i686-linux-gnu-gcc -pthread -shared -Wl,-Bsymbolic-functions -Wl,-Bsymbolic-functions -Wl,-z,relro -fno-strict-aliasing -DNDEBUG -g -fwrapv -Wall -Wstrict-prototypes -D_FORTIFY_SOURCE=2 -g -fstack-protector --param=ssp-buffer-size=4 -Wformat -Werror=format-security c_gorbits.o -o c_gorbits.so"
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
        #rval = os.system("python setup.py build_ext --inplace")
        #assert rval == 0

#        save = sys.argv[1:]
#        sys.argv[1:] = "build_ext --inplace".split()
#        from distutils.core import setup, Extension
#        #from Cython.Build import cythonize
#        
#        setup(
#          name = 'Orbit compute',
#          #ext_modules = cythonize("orbits_wrap.pyx"),
#          ext_modules=[Extension('c_gorbits', ['c_gorbits.c'])],
#          #extra_compile_args = ["-O0"], 
#        )
#        sys.argv[1:] = save

        build()

    import c_gorbits

    s = "0"*n
    s = c_gorbits.get_uniq(s)
    marked = {s:1}
    bdy = [s]
    print "marked:", len(marked)

    while bdy:
        _bdy = []
        for s0 in bdy:
            #print "bdy:", s0
            for s1 in c_gorbits.get_gauge(s0):
                #print "s1:", s1
                s2 = c_gorbits.get_uniq(s1)
                if s2 not in marked:
                    marked[s2] = 1
                    _bdy.append(s2)
        # new boundary
        bdy = _bdy
        print "marked:", len(marked), "bdy:", len(bdy)

    print "marked:", len(marked)

    return

    for s0 in marked:
        obt = c_gorbits.get_orbit(s0)
        for s1 in obt:
            if s1!=s0 and s1 in marked:
                print s0, s1, "both in marked"
        obt.sort(key = lambda s : s.count('1'))
        #print obt[0]


from argv import Argv
argv = Argv()

if __name__ == "__main__":

    main()





