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



def gen_code(l):

    n = l*l

    keys = [(i, j) for i in range(l) for j in range(l)]
    coords = {} 
    for i, j in keys:
        for di in range(-l, l+1):
          for dj in range(-l, l+1):
            coords[i+di, j+dj] = keys.index(((i+di)%l, (j+dj)%l))

    #print keys

    output = open("c_orbits.c", "w")

    print >>output, """

#include "Python.h"

static void
list_append(PyObject *items, char *s)
{
    PyObject *string = PyString_FromString(s);
    PyList_Append(items, string);
    Py_DECREF(string);
}

static PyObject *
get_orbit(PyObject *self, PyObject *args)
{
    const char *src;
    if (!PyArg_ParseTuple(args, "s", &src))
        return NULL;

    assert(strlen(src)==%(n)s);
    //print "orbit", src
    PyObject *items = PyList_New(0);
    assert(items);

    char target[%(n)s+1];
    char target1[%(n)s+1];
    char tmp;
    target[%(n)s] = 0;
    target1[%(n)s] = 0;
    """%{'n':n}

    for di in range(l):
      for dj in range(l):
        print >>output, "    // di, dj = (%d, %d)"%(di, dj)
        for coord in range(n):
            i, j = keys[coord]
            src, tgt = coords[i, j], coords[i+di, j+dj]
            print >>output, "    target[%d] = src[%d];"%(tgt, src)

        for idx in genidx((2,)*l):
          if sum(idx)%2==0:
            print >>output, "    // idx = ", idx
            for tgt in range(n):
                i, j = keys[tgt]
                if idx[j]==0:
                    print >>output, "    target1[%d] = target[%d];"%(tgt, tgt)
                else:
                    print >>output, "    target1[%d] = '0'+'1'-target[%d];"%(tgt, tgt)
            print >>output, "    list_append(items, (target1));"

            # reflection
            for i in range(l//2):
              for j in range(l):
                print >>output, "    tmp = target1[%d];" % (coords[i, j])
                print >>output, "    target1[%d] = target1[%d];" % (coords[i, j], coords[l-i-1, j])
                print >>output, "    target1[%d] = tmp;" % (coords[l-i-1, j])
            print >>output, "    list_append(items, (target1));"

            # reflection
            for i in range(l):
              for j in range(l//2):
                print >>output, "    tmp = target1[%d];" % (coords[i, j])
                print >>output, "    target1[%d] = target1[%d];" % (coords[i, j], coords[i, l-j-1])
                print >>output, "    target1[%d] = tmp;" % (coords[i, l-j-1])
            print >>output, "    list_append(items, (target1));"

            # reflection
            for i in range(l//2):
              for j in range(l):
                print >>output, "    tmp = target1[%d];" % (coords[i, j])
                print >>output, "    target1[%d] = target1[%d];" % (coords[i, j], coords[l-i-1, j])
                print >>output, "    target1[%d] = tmp;" % (coords[l-i-1, j])
            print >>output, "    list_append(items, (target1));"

    print >>output, "    return items;"
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

    for i in range(l):
      for j in range(l):
        print >>output, "    //", (i, j)
        for idx in range(n):
            if coords[i, j] == idx or coords[i, j+1] == idx:
                print >>output, "    target[%d] = '0'+'1'-src[%d];" % (idx, idx)
            else:
                print >>output, "    target[%d] = src[%d];" % (idx, idx)
        print >>output, "    list_append(items, (target));"
        
    print >>output, "    return items;"
    print >>output, "}"

    print >>output, """

static PyMethodDef OrbitsMethods[] = 
{
    {"get_orbit",  get_orbit, METH_VARARGS, "get the orbit"},
    {"get_gauge",  get_gauge, METH_VARARGS, "get the gauge"},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

PyMODINIT_FUNC
initc_orbits(void)
{
    (void) Py_InitModule("c_orbits", OrbitsMethods);
}

    """
    output.flush()
    output.close()


def main():

    l = int(sys.argv[1])
    assert 3<=l<=6
    n = l**2

    if "compile" in sys.argv:
        gen_code(l)
        try:
            os.unlink("c_orbits.so")
        except:
            pass
        rval = os.system("python setup.py build_ext --inplace")
        assert rval == 0

    import c_orbits


    if l==3 and 0:
        s = "100010000"
        #s = "1000100000000000"
        #s = "1000100000000000000000000"
        items = c_orbits.get_orbit(s)
    
        print items
        print len(items), "uniq:", len(set(items))
    
        s = "0"*n
        items = c_orbits.get_gauge(s)
        print items
        print len(items), "uniq:", len(set(items))

    s = "0"*n
    marked = {s:1}
    bdy = [s]
    print "marked:", len(marked)

    while bdy:
        _bdy = []
        for s0 in bdy:
            #print "bdy:", s0
            for s1 in c_orbits.get_gauge(s0):
                #print "s1:", s1
                # is this a new state ?
                if s1 in marked:
                    #print "*"
                    continue
                for s2 in c_orbits.get_orbit(s1):
                    #print "    orbit(s1):", s2
                    if s2 in marked:
                        #print "+"
                        break
                else:
                    #print "add"
                    marked[s2] = 1
                    _bdy.append(s2)
                    #sys.stdout.write('.');sys.stdout.flush()
                    #if len(_bdy)>10000:
                        #return
        # new boundary
        bdy = _bdy
        print "orbits:", len(marked), "bdy:", len(bdy)

    #print "orbits:", marked
    print "orbits:", len(marked)

    for s0 in marked:
        obt = c_orbits.get_orbit(s0)
        for s1 in obt:
            if s1!=s0 and s1 in marked:
                print s0, s1, "both in marked"
        obt.sort(key = lambda s : s.count('1'))
        #print obt[0]



if __name__ == "__main__":

    main()





