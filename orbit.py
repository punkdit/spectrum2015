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

    output = open("c_orbits.pyx", "w")

    print >>output, """
def orbit(char *src):
    assert len(src)==%(n)s
    #print "orbit", src
    items = []

    cdef char target[%(n)s+1]
    cdef char target1[%(n)s+1]
    cdef char tmp
    target[%(n)s] = 0
    target1[%(n)s] = 0
    """%{'n':n}

    for di in range(l):
      for dj in range(l):
        print >>output, "    # di, dj = (%d, %d)"%(di, dj)
        for coord in range(n):
            i, j = keys[coord]
            src, tgt = coords[i, j], coords[i+di, j+dj]
            print >>output, "    target[%d] = src[%d]"%(tgt, src)

        for idx in genidx((2,)*l):
          if sum(idx)%2==0:
            print >>output, "    # idx = ", idx
            for tgt in range(n):
                i, j = keys[tgt]
                if idx[j]==0:
                    print >>output, "    target1[%d] = target[%d]"%(tgt, tgt)
                else:
                    print >>output, "    target1[%d] = ord('0')+ord('1')-target[%d]"%(tgt, tgt)
            print >>output, "    items.append(target1)"

            # reflection
            for i in range(l//2):
              for j in range(l):
                print >>output, "    tmp = target1[%d]" % (coords[i, j])
                print >>output, "    target1[%d] = target1[%d]" % (coords[i, j], coords[l-i-1, j])
                print >>output, "    target1[%d] = tmp" % (coords[l-i-1, j])
            print >>output, "    items.append(target1)"

            # reflection
            for i in range(l):
              for j in range(l//2):
                print >>output, "    tmp = target1[%d]" % (coords[i, j])
                print >>output, "    target1[%d] = target1[%d]" % (coords[i, j], coords[i, l-j-1])
                print >>output, "    target1[%d] = tmp" % (coords[i, l-j-1])
            print >>output, "    items.append(target1)"

            # reflection
            for i in range(l//2):
              for j in range(l):
                print >>output, "    tmp = target1[%d]" % (coords[i, j])
                print >>output, "    target1[%d] = target1[%d]" % (coords[i, j], coords[l-i-1, j])
                print >>output, "    target1[%d] = tmp" % (coords[l-i-1, j])
            print >>output, "    items.append(target1)"

    print >>output, "    return items"


    print >>output, """
def gauge(char *src):
    assert len(src)==%(n)s
    #print "orbit", src
    items = []

    cdef char target[%(n)s+1]
    target[%(n)s] = 0
    """%{'n':n}

    for i in range(l):
      for j in range(l):
        print >>output, "    #", (i, j)
        for idx in range(n):
            if coords[i, j] == idx or coords[i, j+1] == idx:
                print >>output, "    target[%d] = ord('0')+ord('1')-src[%d]" % (idx, idx)
            else:
                print >>output, "    target[%d] = src[%d]" % (idx, idx)
        print >>output, "    items.append(target)"
        
    print >>output, "    return items"
    output.flush()
    output.close()


def main():

    l = int(sys.argv[1])
    assert 3<=l<=6
    n = l**2

    if "compile" in sys.argv:
        gen_code(l)
        os.unlink("c_orbits.so")
        rval = os.system("python setup.py build_ext --inplace")
        assert rval == 0

    import c_orbits


    if l==3 and 0:
        s = "100010000"
        #s = "1000100000000000"
        #s = "1000100000000000000000000"
        items = c_orbits.orbit(s)
    
        print items
        print len(items), "uniq:", len(set(items))
    
        s = "0"*n
        items = c_orbits.gauge(s)
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
            for s1 in c_orbits.gauge(s0):
                #print "s1:", s1
                # is this a new state ?
                if s1 in marked:
                    #print "*"
                    continue
                for s2 in c_orbits.orbit(s1):
                    #print "    orbit(s1):", s2
                    if s2 in marked:
                        #print "+"
                        break
                else:
                    #print "add"
                    marked[s2] = 1
                    _bdy.append(s2)
        # new boundary
        bdy = _bdy
        print "orbits:", len(marked), "bdy:", len(bdy)

    #print "orbits:", marked
    print "orbits:", len(marked)

    for s0 in marked:
        obt = c_orbits.orbit(s0)
        for s1 in obt:
            if s1!=s0 and s1 in marked:
                print s0, s1, "both in marked"
        obt.sort(key = lambda s : s.count('1'))
        #print obt[0]



if __name__ == "__main__":

    main()





