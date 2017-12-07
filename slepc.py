#!/usr/bin/env python

import sys, os

from argv import argv
from solve import shortstr, shortstrx, parse, eq2, dot2, zeros2, array2, identity2


# copied from zorbit.py:

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
    s = ''.join(str(i) for i in v)
    return '0b'+s



def slepc(Gx, Gz, Hx, Hz, Rx, Rz, Pxt, Qx, Pz, Tx, excite=None, **kw):

    name = argv.get("name", "ex3.tmp.c")
    print "slepc: name=%s"%name

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

    offset = argv.get("offset", None)

    mz = len(Gz)
    t = None

    if excite is not None:
        print "excite:", excite
        if type(excite) in (int, long):
            t = Tx[excite]
        else:
            t = dot2(excite, Tx)
        print "t:", shortstr(t)
        Gzt = dot2(Gz, t)
        print "Gzt:", shortstr(Gzt)

    weights = kw.get("weights")
    if weights is not None:
        assert len(weights)==len(Gx)

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

    if weights is None:
        for gx in uniq_gxs:
            s = '+'.join(['pxv']*gxs.count(gx))
            code.append("py[v^%s] += %s;" % (gx, s))
    else:
        gxs = [getnum(dot2(gx, PxtQx)) for gx in Gx]
        for i, gx in enumerate(gxs):
            code.append("py[v^%s] += %s*pxv;" % (gx, weights[i]))

    code.end()
    code.end()

    if name is None:
        return

    s = code.output()

    src = open("ex3.c").read()
    match = '\n#include "body.h"\n'
    assert match in src
    src = src.replace(match, s)
    assert name and name.endswith(".c")
    f = open(name, 'w')
    tag = hash(src)
    print("hash(src):", tag)
    f.write(src)
    f.close()

    import socket
    host = socket.gethostname()
    if host == "bucket":
        cmd = "gcc MATCH.c -O3 -o MATCH -I/home/simon/local/petsc/arch-linux2-c-debug/include -I/home/simon/local/petsc/include/petsc/mpiuni -I/home/simon/local/petsc/include -I/home/simon/local/slepc-3.7.1/include -I/home/simon/local/slepc-3.7.1/arch-linux2-c-debug/include/ -L/home/simon/local/petsc/arch-linux2-c-debug/lib -L/home/simon/local/slepc-3.7.1/arch-linux2-c-debug/lib -lpetsc -lslepc"
    elif host == "hero":
        cmd = "gcc MATCH.c -O3 -o MATCH -I/usr/include/openmpi -I/usr/include/petsc -I/usr/include/slepc -lpetsc -lslepc -lmpi"
    else:
        cmd = "gcc -O3 MATCH.c -I/suphys/sburton/include/ -o MATCH -lpetsc -L$PETSC_DIR/$PETSC_ARCH/lib -L$SLEPC_DIR/$PETSC_ARCH/lib -lslepc"

    cmd = cmd.replace("MATCH.c", name)
    stem = name[:-2]
    cmd = cmd.replace("MATCH", stem)

    rval = os.system(cmd)
    assert rval == 0
    #print("exec:", hash(open(stem).read()))

    nev = argv.get("nev", 1)
    cmd = "./%s -eps_nev %d -eps_ncv %d -eps_largest_real" 

    if argv.plot:
        cmd += " -eps_view_vectors binary:evec.bin "

    cmd = cmd%(stem, nev, max(2*nev, 4))
    eps_tol = argv.get("eps_tol", 1e-4)
    if eps_tol is not None:
        cmd += " -eps_tol %s "%eps_tol

    #cmd += " -eps_type arnoldi -info -eps_monitor -eps_tol 1e-3"
    print cmd
    #rval = os.system(cmd)
    #assert rval == 0
    f = os.popen(cmd)
    s = f.read()
    #print(s)
    lines = s.split('\n')
    vals = []
    for line in lines:
        line = line.strip()
        flds = line.split()
        #print("parse", flds)
        try:
            a, b = flds
            a = float(a)
            b = float(b)
            vals.append(a)
        except:
            pass

    if not argv.plot:
        print("vals:", vals)
        return vals

    assert argv.plot.endswith(".pdf")

    s = open("evec.bin").read()
    sz = 8*2**r

    if len(s)==sz+8:
        s = s[8:]
    elif len(s)==sz+16:
        s = s[16:]
    #elif len(s)==2*(sz+16): # got two vectors
    #    s = s[16:16+sz] # pick the first vector
    elif len(s)%(sz+16) == 0:
        count = len(s)/(sz+16)
#        s = s[16:16+sz] # pick the first vector
        ev_idx = argv.get("ev_idx", 0)
        s = s[16+ev_idx*(16+sz):(ev_idx+1)*(16+sz)]
    else:
        assert 0, "sz=%d but s=%s"%(sz, len(s))

    vec0 = numpy.fromstring(s, dtype=">d")
    assert len(vec0)==2**r

    assert excite is None

    print "graphing..."
    gz, n = Gz.shape
    xdata = []
    lookup = {}
    GzRxt = dot2(Gz, Rx.transpose())
    for i, v in enumerate(genidx((2,)*r)):
        v = array2(v)
        lookup[v.tostring()] = i
        syndrome = dot2(GzRxt, v)
        value = gz - 2*syndrome.sum()
        xdata.append(value)

    pdata = {}
    ndata = {}
    my = 20. # mul y
    EPSILON = argv.get("EPSILON", 1e-6)

    def yfunc(y):
        y = log2(abs(y))
        y = int(round(my*y))
        return y

    for i in range(len(vec0)):
        x = xdata[i] # integer
        y = vec0[i]
        if abs(y) < EPSILON:
            continue
        if y > 0.:
            y = yfunc(y)
            pdata[x, y] = pdata.get((x, y), 0) + 1
        else:
            y = yfunc(y)
            ndata[x, y] = ndata.get((x, y), 0) + 1

    from pyx import graph, canvas, path, trafo, color, deco, text
    
    north = [text.halign.boxcenter, text.valign.top]
    northeast = [text.halign.boxright, text.valign.top]
    northwest = [text.halign.boxleft, text.valign.top]
    south = [text.halign.boxcenter, text.valign.bottom]
    southeast = [text.halign.boxright, text.valign.bottom]
    southwest = [text.halign.boxleft, text.valign.bottom]
    east = [text.halign.boxright, text.valign.middle]
    west = [text.halign.boxleft, text.valign.middle]
    center = [text.halign.boxcenter, text.valign.middle]
    
    c = canvas.canvas()
    sx = 0.4
    sy = 1.4
    tr = trafo.scale(sx, sy)

    green = color.rgb(0.2,0.6,0.2)
    brown = color.rgb(0.8,0.2,0.2)
    grey = color.rgb(0.4,0.4,0.4)
    lgrey = color.rgb(0.8,0.8,0.8)

    W = 2*gz
    H = log2(EPSILON)
    dy = 0.5 * 1.2/my

    X0 = -gz
    Y0 = 0.

    def showp(gx, radius):
        v = dot2(gx, PxtQx)
        syndrome = dot2(GzRxt, v)
        x = gz - 2*syndrome.sum()
        i = lookup[v.tostring()]
        #print syndrome, syndrome.sum(), vec0[i]
        y = (1./my)*yfunc(vec0[i]) + 0.5*dy
        #y = 0.5*dy + log2(abs(vec0[i]))
        c.fill(path.circle(-x*sx, y*sy, radius), [lgrey])

    showp(zeros2(n), 0.8)
    for gx in Gx:
        showp(gx, 0.4)

    for gx0 in Gx:
      for gx1 in Gx:
        gx = (gx0+gx1)%2
        if gx.sum()==0:
            continue
        showp(gx, 0.2)

    for i in range(0, gz+1):
        x, y = X0+2*i, Y0
        c.stroke(path.line(x, y, x, y+H), [tr, grey])
        if i%2 == 0:
            c.text(x*sx, y*sy + 0.2, "%d"%i, south)

    c.stroke(path.line(X0, Y0, X0+1.0*W+3.5, Y0), [tr, deco.earrow(size=0.5)])
    c.stroke(path.line(X0, Y0, X0, Y0+1.0*H-0.5), [tr, deco.earrow(size=0.5)])

    y = 1.0
    i = 0
    while y > EPSILON:

        x = X0*sx
        y1 = sy*(1./my)*yfunc(y)
        c.stroke(path.line(x, y1, x-0.1, y1))

        c.text(x-0.3, y1, "%d"%i, east)

        y /= 2.
        i -= 1

    R = 0.15
    for key, value in pdata.items():
        x, y = key
        y = y/my
        x = -x
        value = 1 + math.log(value)
        r = R*value
        #c.stroke(path.circle(x, y, r))
        #c.stroke(path.line(x, y, x+r, y), [brown, tr])
        c.fill(path.rect(x, y, r, dy), [brown, tr])

    for key, value in ndata.items():
        x, y = key
        y = y/my
        x = -x
        value = 1 + math.log(value)
        r = R*value
        #c.stroke(path.circle(x, y, r))
        #c.stroke(path.line(x-r, y, x, y), [green, tr])
        c.fill(path.rect(x-r, y, r, dy), [green, tr])

    c.writePDFfile(argv.plot)


    if 0:
        print "graph.."
    
        g = graph.graphxy(
            width=16,
            x=graph.axis.linear(reverse=True),
            y=graph.axis.linear())
            #y=graph.axis.log(min=0.8*vec0.min(), max=1.2*vec0.max()))
    
        g.plot(graph.data.values(x=xdata, y=ydata))
        g.writePDFfile(argv.plot)


