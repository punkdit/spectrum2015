#!/usr/bin/env python2.7

import sys
from math import *
from random import *

from pyx import canvas, path, deco, trafo, style, text, color, deformer
from pyx.color import rgb


text.set(mode="latex") 
text.set(docopt="12pt")
text.preamble(r"\usepackage{amsmath,amsfonts,amssymb}")

text.preamble(r"\def\ket #1{|#1\rangle}")


def dopath(points, extra=[], fill=[], closepath=False, smooth=0.0):
    ps = [path.moveto(*points[0])]+[path.lineto(*p) for p in points[1:]]
    if closepath:
        ps.append(path.closepath())
    p = path.path(*ps)
    if fill:
        c.fill(p, [deformer.smoothed(smooth)]+extra+fill)
    c.stroke(p, [deformer.smoothed(smooth)]+extra)


class Turtle(object):
    def __init__(self, x, y, theta):
        self.x = x
        self.y = y
        self.theta = theta
        self.ps = [(x, y)]

    def fwd(self, d):
        self.x += d*sin(self.theta)
        self.y += d*cos(self.theta)
        self.ps.append((self.x, self.y))
        return self

    def reverse(self, d):
        self.fwd(-d)
        return self

    def right(self, dtheta, r=0.):
        theta = self.theta
        self.theta += dtheta
        if r==0.:
            return self
        N = 20
        x, y = self.x, self.y
        x0 = x - r*sin(theta-pi/2)
        y0 = y - r*cos(theta-pi/2)
        for i in range(N):
            theta += (1./(N))*dtheta
            x = x0 - r*sin(theta+pi/2)
            y = y0 - r*cos(theta+pi/2)
            self.ps.append((x, y))
        self.x = x
        self.y = y
        return self

    def left(self, dtheta, r=0.):
        self.right(-dtheta, -r)
        return self

    def stroke(self, deco=[], closepath=False):
        dopath(self.ps, deco, closepath)
        return self


black = rgb(0., 0., 0.) 
blue = rgb(0., 0., 0.8)
lred = rgb(1., 0.4, 0.4)
white = rgb(1., 1., 1.) 

#shade = rgb(0.75, 0.55, 0)

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

def arc(*args):
    return path.path(path.arc(*args))

def arcn(*args):
    return path.path(path.arcn(*args))

###############################################################################


c = canvas.canvas()


w = 1.8
h = 1.5

x, y = 0., 0.

c.text(x, y, r"$\ket{0}$", center)
c.text(x+w, y, r"$\ket{1}$", center)


r = 0.4
theta = 45.0

c.text(x-1.5, y, "$Z:$", center)
c.stroke(arc(x-1.2*r, y, r, theta, 360-theta), [deco.earrow()])
c.text(x-1.2*r, y+1.2*r, "$+1$", south)

c.stroke(arcn(x+w+1.2*r, y, r, 180-theta, 180+theta), [deco.earrow()])
c.text(x+w+1.2*r, y+1.2*r, "$-1$", south)

y -= h

r = 1.8
theta = 20.0

c.text(x, y, r"$\ket{0}$", center)
c.text(x+w, y, r"$\ket{1}$", center)

c.text(x-1.5, y, "$X:$", center)

c.text(x+0.5*w, y+0.5, "$+1$", south)
c.stroke(arc(x+0.5*w, y-0.8*r, r, 90-theta, 90+theta), [deco.earrow()])

c.text(x+0.5*w, y-0.5, "$+1$", north)
c.stroke(arc(x+0.5*w, y+0.8*r, r, 270-theta, 270+theta), [deco.earrow()])



c.writePDFfile("pic-zx.pdf")


###############################################################################


c = canvas.canvas()


w = 2.0
h = 4.0

hk = 0.2 * h
hm = 0.6 * h
hr = 0.4 * h

dm = 0.2

def axis(x, y0, y1, label):

    y = (y0+y1)*0.5
    c.stroke(path.line(x, y0, x, y1))
    #c.stroke(path.line(x, y0-0.5*dm, x, y1+0.5*dm))
    c.stroke(path.line(x-0.1, y0-0.5*dm, x+0.1, y0-0.5*dm))
    c.stroke(path.line(x-0.1, y1+0.5*dm, x+0.1, y1+0.5*dm))

    c.fill(path.circle(x, y, 0.20), [white])
    c.text(x, y, label, center)


x, y = 0., 0.

c.stroke(path.rect(x, y, 2*w, hr))
c.stroke(path.line(x+w, y, x+w, y+hr), st_dotted)

c.text(x+0.5*w, y+0.5*hr, "$R$", center)

axis(x-0.5, y, y+hr, "$r$")

y += hr + dm

c.stroke(path.rect(x, y, w-0.5*dm, hm))
c.text(x+0.5*w, y+0.5*hm, "$S$", center)

c.stroke(path.rect(x+w+0.5*dm, y, w-0.5*dm, hm))
c.text(x+1.5*w, y+0.5*hm, "$T$", center)

axis(x-0.5, y, y+hm, "$m$")

y += hm + dm

c.stroke(path.rect(x, y, 2*w, hk))
c.stroke(path.line(x+w, y, x+w, y+hk), st_dotted)
c.text(x+0.5*w, y+0.5*hk, "$L$", center)

axis(x-0.5, y, y+hk, "$k$")

axis(x+2*w+0.5, 0., hr+hm+hk+2*dm, "$n$")

x, y = 0., 0.

dm2 = 0.5*dm
t = Turtle(x-dm2, y, 0.)
t.fwd(hr+dm+hm).right(pi/2, dm2).fwd(w-dm2).right(pi/2, dm2)
t.fwd(hm).left(pi/2, dm2).fwd(w-dm2).right(pi/2, dm2).fwd(hr)
t.right(pi/2, dm2).fwd(2*w).right(pi/2, dm2)
t.stroke()

c.text(1.6*w, -0.38, "$G$", center)

c.writePDFfile("pic-canonical.pdf")



###############################################################################


c = canvas.canvas()


w = 2.0
h = 4.0

hk = 0.2 * h
hmx = 0.3 * h
hmz = 0.3 * h
hm = hmx+hmz
hr = 0.4 * h

dm = 0.2


def pair(x, y0, y1, left, right, ldots=False, rdots=False):
    p = path.rect(x, y0, w-0.5*dm, y1-y0)
    c.stroke(p, (st_dashed if ldots else []))
    c.text(x+0.5*w, 0.5*(y0+y1), "$%s$"%left, center)

    p = path.rect(x+w+0.5*dm, y0, w-0.5*dm, y1-y0)
    c.stroke(p, (st_dashed if rdots else []))
    c.text(x+1.5*w, 0.5*(y0+y1), "$%s$"%right, center)


x, y = 0., 0.

pair(x, y, y+hr, "R_X", "R_Z")

axis(x-0.5, y, y+hr, "$r$")

y += hr + dm

pair(x, y, y+hmx, "T_X", "S_Z", True, False)
axis(x-0.5, y, y+hmx, "$m_Z$")

y += hmx + dm

pair(x, y, y+hmz, "S_X", "T_Z", False, True)
axis(x-0.5, y, y+hmz, "$m_X$")

y += hmz + dm

pair(x, y, y+hk, "L_X", "L_Z", True, True)

axis(x-0.5, y, y+hk, "$k$")

axis(x+2*w+0.5, 0., hr+hm+hk+3*dm, "$n$")


c.writePDFfile("pic-symplectic.pdf")


###############################################################################



c = canvas.canvas()


w = 1.8
h = 1.0

hk = 1.0 * h
hm = 1.7 * h
hr = 1.0 * h

dm = 0.2


x = 0.
    
def mkop(row, col, op):

    x1 = x + col * w + 0.25*w
    y = row * 1.13*h + 0.5*h

    if 'tilde' in op:
        c.text(x1 + 1.5*0.16*w, y, "$%s$"%op, center)
    else:
        for i, opp in enumerate(op):
            dx = i * 0.16*w
            c.text(x1+dx, y, "$%s$"%opp, center)


for side in [-1, 0, 1]:

    y = 0.

    c.stroke(path.rect(x, y, 2*w, hr))
    c.stroke(path.line(x+w, y, x+w, y+hr), st_dotted)
    
    if side == 0:
        ops = """ZIZI XXII ZZZZ IIIX XXXX ZZZI ZZII XIXI""".strip().split()
    else:
        ops = r"""
            \tilde{X}_1 \tilde{Z}_1 
            \tilde{X}_2 \tilde{Z}_2 
            \tilde{X}_3 \tilde{Z}_3
            \tilde{X}_4 \tilde{Z}_4
            """.strip().split()
        ops = list(reversed(ops))
    
    if side>=0:
        for col in [0,1]:
         for irow, row in enumerate([0., 1., 1.7, 2.75]):
            mkop(row, col, ops[col+2*irow])
    
    y += hr + dm
    
    c.stroke(path.rect(x, y, w-0.5*dm, hm))
    c.stroke(path.rect(x+w+0.5*dm, y, w-0.5*dm, hm))
    
    y += hm + dm
    
    c.stroke(path.rect(x, y, 2*w, hk))
    c.stroke(path.line(x+w, y, x+w, y+hk), st_dotted)
    
    y = 0.0
    
    dm2 = 0.5*dm
    t = Turtle(x-dm2, y, 0.)
    t.fwd(hr+dm+hm).right(pi/2, dm2).fwd(w-dm2).right(pi/2, dm2)
    t.fwd(hm).left(pi/2, dm2).fwd(w-dm2).right(pi/2, dm2).fwd(hr)
    t.right(pi/2, dm2).fwd(2*w).right(pi/2, dm2)
    t.stroke()
    
    if 0:
        c.text(x-0.5, 0.5*hr, "$R$", center)
        c.text(x-0.5, hr+dm+0.5*hm, "$S$", center)
        #c.text(x+2*w+0.5, hr+dm+0.5*hm, "$T$", center)
        c.text(x+2*w+0.5, hr+dm+0.2*hm, "$T$", center)
        c.text(x-0.5, hr+dm+hm+dm+0.5*hk, "$L$", center)
        c.text(x+1.6*w, -0.38, "$G$", center)

    if side==-1:
        c.text(x+0.5, 0.5*hr, "$R$", center)
        c.text(x+0.5, hr+dm+0.5*hm, "$S$", center)
        c.text(x+1*w+0.5, hr+dm+0.5*hm, "$T$", center)
        c.text(x+0.5, hr+dm+hm+dm+0.5*hk, "$L$", center)
        c.text(x+1.6*w, -0.38, "$G$", center)

    if side<=0:
        c.text(x+2*w+0.4, hr+dm+0.5*hm, "$=$", center)

    x += 2.5*w


c.writePDFfile("pic-gauge4.pdf")


###############################################################################



def state(x, pos=True):
    #x0, x1 = x.split()
    x0, x1 = x
    return r"$|%s\rangle %s |%s\rangle$" % (x0, "+" if pos else "-", x1)

def orbit((x, y), desc, pos=True, Z=True, X=True, logop=False):

    r = 0.5
    y -= 0.5*rh

    desc = desc.split(" ")
    if logop:
        for i, s in enumerate(desc):
            s = list(s)
            s[0] = '1' if s[0]=='0' else '0'
            s[2] = '1' if s[2]=='0' else '0'
            s = ''.join(s)
            desc[i] = s

    s1 = desc[:2]
    s2 = desc[2:]

    c.text(x, y, state(s1, pos), center)
    c.text(x, y+rh, state(s2, pos), center)

    if frame == 1:
        return

    delta = 50

    if Z:
        c.stroke(arc(x, y+rh+1.5*r, r, -delta, 180+delta), [deco.earrow()])
        c.stroke(arc(x, y-1.5*r, r, 180-delta, delta), [deco.earrow()])
        c.text(x+1.8*r, y+rh+2.2*r, "$-2$", center)
        c.text(x+1.8*r, y-2.2*r, "$+2$", center)

    delta = 40
    if X:
        c.stroke(arc(x+1.9*r, y+0.5*rh, 1.3*r, -delta, delta), [deco.earrow()])
        c.stroke(arc(x-1.9*r, y+0.5*rh, 1.3*r, 180-delta, 180+delta), [deco.earrow()])
        c.text(x+3.9*r, y+0.5*rh, "$+2$", center)
        c.text(x-3.9*r, y+0.5*rh, "$+2$", center)



rh = 1.0

x, y = 0., 0.
w, h = 5., 1.7


x0 = x-0.3*w
y0 = y+2.5*h


for frame in [0, 1]:
    c = canvas.canvas()
    
    if 0:
        c.stroke(path.line(x0, y0, x+2.3*w, y0), st_dashed)
        
        for r in [-1.0, 1.0]:
            c.stroke(path.line(x0, y0-r, x0, y0+r), [deco.earrow(size=0.5)]+ st_Thick)
        p = path.circle(x0, y0, 0.4)
        extra = [trafo.scale(x=x0, y=y0, sx=2, sy=1.)]
        c.fill(p, [white]+extra)
        c.stroke(p, extra)
        c.text(x0, y0, r"$XIXI$", center)
    
    
    pts = {}
    
    for logop in [0, 1]:
        pts[logop, 0, 0] = x, y+0.5*rh
        pts[logop, 1, 0] = x+1.2*w, y+h +0.5*rh
        pts[logop, 0, 1] = x+0.8*w, y-h +0.5*rh
        pts[logop, 1, 1] = x+2*w, y-0.2*h +0.5*rh
        y += 4*h
    
    
    grey = rgb(0.8, 0.8, 0.8)
    def line(i0, i1):
        x0, y0 = pts[i0]
        x1, y1 = pts[i1]
        c.stroke(path.line(x0, y0, x1, y1), [grey]+st_THICK+st_round)
    
    x0, y0 = pts[0, 0, 0]
    x1, y1 = pts[0, 0, 1]
    #c.stroke(path.line(x0, y0, x1, y1), [grey]+st_THICK+st_round)
    line((0,0,0), (0,0,1))
    for i0 in [0,1]:
     for i1 in [0,1]:
      for i2 in [0,1]:
       for j0 in [0,1]:
        for j1 in [0,1]:
         for j2 in [0,1]:
            d0 = j0-i0
            d1 = j1-i1
            d2 = j2-i2
            if d0<0 or d1<0 or d2<0: continue
            if d0+d1+d2!=1: continue
            line((i0,i1,i2),(j0,j1,j2))
    
    for logop in [0, 1]:
    
        orbit(pts[logop, 0, 0], "0000 1111 1100 0011", logop=logop)
        orbit(pts[logop, 1, 0], "0000 1111 1100 0011", pos=False, X=False, logop=logop)
        orbit(pts[logop, 0, 1], "0001 1110 1101 0010", Z=False, logop=logop)
        orbit(pts[logop, 1, 1], "0001 1110 1101 0010", pos=False, Z=False, X=False, logop=logop)
    
    
    if frame == 0:
        c.writePDFfile("pic-orbit.pdf")
        continue


    def conv(a, b, r=0.5):
        return r*a + (1-r)*b

    def darrow(x0, y0, x1, y1):
        x2, y2 = conv(x0, x1), conv(y0, y1)
        c.stroke(path.line(x2, y2, x1, y1), [deco.earrow(size=0.2)]+st_Thick)
        c.stroke(path.line(x2, y2, x0, y0), [deco.earrow(size=0.2)]+st_Thick)

    def operator(x0, y0, x1, y1, label, dx=0., dy=0., arrow=True):
        x0 += dx; x1 += dx; y0 += dy; y1 += dy;
        x2, y2 = conv(x0, x1), conv(y0, y1)
        x0, x1 = conv(x0, x1, 0.8), conv(x0, x1, 0.2)
        y0, y1 = conv(y0, y1, 0.8), conv(y0, y1, 0.2)
        if arrow:
            darrow(x0, y0, x1, y1)

        p = path.circle(x2, y2, 0.4)
        extra = [trafo.scale(x=x2, y=y2, sx=4, sy=1.)]
        c.fill(p, [white]+extra)
        #c.stroke(p, extra)
        c.text(x2, y2, r"$%s$"%label, center)


    x0, y0 = pts[1, 1, 0]
    x1, y1 = pts[1, 1, 1]
    operator(x0, y0, x1, y1, r"\tilde{X}_3 = IIIX", dx=1.5, dy=1.0)

    x0, y0 = pts[0, 0, 0]
    x1, y1 = pts[1, 0, 0]
    operator(x0, y0, x1, y1, r"\tilde{X}_1 = XIXI", dx=-2.0)

    x0, y0 = pts[1, 0, 0]
    x1, y1 = pts[1, 1, 0]
    operator(x0, y0, x1, y1, r"\tilde{X}_2 = ZZZI", dx=-0.5, dy=1.0)

    x0, y0 = pts[0, 1, 0]
    x1, y1 = pts[0, 1, 0]
    y0 -= 0.6
    y1 += 0.6
    operator(x0, y0, x1, y1, r"\tilde{X}_4 = XXII", dx=3.0, arrow=False)
    darrow(x0+1.5, y0, x1+1.5, y1)

    c.writePDFfile("pic-operators.pdf")



