#!/usr/bin/env python

import sys, os
from math import *

from pyx import canvas, path, deco, trafo, style, text, color, deformer
from pyx.color import rgb

text.set(mode="latex") 
text.set(docopt="12pt")
text.preamble(r"\usepackage{amsmath,amsfonts,amssymb}")



north = [text.halign.boxcenter, text.valign.top]
northeast = [text.halign.boxright, text.valign.top]
northwest = [text.halign.boxleft, text.valign.top]
south = [text.halign.boxcenter, text.valign.bottom]
southeast = [text.halign.boxright, text.valign.bottom]
southwest = [text.halign.boxleft, text.valign.bottom]
east = [text.halign.boxright, text.valign.middle]
west = [text.halign.boxleft, text.valign.middle]
center = [text.halign.boxcenter, text.valign.middle]

black = rgb(0., 0., 0.) 
blue = rgb(0.1, 0.1, 0.9)
lred = rgb(0.8, 0.4, 0.2)
green = rgb(0.0, 0.6, 0.0)
white = rgb(1., 1., 1.) 
#shade = rgb(0.75, 0.55, 0)
grey = rgb(0.75, 0.75, 0.75)
yellow = rgb(1., 1., 0.)

st_dashed = [style.linestyle.dashed]
st_dotted = [style.linestyle.dotted]
st_round = [style.linecap.round]

st_thick = [style.linewidth.thick]
st_Thick = [style.linewidth.Thick]
st_THick = [style.linewidth.THick]
st_THICK = [style.linewidth.THICK]
#print dir(style.linewidth)


c = canvas.canvas()

W, H = 14., 8.

DX = 175.0
DY = 4.0
mx = (W/DX)
my = (H/DY)

ticksize = 0.2
def plot_axes(xtext="qubits: $n$", ytext=r"gap: $\lambda_1-\lambda_2$"):
    c.stroke(path.line(0., 0., 1.05*W, 0.), [deco.earrow(size=0.3)])
    c.stroke(path.line(0., 0., 0., 1.08*H), [deco.earrow(size=0.3)])
    
    for i in range(0, 176, 25):
        x = i*mx
        c.stroke(path.line(x, 0., x, ticksize))
        c.text(x, -2*ticksize, str(i), north)
    
    for i in range(5):
        y = i*my
        c.stroke(path.line(0., y, ticksize, y))
        c.text(-2*ticksize, y, str(i), center)
    
    c.text(W/2., -6*ticksize, xtext, center)
    
    x, y = -5*ticksize, H/2.
    c.text(x, y, ytext, [trafo.rotate(90, 0., 0)]+center)

plot_axes()

#g = graph.graphxy(
#    width=14, height=6,
#    x=graph.axis.linear(min=0., title=r"qubits: $n$"), 
#    y=graph.axis.linear(min=0., title=r"gap: $\lambda_1-\lambda_2$"))

def plot(xs, ys, clr, legend, row):

    n = len(xs)
    xs = [mx*x for x in xs]
    ys = [my*y for y in ys]
    for i in range(n-1):
        x0, y0 = xs[i], ys[i]
        x1, y1 = xs[i+1], ys[i+1]
        c.stroke(path.line(x0, y0, x1, y1), st_Thick + [clr])

    r = 0.13
    for x, y in zip(xs, ys):
        p = path.circle(x, y, r)
        c.stroke(p, st_THick)
        c.fill(p, [clr])

    x = 0.6*W
    y = 0.9*H - 0.7*row
    c.text(x+3*r, y, legend, west)
    p = path.circle(x, y, r)
    c.stroke(p, st_THick)
    c.fill(p, [clr])
    


# XY
# note: n=25 has no stabilizers: 31.851942-31.349618
#xs = [16, 24]
#ys = [20.503324-20.109358, 30.645190-30.383016]
#plot(xs, ys, green, "1D $XY$", 0)

# 3D Compass
xs = [27]
ys = [0.538]
plot(xs, ys, yellow, "3D compass", 2)

# Compass
xs = [16, 25, 36]
ys = [0.644, 0.452, 0.316]
plot(xs, ys, blue, "2D compass", 1)

# 2D XY-model
xs = [36, 64, 100, 144]
ys = [1.93021, 1.28720, 0.90439, 0.63184]
plot(xs, ys, green, "2D XY", 3)

# gauge color 
xs = [15, 65, 175]
ys = [3.241, 1.694, 1.049]
plot(xs, ys, lred, "3D gauge color code", 4)


c.writePDFfile("pic-gap.pdf")


###############################################################################

c = canvas.canvas()

plot_axes(ytext=r"gap")

clrs = {0:white, 8:blue, 12:green, 18:lred, 24:yellow}

r = 0.13

for row, weight in enumerate([8, 12, 18, 24]):
    x = 10*mx + 1*ticksize
    y = 1.5*my - 0.7*row

    p = path.circle(x, y, r)
    c.fill(p, [clrs[weight]])
    c.stroke(p, st_thick)

    c.text(x+2*ticksize, y, "$w(s_Z)=%d$"%weight, west)

#data = {}
##for n in [15, 65, 175]:
#for weight in [8, 12, 18, 24]:
#    data[weight] = []

def plot(n, a, weight, a1):
    x = n*mx
    y = my*(a-a1)
    p = path.circle(x, y, r)
    c.fill(p, [clrs[weight]])
    c.stroke(p, st_thick)
#    data[weight].append((x, y))

def doline(n0, n1, delta0, delta1, clr, extra=[]):
    x0 = n0*mx
    x1 = n1*mx
    y0 = my*delta0
    y1 = my*delta1
    c.stroke(path.line(x0, y0, x1, y1), [clr]+extra)


n0, n1, n2 = 15, 65, 175
a0 = 25.455844  
a1 = 104.076026  
a2 = 267.197576  

doline(n0, n1, a0-22.214755, a1-102.382483, black, st_dashed)
doline(n1, n2, a1-102.382483, a2-266.148188, black, st_dashed)

doline(n0, n1, a0-22.214755, a1-100.429340, blue)

doline(n1, n2, a1-100.429340, a2-263.171190 , blue)
doline(n1, n2, a1-100.429340, a2-263.324858 , blue)
doline(n1, n2, a1-100.429340, a2-263.340832 , blue)

doline(n1, n2, a1-100.585413, a2- 264.269635, green)
doline(n1, n2, a1-100.585413, a2- 264.617135, green)
doline(n1, n2, a1-100.585413, a2- 264.745548, green)

doline(n1, n2, a1-101.602340, a2- 264.269635, green)
doline(n1, n2, a1-101.602340, a2- 264.617135, green)
doline(n1, n2, a1-101.602340, a2- 264.745548, green)

doline(n1, n2, a1-102.382483, a2- 264.843629, lred)
doline(n1, n2, a1-102.382483, a2- 265.413935, lred)
doline(n1, n2, a1-102.382483, a2- 265.754772, lred)

#plot(n0, a0, 0, 16.970563)
plot(n0, a0, 8, 22.214755)

#plot(n1, a1, 0, 99.014097)
plot(n1, a1, 8,  100.429340)
plot(n1, a1, 12, 100.585413)
plot(n1, a1, 12, 101.602340)
plot(n1, a1, 18, 102.382483)

plot(n2, a2, 8,  263.171190 )
plot(n2, a2, 8,  263.324858 )
plot(n2, a2, 8,  263.340832 )
plot(n2, a2, 12, 264.269635 )
#plot(n2, a2, 0, 264.250644)
plot(n2, a2, 12, 264.617135)
plot(n2, a2, 12, 264.745548 )
plot(n2, a2, 18, 264.843629)
plot(n2, a2, 18, 265.413935)
plot(n2, a2, 18, 265.754772)
plot(n2, a2, 24, 266.148188)



c.writePDFfile("pic-gap-stabs.pdf")



