#!/usr/bin/env python
import glob
import sys
import ephem
import math
larray=[]
threeples=[]
PRESTR="BINOCULAR-"

def hours_to_float(c):
    ctoks = c.split(":")
    hours = float(ctoks[0])
    hours += float(ctoks[1])/60.0
    hours += float(ctoks[2])/3600.0
    return (hours)


count = 0
import numpy as np
LOW = -35
HIGH = 70
INCR = 5
LOWTEMP = 20.00
counting = [-1]*((HIGH-LOW)+1)
tracking = [[]]*len(counting)
RES=125

if (len(sys.argv) > 1):
    prestr = sys.argv[1]
else:
    prestr = PRESTR

prelen = len(prestr)
print "prestr %s" % prestr
for fn in glob.glob(prestr+"*-tp.dat"):
    ra = fn[prelen:prelen+5]
    ra = float(ra)
    dec = fn[prelen+6:]
    dec = dec.replace("-tp.dat", "")
    dec = int(dec)
    f = open(fn, "r")
    val = f.readline().strip("\n")
    val = float(val)
    f.close()
    if (dec >= LOW and dec <= HIGH):
        threeples.append([ra, dec, val])
        larray.append(val)
        counting[dec-LOW] += 1
        tracking[dec-LOW].append(ra)

minv = min(larray)

rows = (HIGH-LOW)/INCR
rows +=1
rows *= 2


cols = 24 * 4

pixels  = [[0 for i in range(cols)] for j in range(rows)]

for i in range(rows):
    for j in range(cols):
        pixels[i][j] = 1

count = 0

glongs = 360/INCR
glats = 180/INCR

galactic_pixels = [[0.0 for i in range(glongs)] for j in range(glats)]

for t in threeples:
    randx = t[0]
    randx *= 4.0
    randx = int(randx)
    randx %= cols

    decndx = float(t[1])
    decndx += float(-LOW)
    decndx /= float(INCR)
    decndx = int(decndx)
    decndx *= 2

    val = t[2]
    
    if (pixels[decndx][randx] == 0.0):
        pixels[decndx][randx] = val
    else:
        pixels[decndx][randx] += val
        pixels[decndx][randx] /= 2.0
        
    #pixels[decndx][randx] *= (0.000175/1.25)
    #pixels[decndx][randx] *= LOWTEMP
    count += 1

    p = ephem.Equatorial(str(t[0]), str(t[1]))
    galactic = ephem.Galactic(p)
    
    #print "%s %s" % (galactic.lon, galactic.lat)
    
    glong = math.degrees(galactic.lon)
    glong -= 180.0
    
    glat = math.degrees(galactic.lat)
    
    glndx = int(glong)+180
    glndx /= INCR
    glndx -= 1
    
    gladx = int(glat)+90
    gladx /= INCR
    gladx -= 1
    
    #print "%f %f" % (glong, glat)
    
    try:
        if (galactic_pixels[gladx][glndx] == 0.0):
            galactic_pixels[gladx][glndx] = val
        else:
            galactic_pixels[gladx][glndx] += val
            galactic_pixels[gladx][glndx] /= 2.0
    except:
        print "Hmmmm, %d %d" % (gladx, glndx)
    

#
# Since we have "dummy" rows to make the map a bit more square, we need to
#   interpolate between the rows
#
# This applies only to the equatorial map
#
for i in range(rows-1):
    if ((i % 2) != 0):
        for j in range(cols):
            pixels[i][j] = (pixels[i-1][j] + pixels[i+1][j])/2.0

for j in range(cols):
    pixels[rows-1][j] = pixels[rows-2][j] + pixels[rows-3][j]
    pixels[rows-1][j] /= 2.0

#
# Smoothing for galactic map
#
a = 0.4
b = 1.0-a
for i in range(glats-1):
    v = galactic_pixels[i][0]
    for j in range(glongs):
        cv = galactic_pixels[i][j]
        v = (a)*cv + (b)*v
        galactic_pixels[i][j] = v

for i in range(glongs):
    v = galactic_pixels[0][i]
    for j in range(glats):
        cv = galactic_pixels[j][i]
        v = (a)*cv + (b)*v
        galactic_pixels[j][i] = v


# Smoothing for equatorial map
#
a = 0.55
b = 1.0-a

#
# First dec-wise
#
for i in range(rows-1):
    v = pixels[i][0]
    for j in range(cols):
        cv = pixels[i][j]
        v = (a)*cv + (b)*v
        pixels[i][j] = v

#
# Then ra-wise
#
for i in range(cols):
    v = pixels[0][i]
    for j in range(rows):
        cv = pixels[j][i]
        v = (a)*cv + (b)*v
        pixels[j][i] = v
        
if True:
    import matplotlib.pyplot as plt
    
    if (len(sys.argv) > 2):
        RES=int(sys.argv[2])

    plt.figure(1)
    rpixels = pixels[::-1]
    plt.xticks(list(range(0,cols,8)),list(range(0,24,2)))
    plt.yticks(list(range(0,rows,4)),list(np.arange(HIGH,LOW-2,-10)))
    plt.xlabel("Hours RA")
    plt.ylabel("DEC")
    plt.title("H1 Brightness:  +/- 160km/sec")
    plt.imshow(rpixels,cmap='jet', interpolation="gaussian")
    plt.colorbar(shrink=0.7,orientation="horizontal").set_label("Est. Brightness (K)")
    plt.savefig('21cm.png', dpi=RES)
    
    plt.figure(2)
    rpixels = galactic_pixels[::-1]
    plt.xticks(list(range(0,360/INCR,8)),list(range(-180,180,40)))
    plt.yticks(list(range(0,180/INCR,4)),list(np.arange(90,-90,-20)))
    plt.xlabel("Galactic Longitude")
    plt.ylabel("Galactic Latitude")
    plt.title("H1 Brightness:  +/- 160km/sec")
    plt.imshow(rpixels,cmap='jet', interpolation="gaussian")
    plt.colorbar(shrink=0.7,orientation="horizontal").set_label("Est. Brightness (K)")
    plt.savefig('21cm-galactic.png', dpi=RES)

f = open("missing-pixels.log", "w")
for i in range(0,len(counting)):
    if (counting[i] != cols and counting[i] != -1):
        f.write("%d: %d\n" % (i+LOW, cols-counting[i]))
f.close()
