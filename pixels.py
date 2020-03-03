#!/usr/bin/env python
import glob
import sys
larray=[]
threeples=[]
PRESTR="BINOCULAR-"

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
    
    pixels[decndx][randx] = val
    #pixels[decndx][randx] *= (0.000175/1.25)
    #pixels[decndx][randx] *= LOWTEMP
    count += 1


for i in range(rows-1):
    if ((i % 2) != 0):
        for j in range(cols):
            pixels[i][j] = (pixels[i-1][j] + pixels[i+1][j])/2.0

for j in range(cols):
    pixels[rows-1][j] = pixels[rows-2][j] + pixels[rows-3][j]
    pixels[rows-1][j] /= 2.0
        
if True:
    import matplotlib.pyplot as plt
    
    if (len(sys.argv) > 2):
        RES=int(sys.argv[2])

    rpixels = pixels[::-1]
    plt.xticks(list(range(0,cols,8)),list(range(0,24,2)))
    plt.yticks(list(range(0,rows,4)),list(np.arange(HIGH,LOW-2,-10)))
    plt.xlabel("Hours RA")
    plt.ylabel("DEC")
    plt.title("H1 Brightness:  +/- 160km/sec")
    plt.imshow(rpixels,cmap='jet', interpolation="gaussian")
    plt.colorbar(shrink=0.7,orientation="horizontal").set_label("Est. Brightness (K)")
    plt.savefig('21cm.png', dpi=RES)

f = open("missing-pixels.log", "w")
for i in range(0,len(counting)):
    if (counting[i] != cols and counting[i] != -1):
        f.write("%d: %d\n" % (i+LOW, cols-counting[i]))
f.close()
