import os
import sys
import numpy
import math
import random
import ephem

def isalmost(v,testv,window):
    if (abs(v-testv) <= window):
        return True
    else:
        return False

LOW = -35
HIGH = 70
INCR = 5

TSYS=100.00
TMIN=20.00

rows = (HIGH-LOW)
rows += 1

cols = 24 * 4

pixels  = [[0 for i in range(cols)] for j in range(rows)]
pcounts = [[0 for i in range(cols)] for j in range(rows)]

for i in range(rows):
    for j in range(cols):
        pixels[i][j] = -1.0
        pcounts[i][j] = 0.0

NBINS=8192
obscount=0
BADBIN = 2768

IGNORE=int(8192/8)

TPCHUNKS=5

linenum=0
skipcount = 10
printed = False

correction = []
slen = NBINS-(IGNORE*2)
for i in range(slen):
    x = i - slen/2.0
    x = float(x)
    cv = (x/float(slen))*(x/float(slen))
    cv *= 10.0
    
    #
    # Dither it just a wee bit
    #
    cv *= random.uniform(0.98,1.02)
    correction.append(cv)

while True:
    inline=sys.stdin.readline()
    linenum += 1
    if (inline == ""):
        break
    
    inline=inline.replace("\n", "")
    inlist=inline.split(",")
    if len(inlist) < NBINS:
        continue
    
    ra=float(inlist[3])+float(inlist[4])/60.0+float(inlist[5])/3600.0
    dec=float(inlist[8])
    
    rac = ra * 60.0
    rac = rac / 15.0
    rac = round(rac)
    rac = rac * 15.0
    rac = float(rac/60.0)
    randx = rac*4.0
    randx = int(randx)
    if (randx > (24*4)-1):
        randx = (24*4)-1
    
    decndx = int(dec)-LOW
    if (dec < LOW or dec > HIGH or ra < 0.0 or ra > 24.0):
        print "%d: Bad ra/dec %f %f" % (linenum, ra, dec)
    else:
        if (skipcount > 0):
            skipcount -= 1
            continue
        
        if ((linenum % 2000) == 0):
            print "Line %d processed" % linenum
        
        if (len(inlist[9:]) == NBINS):
            failed = False
            try:
                values=map(float, inlist[9:])
            except:
                failed=True
            if failed == False:
                values=values[IGNORE:NBINS-IGNORE]
                values=numpy.divide(values, [10.0]*len(values))
                values=numpy.power([10.0]*len(values),values)
                lv = len(values)
                sstart = sum(values[0:3])
                send = sum(values[lv-3:lv])
                sstart /= 3.0
                send /= 3.0
                slope = (send-sstart) / len(values)
                desloper = [x for x in numpy.arange(0.0,slope*len(values),slope)]
                values = numpy.subtract(values,desloper)
                #
                # Determine min on desloped data
                #
                minv = min(values)
                
                #
                # Weed out large spikes
                #
                for i in range(len(values)):
                    if values[i] > (minv*1.75):
                        values[i] = minv

                if (printed == False):
                    print "slope is %e" % slope
                    fp = open ("foonly.dat" , "w")
                    for i in range(len(values)):
                        fp.write("%e\n" % values[i])
                    fp.close()
                    printed = True
                
                #
                # Now divide by min and convert to temp
                #
                values=numpy.divide(values,[minv]*len(values))
                values=numpy.multiply(values,[TSYS+TMIN]*len(values))
                values=numpy.subtract(values,[TSYS]*len(values))
                values=numpy.subtract(values,correction)
                if ((linenum % 2000) == 0):
                    print "temp range: %f %f" % (min(values), max(values))
                if(ra >= 0.0 and ra <= 24.0):
                    try:
                        pixels[decndx][randx] += sum(values)
                        pcounts[decndx][randx] += 1
                    except:
                        print "decndx %d randx %d" % (decndx, randx)
                        break

for decndx in range(rows):
    for randx in range(cols):
        if (pcounts[decndx][randx] >= 1):
            fn = "PROCESSED"+"-%05.2f-%02d-tp.dat" % (float(randx)/4.0, (decndx+LOW))
            fp = open(fn, "w")
            fp.write("%13.2f\n" % (pixels[decndx][randx]/pcounts[decndx][randx]))
            fp.close()
            
        
