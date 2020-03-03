import os
import sys
import numpy
import math
import random
import ephem

def dumpit(vals, fn, start, incr):
    axis = start
    try:
        fp = open(fn, "w")
        for x in vals:
            fp.write("%e %e\n" % (axis, x))
            axis += incr
    except:
        print "Having trouble in dumpit"
    fp.close()

def interesting(ra, dec):
    if abs(ra - 19.9) < 0.20 and abs(dec - 40.75) < 5.0:
        return True
    else:
        return False
        

def isalmost(v,testv,window):
    if (abs(v-testv) <= window):
        return True
    else:
        return False

LOW = -35
HIGH = 70
INCR = 5

TSYS=120
TMIN=15.0

rows = (HIGH-LOW)
rows += 1

cols = 24 * 4

pixels  = [[0 for i in range(cols)] for j in range(rows)]
pcounts = [[0 for i in range(cols)] for j in range(rows)]
subbands = [[[0,0,0,0,0] for i in range(cols)] for j in range(rows)]


for i in range(rows):
    for j in range(cols):
        pixels[i][j] = -1.0
        pcounts[i][j] = 0.0

NBINS=8192
obscount=0
BADBIN = 2768
Fc=1420.40575
BW=2.5
START=Fc-(BW/2.0)
PERBIN=BW/NBINS

IGNORE=int(8192/8)

TPCHUNKS=5

linenum=0
skipcount = 10
printed = False

PMEDIAN=7

# 
# Median filter/decimate the input by PMEDIAN
#
ix = []
ox = []
def median(p,pm):
    global ix
    global ox
    out = []
    for i in range(0,len(p),pm):
        cm = p[i:i+pm]
        #
        # Calculate a weighted-median value
        #   where the min value in the window has more effect
        #
        mv = min(cm)*3.0
        mv += numpy.median(cm)
        mv /= 4.0
        for j in range(pm):
            tweakage = random.uniform(mv*0.999,mv*1.001)
            out.append(tweakage)
    return (out)

#
# Compute a fixed correction curve, to compensate for "dishing" effect
#
correction = []
def compute_correction(bins,value):
    localcorr = []
    slen = bins
    for i in range(slen):
        x = i - slen/2.0
        x = float(x)
        cv = (x/float(slen))*(x/float(slen))
        cv *= 4.0
        cv *= value
        
        #
        # Dither it just a wee bit
        #
        cv *= random.uniform(0.98,1.02)
        localcorr.append(cv)
    fp = open("correction.dat", "w")
    for v in localcorr:
        fp.write ("%s\n" % v)
    fp.close()
    return localcorr

slopeinit = False
submap = []
histo = []
badspots = []
fndx = 1
frozen = False
curfile = open(sys.argv[fndx], "r")
print "Processing file: %s" % sys.argv[fndx]
pinrted = False
masky1 = []
masky2 = []
lastminv = 0.0
while True:
    inline=curfile.readline()
    linenum += 1
    if (inline == ""):
        curfile.close()
        fndx += 1
        if (fndx >= len(sys.argv)):
            break
        curfile = open(sys.argv[fndx], "r")
        print "Processing file: %s" % sys.argv[fndx]
        linenum = 0
    
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
        #print "%s:%d: Bad ra/dec %f %f" % (sys.argv[fndx], linenum, ra, dec)
        pass
    else:
        if (skipcount > 0):
            skipcount -= 1
            continue
        
        if (len(inlist[9:]) == NBINS):
            failed = False
            try:
                values=map(float, inlist[9:])
            except:
                failed=True
            if failed == False:
                values=values[IGNORE:NBINS-IGNORE]
                q = numpy.array(values)
                if (len(submap) == 0):
                    bndx = 0
                    for q in range(5):
                        submap.append(bndx)
                        bndx += len(values)/5
                values=numpy.divide(values, [10.0]*len(values))
                values=numpy.power([10.0]*len(values),values)
                slopeinit = True
                lv = len(values)
                sstart = sum(values[0:3])
                send = sum(values[lv-3:lv])
                sstart /= 3.0
                send /= 3.0
                slope = (send-sstart) / len(values)
                if (abs(slope) != 0.0):
                    try:
                        desloper = [x for x in numpy.arange(0.0,slope*len(values),slope)]
                    except:
                        print "Failure in desloper near line %d in file %s.  Slope %f" % (linenum, sys.argv[fndx], slope)
                        print "send %f sstart %f len %d" % (send, sstart, len(values))
                        desloper = [0.0]*len(values)
                else:
                    desloper=[0.0]*len(values)
                values = numpy.subtract(values,desloper)
                values = median(values,PMEDIAN)
                #
                # Determine min/max on desloped data
                #
                minv = min(values)
 
                    
                #
                # Weed out large spikes
                #
                qq = range(len(values))
                if (frozen == True):
                    qq = []
                for i in qq:
                    if values[i] > (minv*1.70):
                        #values[i] = random.uniform(0.99*minv, 1.01*minv)
                        
                        #
                        # Make a running histogram of large spikes
                        #
                        if (len(histo) == 0):
                            histo = [0]*len(values)
                        histo[i] += 1
                        
                        #
                        # Once in a while, do something about this histogram
                        #
                        if ((random.randint(0,65536*512) % 2000) == 0):
                            fp = open ("histo.dat", "w")
                            for h in histo:
                                fp.write("%d\n" % h)
                            fp.close()
                            hmax = max(histo)
                            
                            #
                            # Once we have this many samples, we probably
                            #   have a good idea of the RFI distribution
                            #
                            if (hmax > 8000):
                                frozen = True
                                masky1 = [1.0]*len(values)
                                masky1 = numpy.array(masky1)
                                for ndx in badspots:
                                    for offs in range(10):
                                        masky1[ndx-offs] = 0.0
                                        masky1[ndx+offs] = 0.0
                                
                            havg = sum(histo[0:20])
                            havg += sum(histo[-20:])
                            havg /= 40
                            
                            #
                            # We have a "valid" histogram of spikes
                            #
                            if (hmax > 10*havg):
                                q = numpy.array(histo)
                                hres = numpy.where(q >= (hmax-1))
                                hres = hres[0]
                                for h in hres:
                                    if (h not in badspots):
                                        badspots.append(h)
                #
                # minimize bad spots--RFI spikes
                #
                if (lastminv == 0.0):
                    lastminv = minv
                
                #
                # Update masky2 whenever minv changes by a "significant" amount
                #  
                minratio = minv/lastminv
                if ((minratio < 0.99 or minratio > 1.01) and frozen == True):
                    masky2 = [0.0] * len(values)
                    masky2 = numpy.array(masky2)
                    for ndx in numpy.where(masky1 == 0.0)[0]:
                        masky2[ndx] = minv
                    lastminv = minv
                
                if (frozen == False):                             
                    for ndx in badspots:
                        for offset in range(10):
                            values[ndx-offset] = minv
                            values[ndx] = minv
                            values[ndx+offset] = minv
                else:
                    values = numpy.multiply(values,masky1)
                    values = numpy.add(values,masky2)
                
    
                maxv = max(values)
                #
                # Now divide by min and convert to temperature
                #  Normalized range will be limited to 1.0 to 2.5, typically
                #
                # We then map that through a bit of algebra to an
                #  antenna-temperature estimate, based on knowledge of TSYS
                #  and TMIN of the typical sky at 21cm
                #
                values=numpy.divide(values,[minv]*len(values))
                values=numpy.multiply(values,[TSYS*1.5+TMIN]*len(values))
                values=numpy.subtract(values,[TSYS*1.5]*len(values))
                if (len(correction) == 0):
                    correction = compute_correction(len(values),5.5)
                values=numpy.subtract(values,correction)
                if ((int(time.time()) % 180) == 0):
                    fp = open("foo.dat", "w")
                    for v in values:
                        fp.write ("%f\n" % v)
                    fp.close()
                    fp.close()
                if(ra >= 0.0 and ra <= 24.0):
                    try:
                        #
                        # We compute a kind of weighted-average
                        #  between the average temperature across the
                        #  spectrum, and the peak temperature
                        #
                        temp = (sum(values)/len(values))*2.0
                        temp += max(values)
                        temp /= 3.0
                        
                        pixels[decndx][randx] += temp
                        pcounts[decndx][randx] += 1
                        bndx = 0
                        lmap = len(values)/5
                        for q in range(5):
                            subbands[decndx][randx][q] = sum(values[bndx:bndx+lmap])
                            bndx += lmap
                    except:
                        print "file %s:%d decndx %d randx %d" % (sys.argv[fndx], linenum, decndx, randx)
                        break
                else:
                    print "file %s:%d Ra %f dec %f invalid!" % (sys.argv[fndx], linenum, ra, dec)
                    

for decndx in range(rows):
    for randx in range(cols):
        if (pcounts[decndx][randx] >= 1):
            fn = "PROCESSED"+"-%05.2f-%02d-tp.dat" % (float(randx)/4.0, (decndx+LOW))
            fp = open(fn, "w")
            temp = pixels[decndx][randx]/pcounts[decndx][randx]
            fp.write("%13.2f\n" % temp)
            for q in range(5):
                v = subbands[decndx][randx][q]/pcounts[decndx][randx]
                fp.write ("%13.2f " % v)
            fp.write("\n")
            fp.close()
            
        
