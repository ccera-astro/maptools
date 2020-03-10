import os
import sys
import numpy
import math
import random
import ephem
import time
import copy

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

def iir_filter(v,a):
    b = 1.0-a
    y = v[0]
    
    out = []
    for x in v:
        ov = (x*a)+(y*b)
        y = ov
        out.append(y)
    out = numpy.array(out)
    return out

def to_temp(spec,method):
	if (method == "max"):
		t = max(spec)
	elif (method == "avg"):
		t = sum(spec)/len(spec)
	elif (method == "mix"):
		t = sum(spec)/len(spec)
		t += max(spec)*2.0
		t /= 3.0
	elif (method == "sort"):
		spec.sort()
		t = spec[-1]*2.0
		t += spec[-2]
		t /= 3.0
	else:
		t = min(spec)+max(spec)
		t /= 2.0
	return (t)

#
# Compute a fixed correction curve, to compensate for "dishing" effect
#
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
    return localcorr

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
IGNORE=int(8192/8)

def squash(v,ratio):
    out = []
    for i in range(0,len(v),ratio):
        ov = sum(v[i:i+ratio])/ratio
        out.append(ov)
    out = numpy.array(out)
    return out

def get_next_fn():
    if (len(filelist) <= 0):
        return None
    else:
        fn = filelist.pop()
        return fn

def dump_pixels(pixels,subbands,pcounts,rows,cols):
    minoverall = 1.0e9
    for decndx in range(rows):
        for randx in range(cols):
            if (pcounts[decndx][randx] >= 1):
                v = pixels[decndx][randx]/pcounts[decndx][randx]
                if (v < minoverall):
                    minoverall = v

    for decndx in range(rows):
        for randx in range(cols):
            if (pcounts[decndx][randx] >= 1):
                fn = "PROCESSED"+"-%05.2f-%02d-tp.dat" % (float(randx)/4.0, (decndx+LOW))
                fp = open(fn, "w")
                temp = pixels[decndx][randx]/pcounts[decndx][randx]
                temp -= minoverall
                temp *= 3.0
                temp += TMIN
                fp.write("%13.2f\n" % temp)
                for q in range(5):
                    v = subbands[decndx][randx][q]/pcounts[decndx][randx]
                    v -= minoverall
                    v *= 3.0
                    v += TMIN
                    
                    fp.write ("%13.2f " % v)
                fp.write("\n")
                fp.close()

def process_files(tndx,method):
    global pixels
    global pcounts
    global subbands
    global rows
    global cols
    global filelist
    
    
    linenum=0
    skipcount = 10
    
    submap = []
    histo = []
    badspots = []
    fndx = 1
    frozen = False
    masky1 = []
    masky2 = []
    lastminv = 0.0
    correction_dict = {}

    #
    # Pop first one off the queue
    #
    nextfn = get_next_fn()
    print "Thread %d: Processing %s" % (tndx, nextfn)
    if (nextfn != None):
        curfile = open(nextfn, "r")
    else:
        return True
    minv = None
    maxv = None
    while True:
        inline=curfile.readline()
        linenum += 1
        if (inline == ""):
            curfile.close()
            print "%d: Finished with %s" % (tndx, nextfn)
            linenum = 0
            nextfn = get_next_fn()
            if (nextfn == None):
                break
            else:
                print "Thread %d: Processing %s" % (tndx, nextfn)
                curfile = open(nextfn, "r")
        
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
                    values = squash(values,3)
                    q = numpy.array(values)
                    if (len(submap) == 0):
                        bndx = 0
                        for q in range(5):
                            submap.append(bndx)
                            bndx += len(values)/5
                    values=numpy.divide(values, [10.0]*len(values))
                    values=numpy.power([10.0]*len(values),values)
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
                            print "Failure in desloper near line %d in file %s.  Slope %f" % (linenum, nextfn, slope)
                            print "send %f sstart %f len %d" % (send, sstart, len(values))
                            desloper = [0.0]*len(values)
                    else:
                        desloper=[0.0]*len(values)
                    values = numpy.subtract(values,desloper)
                    values = iir_filter(values,0.2)

                    #
                    # Determine min/max on desloped data
                    #
                    minv = min(values)
   
                    #
                    # Weed out large spikes
                    #
                    qq = numpy.where(values > (minv*1.70))[0]
                    if (frozen == True):
                        qq = []
                    for i in qq:
                        if True:
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
                            if ((random.randint(0,10000000) % 1000) == 0):
                                hmax = max(histo)
                                
                                #
                                # Once we have this many samples, we probably
                                #   have a good idea of the RFI distribution
                                #
                                if (hmax > 4000):
                                    frozen = True
                                    masky1 = [1.0]*len(values)
                                    masky1 = numpy.array(masky1)
                                    for ndx in badspots:
                                        for offs in range(6):
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
                    # Initialize lastminv
                    #
                    if (lastminv == 0.0):
                        lastminv = minv
                    
                    #
                    # Update masky2 whenever minv changes by a "significant" amount
                    #  
                    minratio = minv/lastminv
                    if ((minratio < 0.98 or minratio > 1.02) and frozen == True):
                        masky2 = [0.0] * len(values)
                        masky2 = numpy.array(masky2)
                        for ndx in numpy.where(masky1 == 0.0)[0]:
                            masky2[ndx] = minv
                        lastminv = minv
                    #
                    # We don't yet have frozen RFI info, so the masks are not ready
                    # So we "manually" apply an RFI-mask based on "badspots"
                    #
                    if (frozen == False):                             
                        for ndx in badspots:
                            for offset in range(10):
                                values[ndx-offset] = minv
                                values[ndx] = minv
                                values[ndx+offset] = minv
                    #
                    # We're frozen, so the masks are valid
                    #  Apply those masks. The idea is that numpy is
                    #  faster than manually trolling through one-by-one
                    #
                    else:
                        #
                        # First, turn the RFI zone into zeros, preserving
                        #  the non-RFI zones
                        #
                        values = numpy.multiply(values,masky1)
                        
                        #
                        # Then add in the minv into the RFI zone (which will be zeros)
                        #
                        values = numpy.add(values,masky2)
                    
                    #
                    # Now divide by min and convert to temperature
                    #  Normalized range will be limited to 1.0 to 2.5, typically
                    #
                    # We then map that through a bit of algebra to an
                    #  antenna-temperature estimate, based on knowledge of TSYS
                    #  and TMIN of the typical sky at 21cm
                    #
                    values=numpy.divide(values,[minv]*len(values))
                    values=numpy.multiply(values,[TSYS+TMIN]*len(values))
                    values=numpy.subtract(values,[TSYS]*len(values))
                    
                    #
                    # Figure out which de-dishing correction is best
                    #
                    #
                    # Haven't created a corrections array yet? Do so
                    if (len(correction_dict) == 0):
                        for corr in range(1,30):
                            correction_dict[corr] = compute_correction(len(values), corr)
                    #
                    # Find the average of the "ends"
                    #
                    high = sum(values[-5:])/5.0
                    low = sum(values[0:5])/5.0

                    #
                    # Find lowest approximate temp
                    #
                    mintemp=min(values)

                    #
                    # High and low end must be roughly equal
                    # This crude test seems to work out OK to tell
                    #   how "dished" the passband is after
                    #   slope correction.
                    #
                    if ((high/low) >= 0.92 and (high/low) <= 1.08):
                        #
                        # Average of the "ends"
                        #
                        corr_sel = (high+low)/2.0
                        
                        #
                        # Subtract out the min approx temp.
                        #
                        corr_sel -= (mintemp*1.05)
                        
                        #
                        # Quantize to integer
                        #
                        corr_sel = int(corr_sel)
                        
                        #
                        # See if this value is in the correction set
                        #
                        if (corr_sel in correction_dict):
                            correction = correction_dict[corr_sel]
                        elif corr_sel == 0:
                            correction = [0.0]*len(values)
                        else:
                            correction = correction_dict[29]
                    else:
                        correction = [0.0]*len(values)
                    
                    #
                    # Subtract-out the selected correction value
                    # Hopefully, de-dishing the passband
                    #
                    values=numpy.subtract(values,correction)
                    
                    #
                    # Might be time to write out debug data
                    #
                    if (int(time.time()) % 30 == 0):
                        fp = open("foo.dat", "w")
                        for v in values:
                            fp.write("%f\n" % v)
                        fp.close
                        
                        fp = open ("correction.dat", "w")
                        for v in correction:
                            fp.write("%f\n" % v)
                        fp.close()
                        
                    
                    if(ra >= 0.0 and ra <= 24.0):
                        try:
                            #
                            # We compute a kind of weighted-average
                            #  between the average temperature across the
                            #  spectrum, and the peak temperature
                            #
                            temp = to_temp(values, method)
                            
                            #
                            # Temp not in this range? Probably invalid
                            #
                            if (temp >= 10 and temp <= 150):
                                pixels[decndx][randx] += temp
                                pcounts[decndx][randx] += 1
                            bndx = 0
                            lmap = len(values)/5
                            for q in range(5):
                                srange = values[bndx:bndx+lmap]
                                pv = to_temp(srange,method)
                                if (pv >= 10 and pv <= 150):
                                    subbands[decndx][randx][q] = pv

                                bndx += lmap
                        except:
                            print "file %s:%d decndx %d randx %d" % (nextn, linenum, decndx, randx)
                            break
                    else:
                        print "file %s:%d Ra %f dec %f invalid!" % (nextn, linenum, ra, dec)

filelist = copy.deepcopy(sys.argv[1:])
filelist.reverse()

process_files(0,"avg")
dump_pixels(pixels,subbands,pcounts,rows,cols)
