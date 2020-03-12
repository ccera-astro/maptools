import os
import sys
import numpy
import math
import random
import ephem
import time
import copy
import argparse

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

def dump_pixels(pixels,subbands,pcounts,rows,cols,args):
    pixels = numpy.array(pixels)
    minoverall = 1.0e10
    maxoverall = 0.0
    print numpy.where(pixels < 0.0)[0]
    for decndx in range(rows):
        for randx in range(cols):
            if (pcounts[decndx][randx] >= 1):
                v = pixels[decndx][randx]/float(pcounts[decndx][randx])
                if (v < minoverall):
                    minoverall = v
                if (v > maxoverall):
                    maxoverall = v
    for decndx in range(rows):
        for randx in range(cols):
            if (pcounts[decndx][randx] >= 1):
                fn = args.prefix+"-%05.2f-%02d-tp.dat" % (float(randx)/4.0, (decndx+LOW))
                fp = open(fn, "w")
                temp = pixels[decndx][randx]/float(pcounts[decndx][randx])
                if (args.continuum == False):
                    temp -= minoverall
                    temp *= 3.0
                    temp += args.tmin
                else:
                    temp /= minoverall
                    temp *= (args.tmin+args.tsys)
                    temp -= args.tsys
                fp.write("%13.2f\n" % temp)
                for q in range(5):
                    v = subbands[decndx][randx][q]/float(pcounts[decndx][randx])
                    if (args.continuum == False):
                        v -= minoverall
                        v *= 3.0
                        v += args.tmin
                    else:
                        v /= minoverall
                        v *= (args.tmin+args.tsys)
                        v -= args.tsys
                    
                    fp.write ("%6.2f " % v)
                fp.write("\n")
                fp.close()

def process_files(args):
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
    print "Processing %s" % nextfn
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
            print "Finished with %s" % nextfn
            linenum = 0
            nextfn = get_next_fn()
            if (nextfn == None):
                break
            else:
                print "Processing %s" % nextfn
                curfile = open(nextfn, "r")
        
        inline=inline.replace("\n", "")
        inlist=inline.split(args.delim)
        if len(inlist) < args.nbins:
            continue
        
        q = args.lmststart
        ra=float(inlist[q])+float(inlist[q+1])/60.0+float(inlist[q+2])/3600.0
        dec=float(inlist[args.decstart])
        
        #
        # Adjust for DEC bias
        #
        if (args.decbias != 0):
            dec += args.decbias

        #
        # Adjust RA if necessary
        # Linear model.  Probably good enough for sloppy antennae
        #
        if (args.raslope != 0):
            #
            # Determine inflection point, IOW, the Zenith angle
            #
            adjust = dec-args.latitude
            ra = ra+(adjust*args.raslope)
            if (ra > 24.0):
                ra = ra - 24.0
            elif (ra < 0.0):
                ra = 24.0 + ra

        rac = ra * 60.0
        rac = rac / 15.0
        rac = round(rac)
        rac = rac * 15.0
        rac = float(rac/60.0)
        randx = rac*4.0
        randx = int(randx)
        if (randx > (24*4)-1):
            randx = (24*4)-1
        
        decndx = int(dec)-int(args.declow)
        
        #
        # Detect bogons, skip
        #
        if (decndx < 0):
            print "%s:%d decndx %d" % (nextfn, linenum, decndx)
            continue
        if (dec < args.declow or dec > args.dechigh or ra < 0.0 or ra > 24.0):
            print "%s:%d ra %f dec %f" % (nextfn, linenum, ra, dec)
            continue
        else:
            if (len(inlist[9:]) == args.nbins):
                failed = False
                try:
                    values=map(float, inlist[args.fftstart:])
                except:
                    failed=True
                if failed == False:
                    values=values[args.ignore:args.nbins-args.ignore]
                    
                    #
                    # Reduce the resolution by a factor of 3
                    #
                    # This makes stuff after this point run faster,
                    #   and there's very little down-side
                    #
                    values = squash(values,3)
                    
                    #
                    # Create the sub-band position map
                    #
                    if (len(submap) == 0):
                        bndx = 0
                        for q in range(5):
                            submap.append(bndx)
                            bndx += len(values)/5
                    
                    #
                    # Convert from dB to linear form
                    #
                    values=numpy.divide(values, [10.0]*len(values))
                    values=numpy.power([10.0]*len(values),values)
                    
                    #
                    # Data is out-of-valid range
                    #
                    whr = numpy.where((values > args.highlimit) | (values < args.lowlimit))
                    if (len(whr[0]) > 0):
                        continue
                    
                    
                    #
                    # De-slope
                    #
                    # First, compute the slope
                    #
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
                    
                    #
                    # Apply the de-sloping transform
                    #
                    values = numpy.subtract(values,desloper)
                    
                    #
                    # Little bit of filtering across the spectrum
                    #
                    values = iir_filter(values,0.2)

                    #
                    # Determine min/max on desloped data
                    #
                    minv = min(values)
   
                    #
                    # Weed out large spikes
                    #
                    
                    
                    #
                    # If "frozen", no need for computing RFI locations
                    #   from input data--we already have enough data.
                    #
                    if (frozen == True):
                        qq = []
                        
                    else:
                        qq = numpy.where(values > (minv*1.80))[0]
                    
                    #
                    # For all locations with spikes
                    #
                    for i in qq:
                        if True:
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
                                    #
                                    # For each "bad spot", zero-out
                                    #  a zone around each spot
                                    #
                                    for ndx in badspots:
                                        for offs in range(7):
                                            masky1[ndx-offs] = 0.0
                                            masky1[ndx+offs] = 0.0
                                #
                                # Compute average of the ends of the histogram
                                #    
                                havg = sum(histo[0:20])
                                havg += sum(histo[-20:])
                                havg /= 40
                                
                                #
                                # We have a "valid" histogram of spikes
                                #
                                if (hmax > 10*havg):
                                    q = numpy.array(histo)
                                    
                                    #
                                    # This only finds one "spot" for now
                                    #
                                    hres = numpy.where(q >= (hmax-1))[0]
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
                    #  In this case, about +/- 2%
                    #
                    # This saves a bit of redundant computation
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
                            for offset in range(7):
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
                    # Do some converisons and de-dishing
                    #    
                    if (True):
                        #
                        # Now divide by min and convert to temperature
                        #  Normalized range will be limited to 1.0 to 2.5, typically
                        #
                        # We then map that through a bit of algebra to an
                        #  antenna-temperature estimate, based on knowledge of TSYS
                        #  and TMIN of the typical sky at 21cm
                        #
                        if (args.continuum == False):
                            values=numpy.divide(values,[minv]*len(values))
                            values=numpy.multiply(values,[args.tsys+args.tmin]*len(values))
                            values=numpy.subtract(values,[args.tsys]*len(values))
                        
                        #
                        # Figure out which de-dishing correction is best
                        #
                        #
                        # Haven't created a corrections array yet? Do so
                        if (len(correction_dict) == 0):
                            if (args.continuum == False):
                                for corr in range(1,30):
                                    correction_dict[corr] = compute_correction(len(values), corr)
                            else:
                                correction_dict[0] = compute_correction(len(values), 1.0)
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
                            if (args.continuum == False):
                                if (corr_sel in correction_dict):
                                    correction = correction_dict[corr_sel]
                                elif corr_sel == 0:
                                    correction = [0.0]*len(values)
                                else:
                                    correction = correction_dict[29]
                            else:
                                correction = correction_dict[0]
                        else:
                            correction = [0.0]*len(values)
                        
                        #
                        # Subtract-out the selected correction value
                        # Hopefully, de-dishing the passband
                        #
                        if (args.continuum == False):
                            values=numpy.subtract(values,correction)
                        else:
                            maxc = numpy.max(correction)
                            if (maxc != 0.0):
                                bigtemp = (high+low)/2.0
                                newc = numpy.multiply(correction,(bigtemp-mintemp))
                                correction = newc
                                values = numpy.subtract(values,newc)
                    
                    #
                    # Might be time to write out debug data
                    #
                    if (int(time.time()) % 15 == 0):
                        fp = open("foo.dat", "w")
                        for v in values:
                            fp.write("%e\n" % v)
                        fp.close()
                        
                        fp = open ("correction.dat", "w")
                        for v in correction:
                            fp.write("%e\n" % v)
                        fp.close()
                        
                    
                    #
                    # At this point, highly unlikely that this test will
                    #   ever fail, but, paranoia
                    #
                    if(ra >= 0.0 and ra <= 24.0):
                        try:
                            badtemp = False
                            #
                            # Convert to temperature, based on method
                            #
                            if (args.continuum == False):
                                temp = to_temp(values, args.method)
                                if (temp < (args.tmin*0.5) or temp > (args.tmin*10.0)):
                                    badtemp = True
                            #
                            # Just compute the total across the spectrum
                            #
                            else:
                                temp = sum(values)
                                if (temp < args.lowlimit*len(values) or temp > args.highlimit*len(values)):
                                    badtemp = True
                            
                            #
                            # Temp not in range? Probably invalid
                            #
                            if (badtemp == False):
                                pixels[decndx][randx] += temp
                                pcounts[decndx][randx] += 1
                            else:
                                if ((random.randint(0,int(1.0e6)) % 500) == 0):
                                    print "File: %s:%d bad temp %e" % (nextfn,linenum, temp)
                                    
                            #
                            # Now deal with the sub-bands
                            #
                            if (badtemp == False):
                                bndx = 0
                                lmap = len(values)/5
                                for q in range(5):
                                    srange = values[bndx:bndx+lmap]
                                    if (args.continuum == False):
                                        pv = to_temp(srange,args.method)
                                    else:
                                        pv = sum(values)
                                        subbands[decndx][randx][q] = pv

                                    bndx += lmap
                        except:
                            print "file %s:%d decndx %d randx %d" % (nextfn, linenum, decndx, randx)
                            break
                    else:
                        print "file %s:%d Ra %f dec %f invalid!" % (nextfn, linenum, ra, dec)

LOW = -35
HIGH = 70

TSYS=120
TMIN=15.0

NBINS=8192
IGNORE=int(NBINS/8)

parser = argparse.ArgumentParser()
parser.add_argument("--method", default="mix", help="which temperature estimator to use")
parser.add_argument("--prefix", default="PROCESSED", help="output filename prefix")
parser.add_argument("--continuum", action="store_true", help="enable continuum mode")
parser.add_argument("--tmin", default=TMIN, type=float, help="estimated sky MIN temperature")
parser.add_argument("--tsys", default=TSYS, type=float, help="estimated TSYS")
parser.add_argument("--lowlimit", default=1.0e-10, type=float, help="lower magnitude limit--data error checking")
parser.add_argument("--highlimit", default=1.0e-5, type=float, help="upper magnitude limit--data error checking")
parser.add_argument("--declow", default=LOW, type=int, help="lower DEC limit")
parser.add_argument("--dechigh", default=HIGH, type=int, help="upper DEC limit")
parser.add_argument("--nbins", default=NBINS, type=int, help="Expected number of FFT bins")
parser.add_argument("--ignore", default=IGNORE, type=int, help="Number of bins to trim off the sides")
parser.add_argument("--raslope", default=0.0, type=float, help="RA slope correction (hours/degree-of-DEC)")
parser.add_argument("--decbias", default=0.0, type=float, help="DEC correction")
parser.add_argument("--fftstart", default=9, type=int, help="offset of start of FFT tokens")
parser.add_argument("--lmststart", default=3, type=int, help="offset of start of LMST tokens")
parser.add_argument("--decstart", default=8, type=int, help="offset of DEC token")
parser.add_argument("--delim", default=",", help="token delimiter in input file")
parser.add_argument("--latitude", default=44.9, type=float, help="local geographic latitude (WGS84)")
parser.add_argument("files", nargs="+")
args = parser.parse_args()

#
# Adjust DEC range for bias
#
args.declow += args.decbias
args.dechigh += args.decbias

#
# Initialize the pixel data arrays
#
rows = int(args.dechigh-args.declow)
rows += 1

cols = 24 * 4

pixels  = [[0.0 for i in range(cols)] for j in range(rows)]
pcounts = [[0 for i in range(cols)] for j in range(rows)]
subbands = [[[0,0,0,0,0] for i in range(cols)] for j in range(rows)]


for i in range(rows):
    for j in range(cols):
        pixels[i][j] = 0.0
        pcounts[i][j] = 0.0
        
filelist = args.files
filelist.reverse()

process_files(args)
dump_pixels(pixels,subbands,pcounts,rows,cols,args)
