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

NBINS=8192
observation=[0.0]*NBINS
obscount=0
REDRANGE=170
SPUR=84

TPCHUNKS=5
TSYS=85.0
TMIN=20.0

Fc=1420.40575e6
C=299792
bandwidth=2.5e6
binwidth=bandwidth/NBINS

start_red= -(binwidth/Fc)*C
start_red *= (NBINS/2)
red_incr = (binwidth/Fc)*C

desired_dec=float(sys.argv[1])
ratoks=sys.argv[2].split(":")
desired_ra=float(ratoks[0])+float(ratoks[1])/60.0
prefix=sys.argv[3]

#
# A place to hold an L=PMEDIAN filter
#
PMEDIAN=11
current_filter = [-100.0]*PMEDIAN
def median(p):
    global current_filter

    if (current_filter[0] < -50.0):
        current_filter = [p]*PMEDIAN
    
    #
    # Do the shift
    #
    current_filter = [p]+current_filter[0:PMEDIAN-1]
    
    #
    # Median filter
    #
    return(numpy.median(current_filter))


linenum=0
skipcount = 10
while True:
    inline=sys.stdin.readline()
    linenum += 1
    if (inline == ""):
        break
    
    inline=inline.replace("\n", "")
    inlist=inline.split(",")
    
    ra=float(inlist[3])+float(inlist[4])/60.0+float(inlist[5])/3600.0
    dec=float(inlist[8])
    
    #
    # Deal with bug in data writer
    #
    if (abs(dec-desired_dec) <= 0.5):
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
                values=numpy.divide(values, [10.0]*len(values))
                values=numpy.power([10.0]*len(values),values)
                if(abs(ra-desired_ra) <= 0.120):
                    # Update observation
                    observation=numpy.add(observation, values)
                    obscount += 1

obscount = float(obscount)

if (obscount > 1):

    difference = observation
    difference = difference.tolist()
    difference.reverse()
    
    #
    # We just determine the slope between the limits
    #
    red = start_red
    for v in difference:
        if (isalmost(red,-REDRANGE,red_incr)):
            slope_begin = v
        elif (isalmost(red, REDRANGE,red_incr)):
            slope_end = v
        red += red_incr

    #
    # We have the values at the ends, compute slope
    #
    slope = slope_end-slope_begin
    slope /= (float(REDRANGE) * 2.0)
    slope *= -1.0
    
    #
    # Slope compensation, smoothing
    #
    adjust = 0.0
    red = start_red
    pv = -100.0
    desloped = []
    for v in difference:
        #
        # We only care about data in {-REDRANGE,REDRANGE}
        #
        if (red >= -REDRANGE and red <= REDRANGE):
            
            #
            # There's a spur in the data right around this redshift--make it go away
            #
            if (red >= SPUR-4 and red <= SPUR+4):
                val = random.uniform(0.99*goodval,1.01*goodval)
            else:
                val = v
                goodval = val
            
            if (pv <= -100.0):
                pv = val
            
            #
            # Smooth
            #
            pv = 0.2*val + 0.8*pv
            
            #
            # Deslope
            #
            sv = pv + ((red-start_red)*slope)
            
            #
            # Create a buffer of desloped values
            #
            desloped.append(sv)

        red += red_incr
    
    #
    # Determine the minimum value over the range of redshift we care about
    #  (from the desloped values)
    #
    obsmin = 1.0e15
    for v in desloped:
        if (v < obsmin):
            obsmin = v
        
    #
    # Compute a wee flattening curve
    #
    correction = []
    slen = len(desloped)
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

    #
    # Now, we output the values, somewhat smoothed
    #
    red = -REDRANGE
    sval = desloped[0]
    cndx = 0
    tp = 0.0
    nval = len(desloped)
    vincr = nval / TPCHUNKS
    vincr = int(vincr)
    tpchunks = [0.0]*(TPCHUNKS+1)
    chndx = 0
    for v in desloped:
        sval = median(v)# *0.3 + 0.7*sval
        dval = (sval/obsmin)
        dval *= (TMIN+TSYS)
        dval -= TSYS
        
        #
        # Subtract-out the correction curve
        #
        dval -= correction[cndx]
        
        print "%f %8.4f" % (red, dval)
        tp += dval
        tpchunks[chndx] += dval
        red += red_incr
        cndx += 1
        if ((cndx % vincr) == 0):
			chndx += 1
        

	rac = desired_ra * 60.0
	rac = rac / 15.0
	rac = round(rac)
	rac = rac * 15.0
	rac = float(rac/60.0)
    fn = prefix+"-%05.2f-%02d-tp.dat" % (rac, desired_dec)
    v = 0.0
    vcnt = 0
    try:
        f = open(fn, "r")
        vls = f.readlines()
        f.close()
        v = vls[0].strip("\n")
        v = float(v)
        vcnt = 1
    except:
        pass
    v += tp
    vcnt += 1
    v /= vcnt
    sun = ephem.Sun()
    sun.compute()
    beam = ephem.Equatorial(str(rac), str(desired_dec))
    sep = math.degrees(ephem.separation((sun.ra,sun.dec),(beam.ra,beam.dec)))
    #
    # Avoid applying update when Sun is in the beam
    #  At our beam-size, it's the only object that will
    #  really cause any problem.
    #
    if (sep > 11.0):
		f = open(fn, "w")
		f.write("%9.2f\n" % v)
		for tpc in tpchunks:
			f.write("%9.2f " % tpc)
		f.write("\n")
		f.close()
    
        
