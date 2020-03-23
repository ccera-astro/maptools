#!/usr/bin/env python
import glob
import sys
import ephem
import math
import argparse
import numpy as np

def avg_rgb(r1,r2):
    ret = []
    for i in range(len(r1)):
        ret.append((r1[i]+r2[i])/2.0)
    return ret

def rgb_iir(new,old,a):
    b = 1.0 - a
    ret = []
    for i in range(len(new)):
        ret.append(a*new[i] + b*old[i])
    return ret

def rgb_scale(r,scale):
       visual_correction = [1.1,0.8,1.1]
       ret = []
       for i in range(len(r)):
           v = r[i]/scale[i]
           v *= visual_correction[i]
           ret.append(v)
       return ret

def rgb_make(rv):
    expand = 3
    inter_len = len(rv)*3
    inter_len *= expand
    inter_buf = []
    out_buf = []
    #
    # Expand
    #
    for i in range(len(rv)):
        for j in range(3*expand):
            inter_buf.append(rv[i])
    if (len(inter_buf) != inter_len):
        raise ValueError("Problem with length of inter_buf")
    
    #
    # Smooth inter_buf
    #
    v = inter_buf[0]
    a = 0.5
    b = 1.0 - a
    for i in range(len(inter_buf)):
        v = inter_buf[i]*a + v*b
        inter_buf[i] = v
    
    #
    # Now decimate it back down
    #
    prcnt = len(inter_buf)/7
    prcnt = int(prcnt)
    inter_buf=inter_buf[prcnt:-prcnt]
    for i in range(3):
        lv = len(inter_buf)/3
        v = sum(inter_buf[(i*lv):(i+1)*lv])
        v /= lv
        out_buf.append(v)
    return (out_buf)
     
def gal_index(ra,dec):
    p = ephem.Equatorial(str(ra), str(dec))
    galactic = ephem.Galactic(p)
    glong = math.degrees(galactic.lon)
    glat = math.degrees(galactic.lat)

    glndx = int(glong)+180
    glndx %= 360
    glndx /= args.increment
    glndx -= 1

    gladx = int(glat)+90
    gladx /= args.increment
    gladx -= 1
    return ((glndx,gladx))

parser = argparse.ArgumentParser()

parser.add_argument("--prefix", default="PROCESSED-", help="output filename prefix")
parser.add_argument("--declow", default=-35, type=int, help="lower DEC limit")
parser.add_argument("--dechigh", default=70, type=int, help="upper DEC limit")
parser.add_argument("--resolution", default=125, type=int, help="Output resolution, ppi")
parser.add_argument("--increment", default=5, type=int, help="Sky pixel bin size")
parser.add_argument("--oprefix", default="21cm", help="Output file prefix")

args = parser.parse_args()

#
# Input data are stuffed as they arrive into this array
#
fourples=[]

#
# For tracking "missing" pixels
#
counting = [-1]*((args.dechigh-args.declow)+1)
tracking = [[]]*len(counting)

#
# RGB values are taken from the first four subanded values in the input data
#
rgb_map=[(1,1),(2,2),(3,3)]

#
# First, loop through all the files described by args.prefix
#
# Extract simple-pixel value (first line)
# Extract subbanded pixel value (second line)
#
# Filename gives coordinates:
#
#  <PREFIX>RR.RR-DD
#
#  Where:
#     RR.RR == RA in decimal hours
#     DD = DEC in decimal degrees, with possible negative sign
#
prelen = len(args.prefix)

#
# For scaling RGB appropriately
#
tmax = 0.0
rmax = 0.0
gmax = 0.0
bmax = 0.0

rcount = 0
gcount = 0
bcount = 0

#
# Glob on the file prefix, for each file found, process it
#  into the "fourples" array.
#
for fn in glob.glob(args.prefix+"*-tp.dat"):
    
    #
    # First extract RA value from filename
    #
    ra = fn[prelen:prelen+5]
    ra = float(ra)
    
    #
    # Now the DEC
    #
    dec = fn[prelen+6:]
    dec = dec.replace("-tp.dat", "")
    dec = int(dec)
    
    #
    # Finally open the file, read two lines out of it
    #
    f = open(fn, "r")
    
    #
    # Simple scalar total-power pixel value is first line
    #
    val = f.readline().strip("\n")
    val = float(val)
    
    #
    # Second line has sub-banded values that will form RGB
    #
    rgb_toks = f.readline().strip("\n").split()
    if (len(rgb_toks) == 5):
        rgb_vals = []
        for t in rgb_toks:
            rgb_vals.append(float(t))
    else:
        print "Hmmmm, RGB unexpected tokens %s" % (str(rgb_toks))
    f.close()
    
    #
    # Sanity check input
    #
    if (dec >= args.declow and dec <= args.dechigh):
        #
        # Form RGB values from input sub-bands, based on rgb_map
        #
        if (len(rgb_toks) == 5):
            rgb = rgb_make(rgb_vals)
            if (rgb[0] > rmax):
                rmax = rgb[0]
            if (rgb[1] > gmax):
                gmax = rgb[1]
            if (rgb[2] > bmax):
                bmax = rgb[2]
            
            if (rgb[0] > 35):
                rcount += 1
            if (rgb[1] > 35):
                gcount += 1
            if (rgb[2] > 35):
                bcount += 1
        else:
            rgb=[0.0,0.0,0.0]
        #
        # Stuff values into fourples array
        #
        fourples.append([ra, dec, val,rgb])
        
        #
        # Track missing pixels
        #
        counting[dec-args.declow] += 1
        tracking[dec-args.declow].append(ra)

#
# Now our job is to take the fourples data and place it into
#  The appropriate places in the pixels and galactic_pixels array
#
rows = (args.dechigh-args.declow)/args.increment
rows +=1
rows *= 2
cols = 24 * 4

pixels  = [[0.0 for i in range(cols)] for j in range(rows)]
rgb_pixels = [[[0.0,0.0,0.0] for i in range(cols)] for j in range(rows)]

glongs = 360/args.increment
glats = 180/args.increment

galactic_pixels = [[0.0 for i in range(glongs)] for j in range(glats)]
galactic_rgb_pixels = [[[0.0,0.0,0.0] for i in range(glongs)] for j in range(glats)]


for t in fourples:
    randx = t[0]
    randx *= 4.0
    randx = int(randx)
    randx %= cols

    decndx = float(t[1])
    decndx += float(-args.declow)
    decndx /= float(args.increment)
    decndx = int(decndx)
    decndx *= 2

    val = t[2]

    rgb_pixels[decndx][randx] = rgb_scale(t[3],[rmax,gmax,bmax])
    
    if (pixels[decndx][randx] == 0.0):
        pixels[decndx][randx] = val
    else:
        cv = pixels[decndx][randx]
        cr = val/cv
        #
        # Sanity check--if two values getting plonked in differ by more than
        #  +/- 10%, something is wrong
        #
        if (0.90 <= cr and cr <= 1.10):
            pixels[decndx][randx] += val
            pixels[decndx][randx] /= 2.0
        else:
            print "Hmmmm: cv %f val %f" % (cv, val)

    glndx,gladx = gal_index(t[0],t[1])
    galactic_rgb_pixels[gladx][glndx] = rgb_scale(t[3],[rmax,gmax,bmax])
    
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
# This applies to equatorial maps only (normal and RGB)
#
for i in range(rows-1):
    if ((i % 2) != 0):
        for j in range(cols):
            pixels[i][j] = (pixels[i-1][j] + pixels[i+1][j])/2.0
            rgb_pixels[i][j] = avg_rgb(rgb_pixels[i-1][j],rgb_pixels[i+1][j])

for j in range(cols):
    pixels[rows-1][j] = pixels[rows-2][j] + pixels[rows-3][j]
    pixels[rows-1][j] /= 2.0
    rgb_pixels[rows-1][j] = avg_rgb(rgb_pixels[rows-2][j],rgb_pixels[rows-3][j])

#
# Before we smooth the galactic maps, find dead zones, and make a note...
#
# We cover from -90 to args.declow and args.dechigh to +90
#
# This assumes complete coverage of the resulting "live" zone
#
deadzones = []
for ra in np.arange(0,24,0.1):
    
    #
    # Lower range
    #
    for dec in range(-90,args.declow):
        ndx = gal_index(ra,dec)
        if (ndx not in deadzones):
            deadzones.append(gal_index(ra,dec))
    #
    # Upper range
    #
    for dec in range(args.dechigh+1,90):
        ndx = gal_index(ra,dec)
        if (ndx not in deadzones):
            deadzones.append(gal_index(ra,dec))
        
#
# Smoothing for galactic map
#
a = 0.4
b = 1.0-a
for i in range(glats):
    v = galactic_pixels[i][0]
    ov = galactic_rgb_pixels[i][0]
    for j in range(glongs):
        cv = galactic_pixels[i][j]
        v = (a)*cv + (b)*v
        galactic_pixels[i][j] = v

        ov = rgb_iir(galactic_rgb_pixels[i][j], ov, a)
        galactic_rgb_pixels[i][j] = ov

for i in range(glongs):
    v = galactic_pixels[0][i]
    ov = galactic_rgb_pixels[0][i]
    for j in range(glats):
        cv = galactic_pixels[j][i]
        v = (a)*cv + (b)*v
        galactic_pixels[j][i] = v
        ov = rgb_iir(galactic_rgb_pixels[j][i], ov, a)
        galactic_rgb_pixels[j][i] = ov


if True:
    import random
    for c in deadzones:
		
		#
		# For total-power, just set it to a very-low temperature
		#
        galactic_pixels[c[1]][c[0]] = 3.0

        #
        # For RGB, make it dark gray
        #
        gray = 30.0/255.0
        galactic_rgb_pixels[c[1]][c[0]] = [gray,gray,gray]
    
#
# Smoothing for RGB map
#
if True:
    a = 0.4
    b = 1.0-a
    for i in range(rows-1):
        v = rgb_pixels[i][0]
        for j in range(cols):
            cv = rgb_pixels[i][j]
            v = rgb_iir(cv,v,a)
            rgb_pixels[i][j] = v

    for i in range(cols):
        v = rgb_pixels[0][i]
        for j in range(rows):
            cv = rgb_pixels[j][i]
            v = rgb_iir(cv,v,a)
            rgb_pixels[j][i] = v

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
        v = (a*cv) + (b*v)
        pixels[i][j] = v

#
# Then ra-wise
#
for i in range(cols):
    v = pixels[0][i]
    for j in range(rows):
        cv = pixels[j][i]
        v = (a*cv) + (b*v)
        pixels[j][i] = v

if True:
    import matplotlib.pyplot as plt


    plt.figure(1)
    rpixels = pixels[::-1]
    plt.xticks(list(range(0,cols,8)),list(range(0,24,2)))
    plt.yticks(list(range(0,rows,4)),list(np.arange(args.dechigh,args.declow-2,-10)))
    plt.xlabel("Hours RA")
    plt.ylabel("DEC")
    plt.title("H1 Brightness:  +/- 160km/sec")
    plt.imshow(rpixels,cmap='jet', interpolation="gaussian")
    plt.colorbar(shrink=0.7,orientation="horizontal").set_label("Est. Brightness (K)")
    plt.savefig(args.oprefix+".png", dpi=args.resolution)

    plt.figure(2)
    rpixels = galactic_pixels[::-1]
    xtlist = list(range(0,180,40))
    xtlist = list(range(180,360,40))+xtlist
    plt.xticks(list(range(0,360/args.increment,8)),xtlist)
    plt.yticks(list(range(0,180/args.increment,4)),list(np.arange(90,-90,-20)))
    plt.xlabel("Galactic Longitude")
    plt.ylabel("Galactic Latitude")
    plt.title("H1 Brightness:  +/- 160km/sec")
    plt.imshow(rpixels,cmap='jet', interpolation="gaussian")
    plt.colorbar(shrink=0.7,orientation="horizontal").set_label("Est. Brightness (K)")
    plt.savefig(args.oprefix+"-galactic.png", dpi=args.resolution)
    
    plt.figure(3)
    rpixels = rgb_pixels[::-1]
    plt.xticks(list(range(0,cols,8)),list(range(0,24,2)))
    plt.yticks(list(range(0,rows,4)),list(np.arange(args.dechigh,args.declow,-10)))
    plt.xlabel("Hours RA")
    plt.ylabel("DEC")
    plt.title("H1 Red/Blueshift")
    plt.imshow(rpixels,cmap='jet', interpolation="gaussian")
    if False:
        cbar = plt.colorbar(shrink=0.7,orientation="horizontal")
        cbar.set_label("Redshift km/sec")
        numticks = np.arange(0.0,1.0,0.1)
        tl = 280/(len(numticks)-2)
        tks = range(-140,140+1,tl)
        cbar.set_ticks(np.arange(0.0,1.0,0.1))
        cbar.ax.set_xticklabels(tks)
    plt.savefig(args.oprefix+"-rgb.png", dpi=args.resolution)
    

    plt.figure(4)
    rpixels = galactic_rgb_pixels[::-1]
    xtlist = list(range(0,180,40))
    xtlist = list(range(180,360,40))+xtlist
    plt.xticks(list(range(0,360/args.increment,8)),xtlist)
    plt.yticks(list(range(0,180/args.increment,4)),list(np.arange(90,-90,-20)))
    plt.xlabel("Galactic Longitude")
    plt.ylabel("Galactic Latitude")
    plt.title("H1 Red/Blueshift")
    plt.imshow(rpixels,cmap='jet', interpolation="gaussian")
    #plt.colorbar(shrink=0.7,orientation="horizontal").set_label("Est. Brightness (K)")
    plt.savefig(args.oprefix+"-galactic-rgb.png", dpi=args.resolution)
    print "rmax %f gmax %f bmax %f" % (rmax, gmax, bmax)
    print "rcount %d gcount %d bcount %d" % (rcount, gcount, bcount)

f = open("missing-pixels.log", "w")
for i in range(0,len(counting)):
    if (counting[i] != cols and counting[i] != -1):
        f.write("%d: %d\n" % (i+args.declow, cols-counting[i]))
f.close()
