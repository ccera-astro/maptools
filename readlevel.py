#!/usr/bin/env python2

#
# This simple program reads one or more elevation sensors, made by
#
# Level Developments, Inc, UK
#
# They manifest as HID devices, we use the "hidapi" Python library to
#  interact with the devices.
#
#
import hid
import os
import struct
import time
import math
import sys

LVL_VID=0x0461
LVL_PID=0x0021

#
# Enumerate all the HID devices with the Level Developments VID/PID
#
# Open any devices found
#
dlist = hid.enumerate(LVL_VID, LVL_PID)


if len(dlist) <= 0:
	raise OSError("Unable to locate any elevation sensors")

devs = []
for dv in dlist:
    d = hid.device()
    d.open_path(dv["path"])
    devs.append(d)

buf=bytearray(64)

#
# Gather the serial numbers of each device we found
#
sernums = []
for d in devs:
	
	#
	# Form a "tell me your serial number" command
	#
    buf[0] = 0x00
    buf[1] = 0x00
    buf[2] = 0x01
    
    b = bytearray(4)
    
    #
    # Write that command
    #
    d.write(buf)
    
    #
    # Read the response from the device
    #
    x = d.read(6)
    for i in range(0,len(b)):
        b[i] = x[i+2]
    
    #
    # Unpack the binary
    #
    ser = struct.unpack(">i", b)
    
    #
    # Put it in sernum array
    #
    sernums.append(ser[0])
    
alpha=0.1
beta=1.0-alpha
avgangs=[-90000.0]*len(devs)
wait=10
count = 0
done=False
current_rate=None
while done==False:
	
	#
	# For each device in the list
	#
    for dind in range(0,len(devs)):
		
		#
		# Make up a "give met the angle" command
		#
        buf[0] = 0x00
        buf[1] = 0x00
        buf[2] = 0x05
        
        b = bytearray(4)
        
        #
        # Write the command
        #
        devs[dind].write(buf)
        
        #
        # Read back the response
        #
        x = devs[dind].read(6)
        
        #
        # Load it into a bytearray
        #
        for i in range(0,len(b)):
            b[i] = x[i+2]
        
        #
        # Unpack that bytearray into an integer
        #
        angle = struct.unpack(">i", b)
        
        #
        # Convert into floating-point angle estimate
        #
        ang = float(angle[0])
        ang = ang/1000.0
        
        #
        # Initialize the single-pole IIR filter with first value
        #
        if (avgangs[dind] < -1000):
            avgangs[dind] = ang
        
        avgangs[dind] = (ang*alpha) + (avgangs[dind]*beta)
        lastang = avgangs[0]
        
        print "Serial: %d angle %f" % (sernums[dind], avgangs[dind])
        time.sleep(0.2)
