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
# It also drives a linear-actuator motor at a fixed rate either UP or
#   DOWN to achieve a given dish position.
#
#
import hid
import os
import struct
import time
import math
import sys
from pylibftdi import BitBangDevice
import signal

bbd = None

def readangle(devs,dind):
        
    #
    # Make up a "give me the angle" command
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
    ang *= -1.0
    
    return(ang)

def sighandle(num, stack):
    global bbd
    
    bbd.port = 0
    print "Resetting relay board"
    time.sleep(0.2)
    print "Exiting"
    sys.exit()

signal.signal(signal.SIGINT, sighandle)
signal.signal(signal.SIGHUP, sighandle)
signal.signal(signal.SIGQUIT, sighandle)
signal.signal(signal.SIGSEGV, sighandle)
signal.signal(signal.SIGBUS, sighandle)

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
# Gather the serial numbers of each sensor device we found
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

bbd = None
#
# Create an FTDI BigBangDevice instance, based on serial number passed in
#
bbd = BitBangDevice(device_index=0)

#
# Set direction regiser
#
bbd.direction = 0xFF
bbd.port = 0

#
# For the way we've wired up the linear actuator motor, this should
#   produce an "UP" and "DOWN" motor direction.
#
# Two of the relay channels--the low-order bits just control the +/-
#   that go to the other two channels.  Those channels are wired in
#   a kind of exclusive OR, so that only bits with opposite settings
#   produce any motion at all.
#
DOWN=0b00001101
UP=0b00001110

#
# Desired position is given on the command line
#
position=float(sys.argv[1])

then=time.time()

motor_started = False
done = False
movetimer = 120
timeout = False
oldport = 0
while done == False:
    
    #
    # For each device in the list
    #
    for dind in range(0,len(devs)):
        ang = readangle(devs, dind)
        #
        # Initialize the single-pole IIR filter with first value
        #
        if (avgangs[dind] < -1000):
            avgangs[dind] = ang
            
            #
            # Also initialize correct move-timer estimate
            #
            # We move at about 0.5deg/sec
            #
            movetimer = abs(ang-position)*2.0
            movetimer *= 1.2
        
        avgangs[dind] = (ang*alpha) + (avgangs[dind]*beta)
        lastang = avgangs[0]

        if (motor_started == False):
            if (avgangs[dind] < position):
                bbd.port = UP
            elif (avgangs[dind] > position):
                bbd.port = DOWN
            motor_started = True
        
        else:
            #
            # There will be some amount of overshoot, so we
            #  true to compensate for this by stopping early.
            #
            if abs(position-avgangs[dind]) <= 1.125:
                oldport = bbd.port
                bbd.port = 0
                print "Achieved position: %f" % position
                done = True
        
        #
        # Motion shouldn't take more than movetimer
        #
        if ((time.time() - then) > movetimer):
            bbd.port = 0
            print "Motion timed out"
            timeout = True
            done = True

        time.sleep(0.05)


#
# Let angular momentum damp-out
#
time.sleep(2.0)

#
# Get a fresh reading
#
ang = readangle(devs,0)
time.sleep (0.10)
ang += readangle(devs,0)
time.sleep(0.10)
ang += readangle(devs,0)
ang /= 3.0

#
# If we're out by more than 0.5 degree
# Then run the motor for just a wee bit
#  in the opposite direction.
#
if (abs(ang-position) > 0.5):
	if (ang-position) < 0:
		bbd.port = DOWN
	else:
		bbd.port = UP
	time.sleep(0.5)
	bbd.port = 0

#
# We really want that motor to stop
#
bbd.port = 0
if (timeout == True):
    sys.exit(3)
else:
    sys.exit(0)
    
