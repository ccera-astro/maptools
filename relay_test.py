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
import os
import struct
import time
import math
import sys
from pylibftdi import BitBangDevice
import signal

bbd = None

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

#
# Create an FTDI BigBangDevice instance, based on serial number passed in
#
bbd = BitBangDevice(device_index=0)

#
# Set direction regiser
#
bbd.direction = 0xFF
bbd.port = 0


DOWN=0b00001101
UP=0b00001110

while True:
    bbd.port = DOWN
    time.sleep(1)
    bbd.port = UP
    time.sleep(1)

# We really want that motor to stop
#
bbd.port = 0
if (timeout == True):
    sys.exit(3)
else:
    sys.exit(0)
    
