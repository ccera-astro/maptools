import os
import time
import sys
days=int(sys.argv[1])
t=time.time()-(days*(24*3600))
ltp=time.gmtime(t)

print "%04d%02d%02d" % (ltp.tm_year, ltp.tm_mon, ltp.tm_mday)
