#!/usr/bin/env python

import os, sys
import time

t = TIME = None

while 1:



    try:
        flds = os.stat(sys.argv[1])
    
        t = flds.st_mtime
        if TIME is None:
            TIME = t
    except OSError:
        pass

    if t > TIME:
        rval = os.system("make")
        assert rval == 0
        TIME = t

    time.sleep(1)



