#!/usr/bin/env python

import os, sys
import time

TIME = None

while 1:



    flds = os.stat(sys.argv[1])

    t = flds.st_mtime
    if TIME is None:
        TIME = t

    if t > TIME:
        rval = os.system("make")
        assert rval == 0
        TIME = t

    time.sleep(1)



