#!/usr/bin/env python

import os
import time

start = time.time()
os.system("smoldyn reaction.txt --define DELTA_T=6.67e-5 --define FILE=smoldyn_small_dt.csv")
end = time.time()
duration = end-start
print duration,"s"
