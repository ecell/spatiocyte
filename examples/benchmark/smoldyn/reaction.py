#!/usr/bin/env python

import os
import time

start = time.time()
os.system("smoldyn reaction.txt --define DELTA_T=0.001 --define FILE=smoldyn.csv")
end = time.time()
duration = end-start
print duration,"s"
