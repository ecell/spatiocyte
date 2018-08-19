#!/usr/bin/env python

import os
import time

start = time.time()
os.system("smoldyn reaction_andrews_2018.txt --define DELTA_T=0.001 --define FILE=smoldyn_andrews_2018.csv")
end = time.time()
duration = end-start
print duration,"s"
