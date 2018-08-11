#!/usr/bin/env python

import sys
import math
import os
import time

#T = requested simulation duration
#V = volume in meter cube
#N = number of molecules
#R = voxel radius
#D = diffusion coefficient

def run_single(T, V, N, R, D):
  # V is in m^3, get the cube length, L
  L = math.pow(V, 1.0 / 3.0)
  param = ("%e %e %d %e %e" %(T, L, int(N), R, D))
  print("Now running",int(N),"molecules for",T,"s")
  print param
  start = time.time()
  os.system("./fastbd " + param)
  end = time.time()
  duration = end-start
  print(duration,"s")
  return duration

if __name__ == '__main__':
  T = 0.1
  V = 3e-18
  N = 100
  R = 2.5e-9
  D = 1e-12
  #python run_single.py 1e-5 1e-15 1e+3 2.5e-9 1e-12
  if(len(sys.argv) == 6):
    T = float(sys.argv[1])
    V = float(sys.argv[2])
    N = float(sys.argv[3])
    R = float(sys.argv[4])
    D = float(sys.argv[5])
  run_single(T, V, N, R, D)
