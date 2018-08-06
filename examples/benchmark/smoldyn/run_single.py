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

def run_single(T, V, N, R, D, M):
  # V is in m^3, get the cube length, L
  L = math.pow(V, 1.0 / 3.0)
  DELTA_T = 2.0*R*R/(3.0*D)
  param = ("--define T=%e --define L=%e --define N=%d --define D=%e --define R=%e --define DELTA_T=%e" %(T, L, int(N), D, R*2, DELTA_T))
  print("Now running",int(N),"molecules for",T,"s")
  start = time.time()
  os.system("smoldyn " + M + " " + param)
  end = time.time()
  duration = end-start
  print(duration,"s")
  return duration

if __name__ == '__main__':
  T = 1e-5
  V = 1e-18
  N = 1e+3
  R = 2.5e-9
  D = 1e-12
  M = "diffusion_excluded_volume.txt"
  #python run_single.py 1e-5 1e-15 1e+3 2.5e-9 1e-12 diffusion.txt
  if(len(sys.argv) == 6):
    T = float(sys.argv[1])
    V = float(sys.argv[2])
    N = float(sys.argv[3])
    R = float(sys.argv[4])
    D = float(sys.argv[5])
    M = str(sys.argv[6])
  run_single(T, V, N, R, D, M)
