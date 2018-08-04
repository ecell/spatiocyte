#!/usr/bin/env python

import sys
import math
import os
import pickle
import random

#T = requested simulation duration
#V = volume in liters
#N = number of molecules
#R = voxel radius
#D = diffusion coefficient

def run_single(T, V, N, R, D):
  print('T (sim time) =', T, '; V (liter volume) = ', V, '; N (#molecules) =', N, '; R (voxel radius) =', R, '; D (diffusion coef) =', D)
  filename = ('%f' %random.random())
  # V is in m^3, get the cube length, L
  L = math.pow(V, 1.0 / 3.0)
  param = ("--parameters=\"{'T':%e, 'L':%e, 'N':%e, 'R':%e, 'D':%e, 'filename':'%s'}\"" %(T, L, N, R, D, filename))
  os.system("ecell3-session " + param + " diffusion.py")
  with open(filename, 'rb') as f:
      timing = pickle.load(f)
      f.close()
  os.remove(filename)
  print('total runtime:', timing, '\n')
  return timing

if __name__ == '__main__':
  T = 1e-5
  V = 4e-18
  N = 1e+3
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
