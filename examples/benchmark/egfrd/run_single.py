import numpy
import time
import sys
import gc
from ecell4 import *

def run_single(T, V, N, R, D):
  L = numpy.power(V, 1.0/3.0) # cuboid side length
  with species_attributes():
    A | {'D': D, 'radius': R}
  m = get_model() 
  M = max(int(min(cbrt(N + N), L / (2*R))), 3)
  w = egfrd.EGFRDWorld(Real3(L, L, L), Integer3(M, M, M))
  w.bind_to(m)
  w.add_molecules(Species('A'), N)
  sim = egfrd.EGFRDSimulator(w)
  stirTime = T*0.01
  t = 0.0
  gc.disable
  sim.run(stirTime)
  print("Now running",int(N),"molecules for",T,"s")
  endTime = stirTime+T
  start = time.time()
  sim.run(T)
  end = time.time()
  gc.collect()
  gc.enable()
  duration = end-start
  print("time taken:",duration,"s")
  return duration

if __name__ == '__main__':
  T = 1e-2
  V = 3e-15
  N = 1e+4
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

