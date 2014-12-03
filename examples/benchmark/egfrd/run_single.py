import _gfrd
import model
import gfrdbase
import numpy
import myrandom
import time
import sys
import gc

def run_single(T, V, N, R, D):
  L = numpy.power(V, 1.0/3.0) # cuboid side length
  # matrix_size = 3 # min matrix size
  matrix_size = max(3, int(min(numpy.power(N, 1.0 / 3.0), L / (2 * R)))) 
  m = model.ParticleModel(L)
  A = model.Species('A', D, R)
  m.add_species_type(A)
  m.set_all_repulsive() 
  w = gfrdbase.create_world(m, matrix_size)
  gfrdbase.throw_in_particles(w, A, N)
  nrw = gfrdbase.create_network_rules_wrapper(m)
  sim = _gfrd._EGFRDSimulator(w, nrw, myrandom.rng) 

  stirTime = T*0.01
  t = 0.0
  gc.disable
  while (sim.step(stirTime)):
      pass
  print "Now running",int(N),"molecules for",T,"s"
  endTime = stirTime+T
  start = time.time()
  while (sim.step(endTime)):
      pass
  print "time",sim.t,endTime
  end = time.time()
  gc.collect()
  gc.enable()
  duration = end-start
  print duration,"s"
  return  duration

if __name__ == '__main__':
  myrandom.seed(0)
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

