import _gfrd
import model
import gfrdbase
import numpy
import myrandom
import time
import sys
import gc

T = 1e+2
V = 90.9e-18
R = 10e-9
D = 1e-12
L = numpy.power(V, 1.0/3.0) # cuboid side length
nE = 91
nS = 909
# matrix_size = 3 # min matrix size
matrix_size = max(3, int(min(numpy.power(nE+nS, 1.0 / 3.0), L / (2 * R)))) 
m = model.ParticleModel(L)

E = model.Species('E', D, R)
S = model.Species('S', D, R)
ES = model.Species('ES', D, R)
P = model.Species('P', D, R)
m.add_species_type(E)
m.add_species_type(S)
m.add_species_type(ES)
m.add_species_type(P)

fwd = model.create_binding_reaction_rule(E, S, ES, 0.01e-18)
back = model.create_unbinding_reaction_rule(ES, E, S, 0.1)
prod = model.create_unbinding_reaction_rule(ES, E, P, 0.1)
m.network_rules.add_reaction_rule(fwd)
m.network_rules.add_reaction_rule(back)
m.network_rules.add_reaction_rule(prod)
m.set_all_repulsive() 

w = gfrdbase.create_world(m, matrix_size)
gfrdbase.throw_in_particles(w, E, nE)
gfrdbase.throw_in_particles(w, S, nS)

nrw = gfrdbase.create_network_rules_wrapper(m)
sim = _gfrd._EGFRDSimulator(w, nrw, myrandom.rng) 

start = time.time()
t = 0
f = open('egfrd.csv','w')
while (sim.step(T)):
    t = t+1e-2
    while (sim.step(t)): pass
    f.write(str(sim.t)+','+str(len(w.get_particle_ids(E)))+','+str(len(w.get_particle_ids(S)))+','+str(len(w.get_particle_ids(ES)))+','+str(len(w.get_particle_ids(P)))+'\n')
    f.flush()
print "time",sim.t,T
end = time.time()
duration = end-start
print duration,"s"
