import time
import pickle
import gc

try:
  T
except NameError:
  T = 1e-1
  L = 1e-6
  N = 1e+4
  R = 2.5e-9
  D = 1e-12
  filename = 'tempfile.txt'

stirTime = T*0.05

print "parameters=\"{'T':%0.6e, 'L':%0.6e, 'N':%0.6e, 'R':%0.6e, 'D':%0.6e, 'filename':%s}\"" %(T, L, N, R, D, filename)
sim = theSimulator.createStepper('SpatiocyteStepper', 'SS')
sim.VoxelRadius = R
theSimulator.rootSystem.StepperID = 'SS'
theSimulator.createEntity('Variable', 'Variable:/:GEOMETRY').Value = 0
theSimulator.createEntity('Variable', 'Variable:/:LENGTHX').Value = L
theSimulator.createEntity('Variable', 'Variable:/:LENGTHY').Value = L
theSimulator.createEntity('Variable', 'Variable:/:LENGTHZ').Value = L
#periodic boundaries
theSimulator.createEntity('Variable', 'Variable:/:XYPLANE').Value = 1
theSimulator.createEntity('Variable', 'Variable:/:XZPLANE').Value = 1
theSimulator.createEntity('Variable', 'Variable:/:YZPLANE').Value = 1

theSimulator.createEntity('Variable', 'Variable:/:VACANT')
theSimulator.createEntity('Variable', 'Variable:/:A').Value = N
#log = theSimulator.createEntity('VisualizationLogProcess', 'Process:/:log')
#log.VariableReferenceList = [['_', 'Variable:.:A']] 
pop = theSimulator.createEntity('MoleculePopulateProcess', 'Process:/:pop')
pop.VariableReferenceList = [['_', 'Variable:.:A']] 
dif = theSimulator.createEntity('DiffusionProcess', 'Process:/:diff')
dif.VariableReferenceList = [['_', 'Variable:.:A']]
dif.D = D
dif.WalkPoint = 0
gc.disable
run(stirTime)
print "Now running",int(N),"molecules for",T,"s"
start = time.time()
run(T)
end = time.time()
gc.collect()
gc.enable()
duration = end-start
f = open(filename, 'w')
pickle.dump(duration, f)
#print duration
f.close()

