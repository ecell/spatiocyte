try:
  T
except NameError:
  T = 1000
  V1 = 0.1
  Iterations = 10

filename = "spatiocyte_3D_k%.3f.csv" %V1

import math 
sim = theSimulator
s = sim.createStepper('SpatiocyteStepper', 'SS')
s.VoxelRadius = 0.05
s.SearchVacant = 0
#s.DebugLevel = 0

sim.rootSystem.StepperID = 'SS'
sim.createEntity('Variable', 'Variable:/:GEOMETRY').Value = 0
sim.createEntity('Variable', 'Variable:/:LENGTHX').Value = 10
sim.createEntity('Variable', 'Variable:/:LENGTHY').Value = 10
sim.createEntity('Variable', 'Variable:/:LENGTHZ').Value = 10
sim.createEntity('Variable', 'Variable:/:XYPLANE').Value = 1 #periodic
sim.createEntity('Variable', 'Variable:/:XZPLANE').Value = 1 #periodic
sim.createEntity('Variable', 'Variable:/:YZPLANE').Value = 1 #periodic
sim.createEntity('Variable', 'Variable:/:VACANT')
sim.createEntity('Variable', 'Variable:/:A').Value = 100
sim.createEntity('Variable', 'Variable:/:B').Value = 100

l = sim.createEntity('VisualizationLogProcess', 'Process:/:logger')
l.VariableReferenceList = [['_', 'Variable:/:A']]
l.VariableReferenceList = [['_', 'Variable:/:B']]
l.LogInterval = 1

p = sim.createEntity('MoleculePopulateProcess', 'Process:/:p1')
p.VariableReferenceList = [['_', 'Variable:/:A']]

p = sim.createEntity('MoleculePopulateProcess', 'Process:/:p2')
p.VariableReferenceList = [['_', 'Variable:/:B']]

d = sim.createEntity('DiffusionProcess', 'Process:/:diffuseA')
d.VariableReferenceList = [['_', 'Variable:/:A']]
d.D = 0.08

d = sim.createEntity('DiffusionProcess', 'Process:/:diffuseB')
d.VariableReferenceList = [['_', 'Variable:/:B']]
d.D = 0.08

b = sim.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:r1')
b.VariableReferenceList = [['_', 'Variable:/:A','-1']]
b.VariableReferenceList = [['_', 'Variable:/:B','-1']]
b.VariableReferenceList = [['_', 'Variable:/:B','1']]
b.k = V1

l = sim.createEntity('IteratingLogProcess', 'Process:/:iter')
l.VariableReferenceList = [['_', 'Variable:/:A']]
l.LogInterval = 1
l.LogEnd = T
l.Iterations = Iterations
l.FileName = filename
l.Verbose = 0

run(T+0.01)
