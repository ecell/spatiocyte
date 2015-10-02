try:
  T
except NameError:
  T = 100
  V1 = 0.1
  Iterations = 1000

filename = "spatiocyte_1D_p%.3f.csv" %V1

import math 
sim = theSimulator
s = sim.createStepper('SpatiocyteStepper', 'SS')
s.VoxelRadius = 10e-9 
s.SearchVacant = 0

sim.rootSystem.StepperID = 'SS'
sim.createEntity('Variable', 'Variable:/:GEOMETRY').Value = 6
sim.createEntity('Variable', 'Variable:/:LENGTHX').Value = 8e-6
sim.createEntity('Variable', 'Variable:/:LENGTHY').Value = 0.25e-6
sim.createEntity('Variable', 'Variable:/:YZPLANE').Value = 1
sim.createEntity('Variable', 'Variable:/:VACANT')
sim.createEntity('Variable', 'Variable:/:Vacant').Value = 0
sim.createEntity('Variable', 'Variable:/:A').Value = 1
sim.createEntity('Variable', 'Variable:/:M').Value = 0
sim.createEntity('Variable', 'Variable:/:P').Value = 0

#l = sim.createEntity('VisualizationLogProcess', 'Process:/:logger')
#l.VariableReferenceList = [['_', 'Variable:/:Vacant']]
#l.VariableReferenceList = [['_', 'Variable:/:A']]
#l.VariableReferenceList = [['_', 'Variable:/:M']]
#l.VariableReferenceList = [['_', 'Variable:/:P']]
#l.LogInterval = 0.01

p = sim.createEntity('MoleculePopulateProcess', 'Process:/:p1')
p.VariableReferenceList = [['_', 'Variable:/:A']]
p.OriginX = -1
p.UniformRadiusWidth = 20e-9
p.UniformRadiusXZ = 0.8e-6

d = sim.createEntity('DiffusionProcess', 'Process:/:diffuseA')
d.VariableReferenceList = [['_', 'Variable:/:A']]
d.D = 1e-12

b = sim.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:r1')
b.VariableReferenceList = [['_', 'Variable:/:M','-1']]
b.VariableReferenceList = [['_', 'Variable:/:A','-1']]
b.VariableReferenceList = [['_', 'Variable:/:M','1']]
b.p = V1

b = sim.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:r2')
b.VariableReferenceList = [['_', 'Variable:/:P','-1']]
b.VariableReferenceList = [['_', 'Variable:/:A','-1']]
b.VariableReferenceList = [['_', 'Variable:/:P','1']]
b.p = 1

f = sim.createEntity('FilamentProcess', 'Process:/:filam')
f.VariableReferenceList = [['_', 'Variable:/:Vacant', '-1']]
f.VariableReferenceList = [['_', 'Variable:/:M', '-2']]
f.VariableReferenceList = [['_', 'Variable:/:P', '-3']]
f.VariableReferenceList = [['_', 'Variable:/:A']]
f.LineX = 1
f.LineY = 0
f.LineZ = 0
f.Autofit = 0

l = sim.createEntity('IteratingLogProcess', 'Process:/:iter')
l.VariableReferenceList = [['_', 'Variable:/:A']]
l.LogInterval = 1e-2
l.LogEnd = T
l.Iterations = Iterations
l.FileName = filename

run(T+0.01)
