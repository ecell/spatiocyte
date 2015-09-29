import math 
sim = theSimulator.createStepper('SpatiocyteStepper', 'SS')
sim.VoxelRadius = 10e-9 
sim.SearchVacant = 0

theSimulator.rootSystem.StepperID = 'SS'
theSimulator.createEntity('Variable', 'Variable:/:GEOMETRY').Value = 6
theSimulator.createEntity('Variable', 'Variable:/:LENGTHX').Value = 16e-6
theSimulator.createEntity('Variable', 'Variable:/:LENGTHY').Value = 0.25e-6
theSimulator.createEntity('Variable', 'Variable:/:YZPLANE').Value = 1
theSimulator.createEntity('Variable', 'Variable:/:VACANT')
theSimulator.createEntity('Variable', 'Variable:/:Vacant').Value = 0
theSimulator.createEntity('Variable', 'Variable:/:A').Value = 2
theSimulator.createEntity('Variable', 'Variable:/:C').Value = 1
theSimulator.createEntity('Variable', 'Variable:/:M').Value = 0
theSimulator.createEntity('Variable', 'Variable:/:P').Value = 0
s = theSimulator.createEntity('Variable', 'Variable:/:sA')
s.Value = 0
s.Name = "HD"

#logger = theSimulator.createEntity('VisualizationLogProcess', 'Process:/:logger')
#logger.VariableReferenceList = [['_', 'Variable:/:Vacant']]
#logger.VariableReferenceList = [['_', 'Variable:/:A']]
#logger.VariableReferenceList = [['_', 'Variable:/:C']]
#logger.VariableReferenceList = [['_', 'Variable:/:M']]
#logger.VariableReferenceList = [['_', 'Variable:/:P']]
#logger.LogInterval = 0.01

pop = theSimulator.createEntity('MoleculePopulateProcess', 'Process:/:p1')
pop.VariableReferenceList = [['_', 'Variable:/:A']]
pop.UniformRadiusWidth = 15e-9
pop.UniformRadiusXZ = 0.8e-6

pop = theSimulator.createEntity('MoleculePopulateProcess', 'Process:/:p2')
pop.VariableReferenceList = [['_', 'Variable:/:C']]
pop.UniformRadiusXZ = 10e-9

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffuseA')
diffuser.VariableReferenceList = [['_', 'Variable:/:A']]
diffuser.D = 1e-12

binder = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:r1')
binder.VariableReferenceList = [['_', 'Variable:/:M','-1']]
binder.VariableReferenceList = [['_', 'Variable:/:A','-1']]
binder.VariableReferenceList = [['_', 'Variable:/:M','1']]
binder.VariableReferenceList = [['_', 'Variable:/:sA','1']]
binder.p = 1

binder = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:r2')
binder.VariableReferenceList = [['_', 'Variable:/:P','-1']]
binder.VariableReferenceList = [['_', 'Variable:/:A','-1']]
binder.VariableReferenceList = [['_', 'Variable:/:P','1']]
binder.VariableReferenceList = [['_', 'Variable:/:sA','1']]
binder.p = 1

binder = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:r3')
binder.VariableReferenceList = [['_', 'Variable:/:C','-1']]
binder.VariableReferenceList = [['_', 'Variable:/:A','-1']]
binder.VariableReferenceList = [['_', 'Variable:/:C','1']]
binder.VariableReferenceList = [['_', 'Variable:/:sA','1']]
binder.p = 1

fil = theSimulator.createEntity('FilamentProcess', 'Process:/:filam')
fil.VariableReferenceList = [['_', 'Variable:/:Vacant', '-1']]
fil.VariableReferenceList = [['_', 'Variable:/:M', '-2']]
fil.VariableReferenceList = [['_', 'Variable:/:P', '-3']]
fil.VariableReferenceList = [['_', 'Variable:/:sA']]
fil.VariableReferenceList = [['_', 'Variable:/:A']]
fil.VariableReferenceList = [['_', 'Variable:/:C']]
fil.LineX = 1
fil.LineY = 0
fil.LineZ = 0
fil.Autofit = 0

logger = theSimulator.createEntity('IteratingLogProcess', 'Process:/:iter')
logger.VariableReferenceList = [['_', 'Variable:/:A']]
logger.LogInterval = 1e-2
logger.LogEnd = 10
logger.Iterations = 10000
logger.FileName = "IterateLogX.csv"

run(10.01)
