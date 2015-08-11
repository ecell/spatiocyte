import math 
sim = theSimulator.createStepper('SpatiocyteStepper', 'SS')
sim.VoxelRadius = 10e-9 
sim.SearchVacant = 0

theSimulator.rootSystem.StepperID = 'SS'
theSimulator.createEntity('Variable', 'Variable:/:GEOMETRY').Value = 6
theSimulator.createEntity('Variable', 'Variable:/:LENGTHX').Value = 0.5e-6
theSimulator.createEntity('Variable', 'Variable:/:LENGTHY').Value = 16e-6
theSimulator.createEntity('Variable', 'Variable:/:YZPLANE').Value = 1
theSimulator.createEntity('Variable', 'Variable:/:VACANT')
theSimulator.createEntity('Variable', 'Variable:/:A').Value = 228
s = theSimulator.createEntity('Variable', 'Variable:/:sA')
s.Value = 0
s.Name = "HD"

theSimulator.createEntity('System', 'System:/:Surface').StepperID = 'SS'
theSimulator.createEntity('Variable', 'Variable:/Surface:DIMENSION').Value = 2
theSimulator.createEntity('Variable', 'Variable:/Surface:VACANT')

logger = theSimulator.createEntity('VisualizationLogProcess', 'Process:/:logger')
logger.VariableReferenceList = [['_', 'Variable:/:A']]
logger.VariableReferenceList = [['_', 'Variable:/Surface:VACANT']]
logger.LogInterval = 0.01

pop = theSimulator.createEntity('MoleculePopulateProcess', 'Process:/:pop')
pop.VariableReferenceList = [['_', 'Variable:/:A']]
pop.UniformRadiusWidth = 5e-9
pop.UniformRadiusYZ = 100e-9

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffuseA')
diffuser.VariableReferenceList = [['_', 'Variable:/:A']]
diffuser.D = 1e-12

binder = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:r2')
binder.VariableReferenceList = [['_', 'Variable:/:A','-1']]
binder.VariableReferenceList = [['_', 'Variable:/Surface:VACANT','-1']]
binder.VariableReferenceList = [['_', 'Variable:/:sA','1']]
binder.p = 1

logger = theSimulator.createEntity('IteratingLogProcess', 'Process:/:iter')
logger.VariableReferenceList = [['_', 'Variable:/:A']]
logger.LogInterval = 1e-2
logger.LogEnd = 10
logger.Iterations = 20
logger.FileName = "IterateLogX.csv"

run(10.01)
