import math 
sim = theSimulator.createStepper('SpatiocyteStepper', 'SS')
sim.VoxelRadius = 100e-9 
sim.SearchVacant = 0

theSimulator.rootSystem.StepperID = 'SS'
theSimulator.createEntity('Variable', 'Variable:/:GEOMETRY').Value = 0
theSimulator.createEntity('Variable', 'Variable:/:LENGTHX').Value = 100e-6
theSimulator.createEntity('Variable', 'Variable:/:LENGTHY').Value = 100e-6
theSimulator.createEntity('Variable', 'Variable:/:LENGTHZ').Value = 100e-6
theSimulator.createEntity('Variable', 'Variable:/:VACANT')
theSimulator.createEntity('Variable', 'Variable:/:A').Value = 10000
theSimulator.createEntity('Variable', 'Variable:/:Vacant').Value = 0
s = theSimulator.createEntity('Variable', 'Variable:/:sA')
s.Value = 0
s.Name = "HD"

logger = theSimulator.createEntity('VisualizationLogProcess', 'Process:/:logger')
logger.VariableReferenceList = [['_', 'Variable:/:A']]

populator = theSimulator.createEntity('MoleculePopulateProcess', 'Process:/:pop')
populator.VariableReferenceList = [['_', 'Variable:/:A']]

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffuseA')
diffuser.VariableReferenceList = [['_', 'Variable:/:A']]
diffuser.D = 1e-12

binder = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:reaction1')
binder.VariableReferenceList = [['_', 'Variable:/:A','-1']]
binder.VariableReferenceList = [['_', 'Variable:/:Vacant','-1']]
binder.VariableReferenceList = [['_', 'Variable:/:sA','1']]
binder.k = 5e-6

logger = theSimulator.createEntity('IteratingLogProcess', 'Process:/:iter')
logger.VariableReferenceList = [['_', 'Variable:/:sA']]
logger.LogInterval = 1
logger.LogEnd = 200
logger.Iterations = 250
logger.FileName = "IterateLogXYd.csv"

fil = theSimulator.createEntity('CompartmentProcess', 'Process:/:filam')
fil.VariableReferenceList = [['_', 'Variable:/:Vacant', '-1']]
fil.VariableReferenceList = [['_', 'Variable:/:sA']]
fil.PlaneXY = 1

run(201)
