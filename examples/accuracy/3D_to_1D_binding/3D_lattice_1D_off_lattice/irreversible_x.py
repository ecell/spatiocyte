import math 
sim = theSimulator.createStepper('SpatiocyteStepper', 'SS')
sim.VoxelRadius = 10e-9 
sim.SearchVacant = 0

theSimulator.rootSystem.StepperID = 'SS'
theSimulator.createEntity('Variable', 'Variable:/:GEOMETRY').Value = 0
theSimulator.createEntity('Variable', 'Variable:/:LENGTHX').Value = 0.5e-6
theSimulator.createEntity('Variable', 'Variable:/:LENGTHY').Value = 16e-6
theSimulator.createEntity('Variable', 'Variable:/:LENGTHZ').Value = 16e-6
theSimulator.createEntity('Variable', 'Variable:/:YZPLANE').Value = 1
theSimulator.createEntity('Variable', 'Variable:/:VACANT')
theSimulator.createEntity('Variable', 'Variable:/:Vacant').Value = 0
theSimulator.createEntity('Variable', 'Variable:/:A').Value = 1000
theSimulator.createEntity('Variable', 'Variable:/:Border').Value = 98022
#theSimulator.createEntity('Variable', 'Variable:/:tmp').Value = 0
s = theSimulator.createEntity('Variable', 'Variable:/:sA')
s.Value = 0
s.Name = "HD"

#logger = theSimulator.createEntity('VisualizationLogProcess', 'Process:/:logger')
#logger.VariableReferenceList = [['_', 'Variable:/:Vacant']]
#logger.VariableReferenceList = [['_', 'Variable:/:A']]
#logger.VariableReferenceList = [['_', 'Variable:/:Border']]
#logger.VariableReferenceList = [['_', 'Variable:/:Interface']]
#logger.LogInterval = 0.01

pop = theSimulator.createEntity('MoleculePopulateProcess', 'Process:/:pop')
pop.VariableReferenceList = [['_', 'Variable:/:A']]
pop.UniformRadiusWidth = 5e-9
pop.UniformRadiusYZ = 0.8e-6

pop = theSimulator.createEntity('MoleculePopulateProcess', 'Process:/:popb')
pop.VariableReferenceList = [['_', 'Variable:/:Border']]
pop.UniformRadiusWidth = 20e-9
pop.UniformRadiusYZ = 8e-6

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffuseA')
diffuser.VariableReferenceList = [['_', 'Variable:/:A']]
diffuser.D = 1e-12

binder = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:r')
binder.VariableReferenceList = [['_', 'Variable:/:A','-1']]
binder.VariableReferenceList = [['_', 'Variable:/:Vacant','-1']]
binder.VariableReferenceList = [['_', 'Variable:/:sA','1']]
binder.p = 1

binder = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:r2')
binder.VariableReferenceList = [['_', 'Variable:/:A','-1']]
binder.VariableReferenceList = [['_', 'Variable:/:Border','-1']]
binder.VariableReferenceList = [['_', 'Variable:/:Border','1']]
binder.p = 1

fil = theSimulator.createEntity('FilamentProcess', 'Process:/:filam')
fil.VariableReferenceList = [['_', 'Variable:/:Vacant', '-1']]
fil.VariableReferenceList = [['_', 'Variable:/:sA']]
fil.LineX = 1
fil.LineY = 0
fil.LineZ = 0
fil.Autofit = 0

logger = theSimulator.createEntity('IteratingLogProcess', 'Process:/:iter')
logger.VariableReferenceList = [['_', 'Variable:/:A']]
logger.LogInterval = 1e-2
logger.LogEnd = 9
logger.Iterations = 5
logger.FileName = "IterateLogX.csv"

run(9.1)
