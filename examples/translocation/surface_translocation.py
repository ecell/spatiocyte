import math 
sim = theSimulator.createStepper('SpatiocyteStepper', 'SS')
sim.VoxelRadius = 10e-9 
sim.SearchVacant = 0

theSimulator.rootSystem.StepperID = 'SS'
theSimulator.createEntity('Variable', 'Variable:/:GEOMETRY').Value = 0
theSimulator.createEntity('Variable', 'Variable:/:LENGTHX').Value = 1e-6
theSimulator.createEntity('Variable', 'Variable:/:LENGTHY').Value = 2e-6
theSimulator.createEntity('Variable', 'Variable:/:LENGTHZ').Value = 3e-6
theSimulator.createEntity('Variable', 'Variable:/:VACANT')
theSimulator.createEntity('Variable', 'Variable:/:Vacant').Value = 0
theSimulator.createEntity('Variable', 'Variable:/:A').Value = 1000
theSimulator.createEntity('Variable', 'Variable:/:B').Value = 1000
theSimulator.createEntity('Variable', 'Variable:/:sA').Value = 0
theSimulator.createEntity('Variable', 'Variable:/:sB').Value = 0

logger = theSimulator.createEntity('VisualizationLogProcess', 'Process:/:logger')
logger.VariableReferenceList = [['_', 'Variable:/:A']]
logger.VariableReferenceList = [['_', 'Variable:/:B']]
logger.VariableReferenceList = [['_', 'Variable:/:sA']]
logger.VariableReferenceList = [['_', 'Variable:/:sB']]
logger.LogInterval = 0.0001

populator = theSimulator.createEntity('MoleculePopulateProcess', 'Process:/:pop')
populator.VariableReferenceList = [['_', 'Variable:/:A']]
populator.OriginX = -1
populator.UniformRadiusX = 1

populator = theSimulator.createEntity('MoleculePopulateProcess', 'Process:/:popB')
populator.VariableReferenceList = [['_', 'Variable:/:B']]
populator.OriginX = 1
populator.UniformRadiusX = 1


diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffuseA')
diffuser.VariableReferenceList = [['_', 'Variable:/:A']]
diffuser.D = 1e-12

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffuseB')
diffuser.VariableReferenceList = [['_', 'Variable:/:B']]
diffuser.D = 1e-12

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffusesA')
diffuser.VariableReferenceList = [['_', 'Variable:/:sA']]
diffuser.D = 0

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffusesB')
diffuser.VariableReferenceList = [['_', 'Variable:/:sB']]
diffuser.D = 0

dis = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:r3')
dis.VariableReferenceList = [['_', 'Variable:/:sA','-1']]
dis.VariableReferenceList = [['_', 'Variable:/:B','1']]
dis.k = 100

dis = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:r4')
dis.VariableReferenceList = [['_', 'Variable:/:sB','-1']]
dis.VariableReferenceList = [['_', 'Variable:/:A','1']]
dis.k = 100

binder = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:r1')
binder.VariableReferenceList = [['_', 'Variable:/:A','-1']]
binder.VariableReferenceList = [['_', 'Variable:/:Vacant','-1']]
binder.VariableReferenceList = [['_', 'Variable:/:sA','1']]
binder.p = 1

binder = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:r2')
binder.VariableReferenceList = [['_', 'Variable:/:B','-1']]
binder.VariableReferenceList = [['_', 'Variable:/:Vacant','-1']]
binder.VariableReferenceList = [['_', 'Variable:/:sB','1']]
binder.p = 1

fil = theSimulator.createEntity('CompartmentProcess', 'Process:/:filam')
fil.VariableReferenceList = [['_', 'Variable:/:Vacant', '-1']]
fil.VariableReferenceList = [['_', 'Variable:/:sA']]
fil.VariableReferenceList = [['_', 'Variable:/:sB']]
fil.PlaneXY = 0
fil.PlaneXZ = 0
fil.PlaneYZ = 1
fil.OriginX = -1
fil.BindingDirection = 1
fil.DissociationDirection = 0

run(0.05)
