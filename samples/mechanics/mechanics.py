sim = theSimulator.createStepper('SpatiocyteStepper', 'SS')
sim.VoxelRadius = 12e-9 
sim.SearchVacant = 0

theSimulator.rootSystem.StepperID = 'SS'
theSimulator.createEntity('Variable', 'Variable:/:GEOMETRY').Value = 0
theSimulator.createEntity('Variable', 'Variable:/:LENGTHX').Value = 1.0264e-6
theSimulator.createEntity('Variable', 'Variable:/:LENGTHY').Value = 1.0264e-6
theSimulator.createEntity('Variable', 'Variable:/:LENGTHZ').Value = 1.0264e-6
theSimulator.createEntity('Variable', 'Variable:/:VACANT')
theSimulator.createEntity('Variable', 'Variable:/:XYPLANE').Value = 0
theSimulator.createEntity('Variable', 'Variable:/:XZPLANE').Value = 0
theSimulator.createEntity('Variable', 'Variable:/:YZPLANE').Value = 0

# Create the surface compartment:
theSimulator.createEntity('System', 'System:/:Surface').StepperID = 'SS'
theSimulator.createEntity('Variable', 'Variable:/Surface:DIMENSION').Value = 2
theSimulator.createEntity('Variable', 'Variable:/Surface:VACANT')
theSimulator.createEntity('Variable', 'Variable:/Surface:A').Value = 500
theSimulator.createEntity('Variable', 'Variable:/Surface:As').Value = 0
theSimulator.createEntity('Variable', 'Variable:/:B').Value = 0

logger = theSimulator.createEntity('VisualizationLogProcess', 'Process:/:logger')
logger.VariableReferenceList = [['_', 'Variable:/Surface:A'], ['_', 'Variable:/Surface:As'],['_', 'Variable:/:Vacant'],['_', 'Variable:/:B']]
logger.LogInterval = 0.0001


populator = theSimulator.createEntity('MoleculePopulateProcess', 'Process:/:pop')
populator.VariableReferenceList = [['_', 'Variable:/Surface:A']]
populator.VariableReferenceList = [['_', 'Variable:/Surface:As']]
populator.VariableReferenceList = [['_', 'Variable:/:B']]

binder = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:reaction1')
binder.VariableReferenceList = [['_', 'Variable:/Surface:A','-1']]
binder.VariableReferenceList = [['_', 'Variable:/Surface:A','-1']]
binder.VariableReferenceList = [['_', 'Variable:/Surface:As','1']]
binder.VariableReferenceList = [['_', 'Variable:/Surface:As','1']]
binder.p = 0.0001

binder = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:reaction2')
binder.VariableReferenceList = [['_', 'Variable:/Surface:A','-1']]
binder.VariableReferenceList = [['_', 'Variable:/Surface:As','-1']]
binder.VariableReferenceList = [['_', 'Variable:/Surface:As','1']]
binder.VariableReferenceList = [['_', 'Variable:/Surface:As','1']]
binder.p = 1

react = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:dissocANIOsLip')
react.VariableReferenceList = [['_', 'Variable:/Surface:As', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:A', '1']]
react.Deoligomerize = 6
react.Rates = [243, 81, 27, 9, 3, 1]
#react.k = 100

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffuseA')
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:A']]
diffuser.D = 1e-13

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffuseB')
diffuser.VariableReferenceList = [['_', 'Variable:/:B']]
diffuser.D = 0

mechanics = theSimulator.createEntity('MechanicsProcess', 'Process:/:mechanics')
mechanics.VariableReferenceList = [['_', 'Variable:/:B','0']]
mechanics.FileName = 'mechanics.000'


run(17809)
