sim = theSimulator.createStepper('SpatiocyteStepper', 'SS')
sim.VoxelRadius = 4.4e-9 
sim.SearchVacant = 0

theSimulator.rootSystem.StepperID = 'SS'
theSimulator.createEntity('Variable', 'Variable:/:GEOMETRY').Value = 0
theSimulator.createEntity('Variable', 'Variable:/:LENGTHX').Value = 1e-8
theSimulator.createEntity('Variable', 'Variable:/:LENGTHY').Value = 1e-6
theSimulator.createEntity('Variable', 'Variable:/:LENGTHZ').Value = 1e-6
theSimulator.createEntity('Variable', 'Variable:/:VACANT')
theSimulator.createEntity('Variable', 'Variable:/:XYPLANE').Value = 5
theSimulator.createEntity('Variable', 'Variable:/:XZPLANE').Value = 5
theSimulator.createEntity('Variable', 'Variable:/:YZPLANE').Value = 4

# Create the surface compartment:
theSimulator.createEntity('System', 'System:/:Surface').StepperID = 'SS'
theSimulator.createEntity('Variable', 'Variable:/Surface:DIMENSION').Value = 2
theSimulator.createEntity('Variable', 'Variable:/Surface:VACANT')
theSimulator.createEntity('Variable', 'Variable:/Surface:A').Value = 250
theSimulator.createEntity('Variable', 'Variable:/Surface:A1').Value = 250
theSimulator.createEntity('Variable', 'Variable:/Surface:As').Value = 0
theSimulator.createEntity('Variable', 'Variable:/Surface:A1s').Value = 0

logger = theSimulator.createEntity('VisualizationLogProcess', 'Process:/:logger')
logger.VariableReferenceList = [['_', 'Variable:/Surface:A'], ['_', 'Variable:/Surface:As']]

populator = theSimulator.createEntity('MoleculePopulateProcess', 'Process:/:pop')
populator.VariableReferenceList = [['_', 'Variable:/Surface:A']]
populator.VariableReferenceList = [['_', 'Variable:/Surface:As']]
populator.VariableReferenceList = [['_', 'Variable:/Surface:A1']]
populator.VariableReferenceList = [['_', 'Variable:/Surface:A1s']]

binder = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:reaction1')
binder.VariableReferenceList = [['_', 'Variable:/Surface:A','-1']]
binder.VariableReferenceList = [['_', 'Variable:/Surface:A','-1']]
binder.VariableReferenceList = [['_', 'Variable:/Surface:As','1']]
binder.VariableReferenceList = [['_', 'Variable:/Surface:As','1']]
binder.p = 0.1

binder = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:reaction2')
binder.VariableReferenceList = [['_', 'Variable:/Surface:A','-1']]
binder.VariableReferenceList = [['_', 'Variable:/Surface:As','-1']]
binder.VariableReferenceList = [['_', 'Variable:/Surface:As','1']]
binder.VariableReferenceList = [['_', 'Variable:/Surface:As','1']]
binder.p = 1

binder = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:reaction3')
binder.VariableReferenceList = [['_', 'Variable:/Surface:A1','-1']]
binder.VariableReferenceList = [['_', 'Variable:/Surface:A1','-1']]
binder.VariableReferenceList = [['_', 'Variable:/Surface:A1s','1']]
binder.VariableReferenceList = [['_', 'Variable:/Surface:A1s','1']]
binder.p = 0.1

binder = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:reaction4')
binder.VariableReferenceList = [['_', 'Variable:/Surface:A','-1']]
binder.VariableReferenceList = [['_', 'Variable:/Surface:A1','-1']]
binder.VariableReferenceList = [['_', 'Variable:/Surface:As','1']]
binder.VariableReferenceList = [['_', 'Variable:/Surface:A1s','1']]
binder.p = 0.1

binder = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:reaction5')
binder.VariableReferenceList = [['_', 'Variable:/Surface:A1','-1']]
binder.VariableReferenceList = [['_', 'Variable:/Surface:A1s','-1']]
binder.VariableReferenceList = [['_', 'Variable:/Surface:A1s','1']]
binder.VariableReferenceList = [['_', 'Variable:/Surface:A1s','1']]
binder.p = 1

binder = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:reaction6')
binder.VariableReferenceList = [['_', 'Variable:/Surface:A1','-1']]
binder.VariableReferenceList = [['_', 'Variable:/Surface:As','-1']]
binder.VariableReferenceList = [['_', 'Variable:/Surface:A1s','1']]
binder.VariableReferenceList = [['_', 'Variable:/Surface:As','1']]
binder.p = 1

react = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:dissocANIOsLip0')
react.VariableReferenceList = [['_', 'Variable:/Surface:As', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:A', '1']]
react.Deoligomerize = 6
react.Rates = [243/2.0, 81/2.0, 27/2.0, 9/2.0, 3/2.0, 1/2.0]

react = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:dissocANIOsLip2')
react.VariableReferenceList = [['_', 'Variable:/Surface:A1s', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Surface:A1', '1']]
react.Deoligomerize = 6
react.Rates = [243/2.0, 81/2.0, 27/2.0, 9/2.0, 3/2.0, 1/2.0]
#react.k = 100

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffuseA')
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:A']]
diffuser.D = 1e-13

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffuseA1')
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:A1']]
diffuser.D = 1e-13

coord = theSimulator.createEntity('IteratingLogProcess', 'Process:/:coord')
coord.VariableReferenceList = [['_', 'Variable:/Surface:A']]
coord.VariableReferenceList = [['_', 'Variable:/Surface:As']]
coord.VariableReferenceList = [['_', 'Variable:/Surface:A1']]
coord.VariableReferenceList = [['_', 'Variable:/Surface:A1s']]
coord.LogEnd = 0.99
coord.LogInterval = 0.001

run(1)
