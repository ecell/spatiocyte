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

theSimulator.createEntity('System', 'System:/:Surface').StepperID = 'SS'
theSimulator.createEntity('Variable', 'Variable:/Surface:DIMENSION').Value = 2
theSimulator.createEntity('Variable', 'Variable:/Surface:VACANT')
theSimulator.createEntity('Variable', 'Variable:/Surface:A').Value = 20
theSimulator.createEntity('Variable', 'Variable:/Surface:As').Value = 0
theSimulator.createEntity('Variable', 'Variable:/Surface:GFP').Value = 0

logger = theSimulator.createEntity('VisualizationLogProcess', 'Process:/:logger')
logger.VariableReferenceList = [['_', 'Variable:/Surface:A']]
logger.VariableReferenceList = [['_', 'Variable:/Surface:As']]
logger.VariableReferenceList = [['_', 'Variable:/Surface:GFP']]
logger.LogInterval = 0.01

#Tag 10 molecules of A with GFP, to get A-GFP. A-GFP can transition to As-GFP:
tagger = theSimulator.createEntity('TagProcess', 'Process:/:tagger')
tagger.VariableReferenceList = [['_', 'Variable:/Surface:GFP', '-1' ]]
tagger.VariableReferenceList = [['_', 'Variable:/Surface:A', '10' ]]
tagger.VariableReferenceList = [['_', 'Variable:/Surface:As']]

coord = theSimulator.createEntity('CoordinateLogProcess', 'Process:/:coord')
coord.VariableReferenceList = [['_', 'Variable:/Surface:A']]
coord.VariableReferenceList = [['_', 'Variable:/Surface:As']]
coord.VariableReferenceList = [['_', 'Variable:/Surface:GFP']]
coord.LogInterval = 0.03
coord.LogEnd = 10

micro = theSimulator.createEntity('MicroscopyTrackingProcess', 'Process:/:micro')
micro.VariableReferenceList = [['_', 'Variable:/Surface:GFP', '1']]
micro.VariableReferenceList = [['_', 'Variable:/Surface:GFP', '-1']]
micro.ExposureTime = 0.5

populator = theSimulator.createEntity('MoleculePopulateProcess', 'Process:/:pop')
populator.VariableReferenceList = [['_', 'Variable:/Surface:A']]
populator.VariableReferenceList = [['_', 'Variable:/Surface:As']]

binder = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:reaction1')
binder.VariableReferenceList = [['_', 'Variable:/Surface:A','-1']]
binder.VariableReferenceList = [['_', 'Variable:/Surface:A','-1']]
binder.VariableReferenceList = [['_', 'Variable:/Surface:As','1']]
binder.VariableReferenceList = [['_', 'Variable:/Surface:As','1']]
binder.p = 0.01

binder = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:reaction2')
binder.VariableReferenceList = [['_', 'Variable:/Surface:As','-1']]
binder.VariableReferenceList = [['_', 'Variable:/Surface:A','1']]
binder.k = 0.5

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffuseA')
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:A']]
diffuser.D = 1e-13

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffuseAs')
diffuser.VariableReferenceList = [['_', 'Variable:/Surface:As']]
diffuser.D = 1e-15

run(100)
