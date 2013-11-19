sim = theSimulator.createStepper('SpatiocyteStepper', 'SS')
sim.VoxelRadius = 6e-8 
sim.SearchVacant = 0

theSimulator.rootSystem.StepperID = 'SS'
theSimulator.createEntity('Variable', 'Variable:/:GEOMETRY').Value = 0
theSimulator.createEntity('Variable', 'Variable:/:LENGTHX').Value = 40e-6
theSimulator.createEntity('Variable', 'Variable:/:LENGTHY').Value = 25e-6
theSimulator.createEntity('Variable', 'Variable:/:LENGTHZ').Value = 11e-6
theSimulator.createEntity('Variable', 'Variable:/:VACANT')
theSimulator.createEntity('Variable', 'Variable:/:XYPLANE').Value = 0
theSimulator.createEntity('Variable', 'Variable:/:XZPLANE').Value = 0
theSimulator.createEntity('Variable', 'Variable:/:YZPLANE').Value = 0

# Create the surface compartment:
theSimulator.createEntity('System', 'System:/:Surface').StepperID = 'SS'
theSimulator.createEntity('Variable', 'Variable:/Surface:DIMENSION').Value = 2
theSimulator.createEntity('Variable', 'Variable:/Surface:VACANT')
theSimulator.createEntity('Variable', 'Variable:/:PIP2m').Value = 1233
theSimulator.createEntity('Variable', 'Variable:/:PIP3m').Value = 0
theSimulator.createEntity('Variable', 'Variable:/:PIP3a').Value = 0
theSimulator.createEntity('Variable', 'Variable:/:PTENm').Value = 309
theSimulator.createEntity('Variable', 'Variable:/:PI3Km').Value = 3096

PIP2 = theSimulator.createEntity('Variable', 'Variable:/:PIP2')
PIP2.Value = 8054
PIP2.Name = "HD"

PI3K = theSimulator.createEntity('Variable', 'Variable:/:PI3K')
PI3K.Value = 9264
PI3K.Name = "HD"

PTEN = theSimulator.createEntity('Variable', 'Variable:/:PTEN')
PTEN.Value = 6194
PTEN.Name = "HD"

logger = theSimulator.createEntity('VisualizationLogProcess', 'Process:/:logger')
#logger.VariableReferenceList = [['_', 'Variable:/:Vacant']]
logger.VariableReferenceList = [['_', 'Variable:/:PIP2m']]
logger.VariableReferenceList = [['_', 'Variable:/:PIP3m']]
logger.VariableReferenceList = [['_', 'Variable:/:PIP3a']]
logger.VariableReferenceList = [['_', 'Variable:/:PTENm']]
logger.VariableReferenceList = [['_', 'Variable:/:PI3Km']]
logger.LogInterval = 1.0 

populator = theSimulator.createEntity('MoleculePopulateProcess', 'Process:/:pop')
populator.VariableReferenceList = [['_', 'Variable:/:PIP2m']]
populator.VariableReferenceList = [['_', 'Variable:/:PIP3m']]
populator.VariableReferenceList = [['_', 'Variable:/:PIP3a']]
populator.VariableReferenceList = [['_', 'Variable:/:PTENm']]
populator.VariableReferenceList = [['_', 'Variable:/:PI3Km']]

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffusePIP2')
diffuser.VariableReferenceList = [['_', 'Variable:/:PIP2m']]
diffuser.D = 1e-14

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffusePIP3')
diffuser.VariableReferenceList = [['_', 'Variable:/:PIP3m']]
diffuser.D = 1e-14

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffusePIP3a')
diffuser.VariableReferenceList = [['_', 'Variable:/:PIP3a']]
diffuser.D = 1e-14

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffusePTEN')
diffuser.VariableReferenceList = [['_', 'Variable:/:PTENm']]
diffuser.D = 1e-14

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffusePI3K')
diffuser.VariableReferenceList = [['_', 'Variable:/:PI3Km']]
diffuser.D = 1e-14

react = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:recruitPIP2')
react.VariableReferenceList = [['_', 'Variable:/:PIP2', '-1']]
react.VariableReferenceList = [['_', 'Variable:/:PIP2m', '1']]
react.k = 4e-2

react = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:recruitPTEN')
react.VariableReferenceList = [['_', 'Variable:/:PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/:PIP2m', '-1']]
react.VariableReferenceList = [['_', 'Variable:/:PTENm', '1']]
react.VariableReferenceList = [['_', 'Variable:/:PIP2m', '1']]
react.k = 2e-14

react = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:recruitPI3Ka')
react.VariableReferenceList = [['_', 'Variable:/:PIP3a', '-1']]
react.VariableReferenceList = [['_', 'Variable:/:PI3K', '-1']]
react.VariableReferenceList = [['_', 'Variable:/:PIP3m', '1']]
react.VariableReferenceList = [['_', 'Variable:/:PI3Km', '1']]
react.k = 1e-13

binder = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:dimerPIP3')
binder.VariableReferenceList = [['_', 'Variable:/:PIP3m','-1']]
binder.VariableReferenceList = [['_', 'Variable:/:PIP3m','-1']]
binder.VariableReferenceList = [['_', 'Variable:/:PIP3a','1']]
binder.VariableReferenceList = [['_', 'Variable:/:PIP3a','1']]
binder.p = 0.65

binder = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:PIP2toPIP3')
binder.VariableReferenceList = [['_', 'Variable:/:PIP2m','-1']]
binder.VariableReferenceList = [['_', 'Variable:/:PI3Km','-1']]
binder.VariableReferenceList = [['_', 'Variable:/:PIP3m','1']]
binder.VariableReferenceList = [['_', 'Variable:/:PI3Km','1']]
binder.p = 0.17

binder = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:PIP3toPIP2')
binder.VariableReferenceList = [['_', 'Variable:/:PIP3m','-1']]
binder.VariableReferenceList = [['_', 'Variable:/:PTENm','-1']]
binder.VariableReferenceList = [['_', 'Variable:/:PIP2m','1']]
binder.VariableReferenceList = [['_', 'Variable:/:PTENm','1']]
binder.p = 1

binder = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:PIP3atoPIP2')
binder.VariableReferenceList = [['_', 'Variable:/:PIP3a','-1']]
binder.VariableReferenceList = [['_', 'Variable:/:PTENm','-1']]
binder.VariableReferenceList = [['_', 'Variable:/:PIP2m','1']]
binder.VariableReferenceList = [['_', 'Variable:/:PTENm','1']]
binder.p = 1  

react = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:dissociatePTEN')
react.VariableReferenceList = [['_', 'Variable:/:PTENm', '-1']]
react.VariableReferenceList = [['_', 'Variable:/:PTEN', '1']]
react.k = 0.09

react = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:dissociatePI3K')
react.VariableReferenceList = [['_', 'Variable:/:PI3Km', '-1']]
react.VariableReferenceList = [['_', 'Variable:/:PI3K', '1']]
react.k = 0.02

react = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:dissociatePIP3')
react.VariableReferenceList = [['_', 'Variable:/:PIP3m', '-1']]
react.VariableReferenceList = [['_', 'Variable:/:PIP2', '1']]
react.k = 0.02

react = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:dissociatePIP3a')
react.VariableReferenceList = [['_', 'Variable:/:PIP3a', '-1']]
react.VariableReferenceList = [['_', 'Variable:/:PIP2', '1']]
react.k = 0.02

react = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:dissociatePIP2')
react.VariableReferenceList = [['_', 'Variable:/:PIP2m', '-1']]
react.VariableReferenceList = [['_', 'Variable:/:PIP2', '1']]
react.k = 0.0001

migration = theSimulator.createEntity('MigrationProcess', 'Process:/:migration')
migration.VariableReferenceList = [['_', 'Variable:/:PIP2m','0']]
migration.VariableReferenceList = [['_', 'Variable:/:PIP3m','0']]
migration.VariableReferenceList = [['_', 'Variable:/:PIP3a','0']]
migration.VariableReferenceList = [['_', 'Variable:/:PTENm','0']]
migration.VariableReferenceList = [['_', 'Variable:/:PI3Km','0']]
migration.VariableReferenceList = [['_', 'Variable:/:PIP2','0']]
migration.VariableReferenceList = [['_', 'Variable:/:PI3K','0']]
migration.VariableReferenceList = [['_', 'Variable:/:PTEN','0']]
migration.minhvecX = -0.00155505
migration.minhvecY = -0.00163323 
migration.minhvecZ = 0
migration.maxhvecX =0.00412311 
migration.maxhvecY =0.00130165 
migration.maxhvecZ =0.000361489 
migration.FileName = 'migration.000'


run(17809)
