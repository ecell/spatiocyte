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

theSimulator.createEntity('System', 'System:/:Cell').StepperID = 'SS'
theSimulator.createEntity('Variable', 'Variable:/Cell:PIP2m').Value = 1233
theSimulator.createEntity('Variable', 'Variable:/Cell:PIP3m').Value = 0
theSimulator.createEntity('Variable', 'Variable:/Cell:PIP3a').Value = 0
theSimulator.createEntity('Variable', 'Variable:/Cell:PTENm').Value = 309
theSimulator.createEntity('Variable', 'Variable:/Cell:PI3Km').Value = 3096

PIP2 = theSimulator.createEntity('Variable', 'Variable:/Cell:PIP2')
PIP2.Value = 8054
PIP2.Name = "HD"

PI3K = theSimulator.createEntity('Variable', 'Variable:/Cell:PI3K')
PI3K.Value = 9264
PI3K.Name = "HD"

PTEN = theSimulator.createEntity('Variable', 'Variable:/Cell:PTEN')
PTEN.Value = 6194
PTEN.Name = "HD"

logger = theSimulator.createEntity('VisualizationLogProcess', 'Process:/Cell:logger')
#logger.VariableReferenceList = [['_', 'Variable:/Cell:VACANT']]
logger.VariableReferenceList = [['_', 'Variable:/Cell:PIP2m']]
logger.VariableReferenceList = [['_', 'Variable:/Cell:PIP3m']]
logger.VariableReferenceList = [['_', 'Variable:/Cell:PIP3a']]
logger.VariableReferenceList = [['_', 'Variable:/Cell:PTENm']]
logger.VariableReferenceList = [['_', 'Variable:/Cell:PI3Km']]
logger.LogInterval = 1.0 

populator = theSimulator.createEntity('MoleculePopulateProcess', 'Process:/Cell:pop')
populator.VariableReferenceList = [['_', 'Variable:/Cell:PIP2m']]
populator.VariableReferenceList = [['_', 'Variable:/Cell:PIP3m']]
populator.VariableReferenceList = [['_', 'Variable:/Cell:PIP3a']]
populator.VariableReferenceList = [['_', 'Variable:/Cell:PTENm']]
populator.VariableReferenceList = [['_', 'Variable:/Cell:PI3Km']]

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/Cell:diffusePIP2')
diffuser.VariableReferenceList = [['_', 'Variable:/Cell:PIP2m']]
diffuser.D = 1e-14

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/Cell:diffusePIP3')
diffuser.VariableReferenceList = [['_', 'Variable:/Cell:PIP3m']]
diffuser.D = 1e-14

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/Cell:diffusePIP3a')
diffuser.VariableReferenceList = [['_', 'Variable:/Cell:PIP3a']]
diffuser.D = 1e-14

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/Cell:diffusePTEN')
diffuser.VariableReferenceList = [['_', 'Variable:/Cell:PTENm']]
diffuser.D = 1e-14

diffuser = theSimulator.createEntity('DiffusionProcess', 'Process:/Cell:diffusePI3K')
diffuser.VariableReferenceList = [['_', 'Variable:/Cell:PI3Km']]
diffuser.D = 1e-14

react = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/Cell:recruitPIP2')
react.VariableReferenceList = [['_', 'Variable:/Cell:PIP2', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Cell:PIP2m', '1']]
react.k = 4e-2

react = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/Cell:recruitPTEN')
react.VariableReferenceList = [['_', 'Variable:/Cell:PTEN', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Cell:PIP2m', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Cell:PTENm', '1']]
react.VariableReferenceList = [['_', 'Variable:/Cell:PIP2m', '1']]
react.k = 2e-14

react = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/Cell:recruitPI3Ka')
react.VariableReferenceList = [['_', 'Variable:/Cell:PIP3a', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Cell:PI3K', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Cell:PIP3m', '1']]
react.VariableReferenceList = [['_', 'Variable:/Cell:PI3Km', '1']]
react.k = 1e-13

binder = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/Cell:dimerPIP3')
binder.VariableReferenceList = [['_', 'Variable:/Cell:PIP3m','-1']]
binder.VariableReferenceList = [['_', 'Variable:/Cell:PIP3m','-1']]
binder.VariableReferenceList = [['_', 'Variable:/Cell:PIP3a','1']]
binder.VariableReferenceList = [['_', 'Variable:/Cell:PIP3a','1']]
binder.p = 0.65

binder = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/Cell:PIP2toPIP3')
binder.VariableReferenceList = [['_', 'Variable:/Cell:PIP2m','-1']]
binder.VariableReferenceList = [['_', 'Variable:/Cell:PI3Km','-1']]
binder.VariableReferenceList = [['_', 'Variable:/Cell:PIP3m','1']]
binder.VariableReferenceList = [['_', 'Variable:/Cell:PI3Km','1']]
binder.p = 0.17

binder = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/Cell:PIP3toPIP2')
binder.VariableReferenceList = [['_', 'Variable:/Cell:PIP3m','-1']]
binder.VariableReferenceList = [['_', 'Variable:/Cell:PTENm','-1']]
binder.VariableReferenceList = [['_', 'Variable:/Cell:PIP2m','1']]
binder.VariableReferenceList = [['_', 'Variable:/Cell:PTENm','1']]
binder.p = 1

binder = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/Cell:PIP3atoPIP2')
binder.VariableReferenceList = [['_', 'Variable:/Cell:PIP3a','-1']]
binder.VariableReferenceList = [['_', 'Variable:/Cell:PTENm','-1']]
binder.VariableReferenceList = [['_', 'Variable:/Cell:PIP2m','1']]
binder.VariableReferenceList = [['_', 'Variable:/Cell:PTENm','1']]
binder.p = 1  

react = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/Cell:dissociatePTEN')
react.VariableReferenceList = [['_', 'Variable:/Cell:PTENm', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Cell:PTEN', '1']]
react.k = 0.09

react = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/Cell:dissociatePI3K')
react.VariableReferenceList = [['_', 'Variable:/Cell:PI3Km', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Cell:PI3K', '1']]
react.k = 0.02

react = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/Cell:dissociatePIP3')
react.VariableReferenceList = [['_', 'Variable:/Cell:PIP3m', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Cell:PIP2', '1']]
react.k = 0.02

react = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/Cell:dissociatePIP3a')
react.VariableReferenceList = [['_', 'Variable:/Cell:PIP3a', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Cell:PIP2', '1']]
react.k = 0.02

react = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/Cell:dissociatePIP2')
react.VariableReferenceList = [['_', 'Variable:/Cell:PIP2m', '-1']]
react.VariableReferenceList = [['_', 'Variable:/Cell:PIP2', '1']]
react.k = 0.0001

mig = theSimulator.createEntity('MigrationProcess', 'Process:/Cell:Migration')
mig.VariableReferenceList = [['_', 'Variable:/Cell:PIP2m','0']]
mig.VariableReferenceList = [['_', 'Variable:/Cell:PIP3m','0']]
mig.VariableReferenceList = [['_', 'Variable:/Cell:PIP3a','0']]
mig.VariableReferenceList = [['_', 'Variable:/Cell:PTENm','0']]
mig.VariableReferenceList = [['_', 'Variable:/Cell:PI3Km','0']]
mig.VariableReferenceList = [['_', 'Variable:/Cell:PIP2','0']]
mig.VariableReferenceList = [['_', 'Variable:/Cell:PI3K','0']]
mig.VariableReferenceList = [['_', 'Variable:/Cell:PTEN','0']]
mig.minhvecX = -0.00155505
mig.minhvecY = -0.00163323 
mig.minhvecZ = 0
mig.maxhvecX =0.00412311 
mig.maxhvecY =0.00130165 
mig.maxhvecZ =0.000361489 
mig.FileName = 'migration.000'


run(17809)
