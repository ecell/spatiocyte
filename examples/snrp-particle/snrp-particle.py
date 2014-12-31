#This model performs multi-algorithm simulation using Spatiocyte next reaction
#method with lattice-based particle reaction-diffusion process.
#A is a HD species with molecules assumed to be homogeneously distributed in the
#compartment. B is a individually diffused species by the particle simulator.
#This model executes the reaction A + B -> C using Spatiocyte next reaction #method.

sim = theSimulator.createStepper('SpatiocyteStepper', 'SS')
sim.VoxelRadius = 4.4e-9 

theSimulator.rootSystem.StepperID = 'SS'
theSimulator.createEntity('Variable', 'Variable:/:GEOMETRY').Value = 0
theSimulator.createEntity('Variable', 'Variable:/:LENGTHX').Value = 5e-7
theSimulator.createEntity('Variable', 'Variable:/:LENGTHY').Value = 5e-7
theSimulator.createEntity('Variable', 'Variable:/:LENGTHZ').Value = 5e-7
theSimulator.createEntity('Variable', 'Variable:/:VACANT')
A = theSimulator.createEntity('Variable', 'Variable:/:A')
A.Value = 1000
A.Name = "HD"
theSimulator.createEntity('Variable', 'Variable:/:B').Value = 1500
theSimulator.createEntity('Variable', 'Variable:/:C').Value = 0

pop = theSimulator.createEntity('MoleculePopulateProcess', 'Process:/:pop')
pop.VariableReferenceList = [['_', 'Variable:/:B']]

dif = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffuseB')
dif.VariableReferenceList = [['_', 'Variable:/:B']]
dif.D = 1e-12

dif = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffuseC')
dif.VariableReferenceList = [['_', 'Variable:/:C']]
dif.D = 1e-12

log = theSimulator.createEntity('VisualizationLogProcess', 'Process:/:logger')
log.VariableReferenceList = [['_', 'Variable:/:B']]
log.VariableReferenceList = [['_', 'Variable:/:C']]
log.LogInterval = 5e-5

ite = theSimulator.createEntity('IteratingLogProcess', 'Process:/:logiter')
ite.VariableReferenceList = [['_', 'Variable:/:A']]
ite.VariableReferenceList = [['_', 'Variable:/:B']]
ite.VariableReferenceList = [['_', 'Variable:/:C']]
ite.LogEnd = 0.03
ite.LogInterval = 1e-4

r = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:r')
r.VariableReferenceList = [['_', 'Variable:/:A', '-1']]
r.VariableReferenceList = [['_', 'Variable:/:B', '-1']]
r.VariableReferenceList = [['_', 'Variable:/:C', '1']]
r.k = 1e-20

run(0.031)
