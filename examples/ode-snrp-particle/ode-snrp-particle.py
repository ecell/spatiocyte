# This model performs multi-algorithm simulation using Spatiocyte next reaction,
# mass-action and lattice-based particle reaction-diffusion methods.

# The following reactions are performed simultaneously:
# With mass action ODE solver:
# E + S -> ES
# ES -> E + S
# ES -> E + P

# With Spatiocyte next-reaction method:
# P + A -> B

# With Spatiocyte diffusion-influenced reaction method:
# A + B -> C

# Homogeneously distributed (HD) species:
# E, S, ES, P

# Individually diffused (non-HD) species:
# A, B, C

ss = theSimulator.createStepper('SpatiocyteStepper', 'SS')
ss.VoxelRadius = 4.4e-9 

de = theSimulator.createStepper('ODEStepper', 'DE')
de.MaxStepInterval = 1e-3

theSimulator.rootSystem.StepperID = 'SS'
theSimulator.createEntity('Variable', 'Variable:/:GEOMETRY').Value = 0
theSimulator.createEntity('Variable', 'Variable:/:LENGTHX').Value = 5e-7
theSimulator.createEntity('Variable', 'Variable:/:LENGTHY').Value = 5e-7
theSimulator.createEntity('Variable', 'Variable:/:LENGTHZ').Value = 5e-7

theSimulator.createEntity('Variable', 'Variable:/:VACANT')
theSimulator.createEntity('Variable', 'Variable:/:A').Value = 1500
theSimulator.createEntity('Variable', 'Variable:/:B').Value = 0
theSimulator.createEntity('Variable', 'Variable:/:C').Value = 0

E = theSimulator.createEntity('Variable', 'Variable:/:E')
E.Value = 100
E.Name = 'HD'
S = theSimulator.createEntity('Variable', 'Variable:/:S')
S.Value = 1000
S.Name = 'HD'
ES = theSimulator.createEntity('Variable', 'Variable:/:ES')
ES.Value = 0
ES.Name = 'HD'
P = theSimulator.createEntity('Variable', 'Variable:/:P')
P.Value = 0
P.Name = 'HD'

# E + S -> ES
fwd = theSimulator.createEntity('MassActionProcess', 'Process:/:fwd')
fwd.StepperID = 'DE'
fwd.VariableReferenceList = [['_', 'Variable:.:E','-1']]
fwd.VariableReferenceList = [['_', 'Variable:.:S','-1']]
fwd.VariableReferenceList = [['_', 'Variable:.:ES','1']]
fwd.k = 1e-22

# ES -> E + S
back = theSimulator.createEntity('MassActionProcess', 'Process:/:back')
back.StepperID = 'DE'
back.VariableReferenceList = [['_', 'Variable:.:ES', '-1']]
back.VariableReferenceList = [['_', 'Variable:.:E', '1']]
back.VariableReferenceList = [['_', 'Variable:.:S', '1']]
back.k = 1e-1

# ES -> E + P
prod = theSimulator.createEntity('MassActionProcess', 'Process:/:prod')
prod.StepperID = 'DE'
prod.VariableReferenceList = [['_', 'Variable:.:ES', '-1']]
prod.VariableReferenceList = [['_', 'Variable:.:E', '1']]
prod.VariableReferenceList = [['_', 'Variable:.:P', '1']]
prod.k = 1e-1

# P + A -> B
s = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:snrp')
s.VariableReferenceList = [['_', 'Variable:/:P', '-1']]
s.VariableReferenceList = [['_', 'Variable:/:A', '-1']]
s.VariableReferenceList = [['_', 'Variable:/:B', '1']]
s.k = 5e-24

# A + B -> C
d = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:dirp')
d.VariableReferenceList = [['_', 'Variable:/:A', '-1']]
d.VariableReferenceList = [['_', 'Variable:/:B', '-1']]
d.VariableReferenceList = [['_', 'Variable:/:C', '1']]
d.k = 5e-24

pop = theSimulator.createEntity('MoleculePopulateProcess', 'Process:/:pop')
pop.VariableReferenceList = [['_', 'Variable:/:A']]

dif = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffuseA')
dif.VariableReferenceList = [['_', 'Variable:/:A']]
dif.D = 5e-16

dif = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffuseB')
dif.VariableReferenceList = [['_', 'Variable:/:B']]
dif.D = 5e-16

dif = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffuseC')
dif.VariableReferenceList = [['_', 'Variable:/:C']]
dif.D = 5e-16

log = theSimulator.createEntity('VisualizationLogProcess', 'Process:/:logger')
log.VariableReferenceList = [['_', 'Variable:/:A']]
log.VariableReferenceList = [['_', 'Variable:/:B']]
log.VariableReferenceList = [['_', 'Variable:/:C']]
log.LogInterval = 1e-1

ite = theSimulator.createEntity('IteratingLogProcess', 'Process:/:logiter')
ite.VariableReferenceList = [['_', 'Variable:/:A']]
ite.VariableReferenceList = [['_', 'Variable:/:B']]
ite.VariableReferenceList = [['_', 'Variable:/:C']]
ite.VariableReferenceList = [['_', 'Variable:.:E']]
ite.VariableReferenceList = [['_', 'Variable:.:S']]
ite.VariableReferenceList = [['_', 'Variable:.:ES']]
ite.VariableReferenceList = [['_', 'Variable:.:P']]
ite.LogInterval = 1e-2
ite.LogEnd = 99

run(100)
