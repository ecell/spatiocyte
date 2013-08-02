sim = theSimulator.createStepper('SpatiocyteStepper', 'SS')
sim.VoxelRadius = 2.5e-9
sim.SearchVacant = 0

theSimulator.rootSystem.StepperID = 'SS'
theSimulator.createEntity('Variable', 'Variable:/:GEOMETRY').Value = 0
theSimulator.createEntity('Variable', 'Variable:/:LENGTHX').Value = 5e-8
theSimulator.createEntity('Variable', 'Variable:/:LENGTHY').Value = 5e-8
theSimulator.createEntity('Variable', 'Variable:/:LENGTHZ').Value = 5e-8
theSimulator.createEntity('Variable', 'Variable:/:VACANT')
theSimulator.createEntity('Variable', 'Variable:/:XYPLANE').Value = 0
theSimulator.createEntity('Variable', 'Variable:/:XZPLANE').Value = 0
theSimulator.createEntity('Variable', 'Variable:/:YZPLANE').Value = 0

theSimulator.createEntity('Variable', 'Variable:/:A').Value = 1

#logger = theSimulator.createEntity('VisualizationLogProcess', 'Process:/:logger')
#logger.VariableReferenceList = [['_', 'Variable:/:A']]

populator = theSimulator.createEntity('MoleculePopulateProcess', 'Process:/:pop')
populator.VariableReferenceList = [['_', 'Variable:/:A']]

diffuser = theSimulator.createEntity('PeriodicBoundaryDiffusionProcess', 'Process:/:diffuse')
diffuser.VariableReferenceList = [['_', 'Variable:/:A']]
diffuser.D = 7.5757e-14
diffuser.Origins = 1

iterator = theSimulator.createEntity('IteratingLogProcess', 'Process:/:iterate')
iterator.VariableReferenceList = [['_', 'Variable:/:A']]
iterator.Iterations = 1
iterator.LogEnd = 80000000
iterator.LogStart = 1
iterator.LogInterval = 10
iterator.Diffusion = 1
#iterator.SquaredDisplacement = 1

import time
run(1e-6)
print "Done stirring. Now running..."
start = time.time()
run(80000000.1)
end = time.time()
duration = end-start
print duration


