import time

ss = theSimulator.createStepper('SpatiocyteStepper', 'SS')
ss.VoxelRadius = 10e-9;

de = theSimulator.createStepper('ODEStepper', 'DE')
de.MaxStepInterval = 1e-4

theSimulator.rootSystem.StepperID = 'SS'
theSimulator.createEntity('Variable', 'Variable:/:GEOMETRY').Value = 0
theSimulator.createEntity('Variable', 'Variable:/:LENGTHX').Value = 4.42e-6;
theSimulator.createEntity('Variable', 'Variable:/:LENGTHY').Value = 4.42e-6;
theSimulator.createEntity('Variable', 'Variable:/:LENGTHZ').Value = 4.42e-6;

theSimulator.createEntity('Variable', 'Variable:/:VACANT')
E = theSimulator.createEntity('Variable', 'Variable:/:E')
E.Value = 909
E.Name = 'HD'
S = theSimulator.createEntity('Variable', 'Variable:/:S')
S.Value = 9091
S.Name = 'HD'
ES = theSimulator.createEntity('Variable', 'Variable:/:ES')
ES.Value = 0
ES.Name = 'HD'
P = theSimulator.createEntity('Variable', 'Variable:/:P')
P.Value = 0
P.Name = 'HD'

fwd = theSimulator.createEntity('MassActionProcess', 'Process:/:fwd')
fwd.StepperID = 'DE'
fwd.VariableReferenceList = [['_', 'Variable:.:E','-1']]
fwd.VariableReferenceList = [['_', 'Variable:.:S','-1']]
fwd.VariableReferenceList = [['_', 'Variable:.:ES','1']]
fwd.k = 0.01e-18

back = theSimulator.createEntity('MassActionProcess', 'Process:/:back')
back.StepperID = 'DE'
back.VariableReferenceList = [['_', 'Variable:.:ES', '-1']]
back.VariableReferenceList = [['_', 'Variable:.:E', '1']]
back.VariableReferenceList = [['_', 'Variable:.:S', '1']]
back.k = 1

prod = theSimulator.createEntity('MassActionProcess', 'Process:/:prod')
prod.StepperID = 'DE'
prod.VariableReferenceList = [['_', 'Variable:.:ES', '-1']]
prod.VariableReferenceList = [['_', 'Variable:.:E', '1']]
prod.VariableReferenceList = [['_', 'Variable:.:P', '1']]
prod.k = 1

log = theSimulator.createEntity('IteratingLogProcess', 'Process:/:log')
log.VariableReferenceList = [['_', 'Variable:.:E']]
log.VariableReferenceList = [['_', 'Variable:.:S']]
log.VariableReferenceList = [['_', 'Variable:.:ES']]
log.VariableReferenceList = [['_', 'Variable:.:P']]
log.LogInterval = 0.01
log.LogEnd = 10
log.FileName = "ode_andrews_2018.csv"

run(0.0001)
start = time.time()
run(10)
end = time.time()
print end-start
