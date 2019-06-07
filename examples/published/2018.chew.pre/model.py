import time
import math
import numpy as np

volume = 1.0
molecule_radius = 0.0025
box_l = math.pow(volume, 1.0/3.0)

ka1, kd1, kcat1 = 0.04483455086786913, 1.35, 1.5
ka2, kd2, kcat2 = 0.09299017957780264, 1.73, 15.0
trel = 1e-6
k7 = math.log(2.)/trel
D = 4. # [4, 0.06]
ratios = np.logspace(-1.5,1.5,12)
NKT = 120 # total K
NPP = int(60./(ratios[0]+1))
NKK = 60-NPP
duration = 200

sim = theSimulator.createStepper('SpatiocyteStepper', 'SS')
sim.VoxelRadius = 1.0208582*molecule_radius

theSimulator.rootSystem.StepperID = 'SS'
theSimulator.createEntity('Variable', 'Variable:/:GEOMETRY').Value = 0
theSimulator.createEntity('Variable', 'Variable:/:LENGTHX').Value = box_l
theSimulator.createEntity('Variable', 'Variable:/:LENGTHY').Value = box_l
theSimulator.createEntity('Variable', 'Variable:/:LENGTHZ').Value = box_l

#periodic boundaries
theSimulator.createEntity('Variable', 'Variable:/:XYPLANE').Value = 1
theSimulator.createEntity('Variable', 'Variable:/:XZPLANE').Value = 1
theSimulator.createEntity('Variable', 'Variable:/:YZPLANE').Value = 1

theSimulator.createEntity('Variable', 'Variable:/:VACANT')
theSimulator.createEntity('Variable', 'Variable:/:K').Value = NKT
theSimulator.createEntity('Variable', 'Variable:/:Kp').Value = 0
theSimulator.createEntity('Variable', 'Variable:/:Kpp').Value = 0
theSimulator.createEntity('Variable', 'Variable:/:KK').Value = NKK
theSimulator.createEntity('Variable', 'Variable:/:PP').Value = NPP
theSimulator.createEntity('Variable', 'Variable:/:K_KK').Value = 0
theSimulator.createEntity('Variable', 'Variable:/:Kp_KK').Value = 0
theSimulator.createEntity('Variable', 'Variable:/:Kpp_PP').Value = 0
theSimulator.createEntity('Variable', 'Variable:/:Kp_PP').Value = 0
theSimulator.createEntity('Variable', 'Variable:/:KKa').Value = 0
theSimulator.createEntity('Variable', 'Variable:/:PPa').Value = 0

log = theSimulator.createEntity('VisualizationLogProcess', 'Process:/:log')
log.VariableReferenceList = [['_', 'Variable:.:KK']]
log.VariableReferenceList = [['_', 'Variable:.:Kpp']]
log.VariableReferenceList = [['_', 'Variable:.:PP']]
log.LogInterval = 0.001

pop = theSimulator.createEntity('MoleculePopulateProcess', 'Process:/:pop')
pop.VariableReferenceList = [['_', 'Variable:.:KK']]
pop.VariableReferenceList = [['_', 'Variable:.:PP']]
pop.VariableReferenceList = [['_', 'Variable:.:K']]

dif = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffK')
dif.VariableReferenceList = [['_', 'Variable:.:K']]
dif.D = D

dif = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffKp')
dif.VariableReferenceList = [['_', 'Variable:.:Kp']]
dif.D = D

dif = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffKpp')
dif.VariableReferenceList = [['_', 'Variable:.:Kpp']]
dif.D = D

dif = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffKK')
dif.VariableReferenceList = [['_', 'Variable:.:KK']]
dif.D = D

dif = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffPP')
dif.VariableReferenceList = [['_', 'Variable:.:PP']]
dif.D = D

dif = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffK_KK')
dif.VariableReferenceList = [['_', 'Variable:.:K_KK']]
dif.D = D

dif = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffKp_KK')
dif.VariableReferenceList = [['_', 'Variable:.:Kp_KK']]
dif.D = D

dif = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffKpp_PP')
dif.VariableReferenceList = [['_', 'Variable:.:Kpp_PP']]
dif.D = D

dif = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffKp_PP')
dif.VariableReferenceList = [['_', 'Variable:.:Kp_PP']]
dif.D = D

dif = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffKKa')
dif.VariableReferenceList = [['_', 'Variable:.:KKa']]
dif.D = D

dif = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffPPa')
dif.VariableReferenceList = [['_', 'Variable:.:PPa']]
dif.D = D

r = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:r1')
r.VariableReferenceList = [['_', 'Variable:.:K','-1']]
r.VariableReferenceList = [['_', 'Variable:.:KK','-1']]
r.VariableReferenceList = [['_', 'Variable:.:K_KK','1']]
r.k = ka1

r = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:r2')
r.VariableReferenceList = [['_', 'Variable:.:K_KK', '-1']]
r.VariableReferenceList = [['_', 'Variable:.:K', '1']]
r.VariableReferenceList = [['_', 'Variable:.:KK', '1']]
r.k = kd1

r = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:r3')
r.VariableReferenceList = [['_', 'Variable:.:K_KK', '-1']]
r.VariableReferenceList = [['_', 'Variable:.:Kp', '1']]
r.VariableReferenceList = [['_', 'Variable:.:KKa', '1']]
r.k = kcat1


r = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:r4')
r.VariableReferenceList = [['_', 'Variable:.:Kp','-1']]
r.VariableReferenceList = [['_', 'Variable:.:KK','-1']]
r.VariableReferenceList = [['_', 'Variable:.:Kp_KK','1']]
r.k = ka2

r = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:r5')
r.VariableReferenceList = [['_', 'Variable:.:Kp_KK', '-1']]
r.VariableReferenceList = [['_', 'Variable:.:Kp', '1']]
r.VariableReferenceList = [['_', 'Variable:.:KK', '1']]
r.k = kd2

r = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:r6')
r.VariableReferenceList = [['_', 'Variable:.:Kp_KK', '-1']]
r.VariableReferenceList = [['_', 'Variable:.:Kpp', '1']]
r.VariableReferenceList = [['_', 'Variable:.:KKa', '1']]
r.k = kcat2


r = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:r7')
r.VariableReferenceList = [['_', 'Variable:.:KKa', '-1']]
r.VariableReferenceList = [['_', 'Variable:.:KK', '1']]
r.k = k7


r = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:r8')
r.VariableReferenceList = [['_', 'Variable:.:Kpp','-1']]
r.VariableReferenceList = [['_', 'Variable:.:PP','-1']]
r.VariableReferenceList = [['_', 'Variable:.:Kpp_PP','1']]
r.k = ka1

r = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:r9')
r.VariableReferenceList = [['_', 'Variable:.:Kpp_PP', '-1']]
r.VariableReferenceList = [['_', 'Variable:.:Kpp', '1']]
r.VariableReferenceList = [['_', 'Variable:.:PP', '1']]
r.k = kd1

r = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:r10')
r.VariableReferenceList = [['_', 'Variable:.:Kpp_PP', '-1']]
r.VariableReferenceList = [['_', 'Variable:.:Kp', '1']]
r.VariableReferenceList = [['_', 'Variable:.:PPa', '1']]
r.k = kcat1

r = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:r11')
r.VariableReferenceList = [['_', 'Variable:.:Kp','-1']]
r.VariableReferenceList = [['_', 'Variable:.:PP','-1']]
r.VariableReferenceList = [['_', 'Variable:.:Kp_PP','1']]
r.k = ka2

r = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:r12')
r.VariableReferenceList = [['_', 'Variable:.:Kp_PP', '-1']]
r.VariableReferenceList = [['_', 'Variable:.:Kp', '1']]
r.VariableReferenceList = [['_', 'Variable:.:PP', '1']]
r.k = kd2

r = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:r13')
r.VariableReferenceList = [['_', 'Variable:.:Kp_PP', '-1']]
r.VariableReferenceList = [['_', 'Variable:.:K', '1']]
r.VariableReferenceList = [['_', 'Variable:.:PPa', '1']]
r.k = kcat2


r = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:r14')
r.VariableReferenceList = [['_', 'Variable:.:PPa', '-1']]
r.VariableReferenceList = [['_', 'Variable:.:PP', '1']]
r.k = k7

log = theSimulator.createEntity('IteratingLogProcess', 'Process:/:log2')
log.VariableReferenceList = [['_', 'Variable:.:KK']]
log.VariableReferenceList = [['_', 'Variable:.:Kpp']]
log.VariableReferenceList = [['_', 'Variable:.:PP']]
log.LogInterval = 1.0
log.LogEnd = duration

run(0.01)
start = time.time()
run(duration)
end = time.time()
print "total runtime:", end-start
