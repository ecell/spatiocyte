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
ratios = np.logspace(-1.5,1.5,12)
ratio = ratios[5]
D = 4
NKT = 120 # total K
NPP = int(60./(ratio+1))
NKK = 60-NPP
print("ratio:", ratio, "NKK:",NKK,"NPP:",NPP)
duration = 200

def kon(k):
  kD = 4*3.14*2*molecule_radius*2*D
  return k*kD/(k+kD)
def koff(kd,ka):
  return kon(ka)*kd/ka 

sim = theSimulator.createStepper('SpatiocyteStepper', 'SS')
sim.VoxelRadius = 1.0208582*molecule_radius

theSimulator.rootSystem.StepperID = 'SS'
theSimulator.createEntity('Variable', 'Variable:/:GEOMETRY').Value = 0
theSimulator.createEntity('Variable', 'Variable:/:LENGTHX').Value = box_l
theSimulator.createEntity('Variable', 'Variable:/:LENGTHY').Value = box_l
theSimulator.createEntity('Variable', 'Variable:/:LENGTHZ').Value = box_l

theSimulator.createEntity('Variable', 'Variable:/:VACANT')

v = theSimulator.createEntity('Variable', 'Variable:/:K')
v.Value = NKT
v.Name = "HD"
v = theSimulator.createEntity('Variable', 'Variable:/:Kp')
v.Value = 0
v.Name = "HD"
v = theSimulator.createEntity('Variable', 'Variable:/:Kpp')
v.Value = 0
v.Name = "HD"
v = theSimulator.createEntity('Variable', 'Variable:/:KK')
v.Value = NKK
v.Name = "HD"
v = theSimulator.createEntity('Variable', 'Variable:/:PP')
v.Value = NPP
v.Name = "HD"
v = theSimulator.createEntity('Variable', 'Variable:/:K_KK')
v.Value = 0
v.Name = "HD"
v = theSimulator.createEntity('Variable', 'Variable:/:Kp_KK')
v.Value = 0
v.Name = "HD"
v = theSimulator.createEntity('Variable', 'Variable:/:Kpp_PP')
v.Value = 0
v.Name = "HD"
v = theSimulator.createEntity('Variable', 'Variable:/:Kp_PP')
v.Value = 0
v.Name = "HD"
v = theSimulator.createEntity('Variable', 'Variable:/:KKa')
v.Value = 0
v.Name = "HD"
v = theSimulator.createEntity('Variable', 'Variable:/:PPa')
v.Value = 0
v.Name = "HD"

r = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:r1')
r.VariableReferenceList = [['_', 'Variable:.:K','-1']]
r.VariableReferenceList = [['_', 'Variable:.:KK','-1']]
r.VariableReferenceList = [['_', 'Variable:.:K_KK','1']]
r.k = kon(ka1)

r = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:r2')
r.VariableReferenceList = [['_', 'Variable:.:K_KK', '-1']]
r.VariableReferenceList = [['_', 'Variable:.:K', '1']]
r.VariableReferenceList = [['_', 'Variable:.:KK', '1']]
r.k = koff(kd1, ka1)

r = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:r3')
r.VariableReferenceList = [['_', 'Variable:.:K_KK', '-1']]
r.VariableReferenceList = [['_', 'Variable:.:Kp', '1']]
r.VariableReferenceList = [['_', 'Variable:.:KKa', '1']]
r.k = kcat1


r = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:r4')
r.VariableReferenceList = [['_', 'Variable:.:Kp','-1']]
r.VariableReferenceList = [['_', 'Variable:.:KK','-1']]
r.VariableReferenceList = [['_', 'Variable:.:Kp_KK','1']]
r.k = kon(ka2)

r = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:r5')
r.VariableReferenceList = [['_', 'Variable:.:Kp_KK', '-1']]
r.VariableReferenceList = [['_', 'Variable:.:Kp', '1']]
r.VariableReferenceList = [['_', 'Variable:.:KK', '1']]
r.k = koff(kd2, ka2)

r = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:r6')
r.VariableReferenceList = [['_', 'Variable:.:Kp_KK', '-1']]
r.VariableReferenceList = [['_', 'Variable:.:Kpp', '1']]
r.VariableReferenceList = [['_', 'Variable:.:KKa', '1']]
r.k = kcat2


r = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:r7')
r.VariableReferenceList = [['_', 'Variable:.:KKa', '-1']]
r.VariableReferenceList = [['_', 'Variable:.:KK', '1']]
r.k = k7


r = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:r8')
r.VariableReferenceList = [['_', 'Variable:.:Kpp','-1']]
r.VariableReferenceList = [['_', 'Variable:.:PP','-1']]
r.VariableReferenceList = [['_', 'Variable:.:Kpp_PP','1']]
r.k = kon(ka1)

r = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:r9')
r.VariableReferenceList = [['_', 'Variable:.:Kpp_PP', '-1']]
r.VariableReferenceList = [['_', 'Variable:.:Kpp', '1']]
r.VariableReferenceList = [['_', 'Variable:.:PP', '1']]
r.k = koff(kd1, ka1)

r = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:r10')
r.VariableReferenceList = [['_', 'Variable:.:Kpp_PP', '-1']]
r.VariableReferenceList = [['_', 'Variable:.:Kp', '1']]
r.VariableReferenceList = [['_', 'Variable:.:PPa', '1']]
r.k = kcat1

r = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:r11')
r.VariableReferenceList = [['_', 'Variable:.:Kp','-1']]
r.VariableReferenceList = [['_', 'Variable:.:PP','-1']]
r.VariableReferenceList = [['_', 'Variable:.:Kp_PP','1']]
r.k = kon(ka2)

r = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:r12')
r.VariableReferenceList = [['_', 'Variable:.:Kp_PP', '-1']]
r.VariableReferenceList = [['_', 'Variable:.:Kp', '1']]
r.VariableReferenceList = [['_', 'Variable:.:PP', '1']]
r.k = koff(kd2, ka2)

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
log.Iterations = 1000
log.FileName = "snrp_distributive.csv"

run(0.01)
start = time.time()
run(duration)
end = time.time()
print "total runtime:", end-start
