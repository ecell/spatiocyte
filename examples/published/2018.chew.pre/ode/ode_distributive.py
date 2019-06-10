import os
import sys
import numpy
import csv
import math
from matplotlib import rc
from pylab import *
from collections import OrderedDict
from scipy.integrate import odeint

volume = 1.0
molecule_radius = 0.0025
box_l = math.pow(volume, 1.0/3.0)

ka1, kd1, kcat1 = 0.04483455086786913, 1.35, 1.5
ka2, kd2, kcat2 = 0.09299017957780264, 1.73, 15.0
trel = 1e-6
k7 = math.log(2.)/trel
D = 4.
duration = 200

def kon(k):
  kD = 4*3.14*2*molecule_radius*2*D
  return k*kD/(k+kD)

def koff(kd,ka):
  return kon(ka)*kd/ka 

kon1 = kon(ka1)
koff1 = koff(kd1, ka1)
kon2 = kon(ka2)
koff2 = koff(kd2, ka2)

# (KK + K == K_KK | (kon1, koff1)
#     > Kp + KKa | kcat1)
# (KK + Kp == Kp_KK | (kon2, koff2)
#     > Kpp + KKa | kcat2)
# (KKa > KK | k7)
# (Kpp + PP == Kpp_PP | (kon1, koff1)
#     > Kp + PPa | kcat1)     
# (Kp + PP == Kp_PP | (kon2, koff2)
#     > K + PPa | kcat2)
# (PPa > PP | k7)

# KK = -kon1*KK*K + koff1*K_KK - kon2*KK*Kp + koff2*Kp_KK + k7*KKa
# K =  -kon1*KK*K + koff1*K_KK + kcat2*Kp_PP
# K_KK = kon1*KK*K - koff1*K_KK - kcat1*K_KK
# Kp = kcat1*K_KK - kon2*KK*Kp + koff2*Kp_KK + kcat1*Kpp_PP - kon2*Kp*PP + koff2*Kp_PP
# KKa = kcat1*K_KK + kcat2*Kp_KK - k7*KKa
# Kp_KK = kon2*KK*Kp - koff2*Kp_KK - kcat2*Kp_KK
# Kpp = kcat2*Kp_KK - kon1*Kpp*PP + koff1*Kpp_PP
# PP = -kon1*Kpp*PP + koff1*Kpp_PP - kon2*Kp*PP + koff2*Kp_PP + k7*PPa
# Kpp_PP = kon1*Kpp*PP - koff1*Kpp_PP - kcat1*Kpp_PP
# PPa = kcat1*Kpp_PP + kcat2*Kp_PP - k7*PPa
# Kp_PP = kon2*Kp*PP - koff2*Kp_PP - kcat2*Kp_PP

# x[0] = -kon1*x[0]*x[1] + koff1*x[2] - kon2*x[0]*x[3] + koff2*x[5] + k7*x[4]
# x[1] =  -kon1*x[0]*x[1] + koff1*x[2] + kcat2*x[10]
# x[2] = kon1*x[0]*x[1] - koff1*x[2] - kcat1*x[2]
# x[3] = kcat1*x[2] - kon2*x[0]*x[3] + koff2*x[5] + kcat1*x[8] - kon2*x[3]*x[7] + koff2*x[10]
# x[4] = kcat1*x[2] + kcat2*x[5] - k7*x[4]
# x[5] = kon2*x[0]*x[3] - koff2*x[5] - kcat2*x[5]
# x[6] = kcat2*x[5] - kon1*x[6]*x[7] + koff1*x[8]
# x[7] = -kon1*x[6]*x[7] + koff1*x[8] - kon2*x[3]*x[7] + koff2*x[10] + k7*x[9]
# x[8] = kon1*x[6]*x[7] - koff1*x[8] - kcat1*x[8]
# x[9] = kcat1*x[8] + kcat2*x[10] - k7*x[9]
# x[10] = kon2*x[3]*x[7] - koff2*x[10] - kcat2*x[10]

def f(x, t0):
  return np.array([ 
    -kon1*x[0]*x[1] + koff1*x[2] - kon2*x[0]*x[3] + koff2*x[5] + k7*x[4],
    -kon1*x[0]*x[1] + koff1*x[2] + kcat2*x[10],
    kon1*x[0]*x[1] - koff1*x[2] - kcat1*x[2],
    kcat1*x[2] - kon2*x[0]*x[3] + koff2*x[5] + kcat1*x[8] - kon2*x[3]*x[7] + koff2*x[10],
    kcat1*x[2] + kcat2*x[5] - k7*x[4],
    kon2*x[0]*x[3] - koff2*x[5] - kcat2*x[5],
    kcat2*x[5] - kon1*x[6]*x[7] + koff1*x[8],
    -kon1*x[6]*x[7] + koff1*x[8] - kon2*x[3]*x[7] + koff2*x[10] + k7*x[9],
    kon1*x[6]*x[7] - koff1*x[8] - kcat1*x[8],
    kcat1*x[8] + kcat2*x[10] - k7*x[9],
    kon2*x[3]*x[7] - koff2*x[10] - kcat2*x[10],
    ])

length = 50
NKT = 120 # total K
ratios = np.logspace(-1.5,1.5,length)
x = np.zeros(length)
y = np.zeros(length)
ode_time = np.linspace(0.,duration,100000)
result = []
for i in range(len(ratios)):
  NPP = int(60./(ratios[i]+1)) # initial PP
  NKK = 60-NPP # initial KK 
  x[i] = float(NKK)/NPP
  init_state = np.array([NKK, NKT, 0., 0., 0., 0., 0., NPP, 0., 0., 0.])/volume
  ode_result = odeint(f, y0=init_state, t=ode_time)
  #y[i:] = ode_result[-1:,(0,6,7)]
  y[i] = ode_result[-1:,(6)] # just save Kpp equilibrium values

y = y/NKT
data = np.c_[x, y]
np.savetxt('ode_distributive.csv', data, delimiter=',')
#fig,ax=plt.subplots(1,1,figsize=(12,8))
#ax.semilogx(x,y)
#plt.show()
