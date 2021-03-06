import os
import sys
import numpy
import csv
import math
from matplotlib import rc
import matplotlib.pyplot as plt
from pylab import *
from collections import OrderedDict
from scipy.integrate import odeint

volume = 1000.0
molecule_radius = 0.085
box_l = math.pow(volume, 1.0/3.0)

ka1, kd1, kcat1 = 0.04483455086786913, 1.35, 1.5
ka2, kd2, kcat2 = 0.09299017957780264, 1.73, 15.0
trel = 1e-6
k7 = math.log(2.)/trel
D = 0.06
duration = 200

def kon(k):
  kD = 4*3.14*2*molecule_radius*2*D
  return k*kD/(k+kD)

def koff(kd,ka):
  return kon(ka)*kd/ka 

kon1 = kon(ka1)
koff1 = koff(kd1, ka1)

# (K + KK == K_KK | (kon1, koff1)
#     > Kp_KK | kcat1
#     > Kpp + KK | kcat2)
# (Kpp + PP == Kpp_PP | (kon1, koff1)
#     > Kp + PP | kcat1
#     > K + PP | kcat2)

# K = -kon1*KK*K + koff1*K_KK + kcat2*Kp*PP
# KK = -kon1*KK*K + koff1*K_KK + kcat2*Kp_KK
# K_KK = kon1*KK*K - koff1*K_KK - kcat1*K_KK
# Kp_KK = kcat1*K_KK - kcat2*Kp_KK 
# Kpp = kcat2*Kp_KK - kon1*Kpp*PP + koff1*Kpp_PP
# PP = -kon1*Kpp*PP + koff1*Kpp_PP + kcat1*Kpp_PP
# Kpp_PP = kon1*Kpp*PP - koff1*Kpp_PP - kcat1*Kpp_PP
# Kp = kcat1*Kpp_PP - kcat2*Kp*PP

# x[0] = -kon1*x[1]*x[0] + koff1*x[2] + kcat2*x[7]*x[5]
# x[1] = -kon1*x[1]*x[0] + koff1*x[2] + kcat2*x[3]
# x[2] = kon1*x[1]*x[0] - koff1*x[2] - kcat1*x[2]
# x[3] = kcat1*x[2] - kcat2*x[3] 
# x[4] = kcat2*x[3] - kon1*x[4]*x[5] + koff1*x[6]
# x[5] = -kon1*x[4]*x[5] + koff1*x[6] + kcat1*x[6]
# x[6] = kon1*x[4]*x[5] - koff1*x[6] - kcat1*x[6]
# x[7] = kcat1*x[6] - kcat2*x[7]*x[5]

def f(x, t0):
  return np.array([ 
    -kon1*x[1]*x[0] + koff1*x[2] + kcat2*x[7]*x[5],
    -kon1*x[1]*x[0] + koff1*x[2] + kcat2*x[3],
    kon1*x[1]*x[0] - koff1*x[2] - kcat1*x[2],
    kcat1*x[2] - kcat2*x[3],
    kcat2*x[3] - kon1*x[4]*x[5] + koff1*x[6],
    -kon1*x[4]*x[5] + koff1*x[6] + kcat1*x[6],
    kon1*x[4]*x[5] - koff1*x[6] - kcat1*x[6],
    kcat1*x[6] - kcat2*x[7]*x[5]
    ])

length = 50
NKT = 120*volume # total K
ratios = np.logspace(-1.5,1.5,length)
x = np.zeros(length)
y = np.zeros(length)
ode_time = np.linspace(0.,duration,100000)
result = []
for i in range(len(ratios)):
  NPP = int(60*volume/(ratios[i]+1)) # initial PP
  NKK = 60*volume-NPP # initial KK 
  x[i] = float(NKK)/NPP
  init_state = np.array([NKT, NKK, 0., 0., 0., NPP, 0., 0.]) / volume
  ode_result = odeint(f, y0=init_state, t=ode_time)
  #y[i:] = ode_result[-1:,(1,4,5)]
  y[i] = ode_result[-1:,(4)] # just save Kpp equilibrium values

y = y*volume/NKT
data = np.c_[x, y]
np.savetxt('ode_processive.csv', data, delimiter=',')
#fig,ax=plt.subplots(1,1,figsize=(12,8))
#ax.semilogx(x,y)
#plt.show()

