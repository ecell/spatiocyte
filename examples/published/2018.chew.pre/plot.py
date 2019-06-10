import os
import numpy as np
import math
import matplotlib.pyplot as plt
import glob

duration = 200.
dirs = glob.glob('ratio*_D_4')
ratios = np.logspace(-1.5,1.5,12)
NKT = 120
x = []
y = []
for adir in dirs:
  data = np.loadtxt(adir+'/IterateLog.csv', delimiter=',', skiprows=1)
  x.append((data[:1,1]/data[:1,3])[0])
  y.append(np.mean(data[-int(duration*0.2):,2])/NKT)

x = np.sort(x)
y = np.sort(y)

fig,ax=plt.subplots(1,1,figsize=(12,8))
ax.semilogx(x,y, 'o', markersize=10)

dirs = glob.glob('ratio*_D_0*')
ratios = np.logspace(-1.5,1.5,12)
x = []
y = []
for adir in dirs:
  data = np.loadtxt(adir+'/IterateLog.csv', delimiter=',', skiprows=1)
  x.append((data[:1,1]/data[:1,3])[0])
  y.append(np.mean(data[-int(duration*0.2):,2])/NKT)

x = np.sort(x)
y = np.sort(y)

ax.semilogx(x,y, 'o', markersize=10)

data = np.loadtxt('ode/ode_processive.csv', delimiter=',')
ax.semilogx(data[:,0], data[:,1], '-', markersize=10)

data = np.loadtxt('ode/ode_distributive.csv', delimiter=',')
ax.semilogx(data[:,0], data[:,1], '-', markersize=10)

plt.show()



