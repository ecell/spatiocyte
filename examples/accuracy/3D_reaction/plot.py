import numpy as np
import pylab as P
from scipy.special import erfc,erf
import math

def g(t):
  n = 10000.0
  A = 100e-6*100e-6
  V = 100e-6*100e-6*100e-6
  k = 5e-6
  D = 1e-6
  h = k/D
  return n*A/(h*V)*(np.exp(h*h*D*t)*erfc(h*np.sqrt(D*t))-1.0+2/np.pi*h*np.sqrt(D*t))

def g2(t):
  n = 10000.0
  A = 100e-6*100e-6
  V = 100e-6*100e-6*100e-6
  k = 5e-6
  D = 1e-6
  h = k/D
  return n*A/(h*V)*(math.exp(h*h*D*t)*(1-erf(h*math.sqrt(D*t)))-1.0+2/math.pi*h*math.sqrt(D*t))

labelFontSize = 20
tickFontSize = 20
legendFontSize = 20
lineFontSize = 20

filenames = ['spatiocyte_3D_k0.100.csv']
titles = ['spatiocyte_3D_k0.100.csv']
lines = ['-', '--', '-', '--', '-', '--', '-', '--', '-', '--', '-', '--', '-', '--']
colors = ['k', 'r', 'y', 'm', 'c', 'g', '#6b420c']

P.xticks(fontsize=tickFontSize)
P.yticks(fontsize=tickFontSize)

div = 100.0
for i in range(len(filenames)):
  f = open(filenames[i], 'r')
  legendTitles = f.readline().strip().split(",")
  logInterval = float(legendTitles[0].split("=")[1])
  len_x = float(legendTitles[1].split("=")[1])
  len_y = float(legendTitles[2].split("=")[1])
  len_z = float(legendTitles[3].split("=")[1])
  voxelRadius = float(legendTitles[4].split("=")[1])
  scale = voxelRadius*2
  speciesNames = []
  speciesRadii = []
  for j in range(len(legendTitles)-5):
    speciesNames.append(legendTitles[j+5].split("=")[0])
    speciesRadii.append(float(legendTitles[j+5].split("=")[1]))
  speciesSize = len(speciesNames)

  data = np.genfromtxt(filenames[i], delimiter=',', skip_header=1).T

  colSize = len(data)-1
  for j in range(colSize):
    P.plot(data[0], data[j+1]/div, ls=lines[i], label=titles[i], linewidth=3)


ax = P.gca()
ax.grid(color='b', linestyle='--')
ax.set_xlim(0, 1000)
#ax.yaxis.set_major_locator(MaxNLocator(14))
leg = P.legend(loc=0, labelspacing=0.2, handletextpad=0.2, fancybox=True, fontsize=18)
P.ylabel('Survival Probability', fontsize=20)
P.xlabel('Time (s)', fontsize=20)
P.ylim(ymax=1.02)
P.suptitle('1D Reaction Validation', fontsize=tickFontSize)
P.show()

