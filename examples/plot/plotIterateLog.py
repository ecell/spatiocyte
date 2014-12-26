import numpy as np
import pylab as P

labelFontSize = 14
tickFontSize = 14
legendFontSize = 14
lineFontSize = 14

fileNames = ["IterateLog.csv"]
legendTitles = []
lines = ['-', '-', '-', '-']
colors = ['y', 'r', 'b', 'm', 'c', 'g']

P.xticks(fontsize=tickFontSize)
P.yticks(fontsize=tickFontSize)

f = open(fileNames[0], 'r')
legendTitles = f.readline().strip().split(",")
logInterval = float(legendTitles[0].split("=")[1])
len_x = float(legendTitles[1].split("=")[1])
len_y = float(legendTitles[2].split("=")[1])
len_z = float(legendTitles[3].split("=")[1])
voxelRadius = float(legendTitles[4].split("=")[1])
scale = voxelRadius*2
speciesNames = []
speciesRadii = []
for i in range(len(legendTitles)-5):
  speciesNames.append(legendTitles[i+5].split("=")[0])
  speciesRadii.append(float(legendTitles[i+5].split("=")[1]))
speciesSize = len(speciesNames)

data = np.genfromtxt(fileNames[0], delimiter=',', skip_header=1).T

colSize = len(data)-1
for i in range(colSize):
  P.plot(data[0], data[i+1], ls=lines[i], color=colors[i], label=speciesNames[i], linewidth=1)

ax = P.gca()
ax.grid(color='b', linestyle='--')
#ax.yaxis.set_major_locator(MaxNLocator(14))
leg = P.legend(loc=0, labelspacing=0.2, handletextpad=0.2, fancybox=True)
for t in leg.get_texts():
  t.set_fontsize(legendFontSize)   
frame = leg.get_frame()
frame.set_linewidth(None)
frame.set_facecolor('0.95')
frame.set_edgecolor('0.75')
P.ylabel('# Molecules')
P.xlabel('Time (s)')
P.show()

