import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import os
import pylab as P

labelFontSize = 14
tickFontSize = 14
legendFontSize = 14
lineFontSize = 14

files = []
fileNames = ["CoordinateLog.csv"]
legendTitles = []
lines = ['-', '-', '-', '-']
colors = ['y', 'r', 'b', 'm', 'c', 'g']

max_frames = 600
fig = P.figure(figsize=(10,6))
ax = fig.add_subplot(111, projection='3d')
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width*0.7, box.height])
P.xticks(fontsize=tickFontSize)
P.yticks(fontsize=tickFontSize)

for i in range(len(fileNames)):
  f = open(fileNames[i], 'r')
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
  logCnt = 0
  lineCnt = 0
  markers = []
  frameCnt = 0
  for line in f:
    frameCnt = frameCnt + 1
    if frameCnt > max_frames*speciesSize:
      break
    coords = line.strip().split(",")
    time = float(coords[0])
    x = []
    y = []
    z = []
    for l in range((len(coords)-1)/3):
      x.append(float(coords[l*3+1])*scale)
      y.append(float(coords[l*3+2])*scale)
      z.append(float(coords[l*3+3])*scale)
    ax.plot(x, z, y, linewidth=0, color=colors[lineCnt], marker='.')
    markers.append(P.Rectangle((0, 0), 1, 1, fc=colors[lineCnt]))
    lineCnt = lineCnt + 1
    if lineCnt == speciesSize:
      ax.set_xlabel('X (m)')
      ax.set_ylabel('Y (m)')
      ax.set_zlabel('Z (m)')
      ax.set_title('t = %.2e s' %time)
      ax.set_xlim(0, len_x*scale+voxelRadius)
      ax.set_ylim(0, len_y*scale+voxelRadius)
      ax.set_zlim(0, len_z*scale+voxelRadius)
      ax.grid(False)
      leg = ax.legend(markers, speciesNames, bbox_to_anchor=(1.0,0.95),
          loc='upper left', labelspacing=0.2, handletextpad=0.2, fancybox=True)
      for t in leg.get_texts():
        t.set_fontsize(legendFontSize)   
      frame = leg.get_frame()
      frame.set_linewidth(None)
      frame.set_facecolor('0.95')
      frame.set_edgecolor('0.75')
      fileName = fileNames[i]+'.%03d.png'%logCnt
      print 'Saving frame', fileName
      fig.savefig(fileName)
      markers = []
      ax.cla()
      logCnt = logCnt + 1
      lineCnt = 0

os.system("avconv -i " + fileNames[0] + ".%03d.png -vcodec qtrle " + 
    fileNames[0] + ".mov")



