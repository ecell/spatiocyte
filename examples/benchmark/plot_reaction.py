import os
import sys
import numpy
import csv
import math
from matplotlib import rc
from pylab import *
#uncomment the following to create valid eps (scribus) and svg (inkscape):
#rc('svg', embed_char_paths=True)
#rc('mathtext', fontset=r'stixsans')

#matplotlib.rcParams["mathtext.fontset"] = "stix"
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble'].append(r'\usepackage{amsmath}')
matplotlib.rcParams['text.latex.preamble'].append(r'\usepackage[helvet]{sfmath}')

labelFontSize = 12
legendFontSize = 12
lineFontSize = 12
markersize=8
matplotlib.rcParams.update({'font.size': labelFontSize})

path, file = os.path.split(os.path.abspath(__file__))
path = path+os.sep

fileNames = ['egfrd/egfrd.csv', 'egfrd/egfrd_small_r.csv','smoldyn/smoldyn.csv','smoldyn/smoldyn_small_dt.csv','spatiocyte/spatiocyte.csv','spatiocyte/spatiocyte_small_dt.csv','spatiocyte/ode.csv',]
legendTitles = ['eGFRD ($r=10\ \mathrm{nm},\ T=2914\ \mathrm{s}$)','eGFRD ($r=1\ \mathrm{nm},\ T=2272\ \mathrm{s}$)','Smoldyn ($\Delta t=1\ \mathrm{ms},\ T=21\ \mathrm{s}$)','Smoldyn ($\Delta t=67\ \mathrm{\mu s},\ T=302\ \mathrm{s}$)','Spatiocyte ($\Delta t=1\ \mathrm{ms},\ r=38.73\ \mathrm{nm},\ T=8\ \mathrm{s}$)','Spatiocyte ($\Delta t=67\ \mathrm{\mu s},\ r=10\ \mathrm{nm},\ T=238\ \mathrm{s}$)','Mass Action']
speciesList = ['E','S','ES','P']
lines = ['-','--','-','--','-','--','-']
colors = ['b', 'b', 'g', 'g', 'r', 'r','k']

for f in range(len(fileNames)):
  if (os.path.isfile(path+fileNames[f])):
    deli = ','
    if (f == 2 or f == 3):
      deli = ' '
    data = genfromtxt(path+fileNames[f], delimiter=deli, skip_header=1).T
    colSize = len(data)-1
    for i in range(colSize):
      if (i == 0):
        plot(data[0], data[i+1], ls=lines[f], color=colors[f], label=legendTitles[f], linewidth=1)
      else:
        plot(data[0], data[i+1], ls=lines[f], color=colors[f], linewidth=1)

annotate('ES', xy=(90, 0),  xycoords='data', xytext=(-29, 8), textcoords='offset points', color='k', size=lineFontSize)

annotate('E', xy=(90, 60),  xycoords='data', xytext=(-27, 15), textcoords='offset points', color='k', size=lineFontSize)

annotate('P', xy=(90, 230),  xycoords='data', xytext=(-27, 12), textcoords='offset points', color='k', size=lineFontSize)

annotate('S', xy=(90, 680),  xycoords='data', xytext=(-27, 12), textcoords='offset points', color='k', size=lineFontSize)

annotate(r'E + S $\overset{k_1}{\underset{k_2}\rightleftharpoons}$ ES $\overset{k_3}{\rightarrow}$ E + P', xy=(70, 780),  xycoords='data', xytext=(-29, 10), textcoords='offset points', color='k', size=lineFontSize)

ax = gca()
leg = legend(loc=(0.02,0.25), labelspacing=0.3, handletextpad=0.2)
for t in leg.get_texts():
  t.set_fontsize(legendFontSize)   
xticks(size=labelFontSize)
yticks(size=labelFontSize)
xlabel('time, $t$ (s)',size=labelFontSize)
ylabel('\# Molecules',size=labelFontSize)
xlim(0,100)
ylim(0,1000)
savefig('reaction.eps', format='eps', dpi=1000)
show()

