import os
import sys
import numpy
import csv
import math
from matplotlib import rc
from pylab import *
from collections import OrderedDict
#uncomment the following to create valid eps (scribus) and svg (inkscape):
#rc('svg', embed_char_paths=True)
#rc('mathtext', fontset=r'stixsans')

#matplotlib.rcParams["mathtext.fontset"] = "stix"
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble'].append(r'\usepackage{amsmath}')
matplotlib.rcParams['text.latex.preamble'].append(r'\usepackage[helvet]{sfmath}')

labelFontSize = 23
legendFontSize = 23
lineFontSize = 23
markersize=16
linewidth=3
mew=2
matplotlib.rcParams.update({'font.size': labelFontSize})
matplotlib.rcParams['figure.figsize'] = 9.1,7

path, file = os.path.split(os.path.abspath(__file__))
path = path+os.sep


linestyles = OrderedDict(
    [('solid',               (0, ())),
     ('loosely dotted',      (0, (1, 10))),
     ('dotted',              (0, (1, 5))),
     ('densely dotted',      (0, (1, 1))),

     ('loosely dashed',      (0, (5, 10))),
     ('dashed',              (0, (5, 5))),
     ('densely dashed',      (0, (5, 1))),

     ('loosely dashdotted',  (0, (3, 10, 1, 10))),
     ('dashdotted',          (0, (3, 5, 1, 5))),
     ('densely dashdotted',  (0, (3, 1, 1, 1))),

     ('loosely dashdotdotted', (0, (3, 10, 1, 10, 1, 10))),
     ('dashdotdotted',         (0, (3, 5, 1, 5, 1, 5))),
     ('densely dashdotdotted', (0, (3, 1, 1, 1, 1, 1)))])

fileNames = ['spatiocyte/spatiocyte.csv','spatiocyte/spatiocyte_small_dt.csv','smoldyn/smoldyn.csv','smoldyn/smoldyn_small_dt.csv','egfrd/egfrd.csv', 'egfrd/egfrd_small_r.csv','spatiocyte/ode.csv',]
legendTitles = ['Spatiocyte ($\Delta t=1\ \mathrm{ms},\ r=38.73\ \mathrm{nm},\ T=13\ \mathrm{s}$)','Spatiocyte ($\Delta t=67\ \mathrm{\mu s},\ r=10\ \mathrm{nm},\ T=276\ \mathrm{s}$)','Smoldyn ($\Delta t=1\ \mathrm{ms},\ T=20\ \mathrm{s}$)','Smoldyn ($\Delta t=67\ \mathrm{\mu s},\ T=298\ \mathrm{s}$)','eGFRD ($r=10\ \mathrm{nm},\ T=3246\ \mathrm{s}$)','eGFRD ($r=1\ \mathrm{nm},\ T=2412\ \mathrm{s}$)','Mass Action']
speciesList = ['E','S','ES','P']
lines = ['-','--','-','--','-','--','-']
#lines = [linestyles['solid'],linestyles['dotted'],linestyles['solid'],linestyles['dotted'],linestyles['solid'],linestyles['dotted'],linestyles['solid']]
#lines = [linestyles['densely dashed'],linestyles['densely dotted'],linestyles['densely dashdotdotted'],linestyles['dashdotted'],linestyles['dashed'],linestyles['dotted'],linestyles['solid']]
colors = ['r', 'r', 'g', 'g', 'b', 'b','k']
opacity = [0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 1]
#colors = ['m', 'm', 'g', 'g', 'orange', 'orange','k']

for f in range(len(fileNames)):
  if (os.path.isfile(path+fileNames[f])):
    deli = ','
    if (f == 2 or f == 3):
      deli = ' '
    data = genfromtxt(path+fileNames[f], delimiter=deli, skip_header=1).T
    colSize = len(data)-1
    for i in range(colSize):
      if (i == 0):
        plot(data[0], data[i+1], ls=lines[f], color=colors[f], label=legendTitles[f], linewidth=1.5, alpha=opacity[f])
      else:
        plot(data[0], data[i+1], ls=lines[f], color=colors[f], linewidth=1.5, alpha=opacity[f])

annotate('ES', xy=(95, 0),  xycoords='data', xytext=(-29, 8), textcoords='offset points', color='k', size=lineFontSize)

annotate('E', xy=(95, 50),  xycoords='data', xytext=(-27, 15), textcoords='offset points', color='k', size=lineFontSize)

annotate('P', xy=(95, 235),  xycoords='data', xytext=(-27, 12), textcoords='offset points', color='k', size=lineFontSize)

annotate('S', xy=(95, 650),  xycoords='data', xytext=(-27, 12), textcoords='offset points', color='k', size=lineFontSize)

annotate(r'E + S $\overset{k_1}{\underset{k_2}\rightleftharpoons}$ ES $\overset{k_3}{\rightarrow}$ E + P', xy=(50, 780),  xycoords='data', xytext=(-29, 10), textcoords='offset points', color='k', size=lineFontSize)

ax = gca()
leg = legend(loc=(0.02,0.17), labelspacing=0.12, handlelength=1.0, handletextpad=0.3, frameon=False)
for t in leg.get_texts():
  t.set_fontsize(legendFontSize)   
xticks(size=labelFontSize)
yticks(size=labelFontSize)

ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.tick_params(axis='both',which='both',direction='in',length=10,width=2)
ax.tick_params(axis='both',which='major',length=10,width=2)
for axis in ['top','bottom','left','right']:
     ax.spines[axis].set_linewidth(2)

xlabel('Time, $t$ (s)',size=labelFontSize)
ylabel('Molecules',size=labelFontSize)
xlim(0,100)
ylim(0,1000)
tight_layout(pad=0)
savefig('reaction.pdf', format='pdf', dpi=1000)#, bbox_inches='tight')
show()

