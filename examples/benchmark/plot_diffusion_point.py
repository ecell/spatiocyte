import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
import imp
import scipy.constants
from matplotlib.ticker import ScalarFormatter
import matplotlib.ticker as ticker

labelFontSize = 15
legendFontSize = 14
lineFontSize = 15
markersize=8
matplotlib.rcParams.update({'font.size': labelFontSize})

path, file = os.path.split(os.path.abspath(__file__))
path = path+os.sep
filename = path+'spatiocyte/run_all.py'
imp.load_source('run_all', filename)
from run_all import Nv
X = np.array(Nv)


fig, ax1 = plt.subplots()
ax2 = ax1.twiny()
lines = []

filename = path+'spatiocyte/spatiocyte_point_out.py'
if (os.path.isfile(filename)):
  imp.load_source('spatiocyte_point_out', filename)
  from spatiocyte_point_out import *
  lines += ax1.plot(Nv, spatiocyte_point_data,'r+', fillstyle='none', markersize=markersize, label=r'Spatiocyte ($V=3\ \mathrm{\mu m}^{3}$)')
  ax1.annotate(r'$T\propto N$', xy=(X[9], spatiocyte_point_data[9][0]),  xycoords='data', xytext=(-13, -28), textcoords='offset points', color='r', size=14)

filename = path+'spatiocyte/spatiocyte_point_dilute_out.py'
if (os.path.isfile(filename)):
  imp.load_source('spatiocyte_point_dilute_out', filename)
  from spatiocyte_point_dilute_out import *
  lines += ax1.plot(Nv, spatiocyte_point_dilute_data,'rD', fillstyle='none', markersize=markersize, label=r'Spatiocyte ($V=30\ \mathrm{\mu m}^{3}$)')
  ax1.loglog(X, 0.28*X, 'r--', linewidth=0.5)

filename = path+'fastbd/fastbd_out.py'
if (os.path.isfile(filename)):
  imp.load_source('fastbd_out', filename)
  from fastbd_out import *
  lines += ax1.plot(Nv, fastbd_data,'b1', fillstyle='none', markersize=markersize, label=r'Fast BD ($V=3\ \mathrm{\mu m}^{3}$)')

filename = path+'fastbd/fastbd_dilute_out.py'
if (os.path.isfile(filename)):
  imp.load_source('fastbd_dilute_out', filename)
  from fastbd_dilute_out import *
  lines += ax1.plot(Nv, fastbd_dilute_data,'b^', fillstyle='none', markersize=markersize, label=r'Fast BD ($V=30\ \mathrm{\mu m}^{3}$)')
  #ax1.loglog(X, 0.72*X, 'k--', linewidth=0.5)

filename = path+'smoldyn/smoldyn_out.py'
if (os.path.isfile(filename)):
  imp.load_source('smoldyn_out', filename)
  from smoldyn_out import *
  lines += ax1.plot(Nv, smoldyn_data,'gx', fillstyle='none', markersize=markersize, label=r'Smoldyn ($V=3\ \mathrm{\mu m}^{3}$)')

filename = path+'smoldyn/smoldyn_dilute_out.py'
if (os.path.isfile(filename)):
  imp.load_source('smoldyn_dilute_out', filename)
  from smoldyn_dilute_out import *
  lines += ax1.plot(Nv, smoldyn_dilute_data,'gs', fillstyle='none', markersize=markersize, label=r'Smoldyn ($V=30\ \mathrm{\mu m}^{3}$)')


ax2.axvline(10, color='purple', linewidth=0.5)
ax2.axvline(0.1, color='purple', linewidth=0.5)

def n_to_c(n):
  volume = 30e-18 #in m^3
  volume = volume*1e+3 #in liter
  conc = n*1.0/(scipy.constants.N_A*volume)*1e+6
  return conc

def log_format(y,pos):
  # Find the number of decimal places required
  decimalplaces = int(np.maximum(-np.log10(y),0))     # =0 for numbers >=1
  # Insert that number into a format string
  formatstring = '{{:.{:1d}f}}'.format(decimalplaces)
  # Return the formatted tick label
  return formatstring.format(y)


Y = np.array([60,3600,3600*24,3600*24*30, 3600*24*30*12])
plt.yticks(Y, ['min', 'hour', 'day', 'mon', 'year'])

ax1.set_xlabel('$N$ (# Molecules)', size=labelFontSize)
ax1.set_ylabel('Run time, $T$', size=labelFontSize)
xmin, xmax = ax1.get_xlim()
ax2.set_xlim(n_to_c(xmin), n_to_c(xmax))
ax2.set_xscale("log")
ax2.xaxis.set_major_formatter(ticker.FuncFormatter(log_format))
ax2.set_xlabel(r"Concentration at $V=30\ \mathrm{\mu m}^{3}$ ($\mu$M)")
ax1.yaxis.grid()
ax1.xaxis.grid()
ax1.grid(color='k', linewidth=0.1)

labels = [l.get_label() for l in lines]
leg = ax2.legend(lines, labels, loc='upper left', labelspacing=0.3, handletextpad=0.2)
for t in leg.get_texts():
  t.set_fontsize(legendFontSize) 

plt.savefig('diffusion.eps', format='eps', dpi=1000, bbox_inches='tight')
plt.show()
