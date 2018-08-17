import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
import imp
import scipy.constants
from matplotlib.ticker import ScalarFormatter
import matplotlib.ticker as ticker

labelFontSize = 12
legendFontSize = 12
lineFontSize = 12
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

filename = path+'spatiocyte/spatiocyte_out.py'
if (os.path.isfile(filename)):
  imp.load_source('spatiocyte_out', filename)
  from spatiocyte_out import *
  lines += ax1.plot(Nv, spatiocyte_data,'r+', fillstyle='none', markersize=markersize, label=r'Spatiocyte ($V=3\ \mathrm{\mu m}^{3}$)')
  ax1.loglog(X, 0.7*X, 'r', linewidth=0.5)
  ax1.annotate(r'$T\propto N$', xy=(X[9], spatiocyte_data[9][0]),  xycoords='data', xytext=(-13, -23), textcoords='offset points', color='r', size=14)

filename = path+'spatiocyte/spatiocyte_dilute_out.py'
if (os.path.isfile(filename)):
  imp.load_source('spatiocyte_dilute_out', filename)
  from spatiocyte_dilute_out import *
  lines += ax1.plot(Nv, spatiocyte_dilute_data,'rD', fillstyle='none', markersize=markersize, label=r'Spatiocyte ($V=30\ \mathrm{\mu m}^{3}$)')


filename = path+'smoldyn/back_smoldyn_excluded_volume_out.py'
if (os.path.isfile(filename)):
  imp.load_source('smoldyn_excluded_volume_out', filename)
  from smoldyn_excluded_volume_out import *
  lines += ax1.plot(Nv, smoldyn_excluded_volume_data,'gx', fillstyle='none', markersize=markersize, label=r'Smoldyn ($V=3\ \mathrm{\mu m}^{3}$)')

filename = path+'smoldyn/back_smoldyn_excluded_volume_dilute_out.py'
if (os.path.isfile(filename)):
  imp.load_source('smoldyn_excluded_volume_dilute_out', filename)
  from smoldyn_excluded_volume_dilute_out import *
  lines += ax1.plot(Nv, smoldyn_excluded_volume_dilute_data,'gs', fillstyle='none', markersize=markersize, label=r'Smoldyn ($V=30\ \mathrm{\mu m}^{3}$)')

filename = path+'egfrd/egfrd_dense_out.py'
if (os.path.isfile(filename)):
  imp.load_source('egfrd_dense_out', filename)
  from egfrd_dense_out import *
  lines += ax1.plot(Nv, egfrd_dense_data,'b1', fillstyle='none', markersize=markersize, label=r'eGFRD ($V=3\ \mathrm{\mu m}^{3}$)')
  #ax1.loglog(X, 2e-2*X**(5.0/3.0), 'b--', linewidth=0.5)

filename = path+'egfrd/egfrd_out.py'
if (os.path.isfile(filename)):
  imp.load_source('egfrd_out', filename)
  from egfrd_out import *
  lines += ax1.plot(Nv, egfrd_data,'bo', fillstyle='none', markersize=markersize, label=r'eGFRD ($V=30\ \mathrm{\mu m}^{3}$)')
  #ax1.loglog(X, 4e-3*X**(5.0/3.0), 'b--', linewidth=0.5)



ax2.axvline(10, color='brown', linestyle="--", linewidth=0.8)
ax2.axvline(0.1, color='brown', linestyle="--", linewidth=0.8)

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
ax1.grid(color='lightgray', linewidth=0.1)

labels = [l.get_label() for l in lines]
leg = ax2.legend(lines, labels, loc='upper left', labelspacing=0.1, handletextpad=0.1)
for t in leg.get_texts():
  t.set_fontsize(legendFontSize) 
plt.tight_layout()
plt.savefig('diffusion.pdf', format='pdf', dpi=1000, bbox_inches='tight')
plt.show()
