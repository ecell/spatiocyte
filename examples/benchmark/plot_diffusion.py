#!/usr/bin/env/python

import sys
import os
import numpy
import imp
import scipy.io
from matplotlib.pylab import *
#uncomment the following to create valid eps (scribus) and svg (inkscape):
#rc('svg', embed_char_paths=True)
#rc('mathtext', fontset=r'stixsans')

#matplotlib.rcParams["mathtext.fontset"] = "stix"
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble'].append(r'\usepackage{amsmath}')
matplotlib.rcParams['text.latex.preamble'].append(r'\usepackage{upgreek}')
matplotlib.rcParams['text.latex.preamble'].append(r'\usepackage[helvet]{sfmath}')

labelFontSize = 15
legendFontSize = 14
lineFontSize = 15

def plot_data(N, T, fmt):
    T = numpy.array(T)
    mean = T.mean(1)
    std_err = T.std()/math.sqrt(len(T))
    #errorbar(N, mean, yerr=std_err, fmt=fmt)
    print(N, mean)
    loglog(N, mean, fmt)

axes([.12,.14,.86,.83])
path, file = os.path.split(os.path.abspath(__file__))
path = path+os.sep

filename = path+'spatiocyte/run_all.py'
imp.load_source('run_all', filename)
from run_all import Nv
X = numpy.array(Nv)

filename = path+'egfrd_ecell4/dense_back.py'
if (os.path.isfile(filename)):
  imp.load_source('egfrd_dense_out', filename)
  from egfrd_dense_out import *
  plot(Nv, egfrd_dense_data,'b^', fillstyle='none', label=r'eGFRD ($V=3\ \mathrm{\mu m}^{3}$)')
  loglog(X, 2e-2*X**(5.0/3.0), 'b--')
  annotate(r'$T\propto N^{\mathsf{\frac{5}{3}}}$', xy=(X[2], egfrd_dense_data[2][0]),  xycoords='data', xytext=(-20, 15), textcoords='offset points', color='b', size=14)

filename = path+'egfrd_ecell4/egfrd_out.py'
if (os.path.isfile(filename)):
  imp.load_source('egfrd_out', filename)
  from egfrd_out import *
  plot(Nv, egfrd_data,'bv', fillstyle='none', label=r'eGFRD ($V=3000\ \mathrm{\mu m}^{3}$)')
  loglog(X, 2.05e-4*X**(5.0/3.0), 'b--')
  annotate(r'$T\propto N^{\mathsf{\frac{5}{3}}}$', xy=(X[3], egfrd_data[3][0]),  xycoords='data', xytext=(-15, -25), textcoords='offset points', color='b', size=14)

filename = path+'smoldyn/smoldyn_out.py'
if (os.path.isfile(filename)):
  imp.load_source('smoldyn_out', filename)
  from smoldyn_out import *
  plot(Nv, smoldyn_data,'gs', fillstyle='none', label=r'Smoldyn ($V=3\ \mathrm{\mu m}^{3}$)')
  loglog(X, 0.4*X, 'g--')
  #annotate(r'$T\propto N$', xy=(X[1], smoldyn_data[1][0]),  xycoords='data', xytext=(-25, -23), textcoords='offset points', color='g', size=14)

filename = path+'smoldyn/smoldyn_excluded_volume_out.py'
if (os.path.isfile(filename)):
  imp.load_source('smoldyn_excluded_volume_out', filename)
  from smoldyn_excluded_volume_out import *
  plot(Nv, smoldyn_excluded_volume_data,'go', fillstyle='none', label=r'Smoldyn excluded volume ($V=3\ \mathrm{\mu m}^{3}$)')
  loglog(X, 0.4*X, 'g--')
  #annotate(r'$T\propto N$', xy=(X[1], smoldyn_excluded_volume_data[1][0]),  xycoords='data', xytext=(-25, -23), textcoords='offset points', color='g', size=14)

filename = path+'smoldyn/smoldyn_dillute_out.py'
if (os.path.isfile(filename)):
  imp.load_source('smoldyn_dillute_out', filename)
  from smoldyn_dillute_out import *
  plot(Nv, smoldyn_dillute_data,'gv', label=r'Smoldyn ($V=30\ \mathrm{\mu m}^{3}$)')

filename = path+'spatiocyte/spatiocyte_out.py'
if (os.path.isfile(filename)):
  imp.load_source('spatiocyte_out', filename)
  from spatiocyte_out import *
  plot(Nv, spatiocyte_data,'rx', fillstyle='none', label=r'Spatiocyte ($V=3\ \mathrm{\mu m}^{3}$)')
  loglog(X, 0.8*X, 'r--')
  #annotate(r'$T\propto N$', xy=(X[5], spatiocyte_data[5][0]),  xycoords='data', xytext=(-15, 12), textcoords='offset points', color='r', size=14)

filename = path+'spatiocyte/spatiocyte_dillute_out.py'
if (os.path.isfile(filename)):
  imp.load_source('spatiocyte_dillute_out', filename)
  from spatiocyte_dillute_out import *
  plot(Nv, spatiocyte_dillute_data,'r+', fillstyle='none', label=r'Spatiocyte ($V=30\ \mathrm{\mu m}^{3}$)')


annotate('B', xy=(0, 1),  xycoords='axes fraction', xytext=(-60,-20), textcoords='offset points', color='k', size=30)

leg = legend(loc='upper left', labelspacing=0.3, handletextpad=0.2)
for t in leg.get_texts():
  t.set_fontsize(legendFontSize) 
#legend(loc='upper left', labelspacing=0.2, handletextpad=0.2, fancybox=True)

xlabel('$N$ (\# Molecules)', size=17)
ylabel('Run time, $T$', size=17)

Y = numpy.array([60,3600,3600*24,3600*24*30, 3600*24*30*12])

xlim(X[0]*0.9,X[len(X)-1]*1.1)

xticks(size=17)
yticks(Y, ['min', 'hour', 'day', 'mon', 'year'], size=14)

savefig('diffusion.eps', format='eps', dpi=1000)
show()
