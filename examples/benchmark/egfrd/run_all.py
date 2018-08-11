#!/usr/bin/env python
import sys
import numpy as np
import math

def run_set(outfile, name, V_list, N_list, T_list, R, D, REPEAT):
    import run_single
    outfile.write('%s = [\n' % name)
    for i in range(len(V_list)):
        outfile.write('# T=%g, N=%g, V=%g\n' % 
                      (T_list[i], N_list[i], V_list[i]))
        run_times = []
        est_times = []
        for c in range(REPEAT):
            run_time = run_single.run_single(T_list[i], V_list[i], N_list[i], R, D)
            est_time = run_time * (T / T_list[i])
            run_times.append(run_time)
            est_times.append(est_time)
        outfile.write('# run_times = %s\n' % str(run_times))
        outfile.write('%s,\n' % str(est_times))
        outfile.flush()
    outfile.write(']\n')

T = 10.
Rv = 2.5e-9 #voxel radius
Dv = 1e-12 #diffusion coefficient
Vv = [3e-15, ] * 11 #simulation volume in m^3
Nv = [100,300,1000,3000,10000,30000,100000,300000,1000000,3000000,10000000] #number of molecules
Tr = np.array([333.33, 51.33, 6.831, 1.063, 0.134, 0.02049, 0.002361, 0.0003059, 3.276e-5, 2.533e-6, 1.072e-7])
Tx = np.array([18.62, 17.38, 17.81, 17.66, 16.97, 17.40, 17.90, 21.95, 52.57, 58.873, 344.21])
Tv = 10.0/Tx*Tr
REPEAT = 1

if __name__ == '__main__':
    mode = "egfrd"
    postfix = '_out'
    outfile = open(mode+postfix+'.py','w')
    dataname = mode+'_data'
    run_set(outfile, dataname, Vv, Nv, Tv, Rv, Dv, REPEAT)
    outfile.write('\n\n')

