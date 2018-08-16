#!/usr/bin/env python
import sys
import numpy as np
import math

def run_set(outfile, name, V_list, N_list, T_list, R, D, REPEAT):
    import run_single
    outfile.write('%s = [\n' % name)
    for i in range(len(V_list)):
        if (i > 8):
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
Vv = [3e-17, ] * 11 #simulation volume in m^3
Nv = [100,300,1000,3000,10000,30000,100000,300000,1000000,3000000,10000000] #number of molecules
Tr = np.array([3580.3, 590.6, 76.7, 12.03, 1.579, 0.2355, 0.02638, 0.002638, 0.00.0002638, 1.11732e-5, 1.11732e-5])
Tx = np.array([2016.6, 2078.7, 2142.5, 2305.5, 2627.9, 3254.5, 5096.7, 7078.1, 15430.9, 11874.8, 159672.3])
Tv = 10000./Tx*Tr
REPEAT = 1

if __name__ == '__main__':
    mode = "egfrd"
    postfix = '_out'
    outfile = open(mode+postfix+'.py','w')
    dataname = mode+'_data'
    run_set(outfile, dataname, Vv, Nv, Tv, Rv, Dv, REPEAT)
    outfile.write('\n\n')

