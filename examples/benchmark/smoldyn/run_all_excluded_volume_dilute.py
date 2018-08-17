#!/usr/bin/env python
import sys
import math
import numpy as np

def run_set(outfile, name, V_list, N_list, T_list, R, D, M, REPEAT):
    import run_single
    outfile.write('%s = [\n' % name)
    for i in range(len(V_list)):
        if(i == i):
            outfile.write('# T=%g, N=%g, V=%g\n' % 
                          (T_list[i], N_list[i], V_list[i]))
            run_times = []
            est_times = []
            for c in range(REPEAT):
                run_time = run_single.run_single(T_list[i], V_list[i], N_list[i], R, D, M)
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
Tx = np.array([60.28, 17.85, 5.347, 1.555, 0.2154, 0.05773, 0.01369, 0.002003, 0.0004342, 0.0001159, 5.8783e-5]) 
Tr = np.array([690.4, 560.1, 715.1, 777.7, 708.2, 749.9, 746.0, 675.5, 681.1, 596.2, 1031.7])
Tv = 4000.0/Tr*Tx
Mv = "diffusion_excluded_volume.txt"
REPEAT = 1

if __name__ == '__main__':
    mode = "smoldyn_excluded_volume_dilute"
    postfix = '_out'
    outfile = open(mode+postfix+'.py','w'); 
    dataname = mode+'_data'
    run_set(outfile, dataname, Vv, Nv, Tv, Rv, Dv, Mv, REPEAT);
    outfile.write('\n\n')
