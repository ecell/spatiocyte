#!/usr/bin/env python
import sys
import math
import numpy as np

def run_set(outfile, name, V_list, N_list, T_list, R, D, REPEAT):
    import run_single
    outfile.write('%s = [\n' % name)
    for i in range(len(V_list)):
        if(i == i):
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
Vv = [3e-18, ] * 11 #simulation volume in m^3
Nv = [100,300,1000,3000,10000,30000,100000,300000,1000000,3000000,10000000] #number of molecules
Tr = np.array([1.43, 0.358, 0.0841, 0.0266, 0.00780, 0.00251, 1.386e-3, 4.539e-4, 1.269e-4, 3.757e-5, 1.798e-5])
Tx = np.array([10.27, 7.54, 5.91, 5.61, 5.48, 5.29, 9.7, 9.46, 8.78, 7.89, 11.68])
Tv = 300.0/Tx*Tr
REPEAT = 1

if __name__ == '__main__':
    mode = "fastbd"
    postfix = '_out'
    outfile = open(mode+postfix+'.py','w'); 
    dataname = mode+'_data'
    run_set(outfile, dataname, Vv, Nv, Tv, Rv, Dv, REPEAT);
    outfile.write('\n\n')
