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
Vv = [3e-18, ] * 11 #simulation volume in m^3
Nv = [100,300,1000,3000,10000,30000,100000,300000,1000000,3000000,10000000] #number of molecules
Tx = np.array([max(1e-5, min(T, 20e0 / math.pow(N, 2.0 / 3.0))) for N in Nv]) #duration
Tr = np.array([0.77, 1.25, 1.87, 3.09, 10.00, 17.94, 33.91, 111.42, 230.33, 414.89, 733.01])
Tv = 600.0/Tr*Tx
Mv = "diffusion.txt"
REPEAT = 1

if __name__ == '__main__':
    mode = "smoldyn"
    postfix = '_out'
    outfile = open(mode+postfix+'.py','w'); 
    dataname = mode+'_data'
    run_set(outfile, dataname, Vv, Nv, Tv, Rv, Dv, Mv, REPEAT);
    outfile.write('\n\n')
