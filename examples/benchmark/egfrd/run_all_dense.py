#!/usr/bin/env python
import sys
import numpy
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
Vv = [3e-18, ] * 11 #simulation volume in m^3
Nv = [3000000,10000000] #number of molecules
#Tv = [max(1e-7, min(T, 1e-3 / math.pow(N, 2.0 / 3.0))) for N in Nv] #duration
#Tv = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,3e-7,2e-7,1e-7,1e-7,1e-7]
Tv = [3e-6,1.5e-6]
REPEAT = 1

if __name__ == '__main__':
    mode = "egfrd_dense"
    postfix = '_out'
    outfile = open(mode+postfix+'.py','w')
    dataname = mode+'_data'
    run_set(outfile, dataname, Vv, Nv, Tv, Rv, Dv, REPEAT)
    outfile.write('\n\n')

