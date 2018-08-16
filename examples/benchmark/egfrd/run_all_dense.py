#!/usr/bin/env python
import sys
import numpy as np
import math

def run_set(outfile, name, V_list, N_list, T_list, R, D, REPEAT):
    import run_single
    outfile.write('%s = [\n' % name)
    for i in range(len(V_list)):
        if (i == 10):
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
Tr = np.array([500.85, 77.26, 9.33, 1.14, 0.0921, 0.00684, 0.000332, 1.8443e-5, 5.0517e-6, 4.71252e-7, 4.26392e-8])
Tx = np.array([1345.2, 1365.1, 1427.5, 1397.8, 1466.6, 1629.1, 1747.5, 1575.6, 5711.4, 5722.8, 9006.8])
Tv = 10000.0/Tx*Tr
REPEAT = 1

if __name__ == '__main__':
    mode = "egfrd_dense"
    postfix = '_out'
    outfile = open(mode+postfix+'.py','w')
    dataname = mode+'_data'
    run_set(outfile, dataname, Vv, Nv, Tv, Rv, Dv, REPEAT)
    outfile.write('\n\n')

