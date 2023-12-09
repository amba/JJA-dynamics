#!/usr/bin/env python3

from multiprocessing import Pool, cpu_count
import subprocess
import numpy as np

def run(params):
    subprocess.run(["./ground-state-annealer", *params]) 
    return 0

N = 200
n = 20000000

f_vals = np.arange(0.33, 0.43, 0.01)
print("f_vals = ", f_vals)

runs = []
for f in f_vals:
    runs.append(['-N', str(N), '-f', str(f), '-n', str(n)])

          

if __name__ == '__main__':
    num_cores = cpu_count()
    print("num cores = ", num_cores)
    with Pool(num_cores) as p:
        print(p.map(run, runs))
        
