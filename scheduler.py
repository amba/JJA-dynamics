#!/usr/bin/env python3

from multiprocessing import Pool, cpu_count
import subprocess
import numpy as np

def run(params):
    subprocess.run(["./ground-state-annealer", *params]) 
    return 0

N = 100
n = 20000000

f_vals = (7/20,)
T_start = 2
n_vals = (100000, 500000, 1000000, 5000000, 10000000, 50000000, 100000000) 
print("f_vals = ", f_vals)

runs = []
for f in f_vals:
    for n in n_vals:
        runs.append(['-N', str(N), '-f', str(f), '-n', str(n), '-t', str(T_start)])

          

if __name__ == '__main__':
    num_cores = cpu_count()
    print("num cores = ", num_cores)
    with Pool(num_cores) as p:
        print(p.map(run, runs))
        
