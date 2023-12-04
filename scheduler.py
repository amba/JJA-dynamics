#!/usr/bin/env python3

from multiprocessing import Pool, cpu_count
import subprocess

def run(params):
    subprocess.run(["./ground-state-annealer", *params]) 
    return 0

N = 200
n = 10000000

f_vals = (1/3 + (1/N)**2 *4, 1/3 + (1/N)**2 *8, 1/3 + (1/N)**2 * 16, 1/3 + (1/N)**2 * 32, 1/3 + (1/N)**2 * 64, 1/3 + (1/N)**2 * 128)

runs = []
for f in f_vals:
    runs.append(['-N', str(N), '-f', str(f), '-n', str(n)])

          

if __name__ == '__main__':
    num_cores = cpu_count()
    print("num cores = ", num_cores)
    with Pool(num_cores) as p:
        print(p.map(run, runs))
        
