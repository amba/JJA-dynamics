#!/usr/bin/env python3 

import numpy as np
import sys
import glob
import argparse
import io
import os.path
import sys
import matplotlib.pyplot as plt
import scipy.signal


if np.__version__ < '1.14.1':
    sys.exit("numpy version " + np.__version__ + " is too old")
    
def open_3d_file(file):
    fh = open(file, 'r')
    header = fh.readline().rstrip()
    contents = fh.read().rstrip()
    
    list_of_blocks = contents.split("\n\n")
    num_blocks = len(list_of_blocks)
    arrays = []
    for i, block in enumerate(list_of_blocks):
#        print("reading block %d / %d" % (i, num_blocks))
        arrays.append(np.genfromtxt(io.StringIO(block)))
    first_shape = arrays[0].shape
    for i in range(len(arrays)-1, -1, -1):
        shape = arrays[i].shape
        if shape != first_shape:
            print("block ", i, " with first line", arrays[i][0], " does not match :", shape, " != ", first_shape)
            del arrays[i]
    return np.stack(arrays), header


def save_3d_file(output_file, data, header):
    fh = open(output_file, 'w')
    fh.write(header + "\n")
    for block in data:
        np.savetxt(fh, block, fmt="%.17g", delimiter="\t")
        fh.write("\n")
    fh.close()


parser = argparse.ArgumentParser()

parser.add_argument('-o', '--output', help="basename for output data file")
parser.add_argument('-f', '--force', action="store_true", help="overwrite existing files")

args = parser.parse_args()
output_ix = glob.glob('*_Ix.dat')
output_iy = glob.glob('*_Iy.dat')
print(output_ix, output_iy)
data_x, header = open_3d_file(output_ix[0])
data_y, header = open_3d_file(output_iy[0])
I_x = np.swapaxes(data_x, 0, 1)
I_y = np.swapaxes(data_y, 0, 1)

I_x = I_x[:,:,2]
I_y = I_y[:,:,2]

print("shape I_x = ", I_x.shape)
print("shape I_y = ", I_y.shape)
# first index: i (x/col)
# second index: j (y/row)


Nx = I_x.shape[0] + 1
Ny = I_x.shape[1] 
print("Nx = ", Nx, ", Ny = ", Ny)

def plot_flux():
    plt.clf()
    m = np.zeros((Nx - 1, Ny - 1))
    for i in range(Nx - 1):
        for j in range(Ny - 1):
            m[i,j] = I_x[i, j] - I_x[i,j+1]  + I_y[i+1,j] - I_y[i,j]
            
    m /= 4
    m = np.flip(m, axis=1)
    m = np.swapaxes(m, 0, 1)

    plt.imshow(m, aspect='equal', cmap='gray')
    plt.colorbar(format="%.1f", label='flux')

plot_flux()
plt.show()
