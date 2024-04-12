import os
import sys
import time
import multiprocessing
import numpy as np
import argparse
from tqdm import tqdm
from amr_kitchen import HeaderData

# Argument parser
parser = argparse.ArgumentParser(
        description="A quick taste of the plotfile")

parser.add_argument(
        "plotfile", type=str,
        help="Path of the plotfile to filter")

parser.add_argument(
        "--variable", "-v", type=str,
        help="Variable to taste")

parser.add_argument(
        "--no_log", "-n", action='store_true')

args = parser.parse_args()

# Read the header data
hdr = HeaderData(args.plotfile)
N_FIELDS = len(hdr.fields)
FIELD_INDEX = hdr.fields[args.variable]

# Define the multiprocessing functions

def expand_array3d(arr, factor):
    """
    Data reading utility
    ----
    Expand lower resolution 2D array by [factor]
    to broadcast it to a higher level grid.
    This allows broadcasting lower resolution arrays to a higher 
    AMR level grid without altering the data.
    ----
    """
    return np.repeat(np.repeat(np.repeat(arr, factor, axis=0),
                               factor, axis=1),
                     factor, axis=2)

def indexes_from_header(header):
    h = header.decode('ascii')
    start, stop, _, nfields = h.split()[-4:]
    nfields = int(nfields)
    start = np.array(start.split('(')[-1].replace(')','').split(','), dtype=int)
    stop = np.array(stop.replace('(', '').replace(')', '').split(','), dtype=int)
    return [start, stop]

def readfieldfrombinfile(fname):
    indexes = []
    arrays = []
    with open(fname, 'rb') as bfile:
        while True:
            try:
                h = bfile.readline()
                idx = indexes_from_header(h)
                shape = idx[1] - idx[0] + 1
                tshape = [shape[0], shape[1], shape[2], N_FIELDS]
                bfile.seek(np.prod(shape)*FIELD_INDEX*8, 1)
                arr = np.fromfile(bfile, 'float64', np.prod(shape))
                arrays.append(arr.reshape(shape, order='F'))
                indexes.append(idx)
                remainder = np.prod(tshape)*8 - np.prod(shape)*FIELD_INDEX*8 - np.prod(shape)*8
                bfile.seek(remainder, 1)
            except Exception as e:
                break
    return indexes, arrays

# Covering grid array
data = np.zeros(hdr.grid_sizes[hdr.limit_level], dtype='double')

start = time.time()
for lv in range(hdr.limit_level + 1):
    # Order binary files by size
    binfiles = np.unique(hdr.cells[lv]['files'])
    sizes = np.array([os.path.getsize(f) for f in binfiles])
    read_order = np.flip(np.argsort(sizes))

    if not args.no_log:
        print('Level', lv)
    factor = 2**(hdr.limit_level - lv)
    with multiprocessing.Pool() as pool:
        if not args.no_log:
            prog = tqdm(total=len(binfiles))
        for res in pool.imap_unordered(readfieldfrombinfile, binfiles[read_order]):
            if not args.no_log:
                prog.update(1)
            for idx, arr in zip(res[0], res[1]):
                data[factor * idx[0][0]:(idx[1][0]+1) * factor,
                     factor * idx[0][1]:(idx[1][1]+1) * factor,
                     factor * idx[0][2]:(idx[1][2]+1) * factor] = expand_array3d(arr, factor)

data_sum = np.sum(data)
integral = np.prod(hdr.dx[hdr.limit_level]) * data_sum
print(f"Integral of {args.variable} = {integral}")


