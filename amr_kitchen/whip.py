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
        description="3D covering grids")

parser.add_argument(
        "plotfile", type=str,
        help="Path of the plotfile to filter")

parser.add_argument(
        "--variable", "-v", type=str,
        help="Variable to load")

args = parser.parse_args()

# Read the header data
hdr = HeaderData(args.plotfile)
if hdr.ndims < 3:
    raise ValueError("Use mandoline to create 2D covering grids")

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
    """
    Reads the cell header in the binary file to get the
    shape of the array to read in memory
    -
    This is a bit faster than using the offsets in the cell header
    """
    h = header.decode('ascii')
    start, stop, _, nfields = h.split()[-4:]
    nfields = int(nfields)
    start = np.array(start.split('(')[-1].replace(')','').split(','), dtype=int)
    stop = np.array(stop.replace('(', '').replace(')', '').split(','), dtype=int)
    return [start, stop]

def readfieldfrombinfile(fname):
    """
    Read all arrays for a single field in a binary file
    -
    The index of the field in the plotfile and the number of fields
    are defined as global variables to use multi processing
    """
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
# TODO: This may fail after a very long time on some systems
# TODO: float32 could be used to use less memory for large plotfiles
#       but doubles should be used if the array fits in the RAM because
#       dynamic casting is expensive
try:
    data = np.zeros(hdr.grid_sizes[hdr.limit_level], dtype='double')
except MemoryError:
    raise RunTimeError("Sorry plotfile too big, you can use on-disk arrays")

start = time.time()
# For each AMR Level
for lv in range(hdr.limit_level + 1):
    # Order binary files by size so multiprocessing is a bit faster
    # (Large files take more time so we do them first)
    binfiles = np.unique(hdr.cells[lv]['files'])
    sizes = np.array([os.path.getsize(f) for f in binfiles])
    read_order = np.flip(np.argsort(sizes))

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

np.save(f"{args.variable}_cgrid_{args.plotfile.replace('plt', '')}",
        data)
