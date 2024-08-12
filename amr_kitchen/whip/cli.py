import os
import sys
import argparse
from tqdm import tqdm
import multiprocessing
import numpy as np
from humanize import naturalsize
from amr_kitchen import PlotfileCooker
from amr_kitchen.utils import (expand_array3d,
                               indices_from_header)

def readfieldfrombinfile(args):
    """
    Read all arrays for a single field in a binary file
    """
    N_FIELDS = args['N_FIELDS']
    FIELD_INDEX = args['FIELD_INDEX']
    fname = args['fname']
    indexes = []
    arrays = []
    with open(fname, 'rb') as bfile:
        while True:
            try:
                h = bfile.readline()
                idx = indices_from_header(h)
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

def main():
    # Argument parser
    parser = argparse.ArgumentParser(
            description=("Creates a uniform 3D array of a plotfile field."
                         " Lower level data is broadcasted to the finest"
                         " AMR Level."))

    parser.add_argument(
            "--variable", "-v", type=str,
            help="Variable used to create the uniform grid")
    parser.add_argument(
            "--limit_level", "-l", type=int,
            help="Maximum AMR Level considered")
    parser.add_argument(
            "--outfile", "-o", type=str,
            help="Output file to override the default")
    parser.add_argument(
            "--dtype", "-d", type=str, default="float64",
            help=("Data type used (defaults to float64), but float32 or integer"
                  " types can be used to save space (this is slower)"))
    parser.add_argument(
            "--nochecks", "-y", action="store_true",
            help="Do not prompt if the required memory is acceptable")
    parser.add_argument(
            "plotfile", type=str,
            help="Path of the plotfile used to make the uniform grid")

    args = parser.parse_args()
    """
    Input arguments sanity check
    """
    if args.plotfile is None:
        raise ValueError("Must specify a plotfile to integrate")
    # Creating a PlotfileCooker instance
    pck = PlotfileCooker(args.plotfile)
    if pck.ndims < 3:
        print(("This tool is not supported for plotfiles with ndims < 3."
               " You can use mandoline to create a 2D uniform grid."))
        sys.exit()

    # Field infos
    N_FIELDS = len(pck.fields)
    FIELD_INDEX = pck.fields[args.variable]

    # Covering grid array size
    grid_shape = pck.grid_sizes[pck.limit_level]
    dtype_size = np.array([], dtype=args.dtype).itemsize
    total_size = np.prod(grid_shape) * dtype_size
    human_size = naturalsize(total_size)
    if args.nochecks:
        do_grid = 'y'
    else:
        do_grid = input((f"The uniform grid will require {human_size} in memory."
                          " go ahead with the computation (y/n)? "))
    if do_grid == 'n':
        print("You can use less precise data types to reduce the required memory.")
        sys.exit()
    elif do_grid == 'y':
        data = np.zeros(pck.grid_sizes[pck.limit_level], dtype=args.dtype)
        # For each AMR Level
        for lv in range(pck.limit_level + 1):
            # Order binary files by size so multiprocessing is a bit faster
            # (Large files take more time so we do them first)
            binfiles = np.unique(pck.cells[lv]['files'])
            sizes = np.array([os.path.getsize(f) for f in binfiles])
            read_order = np.flip(np.argsort(sizes))
            print('Level', lv)
            factor = 2**(pck.limit_level - lv)
            mp_inputs = [{'N_FIELDS':N_FIELDS,
                          'FIELD_INDEX':FIELD_INDEX,
                          'fname':binfiles[i]} for i in read_order]
            with multiprocessing.Pool() as pool:
                prog = tqdm(total=len(binfiles))
                for res in pool.imap_unordered(readfieldfrombinfile, mp_inputs):
                    prog.update(1)
                    for idx, arr in zip(res[0], res[1]):
                        data[factor * idx[0][0]:(idx[1][0]+1) * factor,
                             factor * idx[0][1]:(idx[1][1]+1) * factor,
                             factor * idx[0][2]:(idx[1][2]+1) * factor] = expand_array3d(arr, factor)

        if args.outfile is None:
            np.save(f"{args.variable}_ugrid_{args.plotfile.replace('plt', '')}",
                    data)
        else:
            np.save(args.outfile, data)

if __name__ == "__main__":
    main()
