import numpy as np
import multiprocessing
from concurrent.futures import ProcessPoolExecutor as Pool
from functools import partial
from p_tqdm import p_map
import time

def parallel_cook(args):
    """
    Cook a single cell binary file
    Multiprocessing function
    """
    # New offsets
    offsets = np.zeros(args["ncells"], dtype=int)
    # Open the read and write
    print(f"Cooking file {args['bfile']}...")
    now = time.time()
    with open(args['bfile'], 'rb') as bfr, open(args['bfile_w'], 'wb') as bfw:
        for indexes, fst_r, idx in zip(args['indexes'], args['offsets_r'], args['cell_indexes']):
            # Store pointer position to update header
            offsets[idx] = bfw.tell()
            # Go to the data
            bfr.seek(fst_r)
            # Get the header
            header = bfr.readline().decode('ascii')
            # Replace with single var just to be sure
            header_w = header.replace(f'{args["nvars"]}\n', '1\n')
            # Write to binary file
            bfw.write(header_w.encode('ascii'))
             # Read the data
            shape = indexes[1] - indexes[0] + 1
            total_shape = (shape[0], shape[1], shape[2], args['nvars'])
            arr = np.fromfile(bfr, "float64", np.prod(total_shape))
            arr = arr.reshape(total_shape, order="F")
            # Reshape into dicts
            data = {}
            for field in args["fields"]:
                fidx = args["fields"][field]
                data[field] = arr[:, :, :, fidx]
            # Cook the recipe
            arr_out = args['recipe'](data)
            arr_bytes = arr_out.flatten(order="F").tobytes()
            bfw.write(arr_bytes)
    print(f"Done!", np.around(time.time() - now, 2), "s")
    return offsets

def parallel_compute_array(arr, fields=None, recipe=None):
    # Go to the data
    # Reshape into dicts
    data = {}
    for field in fields:
        fidx = fields[field]
        data[field] = arr[:, :, :, fidx]
    # Cook the recipe
    arr_out = recipe(data)
    arr_bytes = arr_out.flatten(order="F").tobytes()
    return arr_bytes


def parallel_cook_byarray(args):
    # New offsets
    offsets = np.zeros(args["ncells"], dtype=int)
    # Read the data
    print(f"Reading data from {args['bfile']} ...")
    read_start = time.time()
    headers_w = []
    arrays = []
    with open(args['bfile'], 'rb') as bfr:
        for indexes, fst_r in zip(args['indexes'], args['offsets_r']):
            # got to array start
            bfr.seek(fst_r)
            # Get the header
            header = bfr.readline().decode('ascii')
            # Replace with single var just to be sure
            header_w = header.replace(f'{args["nvars"]}\n', '1\n')
            # Store for later
            headers_w.append(header_w.encode('ascii'))
            # Read the data
            shape = indexes[1] - indexes[0] + 1
            total_shape = (shape[0], shape[1], shape[2], args['nvars'])
            arr = np.fromfile(bfr, "float64", np.prod(total_shape))
            arr = arr.reshape(total_shape, order="F")
            arrays.append(arr)
    print(f"Read time for file {args['bfile']}", 
          np.around(time.time() - read_start, 2), 's')

    print(f"Cooking file {args['bfile']}...")
    now = time.time()
    with Pool() as pool:
        output = pool.map(partial(parallel_compute_array, 
                           recipe=args['recipe'],
                           fields=args['fields']),
                           arrays)
    print(f"Done!", np.around(time.time() - now, 2), "s")

    with open(args['bfile_w'], 'wb') as bfw:
        for idx, header, arr in zip(args['cell_indexes'], headers_w, output):
            # Store pointer position to update header
            offsets[idx] = bfw.tell()
            # Write to binary file
            bfw.write(header)
            bfw.write(arr)
    return offsets
