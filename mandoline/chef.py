import numpy as np

def parallel_cook(args):
    """
    Cook a single cell binary file
    Multiprocessing function
    """
    # New offsets
    offsets = np.zeros(args["ncells"], dtype=int)
    # Open the read and write
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
    return offsets
