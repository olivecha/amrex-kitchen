import numpy as np

def shape_from_header(h):
    """
    Infer the shape the box and the number of fields
    from the header in a plotfile binary file
    (Only works for 3D plotfiles)
    header: bytestring of the header line
    """
    start, stop, _, nfields = h.split()[-4:]
    nfields = int(nfields)
    start = np.array(start.split('(')[-1].replace(')','').split(','), dtype=int)
    stop = np.array(stop.replace('(', '').replace(')', '').split(','), dtype=int)
    shape = stop - start + 1
    total_shape = [shape[0], shape[1], shape[2], nfields]
    return total_shape
