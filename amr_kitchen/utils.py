import numpy as np

class TastesBadError(Exception):
    """
    Custom exception to tell that something
    is wrong with a plotfile
    """
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

# TODO: replace with more specific errors
# depending on where the problem is
class BadTastingHeadersError(Exception):
    pass

class BadTastingBinariesError(Exception):
    pass                    

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

# AMReX FAB format descriptors for the supported floating point
# data types. These are the exact descriptor strings AMReX writes for
# native plotfiles. The byte size of Real is encoded by the bit width
# at the start of the inner format tuple (64 -> float64, 32 -> float32),
# *not* by the leading number, which can be 8 even for float32 data.
REAL_DESCRIPTORS = {
    8: "FAB ((8, (64 11 52 0 1 12 0 1023)),(8, (8 7 6 5 4 3 2 1)))",
    4: "FAB ((8, (32 8 23 0 1 9 0 127)),(4, (4 3 2 1)))",
}

# Map the byte size of Real to the matching numpy dtype
DTYPE_FROM_SIZE = {8: 'float64', 4: 'float32'}


def dtype_from_header(h):
    """
    Infer the floating point data type and its byte size from the
    FAB header in a plotfile binary file. AMReX encodes the format of
    Real as a tuple whose first number is the total bit width, e.g.
    "FAB ((8, (64 11 52 ..." for float64 and
    "FAB ((8, (32 8 23 ..." for float32. Note the leading number can be
    8 even for float32 data, so the bit width of the inner tuple is the
    reliable field to parse.
    h: string of the header line
    returns: (numpy dtype string, itemsize in bytes)
    """
    # The inner format tuple is the third '('-delimited group, e.g.
    # "FAB ((8, (64 11 52 ..." -> "64 11 52 ..." whose first token is
    # the total number of bits used to store a Real.
    nbits = int(h.split('(')[3].split()[0])
    real_size = nbits // 8
    if real_size not in DTYPE_FROM_SIZE:
        raise ValueError((f"Unsupported Real bit width in FAB header: "
                          f"{nbits} (expected 64 for float64 or 32 for"
                          f" float32)"))
    return DTYPE_FROM_SIZE[real_size], real_size


def shape_from_header(h):
    """
    Infer the shape the box and the number of fields
    from the header in a plotfile binary file
    (Only works for 3D plotfiles)
    h: string of the header line
    """
    start, stop, _, nfields = h.split()[-4:]
    nfields = int(nfields)
    start = np.array(start.split('(')[-1].replace(')','').split(','), dtype=int)
    stop = np.array(stop.replace('(', '').replace(')', '').split(','), dtype=int)
    shape = stop - start + 1
    #total_shape = [shape[0], shape[1], shape[2], nfields]
    total_shape = np.append(shape, nfields)
    return total_shape

def indices_from_header(h):
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
    return [start, stop]

def header_from_indices(start, stop, nfields, real_size=8):
    """
    Creates a binary file header from the box
    global indices and number of fields
    start: start indices
    stop: stop indices
    nfields: number of fields in the plotfile
    real_size: byte size of Real (8 for float64, 4 for float32)
               selects the matching FAB format descriptor
    """
    header_const = REAL_DESCRIPTORS[real_size]
    header_indices = (f"((" + ','.join([str(s) for s in start]) + ')'
                      f" (" + ','.join([str(s) for s in stop]) + ")"
                      f" (" + ','.join(["0" for _ in stop]) + f")) {nfields}\n")
    header = header_const + header_indices
    return header.encode('ascii')

def indexes_and_shape_from_header(header):
    """
    This takes the byte string of a box header in a plotfile binary
    file and infers the indexes of the box and the number of fields
    in the plotfile
    """
    h = header.decode('ascii')
    start, stop, _, nfields = h.split()[-4:]
    nfields = int(nfields)
    start = np.array(start.split('(')[-1].replace(')','').split(','), dtype=int)
    stop = np.array(stop.replace('(', '').replace(')', '').split(','), dtype=int)
    shape = stop - start + 1
    shape = [s for s in shape]
    shape.append(nfields)
    return [start, stop], tuple(shape)


def shapes_from_header_vardims(header, ndim):
    """
    Function to infer the data shape from the 
    binary header which also words for 2D
    plotfiles
    """
    h = header.decode("ascii")
    start, stop, _, nfields = h.split()[-4:]
    nfields = int(nfields)
    start = np.array(start.split('(')[-1].replace(')','').split(','), dtype=int)
    stop = np.array(stop.replace('(', '').replace(')','').split(','), dtype=int)
    shape = stop - start + 1
    total_shape = []
    for i in range(ndim):
        total_shape.append(shape[i])
    total_shape.append(nfields)
    return total_shape

def global2local(indices, refindices, n_ghost=1):
    """
    Convert global box indexes to indexes in a
    ghost cell padded reference box
    """
    i_start = max(refindices[0][0] - n_ghost, 0)
    j_start = max(refindices[0][1] - n_ghost, 0)
    k_start = max(refindices[0][2] - n_ghost, 0)
    return [[indices[0][0] - i_start,
             indices[0][1] - j_start,
             indices[0][2] - k_start],
            [indices[1][0] - i_start,
             indices[1][1] - j_start,
             indices[1][2] - k_start]]

