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
    total_shape = [shape[0], shape[1], shape[2], nfields]
    return total_shape

def indices_from_header(h):
    """
    Infer the global indices of the box
    from the header in a plotfile binary file
    h: string of the header line
    """
    start, stop, _, nfields = h.split()[-4:]
    nfields = int(nfields)
    start = np.array(start.split('(')[-1].replace(')','').split(','), dtype=int)
    stop = np.array(stop.replace('(', '').replace(')', '').split(','), dtype=int)
    return [start, stop]

def header_from_indices(start, stop, nfields):
    """
    Creates a binary file header from the box
    global indices and number of fields
    start: start indices
    stop: stop indices
    nfields: number of fields in the plotfile
    """
    header_const = "FAB ((8, (64 11 52 0 1 12 0 1023)),(8, (8 7 6 5 4 3 2 1)))"
    header_indices = (f"(({start[0]},{start[1]},{start[2]})"
                      f" ({stop[0]},{stop[1]},{stop[2]})"
                      f" (0,0,0)) {nfields}\n")
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

