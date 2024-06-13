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
    header: bytestring of the header line
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
