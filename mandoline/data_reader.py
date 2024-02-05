import numpy as np
from mandoline.utils import expand_array


def slice_box(args):
    """
    Multiprocessing function reading and slicing a AMR cell
    ----
    Input:
    single input argument for multiprocessing

    (index, level, slicedata) = args 
    index: Index of the cell to read (int)
    level: Current AMR Level of the cell (int)
    slicedata: SliceData instance containing the slice information

    Output:

    The output is divided between right and left sides of the slice plane:
    output = [left, right]
    For each side, the default value is None if no data is available.
    This is the case if the slice is at the left of the last points in
    a box.
    If slice data exists the output for a given side is a dict:
    left = {'sx'    : [x_start, x_stop],           # Slice indexes in x-direction
            'sy'    : [y_start, y_stop],           # Slice indexes in y-direction
            'data'  : [arr[x_shape, y_shape], ...] # Arrays containing the data
            'normal': normal_coord,                # Normal coordinate of the data
            'level' : amr_level}                   # AMR level of the data
    """
    # Unpack input
    idx = args['idx']
    Lv = args['Lv']
    pos = args['pos']
    fidxs = args['fidxs']
    limit_level = args['limit_level']
    # Get the cell data from the HeaderData class
    indexes = args['indexes']
    cfile = args['cfile']
    offset = args['offset']
    box = args['box']
    cx = args['cx']
    cy = args['cy']
    cn = args['cn']
    dx = args['dx']
    # Factor between curent grid and covering grid
    factor = 2**(limit_level - Lv)
    # Compute the slice indexes for the global grid
    x_start = indexes[0][cx] * factor
    x_stop = (indexes[1][cx] + 1) * factor
    y_start = indexes[0][cy] * factor
    y_stop = (indexes[1][cy] + 1) * factor
    shape = (indexes[1][0] - indexes[0][0] + 1,
             indexes[1][1] - indexes[0][1] + 1,
             indexes[1][2] - indexes[0][2] + 1)

    # Compute the grid in the normal direction
    normal_grid = np.linspace(box[cn][0] + dx[Lv][cn]/2,
                              box[cn][1] - dx[Lv][cn]/2,
                              shape[cn])
    
    # Size on disk of the data slice
    byte_size = np.prod(shape)
    # A list of arrays for output
    data_arrays = []
    # Open the binary file even if we may be only looking
    # For the grid level
    with open(cfile, "rb") as f:
        # For each field index
        for fidx in fidxs:
            # Try to catch fidx = None
            try:
                # Go to the cell start
                f.seek(offset)
                # Skip the header
                f.readline()
                # Go to the field we want
                f.seek(byte_size*8*fidx, 1)
                # Could be optimized by reading contiguous fields
                # At once especially if all the data is requested
                # Read the data
                arr = np.fromfile(f, "float64", byte_size)
                # Fortran order perhaps a legacy of the early AMReX
                # versions
                arr = arr.reshape(shape, order="F")
                data_arrays.append(arr)
            # If fidx is None (for grid_level) we catch it 
            except TypeError:
                # level is always added to the output
                pass

    # Slice indexes for the data array
    arr_indices = [0, 0, 0]
    arr_indices[cx] = slice(0, shape[cx])
    arr_indices[cy] = slice(0, shape[cy])
    # Create output depending on slice position
    output = [None, None]
    # Case when the plane is between the last point and the left edge
    if pos > normal_grid[shape[cn] - 1]:
        # Slice is at the end of the box in the normal direction
        arr_indices[cn] = shape[cn] - 1
        # THIS IS WHERE THE MANDOLINE HAPPENS YEHAW
        arr_data = [arr[tuple(arr_indices)] for arr in data_arrays]
        # Only write to left interpolation plane
        output[0] = {'sx':[x_start, x_stop], # Slice in limit_level grid
                     'sy':[y_start, y_stop], 
                     'data':[expand_array(arr, factor) for arr in arr_data], # Expanded data
                     'normal':normal_grid[shape[cn] - 1], # normal position for interpolation
                     'level':Lv}
    
    # Case when the plane is between the first point and the right edge
    elif pos < normal_grid[0]:
        # Slice is at the beginning of the box
        arr_indices[cn] = 0
        # Slice the data (MANDOLINE REAL GOURMET)
        arr_data = [arr[tuple(arr_indices)] for arr in data_arrays]
        # Only write to right interpolation plane
        output[1] = {'sx':[x_start, x_stop],
                     'sy':[y_start, y_stop],  # Slice in limit_level grid
                     'data':[expand_array(arr, factor) for arr in arr_data],  # Expanded data
                     'normal':normal_grid[0],  # normal position for interpolation
                     'level':Lv}
        
    # Case when the plane lands on a grid point
    elif np.isclose(pos, normal_grid).any():
        # Find out which point index
        match_idx = np.where(np.isclose(pos, normal_grid))[0][0]
        # Slice directly at the intersect point
        arr_indices[cn] = match_idx
        arr_data = [arr[tuple(arr_indices)] for arr in data_arrays]
        # Put left and right interpolation on the same point
        output[0] = {'sx':[x_start, x_stop], 
                     'sy':[y_start, y_stop],  # Slice in limit_level grid
                     'data':[expand_array(arr, factor) for arr in arr_data],  # Expanded data
                     'normal':normal_grid[match_idx],  # normal position for interpolation
                     'level':Lv}
        # We could only keep it on one side since its the same data and remember which side
        # It is when we do the interpolation, and dont interpolate when the normal coordinate
        # is equal to the slice coordinate. (#TODO)
        output[1] = {'sx':[x_start, x_stop], 
                     'sy':[y_start, y_stop],  # Slice in limit_level grid
                     'data':[expand_array(arr, factor) for arr in arr_data],  # Expanded data
                     'normal':normal_grid[match_idx],  # normal position for interpolation
                     'level':Lv}

    # Case when the plane is between grid points in the box
    else:
        # Slice on both sides
        idx_left = np.where(pos > normal_grid)[0][-1]
        idx_right = np.where(pos < normal_grid)[0][0]
        # Each point to its interpolation plane
        arr_indices[cn] = idx_left
        arr_data = [arr[tuple(arr_indices)] for arr in data_arrays]
        # Maybe we dont need to copy the slice (#TODO)
        output[0] = {'sx':[x_start, x_stop], 
                     'sy':[y_start, y_stop],  # Slice in limit_level grid
                     'data':[expand_array(arr, factor).copy() for arr in arr_data], 
                     'normal':normal_grid[idx_left],  # normal position for interpolation
                     'level':Lv}
        # Right plane
        arr_indices[cn] = idx_right
        arr_data = [arr[tuple(arr_indices)] for arr in data_arrays]
        # Maybe we dont need to copy the slice (#TODO)
        output[1] = {'sx':[x_start, x_stop], 
                     'sy':[y_start, y_stop],  # Slice in limit_level grid
                     'data':[expand_array(arr, factor).copy() for arr in arr_data], 
                     'normal':normal_grid[idx_right],  # normal position for interpolation
                     'level':Lv}

    return output 


def plate():
    pass
