import time
import pickle
import multiprocessing
import numpy as np
import matplotlib.pyplot as plt
from amr_kitchen import PlotfileCooker
from amr_kitchen.utils import shape_from_header

# TODO: rename increment_sum_masked and increment_sum
# so the functions name make sens outside of the context
# they are called
def first_increment(args):
    """
    Increments sum with data from each bfile of a level  
    """
    with open(args["file"], 'rb') as bf:
                bf.seek(args["offset"])
                h = bf.readline()
                shape = shape_from_header(h.decode('ascii'))
                # TODO: only read the integration field its going to be alot faster
                # Also do in the last level function
                # This might do the trick instead of iterating by binary file
                # box_shape = (shape[0], shape[1], shape[2])
                # skip field_pos * 8 bytes per float * points in box
                # bf.seek(np.prod(box_shape)*args["INT_FIELD"]*8, 1)
                # Only read the data from one box
                # data = np.fromfile(bf, 'float64', np.prod(box_shape))
                # data = data.reshape(box_shape, order='F')
                data = np.fromfile(bf, 'float64', np.prod(shape))
                data = data.reshape(shape, order='F')
                # TODO: use lower-case as it is not a global variable anymore
                data = data[..., args["INT_FIELD"]]
                # TODO: maybe rename the key to singular
                return np.sum(data[args["covering_masks"]] * args["dV"])

# TODO: same thing for function name
def second_increment(args):
     """
     Increments sum with data from each bfile of the finest level  
     """
     with open(args["file"], 'rb') as bf:
            bf.seek(args["offset"])
            h = bf.readline()
            shape = shape_from_header(h.decode('ascii'))
            # TODO: also only read the field data
            data = np.fromfile(bf, 'float64', np.prod(shape))
            data = data.reshape(shape, order='F')
            data = data[..., args["INT_FIELD"]]
            return np.sum(data) * args["dV"]

# Move this to amr_kitchen.utils
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


def volume_integration(plotfile, field):
    """
    Prints the volume integral of the chosen field 
    """
    # Creating a PlotfileCooker instance
    pck = PlotfileCooker(plotfile, ghost=True)
    # TODO: if the plotfile has a 'volFrac' field
    # it means there are cells truncated by an embedded boundary
    # you could create two others integration functions where dV
    # is scaled by the local volFrac if there are volFrac values
    # smaller than 1.0 in the box (sum += data[mask] * volfrac * dV)

    INT_FIELD = pck.fields[field]
    # TODO: you can remove this it was for you
    # 1. Find out the upper level mask for each box
    #    That is which points in the box are covered by
    #    Higher level data
    #   - Assume levels intersections dont jump more
    #     than one level (Dont mask lv 2 with lv 4)
    covering_masks = []
    for lv in range(pck.limit_level): # Last level is not masked
        # use box indices in list
        lv_masks = []
        next_lv_factors = pck.grid_sizes[lv + 1] // pck.box_arrays[lv + 1].shape
        for idx, indices in enumerate(pck.cells[lv]['indexes']):
            # Convert to box array indices at lv + 1
            barr_starts = np.array((indices[0] * 2) // next_lv_factors, dtype=int)
            barr_ends = np.array((indices[1] * 2) // next_lv_factors, dtype=int)
            next_level_boxes = pck.box_arrays[lv + 1][barr_starts[0]:barr_ends[0]+1,
                                                    barr_starts[1]:barr_ends[1]+1,
                                                    barr_starts[2]:barr_ends[2]+1]
            # Convert the upper level box slice to lower level bool
            bcast_factor = next_lv_factors[0] // 2
            next_lv_map = expand_array3d(next_level_boxes, bcast_factor)
            mask = np.zeros_like(next_lv_map, dtype=bool)
            # mask is true where the value is -1
            # TODO: (facultatif) here we have a bool array with the number
            # of points in the box for each box even when the box is fully
            # covered by the higher level (not added to the sum, so we read
            # data and then do nothing with it) and fully not covered so we
            # pass a full array of ones to the multiprocessing pipe (also a
            # bit wastefull) it would be nice to add flags when:
            # np.all(next_lv_map == -1) -> no need to mask
            # np.all(next_lv_map != -1) -> no need to read the data
            # And handle them in the integration part
            mask[next_lv_map == -1] = 1
            # -1 means there is no box at this level
            lv_masks.append(mask)
        covering_masks.append(lv_masks)
                
    # TODO: also can be removed as this was for you
    # 2. Do the integration with the masks
    #   No need to use masks for the finest level
    integral = 0
    mp_calls = []
    for lv in range(pck.limit_level):
        #  TODO: (Could be re done by iterating per file (maybe faster))
        #  but check if only reading the integration data is faster first
        dV = np.prod(pck.dx[lv])
        for bid, file, offset in zip(range(len(pck.boxes[lv])),
                                            pck.cells[lv]['files'],
                                            pck.cells[lv]['offsets']):
            mp_call = {"bid":bid,
                       "file":file,
                       "offset":offset,
                       "INT_FIELD":INT_FIELD,
                       "covering_masks":covering_masks[lv][bid],
                       # TODO: I think you can remove 'lv' and 'bid'
                       "lv":lv,
                       "dV":dV,}
            mp_calls.append(mp_call)
        pool = multiprocessing.Pool()
        sums = pool.map(first_increment,
                                    mp_calls)
        # TODO: Replace 'sum' with something else as it is a python keyword
        # (Function from the base library)
        for sum in sums:
            integral+=sum
    mp_calls = []
    for bid, file, offset in zip(range(len(pck.boxes[pck.limit_level])),
                                 pck.cells[pck.limit_level]['files'],
                                 pck.cells[pck.limit_level]['offsets']):

            mp_call = {"file":file,
                       "offset":offset,
                       "INT_FIELD":INT_FIELD,
                       "dV":dV,}
            mp_calls.append(mp_call)
    pool = multiprocessing.Pool()
    sums = pool.map(second_increment,
                                    mp_calls)
    for sum in sums:
        integral+=sum


    # TODO: return the value and print it in the cli
    # Also make the plotfile argument a PlotfileCooker instance
    # and create the object in the cli so we can integrate multiple
    # fields in a script without re-reading the header data each time
    print(f"Volume integral of {field} in plotfile: {integral}")
