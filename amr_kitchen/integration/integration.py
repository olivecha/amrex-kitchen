import time
import pickle
import multiprocessing
import numpy as np
import matplotlib.pyplot as plt
from amr_kitchen import PlotfileCooker
from amr_kitchen.utils import shape_from_header

def first_increment(args):
    """
    Increments sum with data from each bfile of a level  
    """
    with open(args["file"], 'rb') as bf:
                bf.seek(args["offset"])
                h = bf.readline()
                shape = shape_from_header(h.decode('ascii'))
                data = np.fromfile(bf, 'float64', np.prod(shape))
                data = data.reshape(shape, order='F')
                data = data[..., args["INT_FIELD"]]
                return np.sum(data[args["covering_masks"]] * args["dV"])

def second_increment(args):
     """
     Increments sum with data from each bfile of the finest level  
     """
     with open(args["file"], 'rb') as bf:
            bf.seek(args["offset"])
            h = bf.readline()
            shape = shape_from_header(h.decode('ascii'))
            data = np.fromfile(bf, 'float64', np.prod(shape))
            data = data.reshape(shape, order='F')
            data = data[..., args["INT_FIELD"]]
            return np.sum(data) * args["dV"]

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


def volume_integration(plotfile,field):
    """
    Prints the volume integral of the chosen field 
    """
    # Creating a PlotfileCooker instance
    pck = PlotfileCooker(plotfile,ghost=True)

    INT_FIELD = pck.fields[field]
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
            mask[next_lv_map == -1] = 1
            # -1 means there is no box at this level
            lv_masks.append(mask)
        covering_masks.append(lv_masks)
                
    # 2. Do the integration with the masks
    # - Could be faster by not indexing with full masks (all 1 or 0)
    #   Esp if the masks are piped to multiprocessing
    #   No need to use masks for the finest level
    integral = 0
    mp_calls = []
    for lv in range(pck.limit_level):
        # Could be re done by iterating per file (maybe faster)
        dV = np.prod(pck.dx[lv])
        for bid, file, offset in zip(range(len(pck.boxes[lv])),
                                            pck.cells[lv]['files'],
                                            pck.cells[lv]['offsets']):
            mp_call = {"bid":bid,
                       "file":file,
                       "offset":offset,
                       "INT_FIELD":INT_FIELD,
                       "covering_masks":covering_masks[lv][bid],
                       "lv":lv,
                       "dV":dV,}
            mp_calls.append(mp_call)
        pool = multiprocessing.Pool()
        sums = pool.map(first_increment,
                                    mp_calls)
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


    print(f"Volume integral of {field} in plotfile: {integral}")

