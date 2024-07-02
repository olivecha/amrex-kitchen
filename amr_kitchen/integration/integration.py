import time
import pickle
import multiprocessing
import numpy as np
import matplotlib.pyplot as plt
from amr_kitchen import PlotfileCooker
from amr_kitchen.utils import shape_from_header
from amr_kitchen.utils import expand_array3d

def increment_sum_masked(args):
    """
    Increments sum with data from each bfile of a level  
    """
    with open(args["file"], 'rb') as bf:
        bf.seek(args["offset"])
        h = bf.readline()
        shape = shape_from_header(h.decode('ascii'))
        box_shape = (shape[0], shape[1], shape[2])
        # skip field_pos 
        bf.seek(np.prod(box_shape)*args["int_field"]*8, 1)
        # Only read the data from one box
        data = np.fromfile(bf, 'float64', np.prod(box_shape))
        data = data.reshape(box_shape, order='F')
        return np.sum(data[args["covering_mask"]] * args["dV"])

def increment_sum(args):
     """
     Increments sum with data from each bfile of the finest level  
     """
     with open(args["file"], 'rb') as bf:
        bf.seek(args["offset"])
        h = bf.readline()
        shape = shape_from_header(h.decode('ascii'))
        box_shape = (shape[0], shape[1], shape[2])
        # skip field_pos 
        bf.seek(np.prod(box_shape)*args["int_field"]*8, 1)
        # Only read the data from one box
        data = np.fromfile(bf, 'float64', np.prod(box_shape))
        data = data.reshape(box_shape, order='F')
        return np.sum(data) * args["dV"]

def volume_integration(pck, field):
    """
    Prints the volume integral of the chosen field 
    """
    # TODO: if the plotfile has a 'volFrac' field
    # it means there are cells truncated by an embedded boundary
    # you could create two others integration functions where dV
    # is scaled by the local volFrac if there are volFrac values
    # smaller than 1.0 in the box (sum += data[mask] * volfrac * dV)

    INT_FIELD = pck.fields[field]

    covering_masks = []
    for lv in range(pck.limit_level): # Last level is not masked
        # use box indices in list
        lv_masks = []
        next_lv_factors = pck.grid_sizes[lv + 1] // pck.box_arrays[lv + 1].shape
        for idx, indices in enumerate(pck.cells[lv]['indexes']):
            # Convert to box array indices at lv + 1
            barr_starts = np.array((indices[0] * 2) // next_lv_factors, dtype=int)
            barr_ends = np.array((indices[1] * 2) // next_lv_factors, dtype=int)
            #if 2d plotfile
            if pck.ndims == 2:
                next_level_boxes = pck.box_arrays[lv + 1][barr_starts[0]:barr_ends[0]+1,
                                                        barr_starts[1]:barr_ends[1]+1]
            #if 3d plotfile
            else:
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
                       "int_field":INT_FIELD,
                       "covering_mask":covering_masks[lv][bid],
                       # TODO: I think you can remove 'lv' and 'bid'
                       "lv":lv,
                       "dV":dV,}
            mp_calls.append(mp_call)
        pool = multiprocessing.Pool()
        sums = pool.map(increment_sum_masked,
                                    mp_calls)
        for summation in sums:
            integral+=summation
    mp_calls = []
    for bid, file, offset in zip(range(len(pck.boxes[pck.limit_level])),
                                 pck.cells[pck.limit_level]['files'],
                                 pck.cells[pck.limit_level]['offsets']):

            mp_call = {"file":file,
                       "offset":offset,
                       "int_field":INT_FIELD,
                       "dV":dV,}
            mp_calls.append(mp_call)
    pool = multiprocessing.Pool()
    sums = pool.map(increment_sum,
                                    mp_calls)
    for summation in sums:
        integral+=summation

    return integral
    
