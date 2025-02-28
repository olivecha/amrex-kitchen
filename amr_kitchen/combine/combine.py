import os
import time
import multiprocessing
import numpy as np
from amr_kitchen import PlotfileCooker
from amr_kitchen.utils import (shape_from_header,
                               indices_from_header,
                               header_from_indices,)

def parallel_combine_by_binfile(args):
    """
    Read two binary file and write the contents
    to a third one
    """
    # New offsets
    offsets = []
    # Open the three files
    with open(args['bfile_r1'], 'rb') as bf1:
        with open(args['bfile_r2'], 'rb') as bf2:
            with open(args['bfile_w'], 'wb') as bfw:
                while True:
                    try:
                        # Get the read shape and indices
                        h1 = bf1.readline()
                        h1 = h1.decode('ascii')
                        shape1 = shape_from_header(h1)
                        idx1 = indices_from_header(h1)
                        # For both plotfiles
                        h2 = bf2.readline()
                        h2 = h2.decode('ascii')
                        shape2 = shape_from_header(h2)
                        idx2 = indices_from_header(h2)
                    except:
                        break
                    # Define the write binary header
                    hw = header_from_indices(idx1[0],
                           idx1[1],
                           len(args['vidxs1']) + len(args['vidxs2']))
                    # save the current offset
                    offsets.append(bfw.tell())
                    # Write the header and data
                    bfw.write(hw)
                    data1 = np.fromfile(bf1, 'float64', np.prod(shape1))
                    data1 = data1.reshape(shape1, order='F')[..., args['vidxs1']]
                    data2 = np.fromfile(bf2, 'float64', np.prod(shape2))
                    data2 = data2.reshape(shape2, order='F')[..., args['vidxs2']]
                    dataw = np.concatenate([data1.flatten(order='F'),
                                            data2.flatten(order='F')])
                    bfw.write(dataw.tobytes())
    return offsets

def parallel_combine_by_binfile_offsets(args):
    """
    Read two binary file and write the contents
    to a third one
    """
    # New offsets
    offsets = []
    # Open the three files
    with open(args['bfile_r1'], 'rb') as bf1:
        with open(args['bfile_r2'], 'rb') as bf2:
            with open(args['bfile_w'], 'wb') as bfw:
                for off1, off2 in zip(args["offst_r1"],
                                      args["offst_r2"]):
                    # Go to the box in the file
                    bf1.seek(off1)
                    bf2.seek(off2)
                    # Get the read shape and indices
                    h1 = bf1.readline()
                    h1 = h1.decode('ascii')
                    shape1 = shape_from_header(h1)
                    idx1 = indices_from_header(h1)
                    # For both plotfiles
                    h2 = bf2.readline()
                    h2 = h2.decode('ascii')
                    shape2 = shape_from_header(h2)
                    idx2 = indices_from_header(h2)
                    # Define the write binary header
                    hw = header_from_indices(idx1[0],
                           idx1[1],
                           len(args['vidxs1']) + len(args['vidxs2']))
                    # save the current offset
                    offsets.append(bfw.tell())
                    # Write the header and data
                    bfw.write(hw)
                    data1 = np.fromfile(bf1, 'float64', np.prod(shape1))
                    data1 = data1.reshape(shape1, order='F')[..., args['vidxs1']]
                    data2 = np.fromfile(bf2, 'float64', np.prod(shape2))
                    data2 = data2.reshape(shape2, order='F')[..., args['vidxs2']]
                    dataw = np.concatenate([data1.flatten(order='F'),
                                            data2.flatten(order='F')])
                    bfw.write(dataw.tobytes())
    return offsets

def parallel_combine_by_boxes_offsets(args):
    """
    Read two binary file and write the contents
    to a third one
    """
    # New offsets
    offsets = []
    # Open two files
    with open(args['bfile_r1'], 'rb') as bf1:
        with open(args['bfile_w'], 'wb') as bfw:
            for bf_path2, offset2 in zip(args['bfiles_r2'],
                                        args['offsets2']):
                # Go to the box in the files
                h1 = bf1.readline()
                h1 = h1.decode('ascii')
                shape1 = shape_from_header(h1)
                data1 = np.fromfile(bf1, 'float64', np.prod(shape1))
                data1 = data1.reshape(shape1, order='F')[..., args['vidxs1']]
                idx1 = indices_from_header(h1)
                # Define the write binary header
                hw = header_from_indices(idx1[0],
                       idx1[1],
                       len(args['vidxs1']) + len(args['vidxs2']))
                # save the current offset
                offsets.append(bfw.tell())
                # Write the header
                bfw.write(hw)
                # Get the data in the second file
                with open(bf_path2, 'rb') as bf2:
                    bf2.seek(offset2)
                    h2 = bf2.readline()
                    h2 = h2.decode('ascii')
                    shape2 = shape_from_header(h2)
                    idx2 = indices_from_header(h2)
                    data2 = np.fromfile(bf2, 'float64', np.prod(shape2))
                data2 = data2.reshape(shape2, order='F')[..., args['vidxs2']]
                # Write the header and data
                dataw = np.concatenate([data1.flatten(order='F'),
                                        data2.flatten(order='F')])
                bfw.write(dataw.tobytes())
    return offsets


def rewrite_level_header(pck1: PlotfileCooker,
                         pck2: PlotfileCooker,
                         pltout: str,
                         lv: int,
                         nfields: int,
                         mapped_offsets: np.ndarray[int],
                         field_indices1: list[int],
                         field_indices2: list[int]) -> None:
    # New level header path
    cell_header_w = os.path.join(os.getcwd(),
                                 pltout, 
                                 pck1.cell_paths[lv],
                                 "Cell_H")
    # Header we read
    cell_header_r1 = os.path.join(os.getcwd(),
                                  pck1.pfile,
                                  pck1.cell_paths[lv], 
                                  'Cell_H')
    # Header we read
    cell_header_r2 = os.path.join(os.getcwd(),
                                  pck2.pfile,
                                  pck2.cell_paths[lv], 
                                  'Cell_H')
    with open(cell_header_w, 'w') as ch_w:
        with open(cell_header_r1, 'r') as ch_r1:
            with open(cell_header_r2, 'r') as ch_r2:
                # First two lines
                for i in range(2):
                    l = ch_r1.readline()
                    ch_r2.readline()
                    ch_w.write(l)
                # Number of fields
                ch_r1.readline()
                ch_r2.readline()
                ch_w.write(f"{nfields}\n")
                # Mesh stays the same
                while True:
                    l = ch_r1.readline()
                    ch_r2.readline()
                    if "FabOnDisk:" in l:
                        new_l = l.split()[:-1]
                        new_l.append(str(mapped_offsets[0]))
                        ch_w.write(' '.join(new_l) + "\n")
                        break
                    else:
                        ch_w.write(l)
                # Write the cell indexes
                for fst in mapped_offsets[1:]:
                    l = ch_r1.readline()
                    ch_r2.readline()
                    new_l = l.split()[:-1]
                    new_l.append(str(fst))
                    ch_w.write(' '.join(new_l) + "\n")
                # Blank line
                ch_w.write(ch_r1.readline())
                line = ch_r1.readline()
                ch_r2.readline()
                ch_r2.readline()
                ncells, _ = line.split(',')
                ch_w.write(f"{ncells},{nfields}\n")
                # Min vals
                for li in range(int(ncells)):
                    line1 = ch_r1.readline().split(',')[:-1]
                    line2 = ch_r2.readline().split(',')[:-1]
                    min_vals1 = np.array(line1)[field_indices1]
                    min_vals2 = np.array(line2)[field_indices2]
                    min_vals = np.concatenate([min_vals1, min_vals2])
                    ch_w.write(','.join(min_vals) + ',\n')
                # Blank line
                ch_w.write(ch_r1.readline())
                ncells, _ = ch_r1.readline().split(',')
                ch_r2.readline()
                ch_r2.readline()
                ch_w.write(f"{ncells},{nfields}\n")
                # Max vals
                for li in range(int(ncells)):
                    line1 = ch_r1.readline().split(',')[:-1]
                    line2 = ch_r2.readline().split(',')[:-1]
                    max_vals1 = np.array(line1)[field_indices1]
                    max_vals2 = np.array(line2)[field_indices2]
                    max_vals = np.concatenate([max_vals1, max_vals2])
                    ch_w.write(','.join(max_vals) + ',\n')

def validate_combine_input(*args, **kwargs) -> dict:
    """
    Function to validate arbitrary input to the
    combine function. This return the default values
    if no problems are found within the input
    """
    # Store output in a dict
    output = {}
    # Validate plotfile dimensions
    if (args[0].ndims < 3 or
        args[1].ndims < 3):
        # TODO: make work for 2 dimensions
        raise NotImplementedError(
              ("The combine tool only supports 3D plotfiles"
               " for the moment"))
    # Validate plotfile compatibility
    # Do the plotfiles have the same AMR box structure
    if args[0] != args[1]:
        raise ValueError(("The two plotfiles do not have the"
                          " same AMR structure and cannot be"
                          " combined"))
    # Do the plotfile have the same number of levels
    # TODO: this probably fails if plotfile1 != plotfile2 and 
    # could perhaps be removed
    if args[0].limit_level != args[1].limit_level:
        raise ValueError(("The two plotfiles do not have the"
                          " same number of AMR levels and"
                          " cannot be combined"))
    # Fastest and least constrained combination
    # (Keep box order in every binary file)
    # This only works if the box offsets are in the same order
    output["mode"] = "byfile"
    bf_vec = np.vectorize(lambda s:os.path.split(s)[-1])
    for lv in range(args[0].limit_level + 1):
        bfiles_1 = bf_vec(args[0].cells[lv]['files'])
        bfiles_2 = bf_vec(args[1].cells[lv]['files'])
        # Are the boxes of each plotfile in the same binary file
        # (This is not always the case even if plotfile1 == plotfile2)
        if not np.array_equal(bfiles_1,
                              bfiles_2):
            # We need to combine box by box
            output["mode"] = "bybox"
            break # Cannot get worse
        else:
            for bf in bfiles_1:
                offsets_1 = np.array(args[0].cells[lv]['offsets'])[bf == bfiles_1]
                offsets_2 = np.array(args[1].cells[lv]['offsets'])[bf == bfiles_2]
                if not np.array_equal(np.argsort(offsets_1),
                                      np.argsort(offsets_2)):
                    # We need to keep track of the offsets for a given
                    # binary file
                    output["mode"] == "byoffset"
    # Validate the kepts fields do not have any duplicates
    if kwargs["vars1"] is None:
        vars1 = list(args[0].fields.keys())
    else:
        vars1 = []
        for var in kwargs["vars1"]:
            if var in args[0].fields:
                vars1.append(var)
            else:
                print((f"{var} is not a field in {args[0].pfile}"
                        " (skipping)"))
        if len(vars1) == 0:
            raise ValueError((f"No fields to combine from {args[0].pfile}"
                               " from the parsed input"))
    if kwargs["vars2"] is None:
        vars2 = list(args[1].fields.keys())
    else:
        vars2 = []
        for var in kwargs["vars2"]:
            if var in args[1].fields:
                vars2.append(var)
            else:
                print((f"{var} is not a field in {args[1].pfile}"
                        " (skipping)"))
    # Remove duplicates
    vars2 = [var for var in vars2 if var not in vars1]
    if len(vars2) == 0:
        raise ValueError((f"No fields to combine from {args[1].pfile}"
                           " from the parsed input, or that are not in"
                           " in the combined variables from {args[0].pfile}"))
    output["vars1"] = vars1
    output["vars2"] = vars2
    output["cbvars"] = vars1 + vars2
    # Validate the output dir and generate a default value
    if kwargs["pltout"] is None:
        if kwargs["inplace"]:
            # TODO: support inplace combination but warn
            # that it would be dangerous unless a temporary
            # file is used
            raise NotImplementedError(
                  ("The inplace combination of plotfile is"
                   " not supported yet"))
        # Not inplace so default output dir
        else:
            # Concatenate the plotfile directories
            pdir1 = os.path.split(args[0].pfile)[-1]
            pdir2 = os.path.split(args[1].pfile)[-1]
            output["pltout"] = pdir1+pdir2
    else:
        output["pltout"] = kwargs["pltout"]
    return output

def combine(pck1, pck2, pltout=None, 
            vars1=None, vars2=None, inplace=False):
    """
    Function to combine fields between plotfiles
    pck1 : PlotfileCooker instance of the first plotfile
    pck2 : PlotfileCooker instance of the second plotfile
    pckout : Output plotfile dir
    vars1 : Variables to keep in plotfile 1
    vars2 : Variables to keep in plotfile 2
    """
    clean_args = validate_combine_input(pck1, pck2, pltout=pltout,
                                        vars1=vars1, vars2=vars2, inplace=inplace)
    # Read the parsed input args
    pltout = clean_args["pltout"]
    vars1, vars2 = clean_args["vars1"], clean_args["vars2"]
    vidxs1 = [pck1.fields[v] for v in vars1]
    vidxs2 = [pck2.fields[v] for v in vars2]
    cbvars = clean_args["cbvars"]
    cbmode = clean_args["mode"]
    # Number of fields in the output
    nfields = len(cbvars)
    # make the output dir
    pck1.make_dir_tree(pltout)
    # Write the new global header
    pck1.write_global_header_new_fields(pltout, cbvars)
    pool = multiprocessing.Pool()
    # Combine the plotfile using the required mode
    for lv in range(pck1.limit_level + 1):
        lvstart = time.time()
        # Boxes are in the same files in the same order
        if cbmode == "byfile":
            new_offsets = pool.map(parallel_combine_by_binfile,
                                   pck1.by_matched_offsets_output(pck2, lv, pltout,
                                                                  vidxs1=vidxs1,
                                                                  vidxs2=vidxs2))
        # Boxes are in the same files by with different orders
        elif cbmode == "byoffset":
            new_offsets = pool.map(parallel_combine_by_binfile_offsets,
                                   pck1.by_matched_offsets_output(pck2, lv, pltout, 
                                                                  vidxs1=vidxs1,
                                                                  vidxs2=vidxs2))
        # Boxes are in different files (first plotfile structure is kept)
        elif cbmode == "bybox":
            new_offsets = pool.map(parallel_combine_by_boxes_offsets,
                                   pck.by_matched_boxes_output(pck2, lv, pltout,
                                                               vidxs1=vidxs1,
                                                               vidxs2=vidxs2))
        # Reorder the offsets to match the box order
        mapped_offsets = np.empty(len(pck1.boxes[lv]), dtype=int)
        for file_idxs, offsets in zip(pck1.map_bfile_offsets(lv), new_offsets):
            mapped_offsets[file_idxs] = offsets
        # Write the new level header
        rewrite_level_header(pck1, pck2, pltout, lv, nfields, 
                             mapped_offsets, vidxs1, vidxs2)

        print(f"Combined Level {lv} ({time.time() - lvstart:.2f} s)")

