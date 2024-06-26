import os
import time
import multiprocessing
import numpy as np
from amr_kitchen import PlotfileCooker
from amr_kitchen.utils import shape_from_header
from amr_kitchen.utils import indices_from_header
from amr_kitchen.utils import header_from_indices

def parallel_combine_binary_files(args):
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
                           shape1[3] + shape2[3])
                    # save the current offset
                    offsets.append(bfw.tell())
                    # Write the header and data
                    bfw.write(hw)
                    data1 = np.fromfile(bf1, 'float64', np.prod(shape1))
                    data2 = np.fromfile(bf2, 'float64', np.prod(shape2))
                    dataw = np.concatenate([data1, data2])
                    bfw.write(dataw.tobytes())
    return offsets


def combine(pck1, pck2, pltout=None, vars1=None, vars2=None):
    """
    Function to combine fields between plotfiles
    pck1 : PlotfileCooker instance of the first plotfile
    pck2 : PlotfileCooker instance of the second plotfile
    pckout : Output plotfile dir
    vars1 : Variables to keep in plotfile 1
    vars2 : Variables to keep in plotfile 2
    """
    # >>> TODO: Move this to a parse_input function
    # TODO: make work for 2 dimensions
    if pck1.ndims != 3:
        raise NotImplementedError
    # TODO: allow inplace combination
    if pltout is None:
        raise NotImplementedError
    # Check if the plotfiles have the same grid
    assert pck1 == pck2
    # And max AMR Level
    assert pck1.limit_level == pck2.limit_level
    # Also validate each box is in the same file name
    bf_vec = np.vectorize(lambda s:os.path.split(s)[-1])
    for lv in range(pck1.limit_level + 1):
        assert np.array_equal(bf_vec(pck1.cells[lv]['files']),
                              bf_vec(pck2.cells[lv]['files']))
    # TODO: allow selected fields
    if (vars1 is not None) or (vars2 is not None):
        raise NotImplementedError
    # Compute the field list
    # TODO: handle duplicates and assume they have the same
    # value (keep only one of two)
    cb_fields = list(pck1.fields.keys()) + list(pck2.fields.keys())
    nfields = len(cb_fields)
    # <<< End of parse input

    # make the output dir
    pck1.make_dir_tree(pltout)
    # For each AMR level
    for lv in range(pck1.limit_level + 1):
        lvstart = time.time()
        bfiles_paths_1 = np.array(pck1.cells[lv]["files"])
        lv_offsets_1 = np.array(pck1.cells[lv]['offsets'])
        bfiles_paths_2 = np.array(pck2.cells[lv]["files"])
        lv_offsets_2 = np.array(pck2.cells[lv]['offsets'])
        ncells = len(bfiles_paths_1)
        # Should be true but would be bad if not
        assert len(bfiles_paths_1) == len(bfiles_paths_2)
        # All indexes of the boxes at lv
        box_indexes = np.arange(ncells)
        # Multiprocessing args list
        mp_calls = []
        # To find out what goes where after
        box_index_map = []
        # On process per binary file
        for bfile_r1 in np.unique(bfiles_paths_1):
            # This must be the same boolean array 
            # for both plotfiles
            bf_match = bfiles_paths_1 == bfile_r1
            # The other binary file we read
            bfile_r2 = bfiles_paths_2[bf_match][0]
            # Indexes of boxes in the binary file
            bf_indexes = box_indexes[bf_match]
            # Offsets of the boxes in the binaries
            # The order of these can differ between
            # Plotfiles for some reason so we need 
            # to keep track of it
            offsets_bf1 = lv_offsets_1[bf_match]
            offsets_bf2 = lv_offsets_2[bf_match]
            # Store the index of the boxes with the current file
            box_index_map.append(bf_indexes)
            # Path to the combined binary files (for Windows)
            bfile_r1 = os.path.join(os.getcwd(), bfile_r1)
            bfile_r2 = os.path.join(os.getcwd(), bfile_r2)
            # Path to the new binary file
            bfile_w = os.path.join(os.getcwd(),
                                   pltout,
                                   os.path.basename(os.path.split(bfile_r1)[0]),
                                   os.path.basename(bfile_r1))
            mp_call = {"bfile_r1":bfile_r1,
                       "offst_r1":offsets_bf1,
                       "bfile_r2":bfile_r2,
                       "offst_r2":offsets_bf2,
                       "bfile_w":bfile_w}
            mp_calls.append(mp_call)
        # Strain in parallel
        pool = multiprocessing.Pool()
        new_offsets = pool.map(parallel_combine_binary_files,
                               mp_calls)
        # Reorder the offsets to match the box order
        mapped_offsets = np.empty(len(pck1.boxes[lv]), dtype=int)
        for file_idxs, offsets in zip(box_index_map, new_offsets):
            mapped_offsets[file_idxs] = offsets
        # Rewrite the cell headers
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
                        min_vals1 = np.array(line1)
                        min_vals2 = np.array(line2)
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
                        max_vals1 = np.array(line1)
                        max_vals2 = np.array(line2)
                        max_vals = np.concatenate([max_vals1, max_vals2])
                        ch_w.write(','.join(max_vals) + ',\n')
        print(f"Strained Level {lv} ({time.time() - lvstart:.2f} s)")
    # Rewrite the global header
    hfile_path = os.path.join(pltout, "Header")
    with open(hfile_path, 'w') as hfile:
        # Plotfile version
        hfile.write(pck1.version)
        # Number of fields
        hfile.write(f"{nfields}" + '\n')
        # Fields
        for f in pck1.fields:
            hfile.write(f + '\n')
        for f in pck2.fields:
            hfile.write(f + '\n')
        # Number of dimensions
        hfile.write(f"{pck1.ndims}\n")
        # Time
        hfile.write(str(pck1.time) + '\n')
        # Max level
        hfile.write(str(pck1.limit_level) + '\n')
        # Lower bounds
        hfile.write(' '.join([str(f) for f in pck1.geo_low]) + '\n')
        # Upper bounds
        hfile.write(' '.join([str(f) for f in pck1.geo_high]) + '\n')
        # Refinement factors
        factors = pck1.factors[:pck1.limit_level + 1]
        hfile.write(' '.join([str(f) for f in factors]) + '\n')
        # Grid sizes
        # Looks like ((0,0,0) (7,7,7) (0,0,0))
        tuples = []
        for lv in range(pck1.limit_level + 1):
            sizes = ",".join([str(s-1) for s in pck1.grid_sizes[lv]])
            if pck1.ndims == 3:
                tup = f"((0,0,0) ({sizes}) (0,0,0))"
            elif pck1.ndims == 2:
                tup = f"((0,0) ({sizes}) (0,0))"
            tuples.append(tup)
        hfile.write(' '.join(tuples) + '\n')
        # By level step numbers
        step_numbers = pck1.step_numbers[:pck1.limit_level + 1]
        hfile.write(' '.join([str(n) for n in step_numbers]) + '\n')
        # Grid resolutions
        for lv in range(pck1.limit_level + 1):
            hfile.write(' '.join([str(dx) for dx in pck1.dx[lv]]) + '\n')
        # Coordinate system
        hfile.write(str(pck1.sys_coord))
        # Zero for parsing
        hfile.write("0\n")
        # Write the boxes
        for lv in range(pck1.limit_level + 1):
            # Write the level info
            hfile.write(f"{lv} {len(pck1.boxes[lv])} {pck1.time}\n")
            # Write the level step
            hfile.write(f"{pck1.step_numbers[lv]}\n")
            # Write the boxes
            for box in pck1.boxes[lv]:
                for d in range(pck1.ndims):
                    hfile.write(f"{box[d][0]} {box[d][1]}\n")
            # Write the Level path info
            hfile.write(f"Level_{lv}/Cell\n")

