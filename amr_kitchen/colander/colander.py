import os
import time
import multiprocessing
import numpy as np
from amr_kitchen import PlotfileCooker

def parallel_strain_3d(args):
    """
    Cook a single cell binary file
    Multiprocessing function
    """
    # New offsets
    offsets = []
    nkept = len(args["kept_fields"])
    # Open the read and write
    with open(args['bfile_r'], 
              'rb') as bfr, open(args['bfile_w'], 
                                 'wb') as bfw:
        for indexes, fst_r, idx in zip(args['box_indexes'], 
                                       args['offsets_r'], 
                                       args['cell_indexes']):
            # Store pointer position to update header
            offsets.append(bfw.tell())
            # Go to the data
            bfr.seek(fst_r)
            # Get the header
            header = bfr.readline().decode('ascii')
            # Replace with number of vars just to be sure
            header_w = header.replace(f'{args["nvars"]}\n', f'{nkept}\n')
            # Write to binary file
            bfw.write(header_w.encode('ascii'))
             # Read the data
            shape = indexes[1] - indexes[0] + 1
            total_shape = (shape[0], shape[1], shape[2], args['nvars'])
            arr = np.fromfile(bfr, "float64", np.prod(total_shape))
            arr = arr.reshape(total_shape, order="F")
            # Reshape into dicts
            arr_out = arr[:, :, :, args["kept_fields"]]
            arr_bytes = arr_out.flatten(order="F").tobytes()
            bfw.write(arr_bytes)
    return offsets


def parallel_strain_2d(args):
    """
    Cook a single cell binary file
    Multiprocessing function
    """
    # New offsets
    offsets = []
    nkept = len(args["kept_fields"])
    # Open the read and write
    with open(args['bfile_r'], 
              'rb') as bfr, open(args['bfile_w'], 
                                 'wb') as bfw:        
        for indexes, fst_r, idx in zip(args['box_indexes'], 
                                       args['offsets_r'], 
                                       args['cell_indexes']):
            # Store pointer position to update header
            offsets.append(bfw.tell())
            # Go to the data
            bfr.seek(fst_r)
            # Get the header
            header = bfr.readline().decode('ascii')
            # Replace with number of vars just to be sure
            header_w = header.replace(f'{args["nvars"]}\n', f'{nkept}\n')
            # Write to binary file
            bfw.write(header_w.encode('ascii'))
             # Read the data
            shape = indexes[1] - indexes[0] + 1
            total_shape = (shape[0], shape[1], args['nvars'])
            arr = np.fromfile(bfr, "float64", np.prod(total_shape))
            arr = arr.reshape(total_shape, order="F")
            # Reshape into dicts
            arr_out = arr[:, :, args["kept_fields"]]
            arr_bytes = arr_out.flatten(order="F").tobytes()
            bfw.write(arr_bytes)
    return offsets


class Colander(PlotfileCooker):
    """
    Class containing the data used to filter the plotfile
    """

    def __init__(self,
                 plotfile=None,
                 limit_level=None,
                 output=None,
                 variables=None):
        """
        Constructor for the Colander data container
        args : parsed command line arguments
        """
        super().__init__(plotfile,
                         limit_level=limit_level)
        # Define the strainer function
        if self.ndims == 2:
            self.strainer = parallel_strain_2d
        elif self.ndims == 3:
            self.strainer = parallel_strain_3d
        # Store user inputs
        self.outdir = output
        self.kept_fields = []
        # Case for keeping all variables
        if (len(variables) == 1) and (variables[0] == "all"):
            self.kept_names = []
            for f in self.fields:
                self.kept_fields.append(self.fields[f])
                self.kept_names.append(f)
        # Case for arbitrary number of variables
        else:
            self.kept_names = variables
            for v in variables:
                try:
                    self.kept_fields.append(self.fields[v])
                except KeyError:
                    raise ValueError(f"Field {v} not found in file")

    def strain(self):
        """
        Save the strained plotfile
        """
        # Make the output plotfile dir tree 
        self.make_dir_tree(self.outdir)
        # For each level
        for lv in range(self.limit_level + 1):
            lvstart = time.time()
            level_files = np.array(self.cells[lv]["files"])
            ncells = len(level_files)
            cell_header_r = os.path.join(self.pfile,
                                         self.cell_paths[lv], 
                                         'Cell_H')
            # All indexes of the boxes at lv
            box_indexes = np.arange(ncells)

            # Multiprocessing args list
            mp_calls = []
            box_index_map = []
            # On process per binary file
            for bfile_r in np.unique(level_files):
                # Indexes of cells in the binary file
                bf_indexes = box_indexes[level_files == bfile_r]
                # Store the index of the boxes with the current file
                box_index_map.append(bf_indexes)
                # Path to the "old" binary file (for Windows)
                bfile_r = os.path.join(os.getcwd(), bfile_r)
                # Path to the new binary file
                bfile_w = os.path.join(os.getcwd(),self.outdir,
                                       os.path.basename(os.path.split(bfile_r)[0]),
                                       os.path.basename(bfile_r))
                # Indexes of the boxes in the global grid
                box_slices = np.array(self.cells[lv]["indexes"])[bf_indexes]
                # All offsets in the current binary file
                offsets_r = np.array(self.cells[lv]["offsets"])[bf_indexes]
                mp_call = {"bfile_r":bfile_r,
                           "bfile_w":bfile_w,
                           "box_indexes":box_slices,
                           "cell_indexes":bf_indexes,
                           "offsets_r":offsets_r,
                           "nvars":self.nvars,
                           "kept_fields":self.kept_fields,
                           "ncells":ncells}
                mp_calls.append(mp_call)
            # Strain in parallel
            with multiprocessing.Pool() as pool:
                new_offsets = pool.map(self.strainer, mp_calls)
            # Reorder the offsets to match the box order
            mapped_offsets = np.empty(len(self.boxes[lv]), dtype=int)
            for file_idxs, offsets in zip(box_index_map, new_offsets):
                mapped_offsets[file_idxs] = offsets

            # Rewrite the cell headers
            self.update_cell_header(lv, cell_header_r, mapped_offsets)
            print(f"Strained Level {lv} ({time.time() - lvstart:.2f} s)")
        # Rewrite the global header
        self.write_strained_global_header()

    
    def update_cell_header(self, lv, cell_header_r, new_offsets):
        """
        Update the new Cell header
        """
        cell_header_w = os.path.join(self.outdir, 
                                     self.cell_paths[lv],
                                     "Cell_H")

        with open(os.path.join(os.getcwd(),cell_header_w), 'w') as ch_w, open(os.path.join(os.getcwd(),cell_header_r), 'r') as ch_r:
            # First two lines
            for i in range(2):
                l = ch_r.readline()
                ch_w.write(l)
            # Number of fields
            _ = ch_r.readline()
            ch_w.write(f"{len(self.kept_fields)}\n")
            # Mesh stays the same
            while True:
                l = ch_r.readline()
                if "FabOnDisk:" in l:
                    new_l = l.split()[:-1]
                    new_l.append(str(new_offsets[0]))
                    ch_w.write(' '.join(new_l) + "\n")
                    break
                else:
                    ch_w.write(l)
            # Write the cell indexes
            for fst in new_offsets[1:]:
                l = ch_r.readline()
                new_l = l.split()[:-1]
                new_l.append(str(fst))
                ch_w.write(' '.join(new_l) + "\n")
            # Blank line
            ch_w.write(ch_r.readline())
            line = ch_r.readline()
            ncells, _ = line.split(',')
            ch_w.write(f"{ncells},{len(self.kept_fields)}\n")
            for li in range(int(ncells)):
                line = ch_r.readline().split(',')[:-1]
                min_vals = np.array(line)
                ch_w.write(','.join(min_vals[self.kept_fields]) + ',\n')
            # Blank line
            ch_w.write(ch_r.readline())
            ncells, _ = ch_r.readline().split(',')
            ch_w.write(f"{ncells},{len(self.kept_fields)}\n")
            for li in range(int(ncells)):
                line = ch_r.readline().split(',')[:-1]
                max_vals = np.array(line)
                ch_w.write(','.join(max_vals[self.kept_fields]) + ',\n')

    def write_strained_global_header(self):
        """
        Write the plotfile header for the plotfile with filtered
        data
        """
        hfile_path = os.path.join(self.outdir, "Header")
        with open(hfile_path, 'w') as hfile:
            # Plotfile version
            hfile.write(self.version)
            # Number of fields
            hfile.write(str(len(self.kept_fields)) + '\n')
            # Fields
            for v in self.kept_names:
                hfile.write(v + '\n')
            # Number of dimensions
            hfile.write(f"{self.ndims}\n")
            # Time
            hfile.write(str(self.time) + '\n')
            # Max level
            hfile.write(str(self.limit_level) + '\n')
            # Lower bounds
            hfile.write(' '.join([str(f) for f in self.geo_low]) + '\n')
            # Upper bounds
            hfile.write(' '.join([str(f) for f in self.geo_high]) + '\n')
            # Refinement factors
            if self.limit_level > 0:
                factors = self.factors[:self.limit_level + 1]
                hfile.write(' '.join([str(f) for f in factors]) + '\n')
            else:
                hfile.write('\n')
            # Grid sizes
            # Looks like ((0,0,0) (7,7,7) (0,0,0))
            tuples = []
            for lv in range(self.limit_level + 1):
                sizes = ",".join([str(s-1) for s in self.grid_sizes[lv]])
                if self.ndims == 3:
                    tup = f"((0,0,0) ({sizes}) (0,0,0))"
                elif self.ndims == 2:
                    tup = f"((0,0) ({sizes}) (0,0))"
                tuples.append(tup)
            hfile.write(' '.join(tuples) + '\n')
            # By level step numbers
            step_numbers = self.step_numbers[:self.limit_level + 1]
            hfile.write(' '.join([str(n) for n in step_numbers]) + '\n')
            # Grid resolutions
            for lv in range(self.limit_level + 1):
                hfile.write(' '.join([str(dx) for dx in self.dx[lv]]) + '\n')
            # Coordinate system
            hfile.write(str(self.sys_coord))
            # Zero for parsing
            hfile.write("0\n")
            # Write the boxes
            for lv in range(self.limit_level + 1):
                # Write the level info
                hfile.write(f"{lv} {len(self.boxes[lv])} {self.time}\n")
                # Write the level step
                hfile.write(f"{self.step_numbers[lv]}\n")
                # Write the boxes
                for box in self.boxes[lv]:
                    for d in range(self.ndims):
                        hfile.write(f"{box[d][0]} {box[d][1]}\n")
                # Write the Level path info
                hfile.write(f"Level_{lv}/Cell\n")


