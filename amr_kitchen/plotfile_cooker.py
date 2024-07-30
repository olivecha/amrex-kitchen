import os
import shutil
import traceback
import linecache
import numpy as np
from tqdm import tqdm
from amr_kitchen.utils import TastesBadError


class PlotfileCooker(object):

    def __init__(self, plotfile, limit_level=None, header_only=False,
                 validate_mode=False,  maxmins=False, ghost=False):
        """
        Parse the header data and save as attributes
        """
        # TODO: Add a description of what the argument do to
        # the docstring
        self.pfile = plotfile
        filepath = os.path.join(plotfile, 'Header')
        with open(filepath) as hfile:
            self.version = hfile.readline()
            # field names
            self.nvars = int(hfile.readline())
            self.fields = {}
            for i in range(self.nvars):
                self.fields[hfile.readline().replace('\n', '')] = i
            # General data
            self.ndims = int(hfile.readline())
            self.time = float(hfile.readline())
            self.max_level = int(hfile.readline())
            self.geo_low = [float(n) for n in hfile.readline().split()]
            self.geo_high = [float(n) for n in hfile.readline().split()]
            self.factors = [int(n) for n in hfile.readline().split()]
            self.grid_sizes = []
            for block in hfile.readline().split()[1::3]:
                grid_size = np.array(block.replace('(', '').replace(")", '').split(','), dtype=int)
                self.grid_sizes.append(grid_size + 1)
            self.step_numbers = [int(n) for n in hfile.readline().split()]
            # Grid resolutions
            resolutions = []
            for i in range(self.max_level + 1):
                resolutions.append([float(n) for n in hfile.readline().split()])
            self.dx = resolutions
            # Coordinate system
            self.sys_coord = hfile.readline()
            # Sanity check
            assert 0 == int(hfile.readline())
            # Define the max level we read
            if limit_level is None:
                self.limit_level = self.max_level
            else:
                self.limit_level=limit_level
            # Read the box geometry
            try:
                self.box_centers, self.boxes = self.read_boxes(hfile)
            except Exception as e:
                # If the class is created from a Taster class
                if validate_mode:
                    # Get the actual exception string
                    catched_tback = traceback.format_exc()
                    raise TastesBadError((f"PlotfileCooker encountered a fatal"
                                          f" exception while reading the boxes"
                                           " coordinates in the method self.read_boxes."
                                           " This could be due to missing or badly"
                                           " formated box data. The exception message is:"
                                          f" {catched_tback}"))
                else:
                    raise e
        # Read the cell data
        if not header_only:
            try:
                self.cells = self.read_cell_headers(maxmins, validate_mode)
            except Exception as e:
                if validate_mode:
                    catched_tback = traceback.format_exc()
                    raise TastesBadError((f"PlotfileCooker encountered a fatal"
                                          f" exception while reading the binary"
                                           " paths and global grid indices in the level"
                                           " headers, inside the method self.read_cell_headers."
                                           " This could be due to missing or badly"
                                           " formated box data. The exception message is:\n"
                                          f" \n {catched_tback}"))
                else:
                    raise e
        # Gets the number fields in the plt_file
        self.nfields = len(self.fields)
        # Compute the ghost boxes map around each box
        if ghost:
            if self.ndims == 3:
                self.box_arrays, self.barr_indices = self.compute_box_array()
                self.ghost_map = self.compute_ghost_map()
            else:
                raise ValueError(("Ghost boxes are not available for plotfiles with"
                                  " ndims < 3"))

    def __eq__(self, other):
        """
        Overload the '==' operator to use it to test for plotfile
        compatibility. This tests that both plotfiles have the same
        mesh refinement structure but allows different number of fields
        and different binary file distribution
        Example:
        hdr1 = PlotfileCooker(plt1000)
        hdr2 = PlotfileCooker(plt2000)
        hdr1 == hdr2 is True if both plotfiles have the same boxes at
        each AMR level
        """
        # Fail if the maximum AMR level is different
        if self.limit_level != other.limit_level:
            return False
        # Compare boxes
        for lv in range(self.limit_level + 1):
            if not np.allclose(self.boxes[lv], other.boxes[lv]):
                return False
        # Compare cell indexes
        for lv in range(self.limit_level + 1):
            if not np.allclose(self.cells[lv]['indexes'],
                               other.cells[lv]['indexes']):
                return False
        # Compare binary files
        # for lv in range(self.limit_level + 1):
        #     if 
        return True

    def read_boxes(self, hfile):
        """
        Read the AMR boxes geometry in the base header file
        """
        # dicts to store box bounds and centers
        points = []
        boxes = []
        self.npoints = []
        self.cell_paths = []
        # Loop over the grid levels
        for lv in range(self.limit_level + 1):
            # Read level and number of cells
            current_level, n_cells, _ = [n for n in hfile.readline().split()]
            current_level = int(current_level)
            n_cells = int(n_cells)
            # Store the lowest level step number
            if int(current_level) == 0:
                self.step = hfile.readline()
            else:
                hfile.readline()
            # Sanity check
            assert current_level == lv
            # Key for the dict
            self.npoints.append(n_cells)
            lv_points = []
            lv_boxes = []
            for i in range(n_cells):
                point = []
                box = []
                for i in range(self.ndims):
                    lo, hi = [float(n) for n in hfile.readline().split()]
                    box.append([lo, hi])
                    point.append(lo + (hi - lo)/2)
                lv_points.append(point)
                lv_boxes.append(box)
            cell_dir = hfile.readline().split('/')[0]
            self.cell_paths.append(cell_dir)
            points.append(lv_points)
            boxes.append(lv_boxes)
        return points, boxes

    def read_cell_headers(self, maxmins, validate_mode):
        """
        Read the cell header data and the maxs/mins for a given level
        """
        cells = []
        all_maxs = []
        all_mins = []
        for i in range(self.limit_level + 1):
            lvcells = {}
            all_maxs.append({})
            all_mins.append({})
            cfile_path = os.path.join(self.pfile, self.cell_paths[i], "Cell_H")
            with open(cfile_path) as cfile:
                # Skip 2 lines
                cfile.readline()
                cfile.readline()
                # Are we good
                assert int(cfile.readline()) == len(self.fields)
                cfile.readline()
                n_cells = int(cfile.readline().split()[0].replace('(', ''))
                indexes = []
                for _ in range(n_cells):
                    #try:
                    start, stop, _ = cfile.readline().split()
                    #except ValueError as e:
                        #raise TastesBadError("")
                    start = np.array(start.replace('(', '').replace(')', '').split(','), dtype=int)
                    stop = np.array(stop.replace('(', '').replace(')', '').split(','), dtype=int)
                    indexes.append([start, stop])
                lvcells["indexes"] = indexes
                cfile.readline()
                assert n_cells == int(cfile.readline())
                files = []
                offsets = []
                for _ in range(n_cells):
                    _, file, offset = cfile.readline().split()
                    files.append(os.path.join(self.pfile, self.cell_paths[i], file))
                    offsets.append(int(offset))
                if maxmins:
                    lvmaxs = []
                    lvmins = []
                    cfile.readline()
                    cfile.readline()
                    for _ in range(n_cells):
                        mins_str = cfile.readline().split(',')
                        lvmins.append(np.array(mins_str[:-1], dtype=float))
                    cfile.readline()
                    cfile.readline()
                    for _ in range(n_cells):
                        maxs_str = cfile.readline().split(',')
                        lvmaxs.append(np.array(maxs_str[:-1], dtype=float))
            lvcells["files"] = files
            lvcells["offsets"] = offsets
            if maxmins:
                lvcells['mins'] = {}
                lvcells['maxs'] = {}
                for field, minvals, maxvals in zip(self.fields, 
                                                   np.transpose(lvmins),
                                                   np.transpose(lvmaxs)):
                    lvcells['mins'][field] = minvals
                    lvcells['maxs'][field] = maxvals
            cells.append(lvcells)
        return cells

    def compute_box_array(self):
        """
        Compute a Nx * Ny * Nz array defining the
        adjacency of the boxes.
        Nx is equal to the number of cells in the
        x direction divided by the smallest box shape
        """
        # Cell resolution in each direction
        box_shapes = self.unique_box_shapes()
        #box_rez = np.min(box_shapes, axis=0)
        box_rez = np.min(box_shapes)
        box_arrays = []
        box_array_indices = []
        for lv in range(self.limit_level + 1):
            box_array_shape = self.grid_sizes[lv] // box_rez
            box_array = -1 * np.ones(box_array_shape, dtype=int)
            lv_barray_indices = []
            for i, idx in enumerate(self.cells[lv]["indexes"]):
                bidx_lo = idx[0] // box_rez
                bidx_hi = idx[1] // box_rez
                box_array[bidx_lo[0]:bidx_hi[0] + 1,
                          bidx_lo[1]:bidx_hi[1] + 1,
                          bidx_lo[2]:bidx_hi[2] + 1] = i
                
                lv_barray_indices.append([bidx_lo, bidx_hi])
            box_arrays.append(box_array)
            box_array_indices.append(lv_barray_indices)
        return box_arrays, box_array_indices

    def compute_ghost_map(self):
        """
        This computes indices of the boxes adjacent
        to a given box. Indices have shape 3x2 for the
        low and high faces of every dimension. If no box
        is adjacent in a given direction the index is set
        to None
        """
        ghost_map = []
        for lv in range(self.limit_level + 1):
            lv_gmap = []
            barr_shape = self.box_arrays[lv].shape
            for box_index, indices in enumerate(self.barr_indices[lv]):
                gmap = [[[], []], [[], []], [[], []]]
                for coo in range(self.ndims):
                    idx_lo = np.copy(indices)
                    idx_lo[0][coo] = max(idx_lo[0][coo] - 1, 0)
                    for bid in np.unique(self.box_arrays[lv][idx_lo[0][0]:idx_lo[1][0],
                                                             idx_lo[0][1]:idx_lo[1][1],
                                                             idx_lo[0][2]:idx_lo[1][2]]):
                        if bid != box_index:
                            gmap[coo][0].append(bid)

                    idx_hi = np.copy(indices)
                    idx_hi[1] += 1
                    idx_hi[1][coo] = min(idx_hi[1][coo] + 1, barr_shape[coo] - 1)
                    for bid in np.unique(self.box_arrays[lv][idx_hi[0][0]:idx_hi[1][0],
                                                             idx_hi[0][1]:idx_hi[1][1],
                                                             idx_hi[0][2]:idx_hi[1][2]]):
                        if bid != box_index:
                            gmap[coo][1].append(bid)
                lv_gmap.append(gmap)
            ghost_map.append(lv_gmap)
        return ghost_map

    def field_index(self, field):
        """ return the index of a data field """
        for i, f in enumerate(self.fields):
            if f == field:
                return i
        raise ValueError(f"""Field {field} was not found in file. 
                             Available fields in {self.pfile.split('/')[-1]} are:
                             {', '.join(self.fields.keys())} and grid_level""")


    def make_dir_tree(self, outpath, limit_level=None):
        """
        Re-Create the tree structure of the plotfile in :outpath:
        """
        if limit_level is None:
            limit_level = self.limit_level
        os.makedirs(os.path.join(os.getcwd(),outpath), exist_ok=True)
        #shutil.copy(os.path.join(self.pfile, 'Header'),
        #           outpath)
        for pth in self.cell_paths[:limit_level + 1]:
            level_dir = pth
            os.makedirs(os.path.join(os.getcwd(),outpath, level_dir), exist_ok=True)
            #shutil.copy(os.path.join(self.pfile, pth + '_H'),
            #            os.path.join(outpath, level_dir))

    def bybinfile(self, lv):
        """
        Iterate over header data at lv
        by individual binary files
        """
        bfiles = np.array(self.cells[lv]['files'])
        indexes = np.array(self.cells[lv]['indexes'])
        offsets = np.array(self.cells[lv]['offsets'])

        box_indexes = np.arange(len(bfiles))
        for bf in np.unique(bfiles):
            bf_indexes = box_indexes[bfiles == bf]
            yield (bf,
                   offsets[bf_indexes],
                   indexes[bf_indexes],)

    def bybinfile_indexed(self, lv):
        """
        Iterate over header data at lv
        by individual binary files
        """
        bfiles = np.array(self.cells[lv]['files'])
        indexes = np.array(self.cells[lv]['indexes'])
        offsets = np.array(self.cells[lv]['offsets'])

        box_indexes = np.arange(len(bfiles))
        for bf in np.unique(bfiles):
            bf_indexes = box_indexes[bfiles == bf]
            yield (bf,
                   offsets[bf_indexes],
                   indexes[bf_indexes],
                   box_indexes)

    def bybox(self, lv):
        """
        Iterate over header data for evey box
        """
        bfiles = np.array(self.cells[lv]['files'])
        indexes = np.array(self.cells[lv]['indexes'])
        offsets = np.array(self.cells[lv]['offsets'])

        for bf, idx, off in zip(bfiles, indexes, offsets):
            yield {"indexes":idx,
                   "bfile":bf,
                   "off":off}

    def byboxcompared(self, other, lv):
        """
        Generator to iterate over the boxes in two plotfiles for
        a given AMR level: lv
        """
        for bf1, bf2,  off1, off2, idxs in zip(self.cells[lv]['files'],
                                               other.cells[lv]['files'],
                                               self.cells[lv]['offsets'],
                                               other.cells[lv]['offsets'],
                                               self.cells[lv]['indexes']):
            output = {"indexes":idxs,
                      "bfile1":bf1,
                      "bfile2":bf2,
                      "off1":off1,
                      "lv":lv,
                      "off2":off2}
            yield output

    def boxesfromindexes(self, indexes):
        """
        Give a list if indexes with shape n_levels x [n_indexes_at_level]
        Compute the corresponding bounding boxes using the header data
        """
        all_boxes = []
        for lv in range(self.limit_level + 1):
            lv_boxes = []
            xgrid = np.linspace(self.geo_low[0] + self.dx[lv][0]/2, 
                                self.geo_high[0] - self.dx[lv][0]/2,
                                self.grid_sizes[lv][0])
            ygrid = np.linspace(self.geo_low[0] + self.dx[lv][0]/2, 
                                self.geo_high[0] - self.dx[lv][0]/2,
                                self.grid_sizes[lv][0])
            zgrid = np.linspace(self.geo_low[0] + self.dx[lv][0]/2, 
                                self.geo_high[0] - self.dx[lv][0]/2,
                                self.grid_sizes[lv][0])
            hdx = self.dx[lv][0]/2
            hdy = self.dx[lv][1]/2
            hdz = self.dx[lv][2]/2
            for idx in indexes[lv]:
                box_x = [xgrid[idx[0][0]] - hdx, xgrid[idx[1][0]] + hdx]
                box_y = [ygrid[idx[0][1]] - hdy, ygrid[idx[1][1]] + hdy]
                box_z = [zgrid[idx[0][2]] - hdz, zgrid[idx[1][2]] + hdz]
                box = [box_x, box_y, box_z]
                lv_boxes.append(box)
            all_boxes.append(lv_boxes)
        return all_boxes

    def writehdrnewboxes(self, pfdir, boxes, fields):
        """
        Write the global header with new boxes
        """
        if pfdir not in os.listdir():
            os.makedirs(os.getcwd(),pfdir)

        with open(os.path.join(os.getcwd(),pfdir, 'Header'), 'w') as hfile:
            # Plotfile version
            hfile.write(self.version)
            # Number of fields
            hfile.write(f"{len(fields)}\n")
            # Fields
            for f in fields:
                hfile.write(f + '\n')
            # Dimension
            hfile.write(f"{self.ndims}\n")
            # Time is unknown
            hfile.write("0.0\n")
            # Max level
            hfile.write(str(self.limit_level) + '\n')
            # Lower bounds
            lo_str = " ".join([f"{self.geo_low[i]}" for i in range(self.ndims)])
            hfile.write(lo_str + '\n')
            # Upper bounds
            hi_str =  " ".join([f"{self.geo_high[i]}" for i in range(self.ndims)])
            hfile.write(hi_str + '\n')
            # Refinement factors
            factors = self.factors[:self.limit_level]
            hfile.write(' '.join([str(f) for f in factors]) + '\n')
            # Grid sizes
            # Looks like ((0,0,0) (7,7,7) (0,0,0))
            tuples = []
            for lv in range(self.limit_level + 1):
                start = ','.join(['0' for _ in range(self.ndims)])
                cente = ','.join([str(self.grid_sizes[lv][i] - 1) for i in range(self.ndims)])
                end = start
                tup = f"(({start}) ({cente}) ({end}))"
                tuples.append(tup)
            hfile.write(' '.join(tuples) + '\n')
            # By level step numbers (all zero)
            step_numbers = [0 for _ in range(self.limit_level + 1)]
            hfile.write(' '.join([str(n) for n in step_numbers]) + '\n')
            # Grid resolutions
            for lv in range(self.limit_level + 1):
                hfile.write(' '.join([f"{self.dx[lv][i]}" for i in range(self.ndims)]) + '\n')
            # Coordinate system
            hfile.write(str(self.sys_coord))
            # Zero for parsing
            hfile.write("0\n")
            # Write the boxes
            for lv in range(self.limit_level + 1):
                # Write the level info
                hfile.write(f"{lv} {len(boxes[lv])} 0.0\n")
                # Write the level step
                hfile.write(f"0\n")
                # Write the 2D boxes
                for box in boxes[lv]:
                    for i in range(self.ndims):
                        hfile.write(f"{box[i][0]} {box[i][1]}\n")
                # Write the Level path info
                hfile.write(f"Level_{lv}/Cell\n")

    def unique_box_shapes(self):
        """
        Find the unique box shape tuples
        for each level
        """
        shapes = []
        for lv in range(self.limit_level + 1):
            for idx in self.cells[lv]['indexes']:
                shape = idx[1] - idx[0] + 1
                shapes.append(tuple(shape))
        shapes = np.unique(shapes, axis=0)
        shapes = [tuple(shape) for shape in shapes]
        return shapes
