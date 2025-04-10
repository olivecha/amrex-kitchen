import os
import shutil
import traceback
import multiprocessing
import numpy as np
from tqdm import tqdm
from scipy.ndimage import map_coordinates, zoom
from amr_kitchen.utils import TastesBadError, shape_from_header

def mp_read_box_single_field(args):
    with open(args[0], 'rb') as bf:
        bf.seek(args[1])
        shape = shape_from_header(bf.readline().decode('ascii'))
        bf.seek(np.prod(shape[:-1]) * args[2] * 8, 1)
        data = np.fromfile(bf, 'float64', np.prod(shape[:-1]))
    return data.reshape(shape[:-1], order='F')

def mp_read_box_slice_field(args):
    with open(args[0], 'rb') as bf:
        bf.seek(args[1])
        shape = shape_from_header(bf.readline().decode('ascii'))
        start, stop = args[2].indices(shape[-1])[:2]
        slice_size = stop - start
        bf.seek(np.prod(shape[:-1]) * start * 8, 1)
        data = np.fromfile(bf, 'float64', np.prod(shape[:-1]) * slice_size)
    data = data.reshape(np.append(shape[:-1], slice_size), order='F')
    return data[..., args[2]]

def mp_read_box_index_field(args):
    diff = args[2][-1] - args[2][0] + 1
    with open(args[0], 'rb') as bf:
        bf.seek(args[1])
        shape = shape_from_header(bf.readline().decode('ascii'))
        bf.seek(np.prod(shape[:-1]) * args[2][0] * 8, 1)
        data = np.fromfile(bf, 'float64', np.prod(shape[:-1])*diff)
    data = data.reshape(np.append(shape[:-1], diff), order='F')
    return data[..., np.array(args[2]) - args[2][0]]

def mp_read_bfile_single_field(args):
    file_data = []
    with open(args[0], 'rb') as bf:
        while True:
            try:
                shape = shape_from_header(bf.readline().decode('ascii'))
                bf.seek(np.prod(shape[:-1]) * args[1] * 8, 1)
                data = np.fromfile(bf, 'float64', np.prod(shape[:-1]))
                bf.seek(np.prod(shape[:-1]) * (shape[-1] - args[1] - 1) * 8, 1)
                file_data.append(data.reshape(shape[:-1], order='F'))
            except:
                break
    return file_data

def mp_read_bfile_slice_field(args):
    file_data = []
    with open(args[0], 'rb') as bf:
        while True:
            try:
                shape = shape_from_header(bf.readline().decode('ascii'))
                start, stop = args[1].indices(shape[-1])[:2]
                slice_size = stop - start
                bf.seek(np.prod(shape[:-1]) * start * 8, 1)
                data = np.fromfile(bf, 'float64', np.prod(shape[:-1]) * slice_size)
                bf.seek(np.prod(shape[:-1]) * (shape[-1] - slice_size - start) * 8, 1)
                data = data.reshape(np.append(shape[:-1], slice_size), order='F')
                file_data.append(data[..., args[1]])
            except Exception as e:
                print(type(e), e)
                break
    return file_data

def mp_read_bfile_index_field(args):
    diff = args[1][-1] - args[1][0] + 1
    file_data = []
    with open(args[0], 'rb') as bf:
        while True:
            try:
                shape = shape_from_header(bf.readline().decode('ascii'))
                bf.seek(np.prod(shape[:-1]) * args[1][0] * 8, 1)
                data = np.fromfile(bf, 'float64', np.prod(shape[:-1]) * diff)
                bf.seek(np.prod(shape[:-1]) * (shape[-1] - args[1][-1] - 1) * 8, 1)
                data = data.reshape(np.append(shape[:-1], diff), order='F')
                file_data.append(data[..., np.array(args[1]) - args[1][0]])
            except:
                break
    return file_data

class LevelDataIterator(object):

    def __init__(self, fun, bfiles, field_arg):
        pool = multiprocessing.Pool()
        self.iterator = pool.imap(fun,
                                  zip(bfiles,
                                      [field_arg]*len(bfiles)))
        self._data = self.iterator.__next__().__iter__()

    def __iter__(self):
        return self

    def __next__(self):
        try:
            return self._data.__next__()
        except StopIteration:
            self._data = self.iterator.__next__().__iter__()
            return self._data.__next__()

class LevelDataStream(object):

    def __init__(self, bfiles, offsets, field_arg):
        self.bfiles = np.array(bfiles)
        self.offsets = np.array(offsets)
        self.size = len(bfiles)
        self.farg = field_arg
        if isinstance(self.farg, int):
            self.read_fun = mp_read_box_single_field
            self.file_fun = mp_read_bfile_single_field
        elif isinstance(self.farg, slice):
            self.read_fun = mp_read_box_slice_field
            self.file_fun = mp_read_bfile_slice_field
        elif (isinstance(self.farg, list) or
              isinstance(self.farg, np.ndarray)):
            self.farg = np.array(self.farg)
            assert self.farg.ndim == 1, "Field slice indices must be one dimensional"
            self.read_fun = mp_read_box_index_field
            self.file_fun = mp_read_bfile_index_field

    def __getitem__(self, idx):
        if isinstance(idx, int):
            return self.read_fun((self.bfiles[idx],
                                  self.offsets[idx],
                                  self.farg))
        elif isinstance(idx, slice):
            slice_size = len(range(*idx.indices(self.size)))
            pool = multiprocessing.Pool()
            return pool.map(self.read_fun,
                            zip(self.bfiles[idx],
                                self.offsets[idx],
                                [self.farg]*slice_size))
        elif (isinstance(idx, list) or
              isinstance(idx, np.ndarray)):
            if len(idx) == 0:
                return []
            idx = np.array(idx)
            assert idx.ndim == 1, "Box slice indices must be one dimensional"
            pool = multiprocessing.Pool()
            if idx.dtype == int:
                count = len(idx)
            elif idx.dtype == bool:
                count = np.count_nonzero(idx)
            return pool.map(self.read_fun,
                            zip(self.bfiles[idx],
                                self.offsets[idx],
                                [self.farg]*count))
    def __iter__(self):
        return LevelDataIterator(self.file_fun,
                                 np.unique(self.bfiles),
                                 self.farg)
    def iter(self, idx):
        """
        Manual data iterator to support reading data
        on the fly for slices
        """
        if isinstance(idx, int):
            return self.read_fun((self.bfiles[idx],
                                  self.offsets[idx],
                                  self.farg))
        elif isinstance(idx, slice):
            slice_size = len(range(*idx.indices(self.size)))
            pool = multiprocessing.Pool()
            return pool.imap(self.read_fun,
                            zip(self.bfiles[idx],
                                self.offsets[idx],
                                [self.farg]*slice_size))
        elif (isinstance(idx, list) or
              isinstance(idx, np.ndarray)):
            if len(idx) == 0:
                return []
            idx = np.array(idx)
            assert idx.ndim == 1, "Box slice indices must be one dimensional"
            pool = multiprocessing.Pool()
            if idx.dtype == int:
                count = len(idx)
            elif idx.dtype == bool:
                count = np.count_nonzero(idx)
            return pool.imap(self.read_fun,
                            zip(self.bfiles[idx],
                                self.offsets[idx],
                                [self.farg]*count))

class LevelDataSelector(object):

    def __init__(self, fields, cells, field_arg, limit_level, boxes = None, dx = None):
        # Convert key to field index
        if isinstance(field_arg, str):
            field_arg = fields[field_arg]
        # Also for tuples of keys
        elif ((isinstance(field_arg, list) or
               isinstance(field_arg, np.ndarray)) and
               isinstance(field_arg[0], str)):
            field_arg = [fields[fname] for fname in field_arg]
        try:
            _ = np.array(list(fields.keys()))[field_arg]
            self.farg = field_arg
        except IndexError:
            raise IndexError((f"The field indexing [{field_arg}] is not"
                              f" compatible with the number of fields"
                              f" in the plotfile ({len(fields)})"))
        self.cells = cells
        self.fields = fields
        self.limit_level = limit_level
        self.boxes = boxes
        self.dx = dx

    def __getitem__(self, key):
        if key > self.limit_level:
            raise ValueError((f"The maximum AMR level of the plotfile"
                              f" is {self.limit_level}"))
        return LevelDataStream(self.cells[key]['files'],
                               self.cells[key]['offsets'],
                               self.farg)

    def __call__(self, *args):

        point = np.array(args)

        # Input validation
        if len(point) not in [2, 3]:
            raise KeyError("Please enter point with valid 2d (x,y) or 3d (x,y,z) format")

        # Finds boxes containing point at each level
        #box_list_less = []
        #box_list_more = []
        box_matches_exact = {}
        box_matches_inner = {}
        box_matches_outer = {}
        for level in range(self.limit_level + 1):
            # Convert box limits to numpy array
            boxes = np.array(self.boxes[level])
            # grid resolution at level
            dx = self.dx[level]
            # Boxes where the point is contained within the box boundary:
            box_match_exact = (boxes[:, 0, 0] <= point[0]) & \
                              (point[0] <= boxes[:, 0, 1]) & \
                              (boxes[:, 1, 0] <= point[1]) & \
                              (point[1] <= boxes[:, 1, 1]) & \
                              (boxes[:, 2, 0] <= point[2]) & \
                              (point[2] <= boxes[:, 2, 1] )
            box_matches_exact[level] = np.nonzero(box_match_exact)[0]
            # These are boxes where the point is contained within the bounding
            # cells centers
            box_match_inner = (boxes[:, 0, 0] + dx[0]/2 <= point[0]) & \
                              (point[0] <= boxes[:, 0, 1] - dx[0]/2) & \
                              (boxes[:, 1, 0] + dx[1]/2 <= point[1]) & \
                              (point[1] <= boxes[:, 1, 1] - dx[1]/2) & \
                              (boxes[:, 2, 0] + dx[2]/2 <= point[2]) & \
                              (point[2] <= boxes[:, 2, 1] - dx[2]/2)
            box_matches_inner[level] = np.nonzero(box_match_inner)[0]

            # These are boxes neigboring boxes containing the point when the
            # point is between the cell center of the boundary cells and the box
            # limit
            box_match_outer = (boxes[:, 0, 0] - dx[0]/2 <= point[0]) & \
                              (point[0] <= boxes[:, 0, 1] + dx[0]/2) & \
                              (boxes[:, 1, 0] - dx[1]/2 <= point[1]) & \
                              (point[1] <= boxes[:, 1, 1] + dx[1]/2) & \
                              (boxes[:, 2, 0] - dx[2]/2 <= point[2]) & \
                              (point[2] <= boxes[:, 2, 1] + dx[2]/2)
            box_matches_outer[level] = np.nonzero(box_match_outer)[0]

        # Finest matching level for each box bounds condition
        match_lv_exact = [lv for lv in box_matches_exact if len(box_matches_exact[lv]) != 0][-1]
        try:
            match_lv_inner = [lv for lv in box_matches_inner if len(box_matches_inner[lv]) != 0][-1]
        except IndexError:
            match_lv_inner = None
        match_lv_outer = [lv for lv in box_matches_outer if len(box_matches_outer[lv]) != 0][-1]

        # CASE 1: point can be interpolated using data at cell centers from a single box
        if match_lv_inner == match_lv_exact:
            # Point can only be contained in a single box
            assert len(box_matches_inner[match_lv_inner]) == 1
            # Validate we only match a single box for the other conditions
            assert ((match_lv_inner == match_lv_exact) & \
                    (match_lv_inner == match_lv_outer)), "Something went wrong"
            assert len(box_matches_exact[match_lv_exact]) == 1, "Something went wrong"
            assert len(box_matches_outer[match_lv_outer]) == 1, "Something went wrong"
            # Log some info
            match_box_id = int(box_matches_inner[match_lv_inner][0])
            # dx at the interpolation level
            dx = self.dx[match_lv_inner]
            # Converts point to indices at this level
            # Scales back by 0.5, because the domain starts at dx/2 (index = 0)
            point_idx = (point / dx) - 0.5
            # 3D data for a single box at the finest matching level
            data_arrays = self[match_lv_inner][match_box_id]
            # Indices of the box we just read
            box_indices = self.cells[match_lv_inner]['indexes'][match_box_id]
            # Point coordinates in term of local box indices
            point_local = point_idx - box_indices[0]
            # Requesting a single field
            if isinstance(self.farg, int):
                # interpolate the 3D array at the point
                point_data_field = map_coordinates(data_arrays,
                                                   np.transpose([point_local]))
                return point_data_field
            # Requesting multiple fields
            else:
                # Array to store the output
                point_data = []
                # For each field index
                for fid in range(len(self.farg)):
                    # interpolate each 3D array individually
                    point_data_field = map_coordinates(data_arrays[..., fid],
                                                       np.transpose([point_local]))
                    point_data.append(point_data_field[0])
                return np.array(point_data)

        # CASE 2: when the point is between boxes
        else: 
            # Find out which levels we need to load data from
            # This assume there are now adjacent boxes differing for more
            # than one AMR level (normally enforced by AMReX)
            # A level below the level containing the point
            load_lv_low = np.max([match_lv_exact - 1, 0]) 
            # The finest level of the box neighbouring the point
            load_lv_hi = match_lv_outer
            # Use finest level dx
            dx = self.dx[load_lv_hi]
            # Converts point to indices at this level
            point_idx = (point / dx) - 0.5
            # Find out the shape of the array of concatenated boxes
            all_box_indices = []
            all_box_data = {}
            # For every AMR Level were data is needed
            for lv in range(load_lv_low, load_lv_hi + 1):
                # Factor to scale the indices
                idx_factor = 2 ** (load_lv_hi - lv)
                # Interpolation can be done using outer matching boxes
                box_indices = np.array(self.cells[lv]['indexes'])[box_matches_outer[lv]]
                # Convert to slice notation
                box_indices[:, 1, :] += 1
                # Scale by current lv factor
                all_box_indices.append(box_indices * idx_factor)
                all_box_data[lv] = self[lv][box_matches_outer[lv]]
            # Shape of concatenated boxes
            all_box_indices = np.concatenate(all_box_indices)
            indices_lo = np.min(all_box_indices[:, 0, :], axis=0)
            indices_hi = np.max(all_box_indices[:, 1, :], axis=0)
            conc_shape = indices_hi - indices_lo
            # Point coordinate in the concatenated array
            point_local = point_idx - indices_lo
            # Convert field indices to array for iteration
            if isinstance(self.farg, int):
                field_ids = [0]
            else:
                field_ids = range(len(self.farg))
            # Empty array to store interpolated data
            point_data = []
            # For each requested field
            for fid in field_ids:
                # Empty array to store the box data
                conc_array = np.zeros(conc_shape)
                for lv in range(load_lv_low, load_lv_hi + 1):
                    # Factor to scale the indices
                    idx_factor = 2 ** (load_lv_hi - lv)
                    # Interpolation can be done using outer matching boxes
                    box_indices = np.array(self.cells[lv]['indexes'])[box_matches_outer[lv]]
                    # Convert to slice notation
                    box_indices[:, 1, :] += 1
                    # Scale to finest level
                    box_indices *= idx_factor
                    # Rescale to conc_array shape
                    box_indices[:, :, 0] -= indices_lo[0]
                    box_indices[:, :, 1] -= indices_lo[1]
                    box_indices[:, :, 2] -= indices_lo[2]
                    for idx, box in zip(box_indices, all_box_data[lv]):
                        # Current field data
                        data = box[..., fid]
                        if idx_factor > 1:
                            # TODO: compare results with map_coordinates interpolation
                            data = zoom(data, idx_factor)
                        conc_array[idx[0, 0]: idx[1, 0],
                                   idx[0, 1]: idx[1, 1],
                                   idx[0, 2]: idx[1, 2]] = data
                point_value = map_coordinates(conc_array, np.transpose([point_local]))
                point_data.append(point_value[0])
            return point_data


class PlotfileCooker(object):

    def __init__(self,
                 plotfile_path: str,
                 limit_level: int = None,
                 header_only: bool = False,
                 validate_mode: bool = False,
                 maxmins: bool = False,
                 ghost: bool = False):
        """
        Parse the header data and save as attributes
        ___
        plotfile_path: path to the plotfile directory
        limit_level: maximum adaptive mesh refinement level
                     considered when reading the headers
        header_only: only read the main plotile header (plotfile/Header)
                     (This is much faster than reading all the box data)
        validate_mode: do not fail when an error is encountered
                       (This can be used to find out problems in a plotfile)
        maxmins: if True the maximum and mimimum values of each field in the
                 boxes are read (a bit slower)
        ghost: if True the ghost cells around each box are computed by creating
               3D arrays where the value is the index of the box for each level
        """
        self.pfile = plotfile_path
        filepath = os.path.join(self.pfile, 'Header')
        with open(filepath) as hfile:
            self.version = hfile.readline()
            # field names
            self.nvars = int(hfile.readline())
            self.fields = {}
            for i in range(self.nvars):
                field_name = hfile.readline().replace('\n', '')
                if field_name not in self.fields:
                    self.fields[field_name] = i
                else:
                    repeat_number = 2
                    while True:
                        field_name_extra =  f"{field_name}_{repeat_number}" 
                        if field_name_extra not in self.fields:
                            self.fields[field_name_extra] = i
                            break
                        else:
                            repeat_number += 1

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
            elif limit_level <= self.max_level:
                self.limit_level=limit_level
            else:
                raise ValueError((f"The limit level must be less or equal than"
                                  f" the maximum AMR level of the plotfile:"
                                  f" {limit_level} > {self.max_level}"))
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

        # Compute the global 1D grids
        self.grids = self.compute_global_grids()

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

    """
    Methods defining operator overloading
    """

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
        return True

    def __getitem__(self, key):
        """
        Slicing of plotfile data is performed by returning classes for
        level selection, and then AMR box selection that each provide
        their __getitem__ methods

        The first layer defines which fields are included in the data outout.
        multiple indexing modes are supported:

        PlotfileCooker["temp"] # a single field using the field name key
        PlotfileCooker[3] # A single field using the field index
        PlotfileCooker[:3] # Multiple fields using a slice
        PlotfileCooker[[0, 3, 10]] # Multiple fields using a list of indices

        The second layer defines which AMR Level is selected. This indexing
        operator returns an iterator for the data at the selected level.
        Only integers indices are supported:

        ```
        # Temperature data at level 0:
        PlotfileCooker["temp"][0]
        ```

        ```
        # Finest level Y(OH):
        PlotfileCooker["Y(OH)"][PlotfileCooker.limit_level]
        ```

        ```
        # Iterate over every density AMR box data:
        for rho_box in PlotfileCooker["density"][2]:
            # rho_box is a n dimensional array
            # containing density data of a single AMR box
            pass
        ```

        The third layer defines from which AMR box the data is selected.
        integer, slice and array like indices are supported. If the index
        argument is not an integer, multiprocessing is used to read the data.
        Because the shape of the data is not consistent between boxes, a list
        of arrays is returned for non integer slices.
        The box data shape has the format `(shape_x, shape_y, shape_z, fields)`.

        ```
        # The first and last AMR boxes at a given level:
        T_fist = PlotfileCooker["temp"][lv][0]
        T_last = PlotfileCooker["temp"][lv][-1]
        ```

        ```
        # All velocities in the 5th box at the finest level:
        vel_5 = PlotfileCooker[["x_velocity",
                                "y_velocity",
                                "z_velocity"]][-1][5]

        # Index the box data according to field
        ux = vel_5[..., 0] # ux is a 3D array with the AMR box shape
        uy = vel_5[..., 1]
        uz = vel_5[..., 2]
        ```

        ```
        # Every other box (could be any slice):
        half_boxes = PlotfileCooker["field"][lv][::2]
        ```

        ```
        # Specific boxes using a mask
        pck = PlotfileCooker("plotfile", maxmins=True)
        # Box indices where T_max > 1000
        mask = np.nonzero(pck.cells[-1]["maxs"]["temp"] > 1000)[0]
        # Box data from another field
        Z_data = pck["mixture_fraction"][-1][mask]
        # Perform any computation
        Z_mean = np.mean(np.hstack(Z_data))
        ```

        **Warning:** the iterator returned by `PlotfileCooker["field"][lv]`
        loops over the binary files without preserving the AMR box order
        in the plotfile headers as it is about 10x faster. To preserve box
        order use: `for box_data in PlotfileCooker["field"][lv][:]`.
        If the plotfile is large this might use a lot of memory. 
        Instead, the box indices can be used to read boxes one at a time:
        ```
        for i in range(len(PlotfileCooker.boxes[lv])):
            box_data = PlotfileCooker["field"][lv][i]
        ```
        """
        return LevelDataSelector(self.fields, self.cells, key, self.limit_level, self.boxes, self.dx)

    """
    Method for constructing the class from plotfile mesh data
    """

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
                n_fields_valid = cfile.readline()
                assert int(n_fields_valid) == len(self.fields)
                cfile.readline()
                n_cells = int(cfile.readline().split()[0].replace('(', ''))
                indexes = []
                for _ in range(n_cells):
                    start, stop, _ = cfile.readline().split()
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

    def field_index(self, field):
        """ return the index of a data field """
        # TODO: create a class to raise KeyError on __getitem__
        for i, f in enumerate(self.fields):
            if f == field:
                return i
        raise ValueError(f"""Field {field} was not found in file. 
                             Available fields in {self.pfile.split('/')[-1]} are:
                             {', '.join(self.fields.keys())} and grid_level""")

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

    def compute_global_grids(self) -> list[np.ndarray[float]]:
        """
        Compute the plotfile grids for each level and dimension
        grids are defined such as grid[box_indices] = box_coordinate
        """
        grids = []
        for lv in range(self.limit_level + 1): 
            lvgrids = []
            for coord in range(self.ndims):
                lv_dx = self.dx[lv][coord]/2
                grid = np.linspace(self.geo_low[coord] + lv_dx,
                                   self.geo_high[coord] - lv_dx,
                                   self.grid_sizes[lv][coord])
                lvgrids.append(grid)
            grids.append(lvgrids)
        return grids

    def box_points(self, lv: int, bid: int) -> tuple[np.ndarray[float],
                                                np.ndarray[float],
                                                np.ndarray[float]]:
        """
        Return 3 x (shape) array of the box coordinate points in an AMR box
        lv: Level of the box
        bid: index of the box
        """
        indices = self.cells[lv]['indexes'][bid]
        shape = indices[1] - indices[0] + 1
        x_loc = self.grids[lv][0][indices[0][0]:indices[1][0]+1]
        y_loc = self.grids[lv][1][indices[0][1]:indices[1][1]+1]
        z_loc = self.grids[lv][2][indices[0][2]:indices[1][2]+1]

        x_box = np.repeat([np.repeat([x_loc], shape[1], axis=0)], shape[2], axis=0).T
        y_box = np.repeat([np.repeat([y_loc], shape[2], axis=0).T], shape[0], axis=0)
        z_box = np.repeat([np.repeat([z_loc], shape[1], axis=0)], shape[0], axis=0)
        return x_box, y_box, z_box

    """
    Iterators to loop over plotfile data manually
    """

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

    def map_bfile_offsets(self, lv: int) -> list[np.ndarray[int]]:
        """
        Compute the index map of the AMR box offsets
        for each binary file
        """
        bfiles = np.array(self.cells[lv]["files"])
        offsets = np.array(self.cells[lv]["offsets"])
        offsets_map = []
        for bf in np.unique(bfiles):
            mask = bf == bfiles
            box_indices = np.flatnonzero(mask)
            offsets_map.append(box_indices)
        return offsets_map

    def by_binfile_output(self, other, lv, pltout, **kwargs):
        """
        Iterate over the binary files in two PlotfileCooker
        instances with the assumption that the AMR boxes are
        in the same order in the binary data
        """
        for bf1 in np.unique(self.cells[lv]['files']):
            # The other binary file we read
            mask = np.array(self.cells[lv]['files']) == bf1
            other_idx = np.flatnonzero(mask)[0]
            bf2 = other.cells[lv]['files'][other_idx]
            # Path to the combined binary files (for Windows)
            bfile_r1 = os.path.join(os.getcwd(), bf1)
            bfile_r2 = os.path.join(os.getcwd(), bf2)
            # Path to the new binary file
            bfile_w = os.path.join(os.getcwd(),
                                   pltout,
                                   os.path.basename(os.path.split(bfile_r1)[0]),
                                   os.path.basename(bfile_r1))
            mp_call = {"bfile_r1":bfile_r1,
                       "bfile_r2":bfile_r2,
                       "bfile_w":bfile_w}
            for ky in kwargs:
                mp_call[ky] = kwargs[ky]
            yield mp_call

    def by_matched_offsets_output(self, other, lv, pltout, **kwargs):
        """
        Iterate over two PlotfileCooker instances and match
        box data offsets in the second plotfile to the first
        plotfile so that they correspond to the same global
        indices with an added output binary file
        """
        # Map of which boxes are in which binary files
        box_index_map = self.map_bfile_offsets(lv)
        # On process per binary file
        for bf1, box_indices in zip(np.unique(self.cells[lv]['files']),
                                    box_index_map):
            # Other binary files
            bfiles_2 = np.array(other.cells[lv]['files'])[box_indices]
            # Offsets of the boxes in the binaries
            offsets_bf1 = np.array(self.cells[lv]['offsets'])[box_indices]
            offsets_bf2 = np.array(other.cells[lv]['offsets'])[box_indices]
            # Path to the combined binary files (for Windows)
            bfile_r1 = os.path.join(os.getcwd(), bf1)
            bfile_r2 = [os.path.join(os.getcwd(), bf2) for bf2 in bfiles_2]
            # Path to the new binary file
            bfile_w = os.path.join(os.getcwd(),
                                   pltout,
                                   os.path.basename(os.path.split(bfile_r1)[0]),
                                   os.path.basename(bfile_r1))
            mp_call = {"bfile_r1":bfile_r1,
                       "bfile_r2":bfile_r2,
                       "offst_r2":offsets_bf2,
                       "bfile_w":bfile_w}
            for ky in kwargs:
                mp_call[ky] = kwargs[ky]
            yield mp_call

    def by_matched_boxes_output(self, other, lv, pltout, **kwargs):
        """
        Iterate over two PlotfileCooker instances and match
        the boxes in multiple files in the other plotfile to
        boxes in a single file in the current plotfile
        """
        # Map of which boxes are in which binary files
        box_index_map = self.map_bfile_offsets(lv)
        # On process per binary file
        for bf1, box_indices in zip(np.unique(self.cells[lv]['files']),
                                    box_index_map):
            # The other binary file we read
            bf2 = other.cells[lv]['files'][box_indices[0]]
            # Offsets of the boxes in the binaries
            offsets_bf1 = np.array(self.cells[lv]['offsets'])[box_indices]
            offsets_bf2 = np.array(other.cells[lv]['offsets'])[box_indices]
            # Path to the combined binary files (for Windows)
            #bfile_r1 = os.path.join(os.getcwd(), bf1)
            #bfile_r2 = os.path.join(os.getcwd(), bf2)
            bfile_r1 = bf1
            bfile_r2 = bf2
            # Path to the new binary file
            bfile_w = os.path.join(pltout,
                                   f"Level_{lv}",
                                   os.path.split(bfile_r1)[-1])
            mp_call = {"bfile_r1":bfile_r1,
                       "offst_r1":offsets_bf1,
                       "bfile_r2":bfile_r2,
                       "offst_r2":offsets_bf2,
                       "bfile_w":bfile_w}
            for ky in kwargs:
                mp_call[ky] = kwargs[ky]
            yield mp_call

    """
    Methods resolving the box adjacency in the plotfile
    """

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


    """
    Methods to write new plotfiles using existing structure
    """

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

    def write_global_header_new_fields(self, 
                                       plt_path: str, 
                                       field_names: list[str]) -> None:
        """
        Rewrite a plotfile global header with different fields
        ___
        plt_path: path of the plotfile directory
        pck_ref: reference PlotfileCooker instance to retrieve the
                 plotfile information
        field_names: names of the fields to include in the new plotfile
                     header
        """
        hfile_path = os.path.join(plt_path, "Header")
        nfields = len(field_names)
        # Check for duplicates
        if len(field_names) != len(np.unique(field_names)):
            raise ValueError(("Cannot write plotfile header with duplicate"
                              " fields"))
        with open(hfile_path, 'w') as hfile:
            # Plotfile version
            hfile.write(self.version)
            # Number of fields
            hfile.write(f"{nfields}" + '\n')
            # Fields
            for f in field_names:
                hfile.write(f + '\n')
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
            factors = self.factors[:self.limit_level + 1]
            hfile.write(' '.join([str(f) for f in factors]) + '\n')
            # Grid sizes
            # Looks like ((0,0,0) (7,7,7) (0,0,0))
            tuples = []
            for lv in range(self.limit_level + 1):
                sizes = ",".join([str(s - 1) for s in self.grid_sizes[lv]])
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

    def boxesfromindices(self, indexes):
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

