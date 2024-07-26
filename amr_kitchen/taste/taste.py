import os
import sys
import time
import traceback
import multiprocessing
import pickle
import numpy as np
from tqdm import tqdm
from amr_kitchen import PlotfileCooker
from amr_kitchen.utils import TastesBadError
from amr_kitchen.utils import indexes_and_shape_from_header
from amr_kitchen.utils import shapes_from_header_vardims
from amr_kitchen.utils import shape_from_header
from amr_kitchen.utils import header_from_indices

# The prototype of this should also live somewhere as 
# it is faster than how binary files are read in colander.py
def mp_read_binary_data(args):
    """
    Multiprocessing function to read binary file data
    with minimal input by using the ascii encoded
    box headers
    """
    bfilename = args[0]
    bfile_data = []
    with open(args, 'rb') as bfile:
        while True:
            try:
                h = bfile.readline().decode('ascii')
                shape = shape_from_header(h)
                arr = np.fromfile(bfile, 'float64', np.prod(shape))
                arr = arr.reshape(shape, order='F')
                bfile_data.append(arr)
            except Exception as e:
                break
    return bfile_data

class Taster(PlotfileCooker):
    """
    A class to test the validity of AMReX plotfiles
    """
    def __init__(self, plt_file,
                 limit_level=None,
                 binary_headers=True,
                 binary_shape=True,
                 binary_data=False,
                 boxes_coordinates=False,
                 nofail=False,
                 verbose=None):
        """
        Constructor for the plotfile tester
        """
        # Define what will be validated
        self.check_binary_headers = binary_headers
        self.check_binary_shape = binary_shape
        self.check_binary_data = binary_data
        self.check_boxes_coordinates = boxes_coordinates

        if nofail:
            self.fail_on_bad = False
        else:
            self.fail_on_bad = True

        if verbose is None:
            self.v = 1
        else:
            self.v = verbose

        # Assume the plotfile is good
        self.isgood = True
        # The attribute value will change in taste
        # if its bad
        try:
            # Instantiate the parent class (PlotfileCooker)
            if self.check_binary_data:
                super().__init__(plt_file,
                                 limit_level=limit_level,
                                 validate_mode=True,
                                 maxmins=True)
            else:
                super().__init__(plt_file,
                                 limit_level=limit_level,
                                 validate_mode=True,
                                 maxmins=False)
            self.taste()
        # This catches errors in PlotfileCooker
        except Exception as e:
            self.isgood = False
            if self.fail_on_bad:
                raise e
            else:
                tb = traceback.format_exc()
                print(tb)

    def __bool__(self):
        """
        Overiding the bool method of the class
        to return False if the plotfile is Bad
        """
        return self.isgood

    def taste(self):
        """
        Main functions validating the sanity of the plotfile
        ___
        Depending on the input arguments more or less 
        attributes of the plotfile are tested.
        |1 - The number and shape of boxes |
        |2 - The mins and maxs             |
        |3 - The boxes' coordinates        |
        |4 - The NaNs                      |
        """
        # First check that no binary files are missing
        self.taste_plotfile_structure()
        # If flagged check that the boxes bounds match
        # the box indexes in the global grid at each
        # level
        if self.check_boxes_coordinates:
            self.taste_box_coordinates()
        # Validate that the headers in the binary files
        # Match those in the level headers
        if self.check_binary_headers:
            self.taste_binary_headers()
        # Validate that the binary data has the right
        # shape. We can do this without actually reading
        # it by skipping the expected number of bytes in
        # the binary file and validating that we are at
        # the next header
        if self.check_binary_shape:
            if not self.check_binary_headers:
                # TODO: print an actual warning and if
                # the class was instantied from the cli
                # tell which flags should be true together
                if self.v > 0:
                    print(("WARNING the validation of the"
                           " binary data shape assumes that"
                           " the binary headers are valid"))
            self.taste_binary_shape()

        if self.check_binary_data:
            if (not self.check_binary_headers or
                not self.check_binary_shape):
                # TODO: print an actual warning and if
                # the class was instantied from the cli
                # tell which flags should be true together
                if self.v > 0:
                    print(("WARNING the validation of the"
                           " binary data assumes that the"
                           " binary headers and the data shape"
                           " are valid"))
                self.taste_binary_data()

    def taste_plotfile_structure(self):
        """
        Check that no binary file is missing
        """
        for lv in range(self.limit_level + 1):
            lv_files = os.listdir(os.path.join(self.pfile,
                                               self.cell_paths[lv]))
            for bfile_path in np.unique(self.cells[lv]['files']):
                bfile = os.path.split(bfile_path)[-1]
                if bfile not in lv_files:
                    error = f"Missing file {bfile} at Level {lv}"
                    self.raise_error(TastesBadError, error)

    def taste_box_coordinates(self):
        """
        Check that the box coordinates match
        their indexes
        """
        if self.v > 0:
            print(("\nValidating the box coordinates match"
                   " the box indexes in the whole plotfile"
                   " grid..."))

        # for each level
        for lv in range(self.limit_level + 1):
            if self.v > 0:
                print(f"Level {lv} ...")
            # For each dimension
            # Define the global coordinate grid
            grids = []
            for dim in range(self.ndims):
                dim_grid = np.linspace(self.geo_low[dim] + self.dx[lv][dim]/2,
                        self.geo_high[dim] - self.dx[lv][dim]/2,
                        self.grid_sizes[lv][dim])
                grids.append(dim_grid)
            # For each box in the domain
            for i, box in enumerate(self.boxes[lv]):
                # Get the corresponding indexes
                idx = self.cells[lv]["indexes"][i]
                for dim in range(self.ndims):
                    # Trouver les coordonnées de la box avec les indexes
                    box_lo = grids[dim][idx[0][dim]] - self.dx[lv][dim]/2
                    box_hi = grids[dim][idx[1][dim]] + self.dx[lv][dim]/2
                    # not matching lower bounds
                    if ~np.isclose(box_lo, box[dim][0]):
                        # this could be moved to a method
                        error = (f"The lower bound of box {i} at level {lv}"
                                 f" ({box[dim][0]} is not equal to the value"
                                 f" found using the box index {idx[0][dim]} in"
                                 f" the coordinate grid between {grids[dim][0]}"
                                 f" and {grids[dim][-1]} with {len(grids[dim])}"
                                 f" points equal to {box_lo}")
                        self.raise_error(TastesBadError, error)
                    # not matching upper bounds
                    if ~np.isclose(box_hi, box[dim][1]):
                        error = (f"The upper bound of box {i} at level {lv}"
                                 f" ({box[dim][1]} is not equal to the value"
                                 f" found using the box index {idx[1][dim]} in"
                                 f" the coordinate grid between {grids[dim][0]}"
                                 f" and {grids[dim][-1]} with {len(grids[dim])}"
                                 f" points equal to {box_hi}")
                        self.raise_error(TastesBadError, error)
            if self.v > 0:
                print("Done!")

    def taste_binary_headers(self):
        """
        Method to validate that the indices and shape
        of the binary headers match those in the
        level headers
        """
        for lv in range(self.limit_level + 1):
            lv_box_ids = np.arange(len(self.cells[lv]['indexes']))
            for bfile in np.unique(self.cells[lv]['files']):
                # Boxes in the current file
                bfile_mask = np.array(self.cells[lv]['files']) == bfile
                offsets = np.array(self.cells[lv]['offsets'])[bfile_mask]
                indices = np.array(self.cells[lv]['indexes'])[bfile_mask]
                box_ids = lv_box_ids[bfile_mask]
                # Iterate over the sorted offsets so we don't jump 
                # around the file 
                # TODO: check if this is actually faster
                indices = indices[np.argsort(offsets)]
                box_ids = box_ids[np.argsort(offsets)]
                offsets = np.sort(offsets)
                with open(bfile, 'rb') as bf:
                    for ofs, idx, bid in zip(offsets,
                                             indices,
                                             box_ids):
                        # This is pretty fast because we barely
                        # read any data (just one line per header)
                        # so we can test this first before iterating
                        # over the binary data
                        bf.seek(ofs)
                        h = bf.readline()
                        hidx, shape = indexes_and_shape_from_header(h)
                        # Check that the box indices match
                        # If they don't match the shape won't match
                        # so we don't have to check it 
                        if not np.array_equal(idx, hidx):
                            error = (f"The box number {bid} at level {lv}"
                                     f" has the indices {idx[0]} to {idx[1]}"
                                     f" in the level header, and indices"
                                     f" {hidx[0]} to {hidx[1]} in the binary"
                                     f" header in file {bfile}")
                            self.raise_error(TastesBadError, error)
                        # Check that the binary header has the right number of
                        # fields
                        if shape[-1] != len(self.fields):
                            error = (f"The number of field in the binary header"
                                     f" ({shape[-1]}) does not match the number"
                                     f" of fields in the plotfile header"
                                     f" ({len(self.fields)}) for the box number"
                                     f" {bid} at level {lv} in the binary file"
                                     f" {bfile}")
                            self.raise_error(TastesBadError, error)

    def taste_binary_shape(self):
        """
        Check that the number of data bytes between
        binary headers matches what is expected from the
        box indices and number of fields
        """
        if self.v > 0:
            print("Validating the boxes shape in the binary data")
        for lv in range(self.limit_level + 1):
            if self.v > 0:
                print(f"Level {lv} ...")
            lv_box_ids = np.arange(len(self.cells[lv]['indexes']))
            for bfile in np.unique(self.cells[lv]['files']):
                # Boxes in the current file
                bfile_mask = np.array(self.cells[lv]['files']) == bfile
                offsets = np.array(self.cells[lv]['offsets'])[bfile_mask]
                box_ids = lv_box_ids[bfile_mask]
                # Sort the box ids so they match what we read
                # in the file 
                box_ids = box_ids[np.argsort(offsets)]
                # Read the binary file
                with open(bfile, 'rb') as bf:
                    # Iterate with the index so we can infer what the
                    # next binary header is
                    # Read the first header
                    h = bf.readline()
                    for i, bid in enumerate(box_ids[:-1]):
                        # We can assume the header has the right shape
                        # as it was validated before
                        shape = shape_from_header(h.decode('ascii'))
                        # 8 bytes per float64
                        nbytes = np.prod(shape)*8
                        # Skip to the next header
                        # This is pretty fast because we only
                        # read the header line and skip over the
                        # data bytes but we can still validate that
                        # the expected number of points is here
                        bf.seek(nbytes, 1)
                        # Indices of the next header
                        exp_indices = self.cells[lv]['indexes'][box_ids[i+1]]
                        # Byte string of the next header
                        exp_header = header_from_indices(exp_indices[0],
                                                         exp_indices[1],
                                                         len(self.fields))
                        h = bf.readline()
                        if h != exp_header:
                            error = (f"The binary data for box {bid}"
                                     f" at level {lv} does not contain"
                                     f" the appropriate number of bytes"
                                     f" ({nbytes}) for the shape"
                                     f" {shape} in the binary file"
                                     f" {bfile}")
                            self.raise_error(TastesBadError, error)
                    # Special case for the last box
                    shape = shape_from_header(h.decode('ascii'))
                    nbytes = np.prod(shape)*8
                    # Go to the last byte and store the position
                    bf.seek(nbytes, 1)
                    file_size = bf.tell()
                    # Compare with actual end of file
                    if file_size != bf.seek(0, 2):
                        error = (f"The binary data for box {box_ids[-1]}"
                                 f" at level {lv} does not contain"
                                 f" the appropriate number of bytes"
                                 f" ({nbytes}) for the shape"
                                 f" {shape} in the binary file"
                                 f" {bfile}")
                        self.raise_error(TastesBadError, error)
            if self.v > 0:
                print("Done!")

    def taste_binary_data(self):
        """
        It is super long to read the binary data so
        this tests that both the max/mins match those
        in the level headers and warns if there are NaNs
        in the data
        """
        if self.v > 0:
            print("Validating the binary data")
        # for each level
        for lv in range(self.limit_level + 1):
            if self.v > 0:
                print(f"Level {lv} ...")
            lv_boxes_ids = np.arange(len(self.cells[lv]['offsets']))
            # For every binary file
            bfile_data = {}
            for bfile in np.unique(self.cells[lv]['files']):
                # mask of the boxes in the binary file
                bf_mask = np.array(self.cells[lv]['files']) == bfile
                offsets = np.array(self.cells[lv]['offsets'])[bf_mask]
                ofst_sort = np.argsort(offsets)
                # Sort everything by offset to read the file
                # sequentially
                box_ids = lv_boxes_ids[bf_mask][ofst_sort]
                # Maxs and mins for the current binary file
                # Separated between fields
                maxs = {}
                mins = {}
                for f in self.fields:
                        maxs[f] = self.cells[lv]['maxs'][f][bf_mask][ofst_sort]
                        mins[f] = self.cells[lv]['mins'][f][bf_mask][ofst_sort]
                # Divide the data between binary file to access with the
                # multiprocessing output
                bfile_data['maxs'] = maxs
                bfile_data['mins'] = mins
                bfile_data['bids'] = box_ids
            pool = multiprocessing.Pool()
            # Iterate over every binary file
            for bfile, data_out in zip(bfile_data.keys(),
                                       pool.imap(mp_read_binary_data,
                                                 bfile_data.keys())):
                # Loop over every box as the data is read
                # TODO: maybe it would be faster to compute the np.max/np.nanmax
                # in the multiprocessing function and send back only the max/mins
                # as this is a lot of data piped trough the multiprocessing pool
                for idx, data in enumerate(data_out):
                    for f in self.fields:
                        # Get the data of a single field
                        fdata = data[..., self.fields[f]]
                        # Check for NaNs
                        # The min and max catches -inf and inf
                        if (not np.isclose(np.max(fdata), np.nanmax(fdata)) or
                            not np.isclose(np.min(fdata), np.nanmin(fdata))):
                            # Its okay if theres NaNs in the file
                            if self.warn_nans:
                                message = (f"{f} data for box"
                                           f" {bfile_data[bfile]['bids'][idx]}"
                                           f" at level {lv} in the binary"
                                           f" file {bfile} Contains NaNs")
                                if self.v > 0:
                                    print(message)
                        # Check that the min/max match those in the level headers
                        fmin = bfile_data[bfile]["mins"][f][idx]
                        fmax = bfile_data[bfile]["maxs"][f][idx]
                        # TODO: not sure how NaN are handled in the level headers
                        # maybe it would be better to validate against np.min as
                        # we don't want to raise and error if both values are NaN
                        if not np.isclose(fmin, np.nanmin(fdata)):
                            error = (f"The minimum {f} value in the level"
                                     f" header ({fmin}) is different from"
                                     f" the value found in the binary file"
                                     f" {np.nanmin(fdata)} for box"
                                     f" {bfile_data[bfile]['bids'][idx]}"
                                     f" at level {lv} in the binary file"
                                     f" {bfile}")
                            self.raise_error(TastesBadError, error)

                        if not np.isclose(fmax, np.nanmax(fdata)):
                            error = (f"The maximum {f} value in the level"
                                     f" header ({fmax}) is different from"
                                     f" the value found in the binary file"
                                     f" {np.nanmax(fdata)} for box"
                                     f" {bfile_data[bfile]['bids'][idx]}"
                                     f" at level {lv} in the binary file"
                                     f" {bfile}")
                            self.raise_error(TastesBadError, error)

    def raise_error(self, error, message):
        """
        A method to wrap arround raise statements to
        allow printing errors instead if raising to
        continue validating after one error if the
        option is set
        """
        self.isgood = False
        if self.fail_on_bad:
            raise error(message)
        else:
            error_name = str(error).split("'")[1]
            print(f"Encountered {error_name}:\n", message)

