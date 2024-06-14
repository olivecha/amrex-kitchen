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

# The prototype of this should also live somewhere as 
# it is faster than how binary files are read in colander.py
def mp_read_binfile_minmax(args):
    """
    Multiprocessing function to read binary file data
    with minimal input by using the ascii encoded
    box headers
    """
    bfilename = args[0]
    indexes = []
    min_vals = []
    max_vals = []
    with open(args, 'rb') as bfile:
        while True:
            try:
                h = bfile.readline()
                idx, shape = indexes_and_shape_from_header(h)
                arr = np.fromfile(bfile, 'float64', np.prod(shape))
                arr = arr.reshape(shape, order='F')
                indexes.append(idx)
                min_vals.append(np.min(arr, axis=len(shape)-1))
                max_vals.append(np.max(arr, axis=len(shape)-1))
            except Exception as e:
                break
    return indexes, min_vals, max_vals

class Taster(PlotfileCooker):
    """
    A class to test the validity of AMReX plotfiles
    """
    def __init__(self, plt_file, limit_level=None, boxes=True, 
                 maxmin=True, coordinates=True, nan=False, 
                 nofail=False):
        """
        Constructor for the plotfile tester
        """
        self.boxes_bounds = boxes
        self.boxes_maxmin = maxmin
        self.coordinates = coordinates
        self.check_nan = nan
        if nofail:
            self.fail_on_bad = False
        else:
            self.fail_on_bad = True
        
        # Assume the plotfile is good
        self.isgood = True
        # The attribute value will change in taste
        # if its bad
        try:
            # Instantiate the parent class (PlotfileCooker)
            super().__init__(plt_file, 
                             limit_level=limit_level, 
                             validate_mode=True,
                             maxmins=True)
            self.taste()
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
        # Checking if all box coordinates match
        # The indexes in the level headers
        self.taste_plotfile_structure()
        if self.coordinates:
            self.taste_box_coordinates()
        if self.boxes_maxmin:
            self.taste_maxmins_in_binaries()
        if self.check_nan:
            self.taste_for_nans_in_binairies()
        if self.boxes_bounds:
            self.taste_for_box_shapes_in_binairies()
        print("\nDone!")

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
                    raise TastesBadError(f"Missing file {bfile} at Level {lv}")

    def taste_box_coordinates(self):
        """
        Check that the box coordinates match 
        their indexes
        """
        print(("\nValidating the box coordinates match"
               " the box indexes in the whole plotfile"
               " grid..."))

        # for each level
        for lv in range(self.limit_level + 1):
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
              
    def taste_maxmins_in_binaries(self):
        """
        Method to test the max mins in the level
        headers match thoses of the binaries
        """
        print(("\nValidating that the boxes' maxs and mins match"
               " the level header's maxs and mins..."))
        # for each level
        for lv in range(self.limit_level + 1):
            print(f"Level {lv} ...")
            # mins and maxs in the header data 
            mins_header, maxs_header = [], []
            all_maxs_header = self.cells[lv]['maxs']
            all_mins_header = self.cells[lv]['mins']
            for ky in all_mins_header:
                # Ici tu calcule les valeurs min/maxs entre toutes les boxes
                # Pour le niveau courant
                mins_header.append(np.min(all_mins_header[ky]))
                maxs_header.append(np.max(all_maxs_header[ky]))
            
            min_header, max_header = np.min(mins_header), np.max(maxs_header)
            all_mins_cells, all_maxs_cells = [], []
            
            # For each file 
            for cells in self.bybinfile(lv):
                mins_cell, maxs_cell = [], []
                _, all_mins_box, all_maxs_box = mp_read_binfile_minmax(cells[0])

                for box in all_mins_box:
                    mins_cell.append(np.min(box))
                if len(mins_cell) == 0:
                    error = (f"The array is empty")
                    self.raise_error(TastesBadError,error)
                all_mins_cells.append(np.min(mins_cell))

                for box in all_maxs_box:
                    maxs_cell.append(np.max(box))
                if len(maxs_cell) == 0:
                    error = (f"The array is empty")
                    self.raise_error(TastesBadError,error)
                all_maxs_cells.append(np.max(maxs_cell))

            for m in all_mins_cells:
                if m < min_header:
                    error = (f"Level {lv}'s min : {min_header}"
                                 f" is not absolute because {m} < {min_header}"
                                 f" at cell {cells[0]}")
                    self.raise_error(TastesBadError,error)
                        
            for m in all_maxs_cells:
                if m > max_header:
                    error = (f"Level {lv}'s max : {max_header}"
                                 f" is not absolute because {m} > {max_header}"
                                 f" at cell {cells[0]}")
                    self.raise_error(TastesBadError,error)


    def taste_for_nans_in_binairies(self):
        """
        Method to test if there are NaNs in 
        the binary files
        """
        print(("\nValidating the presence of NaNs in the boxes..."))
        # for each level
        for lv in range(self.limit_level + 1):
            print(f"Level {lv} ...")
            # For each file  
            for cells in self.bybinfile(lv):
                # Let's check if there are any NaN in the boxes
                with open(cells[0],"rb") as bfile:
                    # METTRE While True to read all boxes !!
                    # Ici, on en lit juste une
                    h = bfile.readline()
                    tshape = shapes_from_header_vardims(h,self.ndims)
                    arr = np.fromfile(bfile, "float64", np.prod(tshape))
                    arr = arr.reshape(tshape,order="F")
                    
                    # For every field
                    for i in range(tshape[-1]):
                        # HOW TO NOT HARDCODE THIS ?
                        if self.ndims == 2:
                            array = arr[:, :, i]
                        elif self.ndims == 3:
                            array = arr[:, :, :, i]
                        # Validation
                        if np.max(array) != np.nanmax(array):
                            error = (f"There are NaN in level {lv} at cell {cells[0]}")
                            self.raise_error(TastesBadError,error)


    def taste_for_box_shapes_in_binairies(self):
        """
        Method validating that the data in the binary
        files has the shape (nx, ny, nz, n_fields)
        for every box
        """
        print(("\nValidating the shape and number of boxes..."))
        # for each level
        for lv in range(self.limit_level + 1):
            print(f"Level {lv} ...")
            # number of cells and fields according to the level header 
            nbr_cells_header = len(self.cells[lv]['files'])
            # Same nbr of fields at each level 
            nbr_fields_header = self.nfields

            nbr_box = 0 
            # For each file  
            for cells in self.bybinfile(lv):
                # Let's add the number of box in each binary file for the level 
                nbr_box += len(mp_read_binfile_minmax(cells[0])[0])

                # Let's read the binary files
                with open(cells[0], 'rb') as bfile:
                    h = bfile.readline()
                    # Index and shape of indiviual cells
                    idx, shape = indexes_and_shape_from_header(h)
                
                # Shape according to the level header 
                shape_header = tuple(np.append(((idx[1]-idx[0])+1),nbr_fields_header))

                # Shape Validation
                if shape != shape_header:
                    error = (f"The shapes are not the same :"
                              f"{shape} != {shape_header} at cell {cells[0]}")
                    self.raise_error(TastesBadError,error)
                
            # Number of Cell Validation
            if nbr_cells_header != nbr_box:
                error = (f"The number of cells is not the same :"
                              f"{nbr_cells_header} != {nbr_box} at level {lv}")
                self.raise_error(TastesBadError,error)
                

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

