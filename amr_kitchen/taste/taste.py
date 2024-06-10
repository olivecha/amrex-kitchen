import os
import sys
import time
import multiprocessing
import pickle
import numpy as np
from tqdm import tqdm
from amr_kitchen import PlotfileCooker
from amr_kitchen.utils import TastesBadError

# This should live somewhere
def indexes_and_shape_from_header(header):
    """
    This takes the byte string of a box header in a plotfile binary
    file and infers the indexes of the box and the number of fields
    in the plotfile
    """
    h = header.decode('ascii')
    start, stop, _, nfields = h.split()[-4:]
    nfields = int(nfields)
    start = np.array(start.split('(')[-1].replace(')','').split(','), dtype=int)
    stop = np.array(stop.replace('(', '').replace(')', '').split(','), dtype=int)
    shape = stop - start + 1
    shape = [s for s in shape]
    shape.append(nfields)
    return [start, stop], tuple(shape)

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

def shapes_from_header(header,ndim):
    h = header.decode("ascii")
    start, stop, _, nfields = h.split()[-4:]
    nfields = int(nfields)
    start = np.array(start.split('(')[-1].replace(')','').split(','), dtype=int)
    stop = np.array(stop.replace('(', '').replace(')','').split(','), dtype=int)
    shape = stop - start + 1
    total_shape = []
    for i in range(ndim):
        total_shape.append(shape[i])
    total_shape.append(nfields)
    return total_shape


class Taster(PlotfileCooker):
    """
    A class to test the validity of AMReX plotfiles
    """
    def __init__(self, plt_file, limit_level=None, boxes=True, 
                 maxmin=True, coordinates=True, nan=True, 
                 nofail=False):
        """
        Constructor for the plotfile tester
        """
        # Instantiate the parent class (PlotfileCooker)
        super().__init__(plt_file, limit_level=limit_level)
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
        self.taste()

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
        if self.coordinates:
            self.taste_box_coordinates()
        if self.boxes_maxmin:
            self.taste_maxmins_in_binaries()
        if self.check_nan:
            self.taste_for_nans_in_binairies()
        if self.boxes_bounds:
            self.taste_for_box_shapes_in_binairies()
        print("\nDone!")

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
            _, all_mins_header, all_maxs_header = self.read_cell_headers()
            for key in all_mins_header[lv]:
                # Ici tu calcule les valeurs min/maxs entre toutes les boxes
                # Pour le niveau courant
                mins_header.append(np.min(all_mins_header[lv][key]))
                maxs_header.append(np.max(all_maxs_header[lv][key]))
            
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
                    tshape = shapes_from_header(h,self.ndims)
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

#bigfile = os.path.join(os.getcwd(),"test_assets","example_plt_3d")
#file = os.path.join(os.getcwd(),"test_assets","example_plt_3d")
#tst = Taster(file)





#######################################################################################

#cell = os.path.join(os.getcwd(),"test_assets","example_plt_2d","Level_0","Cell_D_00000")
#pck = PlotfileCooker(os.path.join("C:\\","Users","franc","Desktop","Stage-POLY-2024","Data","pltlv4temp_Francois"))
"""def tasting(plt_file,
            boxes=True,
            maxmin=True,
            coordinates=True,
            nan=True,
            nofail=False,):
    # For testing 
    if boxes:
        print("boxes yes")
    else:
        print("boxes no")
    if maxmin:
        print("maxmin yes")
    else:
        print("maxmin no")
    if coordinates:
        print("coordinates yes")
    else:
        print("coordinates no")
    if nan:
        print("nan yes")
    else:
        print("nan no")
    if nofail:
        print("nofail yes")
    else:
        print("nofail no")

    # Create the PlotfileCooker class instance
    pck = PlotfileCooker(plt_file)

    # Check that the number of box coordinates matches the

    # **** My bad but we should always use the name "box" instead of cell ****
    # Faudrait combiner les attributs boxs et cells mais c'est un bon petit 
    # projet de refractoring genre pck.boxes devient pck.boxes['geo']
    # et pck.cells['files'] devient pck.boxes['files']

    # number of files in the cell header data
    for lv in range(pck.limit_level + 1):
        print(f"Tasting level {lv} ...")
        # Checking if all coordonates match
        if coordinates:
            for dim in range(pck.ndims):
                dim_grid = np.linspace(pck.geo_low[dim] + pck.dx[lv][dim]/2,
                        pck.geo_high[dim] - pck.dx[lv][dim]/2,
                        pck.grid_sizes[lv][dim])
                for i in range(len(pck.boxes[lv])):
                    # Indexs to validate (each box has an index)
                    idx = pck.cells[lv]["indexes"][i]
                    # Conrespondant box
                    box = pck.boxes[lv][i]
                    # Trouver les coordonnées de la box avec les indexes
                    box_x_lo = dim_grid[idx[0][0]] - pck.dx[lv][dim]/2 # Pourquoi ça ne marcherait pas avec idx[0][dim]
                    box_x_hi = dim_grid[idx[1][0]] + pck.dx[lv][dim]/2 #                                    idx[1][dim] ??
                    # Validation
                    if not nofail:
                        try:
                            raise TastesBadError((box[0][0],box_x_lo))
                        except TastesBadError as e:
                            if not np.isclose(e.value[0], e.value[1]):
                                print(f"Low coordonates are not the same : {e.value[0]} != {e.value[1]} at ...")
                                raise TastesBadError("Low Coordonate Error")
                        try:
                            raise TastesBadError((box[0][1],box_x_hi))
                        except TastesBadError as e:
                            if not np.isclose(e.value[0], e.value[1]):
                                print(f"High coordonates are not the same : {e.value[0]} != {e.value[1]} at ...")
                                raise TastesBadError("High Coordonate Error")
                    else:
                        if not np.isclose(box[0][0], box_x_lo):
                            print(f"Low coordonates are not the same : {box[0][0]} != {box_x_lo} at ...")
                        if not np.isclose(box[0][1], box_x_hi):
                            print(f"High coordonates are not the same : {box[0][1]} != {box_x_hi} at ...")
                

        # number of cells and fields according to the level header 
        nbr_cells_header = len(pck.cells[lv]['files'])
        # Same nbr of fields at each level 
        nbr_fields_header = pck.nfields

        if maxmin:
            # mins and maxs in the header data 
            mins_header, maxs_header = [], []
            _, all_mins_header, all_maxs_header = pck.read_cell_headers()
            for key in all_mins_header[lv]:
                # Ici tu calcule les valeurs min/maxs entre toutes les boxes
                # Pour le niveau courant
                mins_header.append(np.min(all_mins_header[lv][key]))
                maxs_header.append(np.max(all_maxs_header[lv][key]))
            # Ici tu calcule les valeurs max/mins parmis tout les fields
            # Je pense pas c'est ce que tu voulais faire...
            # Par contre c'est intéressant à avoir pour regarder s'il y a des NotANumber
            # Tu pourrais valider np.nanmax(mins_header) == np.max(mins_header)...
            # Tu pourrais profiler mais je pense c'est plus rapide que np.isnan(data).any()
            # Après je sais pas si le module qui sauvegarde les plotfiles filtre les nans 
            # Pour les max mins des boxes je vais regarder
            min_header, max_header = np.min(mins_header), np.max(maxs_header)
            all_mins_cells, all_maxs_cells = [], []

        nbr_box = 0 
        for cells in pck.bybinfile(lv):
            # Let's check if there are any NaN in the boxes
            if nan:
                with open(cells[0],"rb") as bfile:
                    # METTRE While True to read all boxes !!
                    # Ici, on en lit juste une
                    h = bfile.readline()
                    tshape = shapes_from_header(h,pck.ndims)
                    arr = np.fromfile(bfile, "float64", np.prod(tshape))
                    arr = arr.reshape(tshape,order="F")
                    
                    # For every field
                    for i in range(tshape[-1]):
                        # HOW TO NOT HARDCODE THIS ?
                        if pck.ndims == 2:
                            array = arr[:, :, i]
                        elif pck.ndims == 3:
                            array = arr[:, :, :, i]
                        # Validation
                        if not nofail:
                            try:
                                raise TastesBadError((np.max(array),np.nanmax(array)))
                            except TastesBadError as e:
                                if e.value[0] != e.value[1]:
                                    print(f"There are NaN in level {lv} at cell {cells[0]}")
                                    raise TastesBadError("NaN Error") 
                        else:
                            if np.max(array) != np.nanmax(array):
                                print(f"There are NaN in level {lv} at cell {cells[0]}")
                        
            # Let's add the number of box in each binary file for the level 
            nbr_box += len(mp_read_binfile_minmax(cells[0])[0])

            # Let's collect the mins and maxs of every binary file for the level 
            mins_cell, maxs_cell = [], []
            idx, all_mins_box, all_maxs_box = mp_read_binfile_minmax(cells[0])
            # Ici all_mins_box c'est une liste la longueur du nombre de box dans la binary file
            # Et chaque élément c'est les valeurs mins/maxs pour chaque field
            # all_mins_box[100][3] c'est la valeur min dans la box 100 au niveau actuel du field avec 9
            # l'index 3
            if maxmin:
                for box in all_mins_box:
                    mins_cell.append(np.min(box))
                all_mins_cells.append(np.min(mins_cell))

                for box in all_maxs_box:
                    maxs_cell.append(np.max(box))
                all_maxs_cells.append(np.max(maxs_cell))

            if boxes:
                # Let's read the binary files
                with open(cells[0], 'rb') as bfile:
                    h = bfile.readline()
                    # Index and shape of indiviual cells
                    idx, shape = indexes_and_shape_from_header(h)
                
                # Shape according to the level header 
                shape_header = tuple(np.append(((idx[1]-idx[0])+1),nbr_fields_header))

                if not nofail:
                    # Comparison between shapes
                    try:
                        raise TastesBadError((shape,shape_header))
                    except TastesBadError as e:
                        if e.value[0] != e.value[1]:
                            print(f"The shapes are not the same : {e.value[0]} != {e.value[1]} at cell {cells[0]}")
                            raise TastesBadError("Shape Error")
                else:
                    if shape != shape_header:
                        print(f"The shapes are not the same : {shape} != {shape_header} at cell {cells[0]}")

        if boxes:
            if not nofail:
                # Comparison between number of cells 
                try:
                    raise TastesBadError((nbr_cells_header,nbr_box))
                except TastesBadError as e:
                    if e.value[0] != e.value[1]:
                        print(f"The number of cells is not the same : {e.value[0]} != {e.value[1]} at level {lv}")
                        raise TastesBadError("Box Number Error") 
            else:
                if nbr_cells_header != nbr_box:
                    print(f"The number of cells is not the same : {nbr_cells_header} != {nbr_box} at level {lv}")
        
        if maxmin:
            for m in all_mins_cells:
                if not nofail:
                    try:
                        raise TastesBadError((m,min_header))
                    except TastesBadError as e:
                        if e.value[0] < e.value[1]:
                            print(f"Level {lv}'s min : {e.value[1]} --> not absolute because {e.value[0]} < {e.value[1]} at cell {cells[0]}")
                            raise TastesBadError("Minimum Error") 
                else:
                    if m < min_header:
                        print(f"Level {lv}'s min : {min_header} --> not absolute because {m} < {min_header} at cell {cells[0]}")
                    
            for m in all_maxs_cells:
                if not nofail:
                    try:
                        raise TastesBadError((m,max_header))
                    except TastesBadError as e:
                        if e.value[0] > e.value[1]:
                            print(f"Level {lv}'s max : {e.value[1]} --> not absolute because {e.value[0]} > {e.value[1]} at cell {cells[0]}")
                            raise TastesBadError("Maximum Error") 
                else:
                    if m > max_header:
                        print(f"Level {lv}'s max : {max_header} --> not absolute because {m} > {max_header} at cell {cells[0]}")
        print("Done!")
"""