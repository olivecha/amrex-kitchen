import os
import sys
import time
from pathos.multiprocessing import ProcessingPool as Pool
from importlib.util import spec_from_file_location
from importlib.util import module_from_spec
from inspect import signature
import numpy as np
import cantera as ct
from tqdm import tqdm
from amr_kitchen import PlotfileCooker
from amr_kitchen.utils import shape_from_header

# GLOBALS
SARRAYS = None
PRESSURES = None

def chefs_knife_single_field(args):
    """
    __ Multiprocessing function __
    Cooks a single binary file
    Only compute a single field
    """
    # New offsets
    offsets = []
    mins = []
    maxs = []
    # Open the read and write
    with open(args['bfpath'], 'rb') as bfr:
        with open(args['newbfpath'], 'wb') as bfw:
            while True:
                # Read the header line in the old file
                header = bfr.readline().decode('ascii')
                # Infer the box data shape
                try:
                    datashape = shape_from_header(header)
                # This raises when the line is not a header
                except Exception as e:
                    break
                # Store byte position in new file
                offsets.append(bfw.tell())
                # shape of the box
                boxshape = tuple([datashape[i] for i in range(3)])
                # Compute the new binary header (Always one field)
                header_w = header.replace(f'{datashape[-1]}\n',
                                          '1\n')
                # Write the new header
                bfw.write(header_w.encode('ascii'))
                # Read the data
                arr = np.fromfile(bfr, "float64", np.prod(datashape))
                arr = arr.reshape(datashape, order="F")
                # Isolate the temperature and mass fractions
                Y = arr[:, :, :, args['sp_start']:args['sp_end']]
                T = arr[:, :, :, args['id_temp']]
                P = PRESSURES[boxshape]
                # Clean the mass fractions and temp so the thermo
                # state is not undefined
                T[np.isclose(T, 0)] = 1
                Y[np.isclose(np.sum(Y, axis=3), 0), args['idx_O2']] = 1.0
                sarray = SARRAYS[boxshape]
                sarray.TPY = T, P, Y
                newdata = sarray.__getattribute__(args['recipe'])
                mins.append(np.min(newdata))
                maxs.append(np.max(newdata))
                bfw.write(newdata.flatten(order="F").tobytes())

    return offsets, np.array([mins]).T, np.array([maxs]).T


def chefs_knife_byspecies_field(args):
    """
    __ Multiprocessing function __
    Cooks a single binary file
    Compute a single field for every species selected
    """
    # New offsets
    offsets = []
    mins = []
    maxs = []
    outnfields = len(args['sp_indexes'])
    # Open the read and write
    with open(args['bfpath'], 'rb') as bfr:
        with open(args['newbfpath'], 'wb') as bfw:
            while True:
                # Read the header line in the old file
                header = bfr.readline().decode('ascii')
                # Infer the box data shape
                try:
                    datashape = shape_from_header(header)
                # This raises when the line is not a header
                except Exception as e:
                    break
                # Store byte position in new file
                offsets.append(bfw.tell())
                # shape of the box
                boxshape = tuple([datashape[i] for i in range(3)])
                # Compute the new binary header (Number of species)
                header_w = header.replace(f'{datashape[-1]}\n',
                                          f"{outnfields}\n")
                # Write the new header
                bfw.write(header_w.encode('ascii'))
                # Read the data
                arr = np.fromfile(bfr, "float64", np.prod(datashape))
                arr = arr.reshape(datashape, order="F")
                # Isolate the temperature and mass fractions
                Y = arr[:, :, :, args['sp_start']:args['sp_end']]
                T = arr[:, :, :, args['id_temp']]
                P = PRESSURES[boxshape]
                # Clean the mass fractions and temp so the thermo
                # state is not undefined
                T[np.isclose(T, 0)] = 1
                Y[np.isclose(np.sum(Y, axis=3), 0), args['idx_O2']] = 1.0
                sarray = SARRAYS[boxshape]
                sarray.TPY = T, P, Y
                # Get the by species data in the solution array
                newdata = sarray.__getattribute__(args['recipe'])
                # Index the required species
                newdata = newdata[:, :, :, args['sp_indexes']]
                mins.append(np.min(newdata, axis=(0, 1, 2)))
                maxs.append(np.max(newdata, axis=(0, 1, 2)))
                bfw.write(newdata.flatten(order="F").tobytes())
    return offsets, np.array([mins]), np.array([maxs])

def chefs_knife_byreaction_field(args):
    """
    __ Multiprocessing function __
    Cooks a single binary file
    Compute a single field for every reaction selected
    """
    # New offsets
    offsets = []
    mins = []
    maxs = []
    outnfields = len(args['rx_indexes'])
    # Open the read and write
    with open(args['bfpath'], 'rb') as bfr:
        with open(args['newbfpath'], 'wb') as bfw:
            while True:
                # Read the header line in the old file
                header = bfr.readline().decode('ascii')
                # Infer the box data shape
                try:
                    datashape = shape_from_header(header)
                # This raises when the line is not a header
                except Exception as e:
                    break
                # Store byte position in new file
                offsets.append(bfw.tell())
                # shape of the box
                boxshape = tuple([datashape[i] for i in range(3)])
                # Compute the new binary header (Always one field)
                header_w = header.replace(f'{datashape[-1]}\n', 
                                          f'{outnfields}\n')
                # Write the new header
                bfw.write(header_w.encode('ascii'))
                # Read the data
                arr = np.fromfile(bfr, "float64", np.prod(datashape))
                arr = arr.reshape(datashape, order="F")
                # Isolate the temperature and mass fractions
                Y = arr[:, :, :, args['sp_start']:args['sp_end']]
                T = arr[:, :, :, args['id_temp']]
                P = PRESSURES[boxshape]
                # Clean the mass fractions and temp so the thermo
                # state is not undefined
                T[np.isclose(T, 0)] = 1
                Y[np.isclose(np.sum(Y, axis=3), 0), args['idx_O2']] = 1.0
                sarray = SARRAYS[boxshape]
                sarray.TPY = T, P, Y
                # Get the by species data in the solution array
                newdata = sarray.__getattribute__(args['recipe'])
                # Index the required species
                newdata = newdata[:, :, :, args['rx_indexes']]
                mins.append(np.min(newdata, axis=(0, 1, 2)))
                maxs.append(np.max(newdata, axis=(0, 1, 2)))
                bfw.write(newdata.flatten(order="F").tobytes())
    return offsets, np.array([mins]), np.array([maxs])

def chefs_knife_user_sarray(args):
    """
    __ Multiprocessing function __
    Cooks a single binary file
    Compute a single field for every species selected
    """
    # New offsets
    offsets = []
    mins = []
    maxs = []
    # Open the read and write
    with open(args['bfpath'], 'rb') as bfr:
        with open(args['newbfpath'], 'wb') as bfw:
            while True:
                # Read the header line in the old file
                header = bfr.readline().decode('ascii')
                # Infer the box data shape
                try:
                    datashape = shape_from_header(header)
                # This raises when the line is not a header
                except Exception as e:
                    break
                # Store byte position in new file
                offsets.append(bfw.tell())
                # shape of the box
                boxshape = tuple([datashape[i] for i in range(3)])
                # Compute the new binary header (Always one field)
                header_w = header.replace(f'{datashape[-1]}\n', 
                                          '1\n')
                # Write the new header
                bfw.write(header_w.encode('ascii'))
                # Read the data
                arr = np.fromfile(bfr, "float64", np.prod(datashape))
                arr = arr.reshape(datashape, order="F")
                # Isolate the temperature and mass fractions
                Y = arr[:, :, :, args['sp_start']:args['sp_end']]
                T = arr[:, :, :, args['id_temp']]
                P = PRESSURES[boxshape]
                # Clean the mass fractions and temp so the thermo
                # state is not undefined
                T[np.isclose(T, 0)] = 1
                Y[np.isclose(np.sum(Y, axis=3), 0), args['idx_O2']] = 1.0
                sarray = SARRAYS[boxshape]
                sarray.TPY = T, P, Y
                # Call the user defined function with the field indexes
                # box array and solution array as arguments
                newdata = args['recipe'](args['field_indexes'],
                                         arr,
                                         sarray)
                # Index the required species
                mins.append(np.min(newdata))
                maxs.append(np.max(newdata))
                bfw.write(newdata.flatten(order="F").tobytes())

    return offsets, np.array([mins]).T, np.array([maxs]).T

def chefs_knife_user_pfile(args):
    """
    __ Multiprocessing function __
    Cooks a single binary file
    Compute a single field for every species selected
    Only the plotfile data is available to the user
    recipe
    """
    # New offsets
    offsets = []
    mins = []
    maxs = []
    # Open the read and write
    with open(args['bfpath'], 'rb') as bfr:
        with open(args['newbfpath'], 'wb') as bfw:
            while True:
                # Read the header line in the old file
                header = bfr.readline().decode('ascii')
                # Infer the box data shape
                try:
                    datashape = shape_from_header(header)
                # This raises when the line is not a header
                except Exception as e:
                    break
                # Store byte position in new file
                offsets.append(bfw.tell())
                # shape of the box
                boxshape = tuple([datashape[i] for i in range(3)])
                # Compute the new binary header (Always one field)
                header_w = header.replace(f'{datashape[-1]}\n', 
                                          '1\n')
                # Write the new header
                bfw.write(header_w.encode('ascii'))
                # Read the data
                arr = np.fromfile(bfr, "float64", np.prod(datashape))
                arr = arr.reshape(datashape, order="F")
                # Call the user defined function with the field indexes
                newdata = args['recipe'](args['field_indexes'],
                                         arr)
                # Compute the max/mins
                mins.append(np.min(newdata))
                maxs.append(np.max(newdata))
                # Write to the new binary file
                bfw.write(newdata.flatten(order="F").tobytes())

    return offsets, np.array([mins]).T, np.array([maxs]).T

class Chef(PlotfileCooker):

    cookbook = {'HRR': "heat_release_rate",
                'ENT': "enthalpy_mass",
                'SRi': "net_production_rates",
                'SDi': "mix_diff_coeffs_mass",
                'RRi': "net_rates_of_progress",
                'user': None}

    cookbookhelp = {'HRR': "Heat release rate [W/m^3]",
                    'ENT': "Mass based specific enthalpy [J/kg]",
                    'SRi': "Species net source term [kmol/m^3/s]",
                    'SDi': ("Mixture-averaged diffusion coefficients [m^2/s]"
                            " relating the diffusive mass fluxes to gradients"
                            " in the species mass fractions."),
                    'RRi': "Net rates of progress of indivividual reactions",
                    'user': ("A user defined python function named 'recipe'"
                             " defined in a python file supplied as the recipe"
                             " argument.\n"
                             " The docstring of the function will be used as the"
                             " field name in the output plotfile with a default"
                             " value set to 'user_defined'.\n"
                             "The function must have one of the following call"
                             " signatures:\n"
                             "def recipe(field_indexes, box_array):\n"
                             "field_indexes: A dict mapping the field"
                             " names to their indexes in box_array\n"
                             "box_array: A numpy array with shape (nx, ny, nz, nfields)"
                             " with nx/ny/nz the shape of a box in the plotfile and nfields"
                             " the number of fields in the plotfile."
                             " box_array[:, :, :, field_indexes['temp']] is the nx*ny*nz"
                             " array containing the temperature data for the current box\n\n"
                             "def recipe(field_indexes, box_array, sol_array):\n"
                             " sol_array: a nx*ny*nz Cantera solution array with the same"
                             " state as box_array, which has all the attributes of the"
                             " Cantera solution\n"
                             "In both cases, the user defined function must"
                             " return a nx*ny*nz array\n")}

    cookfields = {'HRR': "HeatRelease",
                  'ENT': "Enthalpy",
                  'SRi': "IRm",
                  'RRi': "R",
                  'SDi': "DI",}

    
    def __init__(self, plotfile=None, recipe=None, outfile=None,
                 species=None, reactions=None, mech=None, 
                 pressure=None, serial=False):
        # Instantiate the PlotfileCooker parent
        super().__init__(plotfile)
        # Only work on 3 dimension data
        if self.ndims != 3:
            raise NotImplementedError(
                  ("Computing derived variables is not supported for"
                   " 2D plotfiles you can do it using the uniform"
                   " grid data created with mandoline"))
        # Output dir for the data
        # TODO: check if the file exists and add numbers at the end
        # of the outdir until the name does not exists
        if outfile is None:
            self.outdir = plotfile + "_ck"
        else:
            self.outdir = outfile
        # Flag for multiprocesssing
        self.serial = serial
        # __ Define default values for some class attributes __
        # Does the recipe needs a Cantera Solution
        self.requires_sol = True
        # Indexes of the output species in the Solution
        self.sp_indexes = []
        # Indexes of the output reactions in the Solution
        self.rx_indexes = []
        # Start of the species massfrac in the plotfile
        self.sp_start = None
        # End of the species mass frac in the plotfile
        self.sp_end = None
        # Id of the temperature field in the plotfile
        self.id_temp = None
        # Index of the O2 species in the Solution
        self.idx_O2 = None

        # Case for user recipe in a .py file
        if recipe.split('.')[-1] == 'py':
            self.recipe = self.import_user_recipe(recipe)
            self.requires_sol, self.outfields = self.get_user_recipe_info()
            if self.requires_sol:
                self.knife = chefs_knife_user_sarray
            else:
                self.knife = chefs_knife_user_pfile
            recipe = 'user'

        # Case for the recipe is a callable defined at runtime
        elif callable(recipe):
            self.recipe = recipe
            self.requires_sol, self.outfileds = self.get_user_recipe_info()
            # Define the Cantera solution
            if self.requires_sol:
                self.knife = chefs_knife_user_sarray
            else:
                self.knife = chefs_knife_user_pfile
            recipe = 'user'

        # Case for a predetermined recipe
        elif recipe in self.cookbook:
            self.recipe = self.cookbook[recipe]
            # Define the output fields and if in species mode
            if species is not None:
                self.requires_sol = True
                # Allow 'all' argument to be used
                if species == 'all' or species[0] == 'all':
                    species = [sp.name for sp in self.gas.species()]
                prefix = self.cookfields[recipe]
                self.outfields = [f"{prefix}({sp})" for sp in species]
                self.knife = chefs_knife_byspecies_field
            elif reactions is not None:
                self.requires_sol = True
                self.rx_indexes = reactions
                prefix = self.cookfields[recipe]
                self.outfields = [f"{prefix}{idx}" for idx in self.rx_indexes]
                self.knife = chefs_knife_byreaction_field
                # TODO: allow reactions to be an element and find
                # all the reactions involving the element
                # something like: (but more input validation)
                # if len(reactions) == 0 and type(reactions[0]) == str:
                #     element = reactions[0]
                # or type(reactions) == str:
                #     element = reactions
                # self.rx_indexes = [i for i, rx in gas.reactions() if element in rx]
            else:
                self.requires_sol = True
                self.outfields = [self.cookfields[recipe]]
                self.knife = chefs_knife_single_field

        # Error if no recipe can be found
        else:
            raise ValueError((f"The recipe {recipe} is not in the cookbook\n"
                              "Available recipes are:\n"
                              f"{[self.cookbookhelp[recipe] for recipe in self.cookbook]}"))

        # Define the Cantera Solution if needed
        if self.requires_sol:
            # Define the data from the Cantera solution
            self.gas, self.P, self.sp_start, self.sp_end = self.sarray_input(mech, pressure)
            # Set the global pressure and solution arrays
            self.set_global_sarrays()
            # Temperature ID to recreate the state
            self.id_temp = self.fields['temp']
            # Index of O2 to set Y(O2) = 1 when sum(Y) = 0
            self.idx_O2 = self.gas.species_index('O2')
            # Index of species in the output
            if (recipe != 'user' and species is not None):
                self.sp_indexes = [self.gas.species_index(sp) for sp in species]

    def cook(self):
        """
        Function performing the computing and output
        of the new field data
        """
        # Create the empty output plotfile and dirs
        self.make_dir_tree(self.outdir)
        # Write the global header with self.outfields
        self.write_global_header()
        # Iterate over the existing binary files and write 
        for lv in range(self.limit_level + 1):
            level_files = np.array(self.cells[lv]["files"])
            level_offsets = np.array(self.cells[lv]["offsets"])
            ncells = len(level_files)
            # All indexes of the boxes at lv
            box_indexes = np.arange(ncells)
            box_index_map = []
            mp_calls = []
            for bfpath in np.unique(level_files):
                # Boxes in the current binary file
                bf_mask = level_files == bfpath
                bf_indexes = box_indexes[bf_mask]
                # Sort the box map by offsets values
                # So they match the read order
                bf_offsets_r = level_offsets[bf_mask]
                box_index_map.append(bf_indexes[np.argsort(bf_offsets_r)])
                # New binary file
                newbfpath = os.path.join(self.outdir, 
                                         self.cell_paths[lv], 
                                         os.path.split(bfpath)[-1])
                call = {"recipe":self.recipe,
                        "bfpath":bfpath,
                        "newbfpath":newbfpath,
                        "sp_indexes":self.sp_indexes,
                        "rx_indexes":self.rx_indexes,
                        "sp_start":self.sp_start,
                        "sp_end":self.sp_end,
                        "field_indexes":self.fields,
                        "id_temp":self.id_temp,
                        "idx_O2":self.idx_O2}
                mp_calls.append(call)

            if self.serial:
                print(f"Cooking level {lv} in serial")
                output = []
                for args in tqdm(mp_calls):
                    output.append(self.knife(args))
            else:
                print(f"Cooking level {lv} in parallel")
                pool = Pool()
                output = tqdm(pool.imap(self.knife, mp_calls),
                              total=len(mp_calls))
            #Reorder the offsets to match the box order
            mapped_offsets = np.empty(len(self.boxes[lv]), dtype=int)
            mapped_mins = np.empty((len(self.boxes[lv]), 
                                    len(self.outfields)), dtype=float)
            mapped_maxs = np.empty((len(self.boxes[lv]), 
                                    len(self.outfields)), dtype=float)
            for file_idxs, bfile_result in zip(box_index_map, output):
                mapped_offsets[file_idxs] = bfile_result[0]
                mapped_mins[file_idxs, :] = bfile_result[1]
                mapped_maxs[file_idxs, :] = bfile_result[2]
            # write the level headers
            self.update_cell_header(lv, 
                                    os.path.join(self.pfile,
                                                 self.cell_paths[lv],
                                                 "Cell_H"),
                                    mapped_offsets,
                                    mapped_mins,
                                    mapped_maxs)

    def import_user_recipe(self, recipe_file):
        """
        If the path to a .py file is supplied as the recipe
        this method uses some importlib sys.path manipulation
        gymnastics to import it and save it as self.recipe 
        even if it is clearly not a proper python module
        """
        # Get the absolute path so it is importable
        recipe_file = os.path.abspath(recipe_file)
        # Get the .py location and add it to $PYTHONPATH
        recipe_location = os.path.join(*os.path.split(recipe_file)[:-1])
        sys.path.append(recipe_location)
        # Get the module spec from the python file
        mod_spec = spec_from_file_location(name='user_recipe',
                                           location=recipe_file)
        # Define the module placeholder
        recipe_module = module_from_spec(mod_spec)
        # Import the module
        mod_spec.loader.exec_module(recipe_module)
        # Get the recipe from the module
        recipe_fun = recipe_module.recipe
        return recipe_fun

    def get_user_recipe_info(self):
        """
        Find out the output field name and
        if the recipe requires a SolutionArray
        """
        # Find out how many arguments the function takes
        recipe_signature = str(signature(self.recipe))
        recipe_nargs = len(recipe_signature.split(','))
        # Case if the function only takes two arguments
        # recipe(field_indexes, box_array)
        if recipe_nargs == 2:
            solution_mode = False
        # Case if the function takes three arguments
        # recipe(field_indexes, box_array, sol_array)
        elif recipe_nargs == 3:
            solution_mode = True
        # Can't do much more input validation than checking the number
        # of positional arguments
        else:
            raise ValueError((f"User supplied recipe in {self.recipe} "
                               " with signature:\n"
                              f"recipe{recipe_signature}\n"
                              f"takes the wrong number of arguments ({recipe_nargs})"
                               "the recipe function must take either"
                               " 2 or 3 arguments"))
        # Get the recipe function docstring and use it as a field name
        # in the new plotfile
        try:
            outfields = [self.recipe.__doc__.split()[0]]
        # Default value is "user_defined"
        except AttributeError:
            outfields = ['user_defined']
        return solution_mode, outfields

    def sarray_input(self, mech, pressure):
        """
        Method to compute the input needed to create
        the SolutionArray for each box
        """
        # Create the gas instance
        gas = ct.Solution(mech)
        # Check that the plotfile has all the species
        for sp in gas.species():
            if f"Y({sp.name})" not in self.fields:
                raise ValueError((f"plotfile {plotfile} is missing field"
                                  f" Y({sp.name}) to set thermochemical"
                                  f" state with kinetics {mech}"))
        # Check that the plotfile has temperature
        if "temp" not in self.fields:
            raise ValueError(f"plotfile {plotfile} is missing 'temp' field")
        # We need pressure (could default to 1 atm)
        if pressure is None:
            raise ValueError(("Must specify simulation pressure to recreate the"
                              " thermochemical state"))
        else:
            P = pressure * ct.one_atm
        #Assume species are contiguous but verify
        species_start = self.fields[f"Y({gas.species()[0].name})"]
        species_end = species_start + len(gas.species())
        field_list = list(self.fields.keys())
        for i, sp in zip(range(species_start, species_end),
                         gas.species()):
            pltspname = field_list[i][2:-1]
            if pltspname != sp.name:
                raise ValueError(("The order of the species in the Cantera"
                                  " mechanism must match the order of the"
                                  " species in the plotfile.\n"
                                 f" Species at index {i} is {sp.name} in the"
                                 f" mechanism and {pltspname} in the plotfile"))
        # TODO: find the species indexes in the state array
        # We could create a index map if all the fields are in
        # the plotfile but not in the same order
        return gas, P, species_start, species_end

    def set_global_sarrays(self):
        """
        Set the dictionnary maps for pressure
        and placeholder solution arrays for every
        possible box shape
        """
        global SARRAYS
        global PRESSURES
        SARRAYS = {}
        PRESSURES = {}
        for shape in self.unique_box_shapes():
            SARRAYS[shape] = ct.SolutionArray(self.gas, shape)
            PRESSURES[shape] = self.P * np.ones(shape)

    def write_global_header(self):
        """
        Write the plotfile header for the plotfile with filtered
        data
        """
        hfile_path = os.path.join(self.outdir, "Header")
        with open(hfile_path, 'w') as hfile:
            # Plotfile version
            hfile.write(self.version)
            # Number of fields
            hfile.write(str(len(self.outfields)) + '\n')
            # Fields
            for v in self.outfields:
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
                tup = f"((0,0,0) ({sizes}) (0,0,0))"
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
    
    def update_cell_header(self, lv, cell_header_r, new_offsets,
                           new_mins, new_maxs):
        """
        Update the new Cell header
        """
        cell_header_w = os.path.join(self.outdir, 
                                     self.cell_paths[lv],
                                     "Cell_H")

        with open(os.path.join(os.getcwd(), cell_header_w), 'w') as ch_w:
            with open(os.path.join(os.getcwd(),cell_header_r), 'r') as ch_r:
                # First two lines
                for i in range(2):
                    l = ch_r.readline()
                    ch_w.write(l)
                # Number of fields
                _ = ch_r.readline()
                ch_w.write(f"{len(self.outfields)}\n")
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
                ch_w.write(f"{ncells},{len(self.outfields)}\n")
                for min_vals in new_mins:
                    min_vals = [str(m) for m in min_vals]
                    ch_w.write(','.join(min_vals) + ',\n')
                # Blank line
                ch_w.write('\n')
                ch_w.write(f"{ncells},{len(self.outfields)}\n")
                for max_vals in new_maxs:
                    max_vals = [str(m) for m in max_vals]
                    ch_w.write(','.join(max_vals) + ',\n')
