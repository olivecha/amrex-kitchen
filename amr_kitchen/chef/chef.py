import time
import multiprocessing
import numpy as np
import cantera as ct
from amr_kitchen import PlofileCooker


def multiprocessing_cook(args):
    """
    __ Multiprocessing function __
    Cook a single cell binary file
    """
    # New offsets
    offsets = []
    gas_attr = args['gas_attr']
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

class Chef(PlotfileCooker)

    cookbook = {'HRR': "heat_release_rate",
                'ENT': "enthalpy_mass",
                'SRi': "net_production_rates",
                'SDi': "mix_diff_coeffs_mass"}
    cookfields = {'HRR': "HeatRelease",
                  'ENT': "Enthalpy",
                  'SRi': "IR",
                  'SDi': "DI",}
                  

    
    def __init__(self, plotfile=None, recipe=None, 
                 species=None, mech=None):

        # Instantiate the PlotfileCooker parent
        super().__init__(plotfile)
        # Output dir for the data
        self.outdir = plotfile + "_c"
        # Store the recipe
        self.recipe=recipe
        # Only work on 3 dimension data
        if self.ndims != 3:
            raise NotImplementedError(
                  ("Computing derived variables is not supported for"
                   " 2D plotfiles you can do it using the uniform"
                   " grid data created with mandoline"))
        # Create the Cantera Solution object
        self.gas = ct.Solution(mech)
        # Validate the state is available in the plotfile

        # Define the output fields and if in species mode
        if species is not None:
            self.sp_mode = True
            self.sp_indexes = [self.gas.species_index(sp) for sp in species]
            prefix = self.cookfields[self.recipe]
            self.outfields = [f"{prefix}({sp})" for sp in species]
        else:
            self.sp_mode = False
            self.sp_indexes = []
            self.outfields = [self.cookfields[self.recipe]]

    def cook(self):
        """
        Function performing the computing and output
        of the new field data
        """

        







