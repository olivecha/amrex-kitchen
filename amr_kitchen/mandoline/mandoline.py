import os
import time
import multiprocessing
import shutil
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from amr_kitchen import HeaderData
from .utils import expand_array
from .blades import slice_box, plate_box


class Mandoline(HeaderData):
    """
    Class containing the slicing data to supply to the
    multiprocessing function
    """
    coordnames = {0:'x',
                  1:'y',
                  2:'z',
                  3:'2D'}

    def __init__(self, plotfile, fields=None, limit_level=None,
                 serial=False, verbose=None):
        """
        Constructor for the mandoline object
        ----
        plotfile:    str path to the plotfile to slice

        normal:       int slice normal coordinate 0:x, 1:y, 2:z

        fields:      Field(s) for which a slice is computed

        pos:         slice position in plotfile coordinates, defaults to
                     the center of the domain

        limit_level: Maximum data level used to create the slice, defaults
                     to the maximum level in the file 

        serial:      Wheter or not to run the slicing function in parallel

        verbose:     Verbosity level for the informational output

        """
        # Verbosity level
        if verbose is None:
            self.v = 1
        else:
            self.v = verbose
        # Do stuff in serial
        self.serial = serial
        # Read the header data into the HeaderReader class
        header_start = time.time()
        super().__init__(plotfile, limit_level=limit_level)
        header_time = time.time() - header_start
        self.header_load_info(header_time)
        # Define the fields in the slice
        self.slicefields, self.fidxs, self.do_grid = self.parse_input_fields(fields)
        # Number of fields
        self.nfidxs = len([i for i in self.fidxs if i is not None])
        # Slicing normal and position
        self.cn = None
        self.cx = None
        self.cy = None
        self.pos = None

    def slice(self, normal=None, pos=None, outfile=None, fformat=None, **pltkwargs):
        """
        Slicing function of the mandoline
        """
        if self.ndims == 2:
            return self.plate(outfile, fformat)
        # Define normal and position
        if normal is None:
            normal = self.cn
        if pos is None:
            pos = self.pos
        self.cn, self.cx, self.cy, self.pos = self.define_slicing_coordinates(normal, pos)
        # Default output is user friendly image
        if fformat is None:
            fformat = "image"

        # Object to store the slices
        plane_data = []
        # For a given level
        for Lv in range(self.limit_level + 1):
            # Box reading timer
            read_start = time.time()
            # Multiprocessing inputs
            pool_inputs = self.compute_mpinput_3d(Lv)
            # Read the data in parallel (or serial)
            # The box reader returns the slice at both sides of the slicing
            # Plane if available (i.e. not a boundary, or between boxes). 
            # This allows handling the case when the slice is between boxes
            # This function interpolate between the two planes and returns
            if self.serial:
                plane_data.append(list(map(slice_box, pool_inputs)))
            else:
                with multiprocessing.Pool() as pool:
                    plane_data.append(pool.map(slice_box, pool_inputs))

            if self.v > 0:
                print(f"Time to read Lv {Lv}:", 
                      np.around(time.time() - read_start, 2))

        # Interpolation timer
        interp_start = time.time()
        # Process the sliced box to create uniform grids
        if fformat in ["array", "return", "image"]:
            all_data = self.reducemp_data_ortho(plane_data)
            if fformat == "return":
                output = self.format_array_output(all_data)
                return output
            if fformat == "array":
                if outfile is None:
                    outfile = self.default_output_path()
                outfile += ".npz"
                output = self.format_array_output(all_data)
                np.savez_compressed(outfile, **output)
                return
            if fformat == "image":
                self.plot_slice(all_data, outfile, **pltkwargs)


        # For plotfiles we keep the data needed to reconstruct the multilevel grid
        elif fformat == "plotfile":
            output_start = time.time()
            # Interpolate for each level
            all_data_bylevel, indexes, headers  = self.interpolate_bylevel(plane_data)
            # Define the plotfile name
            if outfile is None:
                outfile = self.default_output_path()
            # Create the plotfile dir
            try:
                os.mkdir(outfile)
            except FileExistsError:
                shutil.rmtree(outfile)
                os.mkdir(outfile)
            # Rewrite the header
            with open(os.path.join(outfile, "Header"), "w") as hfile:
                self.write_2d_slice_global_header(hfile,
                                                  self.fields_in_slice(),
                                                  indexes)
            # Write the level data
            for lv in range(self.limit_level + 1):
                self.write_cell_data_at_level(outfile,
                                              lv,
                                              all_data_bylevel[lv],
                                              indexes[lv])
            if self.v > 0:
                print("Time to save AMReX plotfile: ", 
                      np.around(time.time() - output_start, 2))

    def thick_slice(self, normal=None, pos=None, 
                    outfile=None, thickness=None):
        """
        Non interpolated thick slices saved as AMReX plotfiles
        """
        # Alterate the position to contain uniform points
        # At the lowest level in the given thickness
        # Define the lowest Level grid in the normal
        # direction:
        self.cn, self.cx, self.cy, self.pos = self.define_slicing_coordinates(normal, pos)
        geo_lo = self.geo_low[self.cn] + self.dx[0][self.cn]
        geo_hi = self.geo_high[self.cn] - self.dx[0][self.cn]
        ngrid = np.arange(geo_lo, geo_hi, 
                          self.grid_sizes[0][self.cn])
        print(ngrid)


    def plate(self, outfile=None, fformat=None, **pltkwargs):
        """
        Alternate function to "slice" 2D plotfiles
        This is a conversion to uniform covering grid
        
        output: output path name

        fformat: output file format "image", "array", "plotfile"

        if fformat = "return" the array output is returned without
        saving it to the output path
        """
        #TODO: Use colander to remove fields from 2D plotfiles
        if fformat == "plotfile":
            raise NotImplementedError("TODO (you can use amrex-kitchen/"
                                      "colander")
        # Define the output
        if outfile is None:
            outfile = self.default_output_path()
        # Multiprocessing
        if not self.serial:
            pool = multiprocessing.Pool()
        # The slice is just the header data
        # Object to store the slices
        plane_data = []
        # For a given level
        for Lv in range(self.limit_level + 1):
            # Box reading timer
            read_start = time.time()
            # Compute the pool inputs
            pool_inputs = self.compute_mpinput_2d(Lv)
            # Read the data in parallel or serial
            if self.serial:
                plane_data.append(list(map(plate_box, pool_inputs)))
            else:
                plane_data.append(pool.map(plate_box, pool_inputs))

            if self.v > 0:
                print(f"Time to read Lv {Lv}:", 
                      np.around(time.time() - read_start, 2))

        # Broadcast the data to empty arrays
        all_data = [self.limit_level_arr() for _ in range(self.nfidxs)]

        if self.do_grid:
            grid_level = self.limit_level_arr()

        # Parse the multiprocessing output
        # Do levels sequentially to update with finer data
        for Lv in range(self.limit_level + 1):
            for out in plane_data[Lv]:
                # Add the slices if they are defined
                xa, xo = out['sx']  # x slice
                ya, yo = out['sy']  # y slice
                # add the field data
                for i in range(self.nfidxs):
                    all_data[i][xa:xo, ya:yo] = out['data'][i]
                # add the normal coordinate
                # broadcast the grid level to the grid if needed
                if self.do_grid:
                    grid_level[xa:xo, ya:yo] = out['level']

        if self.do_grid:
            all_data.append(grid_level)

        for i, data in enumerate(all_data):
            all_data[i] = data.T

        if fformat == "return":
            output = self.format_array_output(all_data)
            return output

        elif fformat == "array":
            output = self.format_array_output(all_data)
            outfile += ".npz"
            np.savez_compressed(outfile, **output)
            return

        elif fformat == "image":
            self.plot_slice(all_data, outfile, **pltkwargs)


    def infer_figure_size(self):
        """
        Determine a good figure size given 
        the dataset aspect ratio
        """
        # Aspect ratio
        Lx = self.geo_low[self.cx] - self.geo_high[self.cx]
        Ly = self.geo_low[self.cy] - self.geo_high[self.cy]
        aspect = Lx / Ly
        # Horizontal
        if aspect > 1:
            size_x = 7
            size_y = int(size_x / aspect)
            figsize = (size_x + 1, size_y)
        elif aspect < 1:
            size_y = 7
            size_x = int(size_y * aspect)
            figsize = (size_x + 1, size_y)
        else:
            figsize=(8, 6)
        return figsize

    def plot_slice(self, all_data, outfile, 
                   uselog=False, cmap=None, vmin=None, vmax=None):
        """
        Plot the slice data using matplotlib
        """

        # Get the formatted data
        pltdata = self.format_array_output(all_data)

        # Use log scale ?
        if uselog:
            norm = matplotlib.colors.LogNorm()
        else:
            norm = None
        # User supplied colormap
        if cmap in plt.colormaps():
            pass
        else:
            cmap = 'jet'
        figsize = self.infer_figure_size()

        # A figure per field
        plot_names = [f for f in pltdata if f not in ['x', 'y']]
        for name in plot_names:

            if self.v > 0:
                plot_start = time.time()
                print(f"Plotting {name} with vmax = {vmax}")
            # Pretty large figure
            fig = plt.figure(figsize=figsize)

            # Set face color to colormap minimum so masked
            # Values look good
            plt.gca().set_facecolor(plt.cm.get_cmap(cmap)(0)) 

            # Plot the slice
            plt.pcolormesh(pltdata['x'],
                           pltdata['y'],
                           pltdata[name],
                           vmin=vmin,
                           vmax=vmax,
                           norm=norm,
                           cmap=cmap)

            # Add a colorbar  
            plt.colorbar(label=f'{name}')
            ax = plt.gca()
            ax.set_aspect('equal')
            ax.set_xlabel(f"{self.coordnames[self.cx]} [m]")
            ax.set_ylabel(f"{self.coordnames[self.cy]} [m]")

            if self.v > 0:
                print(f"Done! ({np.around(time.time() - plot_start, 2)} s)")

            if self.v > 0:
                save_start = time.time()
                print(f"Saving {name} plot...")
            
            # save and close
            if outfile is None:
                outfile = self.default_output_path(fieldname=name)
            fig.savefig(outfile, dpi=500)
            plt.close(fig)

            if self.v > 0:
                print(f"Done! ({np.around(time.time() - save_start, 2)} s)")


    def compute_mpinput_2d(self, lv):
        """
        Compute the multiprocessing input for 2D plotfiles
        For the current AMR Level (lv)
        """
        # Multiprocessing inputs
        pool_inputs = []
        for indexes, cfile, offset, box in zip(self.cells[lv]['indexes'],
                                               self.cells[lv]['files'],
                                               self.cells[lv]['offsets'],
                                               self.boxes[lv]):
            # Everything needed by the slice reader
            p_in  = {'cx':0,
                     'cy':1,
                     'dx':self.dx,
                     'limit_level':self.limit_level,
                     'fidxs':self.fidxs,
                     'Lv':lv,
                     'indexes':indexes,
                     'cfile':cfile,
                     'offset':offset,
                     'box':box}
            pool_inputs.append(p_in) # Add to inputs
        return pool_inputs

    def compute_mpinput_3d(self, lv):
        """
        Find the intersecting boxes and add the to the
        multiprocessing input for a given level
        """
        pool_inputs = []
        # For each box in that level
        for idx, box in enumerate(self.boxes[lv]):
            # Check if the box intersect the slicing plane
            if (box[self.cn][0] <= self.pos and 
                box[self.cn][1] >= self.pos):
                # ----
                # Here, intersecting boxes at each level
                # Are added to the slice, even if there are
                # Higher level boxes covering the lower level
                # Box region in the plane. This is a little
                # wastefull as more data is readed than needed
                # to create the plane, but finding if the box is
                # Completely covered is a nightmare (and slower).
                # -----
                # Everything needed by the slice reader
                p_in  = {'cx':self.cx,
                         'cy':self.cy,
                         'cn':self.cn,
                         'dx':self.dx,
                         'pos':self.pos,
                         'limit_level':self.limit_level,
                         'fidxs':self.fidxs,
                         'Lv':lv,
                         'bidx':idx,
                         'indexes':self.cells[lv]['indexes'][idx],
                         'cfile':self.cells[lv]['files'][idx],
                         'offset':self.cells[lv]['offsets'][idx],
                         'box':box}
                pool_inputs.append(p_in) # Add to inputs
        return pool_inputs

    def define_slicing_coordinates(self, normal=None, pos=None):
        """
        Parse the normal and position input to define the slicing
        coordinates (normal in in slice plane)
        """
        # 2D case
        if self.ndims == 2:
            cx, cy = 0, 1
            cn = 3
            pos = 0.
        # 3D
        else:
            # Define the normale and slice plane coordinates
            # If cn = 0, then cx = 1 and cy = 2
            if normal is None:
                cn = 0
            else:
                cn = normal
            cx, cy = [i for i in range(3) if i != cn]

            # define and store the slice position
            if pos is None:
                # Domain center
                pos = (self.geo_high[cn] - self.geo_low[cn])/2
            else:
                if (pos < self.geo_low[cn] or
                    pos > self.geo_high[cn]):
                    coord = self.coordnames[cn]
                    low = self.geo_low[cn]
                    high = self.geo_high[cn]
                    raise ValueError(("Slicing plane at "
                                     f"{coord} = {pos} "
                                      "is outside the domain bounds: "
                                     f"{low} < {coord} < {high}"))

        return cn, cx, cy, pos

    def parse_input_fields(self, fields):
        """
        Parse the fields parsed argument input and define
        the indexes and names of the fields to slice
        """
        # Find out what to slice
        # Default field (same as AMReX fsnapshot.cpp)
        if fields is None:
            fields = ['density']
        # String input (for outside the cli)
        if type(fields) == str:
            fields = [fields]
        # Special case for all fields
        if 'all' in fields:
            fields = ['all']
            fidxs = list(range(self.nvars))
            # With grid level
            fidxs.append(None)
        # Iterate over what we are sure is a list
        else:
            fidxs = []
            for f in fields:
                # Special case for grid_data
                if f == 'grid_level':
                    fidxs.append(None)
                # Special case for all fields
                else:
                    # All field indexes
                    fidxs.append(self.field_index(f))
        # Flag for grid level
        if None in fidxs:
            do_grid_level = True
        else:
            do_grid_level = False

        return fields, fidxs, do_grid_level

    def format_array_output(self, all_data):
        """
        Format the interpolated data to export as array
        """
        # Make a dict with output
        x_grid, y_grid = self.slice_plane_coordinates()
        output = {'x': x_grid,
                  'y': y_grid,}
        # Store in output
        for i, name in enumerate(self.fields_in_slice()):
            output[name] = all_data[i]

        # Works if only grid_level
        if self.do_grid:
            output['grid_level'] = all_data[-1]

        return output


    def fields_in_slice(self):
        """
        Get the field names from the computed field indexes
        """
        all_names = [name for name in self.fields]
        field_names = [all_names[idx] for idx in self.fidxs if idx is not None]
        return field_names

    def default_output_path(self, fieldname=None):
        """
        Define the default output path from the plotfile
        Supplied in the command line arguments
        """
        # Split the input file
        outroot, plotname = os.path.split(self.pfile)
        plotnum = "_" + plotname.replace("plt", '').replace('_', '')
        if self.ndims == 3:
            # Slice fraction in domain
            width = self.geo_high[self.cn] - self.geo_low[self.cn]
            pos_frac = self.pos / width
            # Looks like Sx05000density45000
            slicename = f"S{self.coordnames[self.cn]}"
            pos_string = f"{pos_frac:.4f}".replace('.', '')
            slicename = slicename + pos_string
        elif self.ndims == 2:
            # Looks like S2Ddensity45000
            slicename = f"S2D"
        # Single field slice
        if fieldname is not None:
            fieldstring = fieldname
            fieldstring = fieldstring.replace('(', '').replace(')', '')
            fieldstring = fieldstring[:7]
            slicename = slicename + fieldstring
        elif len(self.slicefields) == 1:
            fieldstring = self.slicefields[0]
            fieldstring = fieldstring.replace('(', '').replace(')', '')
            fieldstring = fieldstring[:7]
            slicename = slicename + fieldstring
        # Multi field slice
        return os.path.join(outroot, slicename + plotnum)

    def header_load_info(self, ltime):
        """
        Print message acoording to the time it took
        to read the header
        """
        if self.v > 0:
            print(f"Time to read header files: {ltime:.2f} s")
            if self.serial:
                pass
            elif ltime < 0.1:
                print("This plotfile seems small you may consider "
                      "reading it in serial")
            else:
                pass

    def limit_level_arr(self):
        """
        Return an empty numpy array with the dimensions of the
        limit_level grid
        """
        if (self.cx is None or
            self.cy is None):
            self.cn, self.cx, self.cy, self.pos = self.define_slicing_coordinates()

        shape = self.grid_sizes[self.limit_level][[self.cx, self.cy]]
        arr = np.empty(shape)
        return arr

    def reducemp_data_ortho(self, plane_data):
        """
        Parse the multiprocessing output and reduce it to
        Left and Right arrays for both sides of the plane
        (For orthogonal planes)
        """
        # Array for the "left" side of the plane (could be down whatever)
        left = {'data':[self.limit_level_arr() for _ in range(self.nfidxs)],
                'normal':self.limit_level_arr()}
        # Array for the "right" or up side of the plane
        # The only convention is that "left" < slice_coordinate < "right"
        right = {'data':[self.limit_level_arr() for _ in range(self.nfidxs)],
                'normal':self.limit_level_arr()}

        if self.do_grid:
            grid_level = {'left':self.limit_level_arr(),
                          'right':self.limit_level_arr()}
        else:
            grid_level = None

        # Parse the multiprocessing output
        # Do levels sequentially to update with finer data
        for Lv in range(self.limit_level + 1):
            first_grid_pt = self.geo_low[self.cn] + self.dx[Lv][self.cn]/2
            last_grid_pt = self.geo_high[self.cn] - self.dx[Lv][self.cn]/2
            for output in plane_data[Lv]:
                # Add the slices if they are defined
                if output[0] is not None:
                    out = output[0]
                    xa, xo = out['sx']  # x slice
                    ya, yo = out['sy']  # y slice
                    # add the field data
                    for i, arr in enumerate(left['data']):
                        left['data'][i][xa:xo, ya:yo] = out['data'][i]
                    # add the normal coordinate
                    left['normal'][xa:xo, ya:yo] = out['normal']
                    # broadcast the grid level to the grid if needed
                    if self.do_grid:
                        grid_level['left'][xa:xo, ya:yo] = out['level']
                    # Case when the slice plane is after the last grid
                    if np.isclose(out['normal'], last_grid_pt):
                        # add the field data
                        for i, arr in enumerate(right['data']):
                            right['data'][i][xa:xo, ya:yo] = out['data'][i]
                        # add the normal coordinate
                        right['normal'][xa:xo, ya:yo] = out['normal']
                        
                # Same for the right side
                if output[1] is not None:
                    out = output[1]
                    xa, xo = out['sx']  # x slice
                    ya, yo = out['sy']  # y slice
                    for i, arr in enumerate(left['data']):
                        right['data'][i][xa:xo, ya:yo] = out['data'][i]
                    right['normal'][xa:xo, ya:yo] = out['normal']
                    # broadcast the grid level to the grid if needed
                    if self.do_grid:
                        grid_level['right'][xa:xo, ya:yo] = out['level']
                    # Case when the slice plane is before the first grid
                    if np.isclose(out['normal'], first_grid_pt):
                        # add the field data
                        for i, arr in enumerate(left['data']):
                            left['data'][i][xa:xo, ya:yo] = out['data'][i]
                        # add the normal coordinate
                        left['normal'][xa:xo, ya:yo] = out['normal']
        # Do the linear interpolation if normals are not the same
        # Empty arrays for the final data
        all_data = []
        # Boolean slicing array for when the points are the same
        bint = ~np.isclose(left['normal'], right['normal'])
        # Iterate with the number of fields
        for i in range(self.nfidxs):
            data = self.limit_level_arr()
            # Linear interpolation
            term1 = left['data'][i][bint] * (right['normal'][bint] - self.pos) 
            term2 = right['data'][i][bint] * (self.pos - left['normal'][bint])
            term3 = right['normal'][bint] - left['normal'][bint]
            data[bint] =  (term1 + term2) / term3
            # Could be either
            data[~bint] = right['data'][i][~bint]
            # For some reason
            all_data.append(data.T)

        if self.do_grid:
            # Concatenate both sides
            all_levels = np.stack([grid_level['right'], grid_level['left']])
            # I guess the min is better for debuging
            # All things considered they should be pretty similar
            all_data.append(np.min(all_levels, axis=0).T)
        return all_data

    def interpolate_bylevel(self, plane_data):
        """
        Interpolate sliced output but keep the level separated
        to be able to save the data to an amrex plotfile
        """
        # Base structure to store interpolated data
        lvarr = {'data':[self.limit_level_arr() for _ in range(self.nfidxs)],
                 'normal':self.limit_level_arr()}
        # An array for each level (left side)
        left = [lvarr for _ in range(self.limit_level + 1)]
        
        # Array for the "right" or up side of the plane
        # The only convention is that "left" < slice_coordinate < "right"
        right = [lvarr for _ in range(self.limit_level + 1)]
        # Keep track of the box headers
        box_headers = []
        box_indexes = []

        # Parse the boxes in the output
        for lv in range(self.limit_level + 1):
            lv_box_headers = {}
            lv_box_indexes = []
            for output in plane_data[lv]:
                lv_box_headers[output[3]] = output[2]
                lv_box_indexes.append(output[3])
                # Add the slices if they are defined
                if output[0] is not None:
                    out = output[0]
                    xa, xo = out['sx']  # x slice
                    ya, yo = out['sy']  # y slice
                    # add the field data
                    for i in range(self.nfidxs):
                        left[lv]['data'][i][xa:xo, ya:yo] = out['data'][i]
                    # add the normal coordinate
                    left[lv]['normal'][xa:xo, ya:yo] = out['normal']

                # Same for the right side
                if output[1] is not None:
                    out = output[1]
                    xa, xo = out['sx']  # x slice
                    ya, yo = out['sy']  # y slice
                    # every field
                    for i in range(self.nfidxs):
                        right[lv]['data'][i][xa:xo, ya:yo] = out['data'][i]
                    # normal
                    right[lv]['normal'][xa:xo, ya:yo] = out['normal']
                
            box_indexes.append(lv_box_indexes)
            box_headers.append(lv_box_headers)

        # Ouput list
        all_data = []
        all_indexes = []
        # Interpolate for each level
        for lv in range(self.limit_level + 1):
            level_data = []
            # Boolean array where the normals are equal
            bint = ~np.isclose(left[lv]['normal'], right[lv]['normal'])
            # Interpolate each field
            for i in range(self.nfidxs):
                # Empty array
                data = self.limit_level_arr()
                # Linear interpolation
                term1 = left[lv]['data'][i][bint]\
                        * (right[lv]['normal'][bint] - self.pos) 
                term2 = right[lv]['data'][i][bint]\
                        * (self.pos - left[lv]['normal'][bint])
                term3 = right[lv]['normal'][bint] - left[lv]['normal'][bint]
                data[bint] =  (term1 + term2) / term3
                # Non interpolated points
                data[~bint] = right[lv]['data'][i][~bint]
                level_data.append(data)
            all_data.append(level_data)

        return all_data, box_indexes, box_headers

    def slice_plane_coordinates(self):
        """
        Return the x and y grid corresponding to the slice
        array coordinates to plot the data after output
        """
        # Define coord sys
        if (self.cx is None or
            self.cy is None):
            self.cn, self.cx, self.cy, self.pos = self.define_slicing_coordinates()
        # grids a from x_lo + dx/2 to y_hi - dx/2 for cell centered data
        x_grid = np.linspace(self.geo_low[self.cx]\
                             + self.dx[self.limit_level][self.cx]/2,
                             self.geo_high[self.cx]\
                             - self.dx[self.limit_level][self.cx]/2,
                             self.grid_sizes[self.limit_level][self.cx])

        y_grid = np.linspace(self.geo_low[self.cy]\
                             + self.dx[self.limit_level][self.cy]/2,
                             self.geo_high[self.cy]\
                             - self.dx[self.limit_level][self.cy]/2,
                             self.grid_sizes[self.limit_level][self.cy])

        return x_grid, y_grid

    def write_2d_slice_global_header(self, fobj, fnames, indexes):
        """
        Write the plotfile header for the plotfile corresponding
        to the 2D slice
        """
        # Plotfile version
        fobj.write(self.version)
        # Number of fields
        fobj.write(str(self.nfidxs) + '\n')
        # Fields
        for f in fnames:
            fobj.write(f + '\n')
        # Always 2D
        fobj.write("2\n")
        # Time
        fobj.write(str(self.time) + '\n')
        # Max level
        fobj.write(str(self.limit_level) + '\n')
        # Lower bounds
        fobj.write(f"{self.geo_low[self.cx]} {self.geo_low[self.cy]}\n")
        # Upper bounds
        fobj.write(f"{self.geo_high[self.cx]} {self.geo_high[self.cy]}\n")
        # Refinement factors
        factors = self.factors[:self.limit_level]
        fobj.write(' '.join([str(f) for f in factors]) + '\n')
        # Grid sizes
        # Looks like ((0,0,0) (7,7,7) (0,0,0))
        tuples = []
        for lv in range(self.limit_level + 1):
            tup = (f"((0,0) ({self.grid_sizes[lv][self.cx] - 1},"
                          f"{self.grid_sizes[lv][self.cy] - 1})"
                           " (0,0))")
            tuples.append(tup)
        fobj.write(' '.join(tuples) + '\n')
        # By level step numbers
        step_numbers = self.step_numbers[:self.limit_level + 1]
        fobj.write(' '.join([str(n) for n in step_numbers]) + '\n')
        # Grid resolutions
        for lv in range(self.limit_level + 1):
            fobj.write(f"{self.dx[lv][self.cx]} {self.dx[lv][self.cy]}\n")
        # Coordinate system
        fobj.write(str(self.sys_coord))
        # Zero for parsing
        fobj.write("0\n")
        # Write the boxes
        for lv in range(self.limit_level + 1):
            # Write the level info
            fobj.write(f"{lv} {len(indexes[lv])} {self.time}\n")
            # Write the level step
            fobj.write(f"{self.step_numbers[lv]}\n")
            # Write the 2D boxes
            for idx in indexes[lv]:
                box = self.boxes[lv][idx]
                fobj.write(f"{box[self.cx][0]} {box[self.cx][1]}\n")
                fobj.write(f"{box[self.cy][0]} {box[self.cy][1]}\n")
            # Write the Level path info
            fobj.write(f"Level_{lv}/Cell\n")

    def write_cell_data_at_level(self, outfile, lv, lvdata,  indexes):
        """
        Write the cell data and cell header at a given level
        outfile: plotfile directory
        lv: current amrlevel
        lvdata: list of 2D arrays for each field in the dataset at the level
        indexes: box indexes for the current level
        """
        # Make the directory
        os.mkdir(os.path.join(outfile, f"Level_{lv}"))
        # Factor between lv and limit level data
        factor = 2**(self.limit_level - lv)
        # Determine how many cell files to use
        # I dont know what amrex does lets target 1 MB per cell file
        total_size = 0
        cell_indexes = []
        for idx in indexes:
            cidx = self.cells[lv]["indexes"][idx]
            cell_indexes.append(cidx)
            size_x = cidx[1][self.cx] - cidx[0][self.cx] + 1
            size_y = cidx[1][self.cy] - cidx[0][self.cy] + 1
            total_size += size_x * size_y * self.nfidxs * 8

        # Divide by 1 MB and add one so there is a file
        nfiles = total_size // int(1e6) + 1
        chunk_size = len(cell_indexes) // nfiles

        # Filenames
        fnames = [f"Cell_D_{n:05d}" for n in range(nfiles + 1)]
        # Store offsets for each file
        offsets = []
        field_max_vals = []
        field_min_vals = []
        # For each chunk
        for cfile, i in zip(fnames, range(0, len(cell_indexes), chunk_size)):
            with open(os.path.join(outfile, f"Level_{lv}", cfile), "wb") as bfile:
                curr_offsets = []
                subcells_indexes = cell_indexes[i:i+chunk_size]
                for idxs in subcells_indexes:
                    offset = bfile.tell()
                    curr_offsets.append(offset)
                    curr_data = []
                    curr_field_min = []
                    curr_field_max = []
                    # Compute the slice indexes for the global grid
                    x_start = idxs[0][self.cx] * factor
                    x_stop = (idxs[1][self.cx] + 1) * factor
                    y_start = idxs[0][self.cy] * factor
                    y_stop = (idxs[1][self.cy] + 1) * factor
                    # For each field
                    for arr in lvdata:
                        # Read the data
                        data = arr[x_start:x_stop, y_start:y_stop]
                        # Downsize to current level
                        data = data[::factor, ::factor]
                        # Store the min/max values in the cell
                        curr_field_min.append(np.min(data))
                        curr_field_max.append(np.max(data))
                        # Append the cell data to be written in binary
                        curr_data.append(data.flatten(order="F"))
                    # Append the min max values for each cell
                    field_min_vals.append(curr_field_min)
                    field_max_vals.append(curr_field_max)
                    # Redefine the header
                    new_header = ("FAB ((8, (64 11 52 0 1 12 0 1023)),"
                                  "(8, (8 7 6 5 4 3 2 1)))")
                    new_header += (f"(({idxs[0][self.cx]},{idxs[0][self.cy]}) "
                                   f"({idxs[1][self.cx]},{idxs[1][self.cy]}) ")
                    new_header +=  f"(0,0)) {self.nfidxs}\n"
                    bfile.write(new_header.encode("ascii"))
                    bfile.write(np.hstack(curr_data).tobytes())
                offsets.append(curr_offsets)

        # Write the cell header
        with open(os.path.join(outfile, f"Level_{lv}", "Cell_H"), "w") as hfile:
            # Always ones
            hfile.write("1\n")
            hfile.write("1\n")
            # Number of fields
            hfile.write(f"{self.nfidxs}\n")
            # Always 0
            hfile.write("0\n")
            # Number of cells
            hfile.write(f"({len(indexes)} 0\n")
            # All cell indexes
            for cidxs in cell_indexes:
                hfile.write(f"(({cidxs[0][self.cx]},{cidxs[0][self.cy]}) "
                             f"({cidxs[1][self.cx]},{cidxs[1][self.cy]}) "
                              "(0,0))\n")
            hfile.write(")\n")
            hfile.write(f"{len(indexes)}\n")
            for bfname, bfoffsets in zip(fnames, offsets):
                for offset in bfoffsets:
                    hfile.write(f"FabOnDisk: {bfname} {offset}\n")
            # Empty line
            hfile.write('\n')
            # Number of cells number of fields
            hfile.write(f"{len(indexes)},{self.nfidxs}\n")
            # For each cell write the minimum value
            for min_arr in field_min_vals:
                # Scientific notation with 16 digits
                line = ','.join([f"{num:.16e}" for num in min_arr])
                hfile.write(line + ',\n')
            # Empty line
            hfile.write('\n')
            # Number of cells number of fields
            hfile.write(f"{len(indexes)},{self.nfidxs}\n")
            # For each cell write the maximum value
            for max_arr in field_max_vals:
                # Scientific notation with 16 digits
                line = ','.join([f"{num:.16e}" for num in max_arr])
                hfile.write(line + ',\n')
