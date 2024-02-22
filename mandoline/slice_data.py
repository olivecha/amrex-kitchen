import os
import numpy as np
from mandoline import HeaderData


class SliceData(HeaderData):
    """
    Class containing the slicing data to supply to the
    multiprocessing function
    """

    def __init__(self, plotfile, fields=None, normal=None, pos=None, limit_level=None):
        """
        Constructor for the mandoline object
        ----
        plotfile:    str path to the plotfile to slice

        coord:       int slice normal coordinate 0:x, 1:y, 2:z

        fields:      Field(s) for which a slice is computed

        pos:         slice position in plotfile coordinates, defaults to
                     the center of the domain

        limit_level: Maximum data level used to create the slice, defaults
                     to the maximum level in the file 

        """
        # Read the header data from the parent class
        super().__init__(plotfile, limit_level=limit_level)
        # 2D case
        if self.ndims == 2:
            self.cx, self.cy = 0, 1
            self.cn = '2D'
            self.pos = 0.
        # 3D
        else:
            # Define the normale and slice plane coordinates
            # If cn = 0, then cx = 1 and cy = 2
            if normal is None:
                self.cn = 0
            else:
                self.cn = normal
            self.cx, self.cy = [i for i in range(3) if i != self.cn]

            # define and store the slice position
            if pos is None:
                # Domain center
                self.pos = (self.geo_high[self.cn] - self.geo_low[self.cn])/2
            else:
                self.pos = pos

        # Find out what to slice
        # Default field (same as AMReX fsnapshot.cpp)
        if fields is None:
            fields = ['density']
        # String input (for outside the cli)
        if type(fields) == str:
            fields = [fields]
        # Special case for all fields
        if 'all' in fields:
            self.fidxs = list(range(self.nvars))
            # With grid level
            self.fidxs.append(None)
        # Iterate over what we are sure is a list
        else:
            self.fidxs = []
            for f in fields:
                # Special case for grid_data
                if f == 'grid_level':
                    self.fidxs.append(None)
                # Special case for all fields
                else:
                    # All field indexes
                    self.fidxs.append(self.field_index(f))
        self.nfidxs = len([i for i in self.fidxs if i is not None])
        if None in self.fidxs:
            self.do_grid_level = True
        else:
            self.do_grid_level = False

    def limit_level_arr(self):
        """
        Return an empty numpy array with the dimensions of the
        limit_level grid
        """
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

        if self.do_grid_level:
            grid_level = {'left':self.limit_level_arr(),
                          'right':self.limit_level_arr()}
        else:
            grid_level = None

        # Parse the multiprocessing output
        # Do levels sequentially to update with finer data
        for Lv in plane_data:
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
                    if self.do_grid_level:
                        grid_level['left'][xa:xo, ya:yo] = out['level']
                # Same for the right side
                if output[1] is not None:
                    out = output[1]
                    xa, xo = out['sx']  # x slice
                    ya, yo = out['sy']  # y slice
                    for i, arr in enumerate(left['data']):
                        right['data'][i][xa:xo, ya:yo] = out['data'][i]
                    right['normal'][xa:xo, ya:yo] = out['normal']
                    # broadcast the grid level to the grid if needed
                    if self.do_grid_level:
                        grid_level['right'][xa:xo, ya:yo] = out['level']
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

        if self.do_grid_level:
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
        for lv in plane_data:
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
        factors = self.factors[:self.limit_level + 1]
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
            cidx = self.cells[f"Lv_{lv}"]["indexes"][idx]
            cell_indexes.append(cidx)
            size_x = cidx[1][self.cx] - cidx[0][self.cx] + 1
            size_y = cidx[1][self.cy] - cidx[0][self.cy] + 1
            total_size += size_x * size_y * self.nfidxs * 8

        # Divide by 1 MB and add one so there is a file
        nfiles = total_size // int(1e6) + 1
        chunk_size = len(cell_indexes) // nfiles

        # Filenames
        fnames = [f"Cell_D_{n:05d}" for n in range(nfiles)]
        # Store offsets for each file
        offsets = []
        # For each chunk
        for cfile, i in zip(fnames, range(0, chunk_size, nfiles)):
            with open(os.path.join(outfile, f"Level_{lv}", cfile), "wb") as bfile:
                curr_offsets = []
                subcells_indexes = cell_indexes[i:i+chunk_size]
                for idxs in subcells_indexes:
                    offset = bfile.tell()
                    curr_offsets.append(offset)
                    curr_data = []
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
                        curr_data.append(data.flatten(order="F"))
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

                






                




        
            
        lv_shape = self.grid_sizes[lv]
        size_x = lv_shape[self.cx]
        size_y = lv_shape[self.cy]
        nfields = self.nfidxs
        lv_size = size_x * size_y * nfields * 8
        
        


