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

