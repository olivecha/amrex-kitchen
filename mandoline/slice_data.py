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
        # Define the normale and slice plane coordinates
        if normal is None:
            self.cn = 0
        else:
            self.cn = normal
        # If cn = 0, then cx = 1 and cy = 2
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

