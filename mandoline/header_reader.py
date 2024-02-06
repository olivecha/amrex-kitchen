import os
import shutil
import linecache
import numpy as np
from tqdm import tqdm


class HeaderData(object):

    def __init__(self, plotfile, limit_level=None):
        """
        Parse the header data and save as attributes
        """
        self.pfile = plotfile
        filepath = os.path.join(plotfile, 'Header')
        with open(filepath) as hfile:
            _ = hfile.readline()
            # field names
            self.nvars = int(hfile.readline())
            self.fields = {}
            for i in range(self.nvars):
                self.fields[hfile.readline().replace('\n', '')] = i
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
            step_numbers = [int(n) for n in hfile.readline().split()]
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
            else:
                self.limit_level=limit_level
            # Read the box geometry
            self.box_centers, self.boxes = self.read_boxes(hfile)
        # Read the cell data
        self.cells = self.read_cell_headers()
            
    def read_boxes(self, hfile):
        """
        Read the AMR boxes geometry in the base header file
        """
        # dicts to store box bounds and centers
        points = {}
        boxes = {}
        self.npoints = {}
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
            level_key = f"Lv_{current_level}"
            self.npoints[level_key] = n_cells
            points[level_key] = []
            boxes[level_key] = []
            for i in range(n_cells):
                point = []
                box = []
                for i in range(self.ndims):
                    lo, hi = [float(n) for n in hfile.readline().split()]
                    box.append([lo, hi])
                    point.append(lo + (hi - lo)/2)
                points[level_key].append(point)
                boxes[level_key].append(box)
            self.cell_paths.append(hfile.readline().replace('\n', ''))
        for ky in points:
            points[ky] = np.array(points[ky])
            boxes[ky] = np.array(boxes[ky])
        return points, boxes

    def read_cell_headers(self):
        """
        Read the cell header data for a given level
        """
        cells = {}
        for i in range(self.limit_level + 1):
            cells[f"Lv_{i}"] = {}
            cfile_path = os.path.join(self.pfile, self.cell_paths[i] + "_H")
            with open(cfile_path) as cfile:
                # Skip 2 lines
                cfile.readline()
                cfile.readline()
                # Are we good
                assert int(cfile.readline()) == len(self.fields)
                cfile.readline()
                n_cells = int(cfile.readline().split()[0].replace('(', ''))
                indexes = []
                for _ in range(n_cells):
                    start, stop, _ = cfile.readline().split()
                    start = np.array(start.replace('(', '').replace(')', '').split(','), dtype=int)
                    stop = np.array(stop.replace('(', '').replace(')', '').split(','), dtype=int)
                    indexes.append([start, stop])
                cells[f"Lv_{i}"]["indexes"] = indexes
                cfile.readline()
                assert n_cells == int(cfile.readline())
                files = []
                offsets = []
                for _ in range(n_cells):
                    _, file, offset = cfile.readline().split()
                    files.append(os.path.join(self.pfile, self.cell_paths[i].replace('Cell', ''), file))
                    offsets.append(int(offset))
                cells[f"Lv_{i}"]["files"] = files
                cells[f"Lv_{i}"]["offsets"] = offsets
        return cells


    def field_index(self, field):
        """ return the index of a data field """
        for i, f in enumerate(self.fields):
            if f == field:
                return i
        raise ValueError(f"""Field {field} was not found in file. 
                             Available fields in {self.pfile.split('/')[-1]} are:
                             {', '.join(self.fields.keys())} and grid_level""")

    
    def make_dir_tree(self, outpath, limit_level=None):
        """
        Re-Create the tree structure of the plotfile in :outpath:
        """
        os.makedirs(outpath, exist_ok=True)
        #shutil.copy(os.path.join(self.pfile, 'Header'),
        #           outpath)
        for pth in self.cell_paths:
            level_dir = pth.split('/')[0]
            os.makedirs(os.path.join(outpath, level_dir), exist_ok=True)
            #shutil.copy(os.path.join(self.pfile, pth + '_H'),
            #            os.path.join(outpath, level_dir))
            
