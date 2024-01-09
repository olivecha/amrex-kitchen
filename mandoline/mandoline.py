import os
import linecache
import numpy as np
from tqdm import tqdm


class HeaderData(object):

    def __init__(self, plotfile, max_level=None):
        """
        Parse the header data and save as attributes
        """
        self.pfile = plotfile
        filepath = os.path.join(plotfile, 'Header')
        with open(filepath) as hfile:
            _ = hfile.readline()
            # field names
            self.nvars = int(hfile.readline())
            self.fields = []
            for _ in range(self.nvars):
                self.fields.append(hfile.readline().replace('\n', ''))
            # General data
            self.ndims = int(hfile.readline())
            self.time = float(hfile.readline())
            self.max_level = int(hfile.readline())
            self.geo_low = [float(n) for n in hfile.readline().split()]
            self.geo_high = [float(n) for n in hfile.readline().split()]
            self.factors = [int(n) for n in hfile.readline().split()]
            self.block_indexes = [b for b in hfile.readline().split()]
            self.step_numbers = [int(n) for n in hfile.readline().split()]
            # Grid resolutions
            resolutions = []
            for i in range(self.max_level + 1):
                resolutions.append([float(n) for n in hfile.readline().split()])
            self.resolutions = resolutions
            # Skip 2 lines
            hfile.readline()
            hfile.readline()
            # dicts to store box bounds and centers
            points = {}
            boxes = {}
            self.npoints = {}
            # Loop over the grid levels
            if max_level is None:
                max_level = self.max_level
            for lv in range(max_level):
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
                for i in tqdm(range(n_cells)):
                    point = []
                    box = []
                    for i in range(self.ndims):
                        lo, hi = [float(n) for n in hfile.readline().split()]
                        box.append([lo, hi])
                        point.append(lo + (hi - lo)/2)
                    points[level_key].append(point)
                    boxes[level_key].append(box)
                print("Loaded", hfile.readline().split('/')[0].replace('_', ' '))
            for ky in points:
                points[ky] = np.array(points[ky])
                boxes[ky] = np.array(boxes[ky])
            self.points = points
            self.boxes = boxes

    def load_cell_header(self, filepath, max_level):
        """
        Read the cell header data for a given level
        """
        with open(filepath) as cfile:
            level_string = f"Lv_{level}"
            target = f"{self.npoints[level_string]},{self.nvars}\n"
            while True:
                line = cfile.readline()
                if line == target:
                    break
            raw_data = []
            for i in tqdm(range(self.npoints[level_string])):
                raw_data.append([float(n) for n in cfile.readline().split(',')[:-1]])
            return np.array(raw_data)

    def field_index(self, field):
        """ return the index of a data field """
        for i, f in enumerate(self.fields):
            if f == field:
                return i
        raise ValueError("Field not found")

    def slice_indexes(self, coord, value=None, max_level=None):
        """
        find boxes indexes in a slice
        """
        print("Combuting box indexes in slice...")
        if max_level is None:
            max_level = self.max_level
        if value is None:
            all_bx = np.vstack([self.boxes[ky] for ky in self.boxes])
            coord_min = np.min(all_bx[:, coord, 0])
            coord_max = np.max(all_bx[:, coord, 1])
            value = coord_min + (coord_max - coord_min) / 2
            value += 1e-6*coord_max # So we are not between boxes

        indexes = {}
        for lv in range(max_level):
            lv_key = f"Lv_{lv}"
            indexes[lv_key] = []
            for i, box in enumerate(self.boxes[lv_key]):
                if (box[coord][0] < value) and (box[coord][1] > value):
                    indexes[lv_key].append(i)

        xi, yi = [i for i in range(3) if i != coord]

        # Remove indexes where higher resolution points exist
        for lv in range(max_level-1):
            lv_key = f"Lv_{lv}"
            nlv_key = f"Lv_{lv+1}"
            boxes = self.boxes[lv_key][indexes[lv_key]]
            points = self.points[nlv_key][indexes[nlv_key]]
            new_indexes = []
            for idx, b in zip(indexes[lv_key], boxes):
                flag = True
                for p in points:
                    if (p[xi] >= b[xi][0]) and (p[xi] <= b[xi][1]) and (p[yi] >= b[yi][0]) and (p[yi] <= b[yi][1]):
                        flag = False
                        #print(f"point ({p[1]}, {p[2]}) is in box [[{b[1]}], [{b[2]}]]")
                        break
                if flag:
                    new_indexes.append(idx)
            indexes[lv_key] = new_indexes
        print("Done!")
        return indexes

    def get_start_line(self, filepath, level):
        """
        Get the line number before the cell centered data starts in Cell_H
        """
        with open(filepath) as cfile:
            level_string = f"Lv_{level}"
            target = f"{self.npoints[level_string]},{self.nvars}\n"
            counter = 0
            while True:
                line = cfile.readline()
                counter += 1
                if line == target:
                    break
        return counter

    def load_field_from_indexes(self, indexes, field):
        """
        Load only cell center data with the corresponding indexes
        """
        if field == 'grid_level':
            return self.grid_level_data(indexes)

        field_data = {}
        for lv, lv_ky in enumerate(indexes):
            filepath = os.path.join(self.pfile, f"Level_{lv}", "Cell_H")
            start = self.get_start_line(filepath, lv)
            idx = self.field_index(field)
            field_data[lv_ky] = []
            for i in tqdm(indexes[lv_ky]):
                line = linecache.getline(filepath, i + start + 1)
                value = float(line.split(',')[idx])
                field_data[lv_ky].append(value)
            print(f"Loaded {field} data at Level {lv}")
        return np.hstack([field_data[ky] for ky in field_data])

    def grid_level_data(self, indexes):
        """
        Compute the grid level from slice indexes
        """
        grid_data = []
        for i, ky in enumerate(indexes):
            grid_data.append(np.repeat(i, len(indexes[ky])))
        return np.hstack(grid_data)

    def points_from_indexes(self, indexes, fix_corners=True):
        """
        Return a dict with ndim x nindexes arrays at each level
        """
        points = {ky:self.points[ky][indexes[ky]].T for ky in indexes}
        if fix_corners:
            return self.fix_corners(points)
        else:
            return np.hstack([points[ky] for ky in points])


    def fix_corners(self, points):
        """
        Return a ndim x nindexes array from indexes dict
        Attemps to fix corners by moving the lower level points at
        the coordinates of the higher level points
        """
        # This is bad but it works
        max_level = [ky for ky in points][-1]
        X, Y, Z = np.hstack([points[ky] for ky in points])
        delta_vals = [np.max(C) - np.min(C) for C in [X, Y, Z]]
        slice_coord = np.argmin(delta_vals)
        coord_x, coord_y = [i for i in range(3) if i != slice_coord]
        print(coord_x, coord_y)
        x_min = np.min(points[max_level][coord_x])
        x_max = np.max(points[max_level][coord_x])
        y_min = np.min(points[max_level][coord_y])
        y_max = np.max(points[max_level][coord_y])
        for ky in points:
            lpts = points[ky]
            lpts[coord_x][lpts[coord_x] == np.min(lpts[coord_x])] = min(x_min, np.min(lpts[coord_x]))
            lpts[coord_x][lpts[coord_x] == np.max(lpts[coord_x])] = max(x_max, np.max(lpts[coord_x]))
            lpts[coord_y][lpts[coord_y] == np.min(lpts[coord_y])] = min(y_min, np.min(lpts[coord_y]))
            lpts[coord_y][lpts[coord_y] == np.max(lpts[coord_y])] = max(y_max, np.max(lpts[coord_y]))
            points[ky] = lpts

        return np.hstack([points[ky] for ky in points])

