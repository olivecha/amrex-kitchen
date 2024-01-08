import numpy as np
import linecache


class HeaderData(object):
    
    def __init__(self, filepath, max_level=0):
        """
        Parse the header data and save as attributes
        """
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
            # Read the cell centers
            current_level, n_cells, _ = [n for n in hfile.readline().split()]
            self.step = hfile.readline()
            current_level = int(current_level)
            level_key = f"Lv_{current_level}"
            n_cells = int(n_cells)
            points = {}
            boxes = {}
            self.npoints = {f'Lv_{current_level}':n_cells}
            while True:
                points[level_key] = []
                boxes[level_key] = []
                while True:
                    try:
                        point = []
                        box = []
                        for i in range(self.ndims):
                            lo, hi = [float(n) for n in hfile.readline().split()]
                            box.append([lo, hi])
                            point.append(lo + (hi - lo)/2)
                        points[level_key].append(point)
                        boxes[level_key].append(box)
                    except:
                        break
                try:
                    current_level, n_cells, _ = [n for n in hfile.readline().split()]
                    current_level = int(current_level)
                    level_key = f"Lv_{current_level}"
                    n_cells = int(n_cells)
                    self.npoints[f'Lv_{current_level}'] = n_cells
                    hfile.readline()
                    if current_level > max_level:
                        break
                except:
                    break
            for ky in points:
                points[ky] = np.array(points[ky])
                boxes[ky] = np.array(boxes[ky])
            self.points = points
            self.boxes = boxes
            
    def load_cell_header(self, filepath, level):
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
            for i in range(self.npoints[level_string]):
                raw_data.append([float(n) for n in cfile.readline().split(',')[:-1]])
            return np.array(raw_data)
        
    def field_index(self, field):
        """ return the index of a data field """
        for i, f in enumerate(self.fields):
            if f == field:
                return i
        raise ValueError("Field not found")
        
    def slice_indexes(self, coord, value, max_level=0):
        """
        find boxes indexes in a slice
        """
        indexes = {}
        for lv in range(max_level +1):
            lv_key = f"Lv_{lv}"
            indexes[lv_key] = []
            for i, box in enumerate(self.boxes[lv_key]):
                if (box[coord][0] < value) and (box[coord][1] > value):
                    indexes[lv_key].append(i)
                    
        xi, yi = [i for i in range(3) if i != coord]
        
        # Remove indexes where higher resolution points exist
        for lv in range(max_level):
            print(f"Level {lv}")
            lv_key = f"Lv_{lv}"
            nlv_key = f"Lv_{lv+1}"
            boxes = self.boxes[lv_key][indexes[lv_key]]
            points = self.points[nlv_key][indexes[nlv_key]]
            new_indexes = []
            for idx, b in zip(indexes[lv_key], boxes):
                flag = True
                for p in points:
                    if (p[xi] > b[xi][0]) and (p[xi] < b[xi][1]) and (p[yi] > b[yi][0]) and (p[yi] < b[yi][1]):
                        flag = False
                        break
                if flag:
                    new_indexes.append(idx)
            indexes[lv_key] = new_indexes
        print(f"Level {max_level}")
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
    
    def load_field_from_indexes(self, filepath, field, indexes, level):
        """
        Load only cell center data with the corresponding indexes
        """
        start = self.get_start_line(filepath, level)
        idx = self.field_index(field)
        raw_data = []
        for i in indexes:
            line = linecache.getline(filepath, i + start + 1)
            value = float(line.split(',')[idx])
            raw_data.append(value)
        return raw_data
    
