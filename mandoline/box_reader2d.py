import time
now = time.time()
import os
import sys
import pickle
import linecache
import multiprocessing
import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import matplotlib.colors as colors

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
            self.block_indexes = [b for b in hfile.readline().split()]
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
            # Read the boxes
            self.box_centers, self.boxes = self.read_boxes(hfile, limit_level)
            self.cells = self.read_cell_headers(limit_level)
            
    def read_boxes(self, hfile, limit_level):
        """
        Read the AMR boxes geometry in the base header file
        """
        # dicts to store box bounds and centers
        points = {}
        boxes = {}
        self.npoints = {}
        self.cell_paths = []
        # Loop over the grid levels
        if limit_level is None:
            limit_level = self.max_level
        for lv in range(limit_level + 1):
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

    def read_cell_headers(self, limit_level):
        """
        Read the cell header data for a given level
        """
        if limit_level is None:
            limit_level = self.max_level
        cells = {}
        for i in range(limit_level + 1):
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
        raise ValueError("Field not found")

        
def read_bin_parallel(args):
    """
    Multiprocessing compatible function to read data for one cell
    args = box, cfile, idxs, offset, Lv 
    box : coordinates of the box [[x_lo, x_hi],[y_lo, y_hi],[z_lo, z_hi]]
    cfile : Binary file containing the cell data
    idxs : indexes of the cell in the level grid
    offset : offset in bytes of the cell in the binary file
    Lv : Current AMR Level in the iteration process
    """
    box, cfile, idxs, offset, Lv = args
    # Read the field data
    shape = (idxs[1][0] - idxs[0][0] + 1,
             idxs[1][1] - idxs[0][1] + 1,)

    x_grid = np.linspace(box[0][0] + hdr.dx[Lv][0],
                         box[0][1] - hdr.dx[Lv][0],
                         shape[0]) 

    y_grid = np.linspace(box[1][0] + hdr.dx[Lv][1],
                         box[1][1] - hdr.dx[Lv][1],
                         shape[1]) 
    x_points, y_points = np.meshgrid(x_grid, y_grid)
    x_points = x_points.flatten(order="C")
    y_points = y_points.flatten(order="C")

    byte_size = np.product(shape)
    with open(cfile, "rb") as f:
        f.seek(offset)
        f.readline()
        f.seek(byte_size*8*fidx, 1)
        arr = np.fromfile(f, "float64", byte_size)
    return x_points, y_points, arr

field = "temp"
limit_level = None
print("Problem setup time:", np.around(time.time() - now, 2), "s")

now = time.time()
filename = sys.argv[1]
hdr = HeaderData(filename, limit_level=limit_level)
limit_level = hdr.max_level
print("Header reading time:", np.around(time.time() - now, 2), "s")
now = time.time()

"""
Remove indexes where higher resolution points exist
"""
# Current not normal indexes
cx, cy = [0, 1]
            
# Store upper level boxes that are in the intersecting box
intersect_boxes = {}
contained_boxes = {}
match_indexes = {}
for Lv in range(limit_level):
    # All boxes in the plane at Lv
    boxes = hdr.boxes[f"Lv_{Lv}"]
    # All boxes in the plane at Lv + 1
    boxes_upper = hdr.boxes[f"Lv_{Lv+1}"]
    # Indexes of the not completly intersecting boxes at Lv
    new_indexes = []
    # Current partially intersecting boxes
    intersect_boxes[Lv] = []
    contained_boxes[Lv] = []
    for idx, b in zip(range(len(boxes)), boxes):
        box_counter = 0
        flagged_boxes = []
        for bu in boxes_upper:
            # Smaller box is within the larger box
            if (bu[0][0] >= b[0][0] and
                bu[0][1] <= b[0][1] and
                bu[1][0] >= b[1][0] and
                bu[1][1] <= b[1][1]):
                box_counter += 1
                flagged_boxes.append(bu)
        # If there is less than ref_factor*ref_factor we keep the box
        if box_counter < hdr.factors[Lv]**2:
            new_indexes.append(idx)
            if box_counter > 0:
                # Intersecting lower level box
                intersect_boxes[Lv].append(b)
                # Upper level boxes in the intersect
                contained_boxes[Lv].append(flagged_boxes)
    match_indexes[Lv] = new_indexes
match_indexes[limit_level] = np.arange(len(hdr.boxes[f"Lv_{limit_level}"]))

match_cells = {}
for Lv in range(limit_level + 1):
    match_cells[Lv] = {}
    match_cells[Lv]["files"] = np.array(hdr.cells[f"Lv_{Lv}"]["files"])[match_indexes[Lv]]
    match_cells[Lv]["indexes"] = np.array(hdr.cells[f"Lv_{Lv}"]["indexes"])[match_indexes[Lv]]
    match_cells[Lv]["offsets"] = np.array(hdr.cells[f"Lv_{Lv}"]["offsets"])[match_indexes[Lv]]

# Compute stuff from input data
fidx = hdr.field_index(field)
# This should be moved to HeaderData
grid_sizes = [1 + np.array(b.replace('(', '').replace(")", '').split(','), dtype=int) for b in hdr.block_indexes[1::3]]

"""
Read the binary cell files into a dict of the plane data
"""
plane_points = {}
for Lv in range(limit_level + 1):
    now = time.time()
    plane_points[Lv] = {'x':[],
                        'y':[],
                        'v':[],}

    # One function call per box
    call_args = []
    for  box, cfile, idxs, offset in zip(hdr.boxes[f"Lv_{Lv}"][match_indexes[Lv]], 
                                         match_cells[Lv]["files"],
                                         match_cells[Lv]["indexes"],
                                         match_cells[Lv]["offsets"]):
        arg = (box, cfile, idxs, offset, Lv)
        call_args.append(arg)
    # Read the cell files and interpolate in parallel
    with multiprocessing.Pool() as pool:
        out = pool.map(read_bin_parallel, call_args)
    # Restructure the output
    for tup in out:
        for i, ky in enumerate(plane_points[Lv]):
            plane_points[Lv][ky].append(tup[i])
    # Concatenate the cells together
    for ky in plane_points[Lv]:
        plane_points[Lv][ky] = np.hstack(plane_points[Lv][ky])
    
    print("Time to read Lv", Lv,":", np.around(time.time() - now, 2), "s")

"""
Remove lower level points intersecting with upper level boxes
"""
now = time.time()
for Lv in range(limit_level):
    points_mask = np.zeros_like(plane_points[Lv]["x"], dtype=bool)
    for i, bi in enumerate(intersect_boxes[Lv]):
        # Points within lower level boxes in x direction
        inlower_x = np.logical_and(plane_points[Lv]["x"] > bi[0][0],
                                   plane_points[Lv]["x"] < bi[0][1])
        # Points within lower level boxes in y direction
        inlower_y = np.logical_and(plane_points[Lv]["y"] > bi[1][0],
                                   plane_points[Lv]["y"] < bi[1][1])
        # Points within flagged lower level boxes
        inlower = np.logical_and(inlower_x,
                                 inlower_y)
        # Bool array for the current lower level box
        inupper_all = np.zeros_like(inlower)
        # For each upper level box in the current lower level box
        for bu_i in contained_boxes[Lv][i]:
            # Contained in x
            inupper_x = np.logical_and(plane_points[Lv]["x"] > bu_i[0][0],
                                       plane_points[Lv]["x"] < bu_i[0][1])
            # Contained in y
            inupper_y = np.logical_and(plane_points[Lv]["y"] > bu_i[1][0],
                                       plane_points[Lv]["y"] < bu_i[1][1])
            # Contained in the upper level box
            inupper = np.logical_and(inupper_x, inupper_y)
            # Update the upper level bool array
            inupper_all = np.logical_or(inupper_all, inupper)
        # Mask points that are in both lower and upper level inter-contained boxes
        box_mask =  np.logical_and(inlower, inupper_all)
        # Update the global bool array
        points_mask = np.logical_or(points_mask, box_mask)
    # Keep only the non intersecing points
    for ky in plane_points[Lv]:
        plane_points[Lv][ky] = plane_points[Lv][ky][~points_mask] 
print("Time to remove redundant points:", np.around(time.time() - now, 2), "s")

#with open("grid_test.pkl", 'wb') as pfile:
#    pickle.dump(plane_points, pfile)


"""
Create a covering grid, this is equivalent as plotting the triangulation
but we can control the resolution
"""
now = time.time()
# Remove level from plane points
data = {}
for ky in ["x", "y", "v"]:
    data[ky] = np.hstack([plane_points[Lv][ky] for Lv in range(limit_level +1)])

vec_x = np.linspace(np.min(data["x"]),
                    np.max(data["x"]),
                    1000)

vec_y = np.linspace(np.min(data["y"]),
                    np.max(data["y"]),
                    1000)

plt_gridx, plt_gridy = np.meshgrid(vec_x, vec_y)

data_out = griddata([(xi, yi) for xi, yi in zip(data["x"], data["y"])],
                    data["v"],
                    [(xi, yi) for xi, yi in zip(plt_gridx.flatten(), plt_gridy.flatten())])
print("Covering grid time:", np.around(time.time() - now, 2), "s")

now = time.time()
VRT = data_out.reshape(1000, 1000)
VRT[np.isnan(VRT)] = 1
VRT[VRT <= 1] = 1
fig, ax = plt.subplots(figsize=(8, 8))
plt.sca(ax)
plt.pcolormesh(vec_x,
               vec_y,
               VRT,
               vmax=45000,
               norm=colors.LogNorm(),
               cmap="magma")
ax.set_aspect("equal")
print("Plotting time:", np.around(time.time() - now), "s")

now = time.time()
plt.gcf().savefig(f"slice_x_{filename}.png", dpi=500)
print("Saving figure time:", np.around(time.time() - now, 2), "s")

