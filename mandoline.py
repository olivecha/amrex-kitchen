import os
import sys
import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from mandoline_utils import HeaderData
matplotlib.use('Agg')

coords_dict = {0:'x',
               1:'y',
               2:'z'}

plotfile = sys.argv[1]
field = sys.argv[2]
coord = int(sys.argv[3])
pos = float(sys.argv[4])
max_level = int(sys.argv[5])

header_pth = os.path.join(plotfile, 'Header')
print("Reading header file...")
now = time.time()
hdr = HeaderData(header_pth, max_level=max_level)
print(f"Took {time.time() - now} s") 
datasets = {}

print("Computing which boxes are in the slice...")
now = time.time()
slice_idx = hdr.slice_indexes(coord, pos, max_level)
print(f"Took {time.time() - now} s")

all_points = np.vstack([hdr.points[ky][slice_idx[ky]] for ky in slice_idx])

field_data = []
print(f"Reading {field} data...")
now = time.time()
for i in range(max_level + 1):
    lv_ky = f"Lv_{i}"
    print(f"Level {i}")
    cell_pth = os.path.join(plotfile, f'Level_{i}', 'state_H')
    ds = hdr.load_field_from_indexes(cell_pth, field, slice_idx[lv_ky], i)
    field_data.append(ds)
print(f"Took {time.time() - now} s")
print(f"Plotting the slice...")

all_values = np.hstack(field_data)
    
x_coord, y_coord = [i for i in range(3) if i != coord]

fig = plt.figure(figsize=(8, 6))

plt.tricontourf(all_points.T[x_coord],
                all_points.T[y_coord],
                all_values,
                100)

plt.colorbar(label='field')
ax = plt.gca()
ax.set_xlabel(f"{coords_dict[x_coord]} [m]")
ax.set_ylabel(f"{coords_dict[y_coord]} [m]")

fig.savefig(f"{plotfile}_{coords_dict[coord]}_{field}", dpi=300)
plt.close(fig)
print("Done!")
