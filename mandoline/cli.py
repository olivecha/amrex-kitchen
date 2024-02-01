import os
import sys
import time
import argparse
import multiprocessing
import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from mandoline import HeaderData
matplotlib.use('Agg')

coords_dict = {0:'x',
               1:'y',
               2:'z'}

def scale_array(arr, factor):
    """
    Expand lower resolution array by [factor]
    to broadcast it to a higher level grid
    """
    exp = np.repeat(arr, factor).reshape(arr.shape[0], arr.shape[1]*factor)
    exp = np.repeat(exp, factor, axis=0).reshape(arr.shape[0]*factor, arr.shape[1]*factor)
    return exp
        

def read_box_parallel(args):
    """
    box : coordinates of the box [[x_lo, x_hi],[y_lo, y_hi],[z_lo, z_hi]]
    cfile : Binary file containing the cell data
    idxs : indexes of the cell in the level grid
    offset : offset in bytes of the cell in the binary file
    Lv : Current AMR Level in the iteration process
    """
    """
    Multiprocessing compatible function to read data for one cell
    args = (idx, Lv, box)
    idx: (int)            index of the cell to read
    Lv:  (int)            AMR Level of the cell
    box: (2 x ndim array) Bounds of the box to read
    Read the data from an AMR box and return the slices
    Using a defined slice coordinate (pos):
    output: [[[x_slice_start, x_slice_end], # Left side
              [y_slice_start, y_slice_end], # Left side
              data_array, # shape = box_shape[x, y] * (Lv - max_lvl)**2
              normal,]  # Normal coordinate of the read data
             [[x_slice_start, x_slice_end], # Right side
              [y_slice_start, y_slice_end], # Right side
              data_array, # shape = box_shape[x, y] * (Lv - max_lvl)**2
              normal,]]  # Normal coordinate of the read data
    output[0] = None if there is no data at the left
    output[1] = None if there is no data at the right
    (Multiprocessing functions I/O are cursed)
    """
    idx, Lv, box = args
    factor = 2**(limit_level - Lv)
     # Get the box indexes and compute shape
    indexes = hdr.cells[f"Lv_{Lv}"]["indexes"][idx]
    cfile = hdr.cells[f"Lv_{Lv}"]["files"][idx]
    offset = hdr.cells[f"Lv_{Lv}"]["offsets"][idx]
    # Compute the slice indexes for the global grid
    x_start = indexes[0][cx] * factor
    x_stop = (indexes[1][cx] + 1) * factor
    y_start = indexes[0][cy] * factor
    y_stop = (indexes[1][cy] + 1) * factor
    shape = (indexes[1][0] - indexes[0][0] + 1,
             indexes[1][1] - indexes[0][1] + 1,
             indexes[1][2] - indexes[0][2] + 1)
    # Shape in the limit level grid
    grid_shape = np.array(shape)*factor

    # Compute the grid in the normal direction
    normal_grid = np.linspace(box[coord][0] + hdr.dx[Lv][coord]/2,
                              box[coord][1] - hdr.dx[Lv][coord]/2,
                              shape[coord])
    
    # Read the field in the binary file
    byte_size = np.product(shape)
    try:
        with open(cfile, "rb") as f:
            # Go to the cell start
            f.seek(offset)
            # Skip the header
            f.readline()
            # Go to the field we want
            f.seek(byte_size*8*fidx, 1)
            arr = np.fromfile(f, "float64", byte_size,)
            arr = arr.reshape(shape, order="F")
    except FileNotFoundError:
        #print(f'Could not find file {cfile}')
        arr = Lv * np.ones(shape)

    # Slice indexes for the data array
    arr_indices = [0, 0, 0]
    arr_indices[cx] = slice(0, shape[cx])
    arr_indices[cy] = slice(0, shape[cy])
    # Create output depending on slice position
    output = [None, None]
    # Case when the plane is between the last point and the left edge
    if pos > normal_grid[shape[coord] - 1]:
        # Get the slice in the field data box array
        arr_indices[coord] = shape[coord] - 1
        arr_data = arr[tuple(arr_indices)]
        # Only write to left interpolation plane
        output[0] = [[x_start, x_stop], [y_start, y_stop], # Slice in limit_level grid
                     scale_array(arr_data, factor), # Expanded data
                     normal_grid[shape[coord] - 1]] # normal position for interpolation
    
    # Case when the plane is between the first point and the right edge
    elif pos < normal_grid[0]:
        arr_indices[coord] = 0
        arr_data = arr[tuple(arr_indices)]
        # Only write to right interpolation plane
        output[1] = [[x_start, x_stop], [y_start, y_stop],  # Slice in limit_level grid
                     scale_array(arr_data, factor),  # Expanded data
                     normal_grid[0]]  # normal position for interpolation
        
    # Case when the plane lands on a grid point
    elif np.isclose(pos, normal_grid).any():
        match_idx = np.where(np.isclose(pos, normal_grid))[0][0]
        arr_indices[coord] = match_idx
        arr_data = arr[tuple(arr_indices)]
        # Put left and right interpolation on the same point
        output[0] = [[x_start, x_stop], [y_start, y_stop],  # Slice in limit_level grid
                     scale_array(arr_data, factor),  # Expanded data
                     normal_grid[match_idx]]  # normal position for interpolation
        output[1] = [[x_start, x_stop], [y_start, y_stop],  # Slice in limit_level grid
                     scale_array(arr_data, factor),  # Expanded data
                     normal_grid[match_idx]]  # normal position for interpolation

    # Case when the plane is between grid points in the box
    else:
        idx_left = np.where(pos > normal_grid)[0][-1]
        idx_right = np.where(pos < normal_grid)[0][0]
        # Each point to its interpolation plane
        arr_indices[coord] = idx_left
        arr_data = arr[tuple(arr_indices)]
        output[0] = [[x_start, x_stop], [y_start, y_stop],  # Slice in limit_level grid
                     scale_array(arr_data, factor).copy(),  # Expanded data
                     normal_grid[idx_left]]  # normal position for interpolation
        # Right plane
        arr_indices[coord] = idx_right
        arr_data = arr[tuple(arr_indices)]
        output[1] = [[x_start, x_stop], [y_start, y_stop],  # Slice in limit_level grid
                     scale_array(arr_data, factor).copy(),  # Expanded data
                     normal_grid[idx_right]]  # normal position for interpolation

    return output 


def main():
    """
    Main function running the mandoline Command Line Tool
    """
    # Must define this inside of main
    """
    Argument parser
    """

    parser = argparse.ArgumentParser(
            description="Fast approximate slices of AMReX plotfiles")

    parser.add_argument(
            "plotfile", type=str,
            help="Path of the plotfile to slice")
    parser.add_argument(
            "--normal", "-n", type=int,
            help="Coordinate normal to the slice x:0, y:1, z:2")
    parser.add_argument(
            "--position", "-p", type=float,
            help="position of the slice, defaults to mid plane")
    parser.add_argument(
            "--variable", "-v", type=str,
            help="variable name, defaults to \"density\"")
    parser.add_argument(
            "--max_level", "-L", type=int,
            help="Maximum AMR level loaded, defaults to finest level")
    parser.add_argument(
            "--format", "-f", type=str,
            help=("Either image or array "
                  "image: creates and saves an image using matplotlib "
                  "array: creates a numpy array with a 500x500 uniform grid "
                  "and saves it as a pickle file of a python dict"))
    parser.add_argument(
            "--output", "-o", type=str,
            help="File name used to override the default value")
    parser.add_argument(
            "--colormap", "-c", type=str,
            help="A named matplotlib colormap, defaults to jet")
    parser.add_argument(
            "--minimum", "-m", type=float,
            help="Minimum value used in the colormap")
    parser.add_argument(
            "--maximum", "-M", type=float,
            help="Maximum value used in the colormap")
    parser.add_argument(
            "--log", "-l", action='store_true',
            help="Flag to use log scale in the plotted image")
    parser.add_argument(
            "--verbose", "-V", type=int,
            help="Verbosity level, defaults to 1")
    
    args = parser.parse_args()

    """
    Define default arguments
    """

    # Default verbosity
    if args.verbose is None:
        verbose = 0
    else:
        verbose = args.verbose

    # Default output file
    if args.output is None:
        outfile = f"{args.plotfile}_{coords_dict[args.normal]}_{args.variable}"
    else:
        outfile = args.output

    # Default variable is density as in AMReX
    if args.variable is None:
        variable = 'density'
    else:
        variable = args.variable

    # Default to slice in the x direction
    global coord
    if args.normal is None:
        coord = 0
    else:
        coord = args.normal

    """
    Read the header files and define field and coords
    """
    # Header timer 
    header_start = time.time()

    # Read the plotfile headers (defaults to max_level)
    global hdr
    hdr = HeaderData(args.plotfile, limit_level=args.max_level)
    
    # Define the limit_level using the header data
    global limit_level
    if args.max_level is None:
        limit_level = hdr.max_level
    else:
        limit_level = args.max_level
    print(limit_level)

    # Slice position defaults to domain center
    global pos
    if args.position is None:
        pos = (hdr.geo_high[coord] - hdr.geo_low[coord])/2
    else:
        pos = args.position

    # Find the field index
    global fidx
    fidx = hdr.field_index(variable)

    # Define x and y coordinates in the slice plane
    global cx
    global cy
    cx, cy = [i for i in range(3) if i != args.normal]

    # Timing info
    if verbose > 0:
        print("Time to read header files:",
              np.around(time.time() - header_start, 2))

    """
    Reading the boxes intersecting with the plane
    """
    # Object to store the slices
    plane_data = {}
    # For a given level
    for Lv in range(limit_level + 1):
        # Box reading timer
        read_start = time.time()
        # Divide by level so we can stack later
        plane_data[Lv] = []
        # Multiprocessing inputs
        pool_inputs = []
        # For each box in that level
        for idx, box in enumerate(hdr.boxes[f"Lv_{Lv}"]):
            # Check if the box intersect the slicing plane
            if (box[coord][0] <= pos and 
                box[coord][1] >= pos):
                # Here, intersecting boxes at each level
                # Are added to the slice, even if there are
                # Higher level boxes covering the lower level
                # Box region in the plane. This is a little
                # wastefull as more data is readed than needed
                # to create the plane, but finding if the box is
                # Completely covered is a nightmare (and slower).
                pool_inputs.append((idx, Lv, box)) # Add to inputs

        # Read the data in parallel
        with multiprocessing.Pool() as pool:
            plane_data[Lv] = pool.map(read_box_parallel, pool_inputs)
        if verbose > 0:
            print(f"Time to read Lv {Lv}:", 
                  np.around(time.time() - read_start, 2))
    """
    Stacking and Interpolation
    """
    # Interpolation timer
    interp_start = time.time()
    # The box reader returns the slice at both sides of the slicing
    # Plane if available (i.e. not a boundary, or between boxes). 
    # This allows handling the case when the slice is between boxes
    # Array for the "left" side of the plane (could be down whatever)
    left = {'data':np.empty(hdr.grid_sizes[limit_level][[cx, cy]]),
            'normal':np.empty(hdr.grid_sizes[limit_level][[cx, cy]])}
    # Array for the "right" or up side of the plane
    # The only convention is that "left" < slice_coordinate < "right"
    right = {'data':np.empty(hdr.grid_sizes[limit_level][[cx, cy]]),
             'normal':np.empty(hdr.grid_sizes[limit_level][[cx, cy]])}

    # Parse the multiprocessing output
    # Do levels sequentially to update with finer data
    for Lv in plane_data:
        for output in plane_data[Lv]:
            # Add the slices if they are defined
            if output[0] is not None:
                out = output[0]
                # x slice
                xa, xo = out[0][0], out[0][1]
                # y slice
                ya, yo = out[1][0], out[1][1]
                # add the field data
                left['data'][xa:xo, ya:yo] = out[2]
                # add the normal coordinate
                left['normal'][xa:xo, ya:yo] = out[3]
            # Same for the right side
            if output[1] is not None:
                out = output[1]
                xa, xo = out[0][0], out[0][1]
                ya, yo = out[1][0], out[1][1]
                right['data'][xa:xo, ya:yo] = out[2]
                right['normal'][xa:xo, ya:yo] =  out[3]

    # Do the linear interpolation if normals are not the same
    # Empty array for the final data
    data = np.empty(hdr.grid_sizes[limit_level][[cx, cy]])
    # This handles the case when the slice is on a point
    bint = ~np.isclose(left['normal'], right['normal'])
    # Linear interpolation
    itrp = (left['data'][bint] * (right['normal'][bint] - pos) + right['data'][bint] * (pos - left['normal'][bint]))
    # Divide by the difference in normal corrdinates
    data[bint] =  itrp/(right['normal'][bint] - left['normal'][bint])
    # Could be either
    data[~bint] = right['data'][~bint]
    # For some reason
    data = data.T

    if verbose > 0:
        print("Interpolation time:", 
              np.around(time.time() - interp_start, 2), "s")

    """
    Output the slice
    """
    output_start = time.time()
    # Define the limit level coodinates vectors
    x_grid = np.linspace(hdr.geo_low[cx]  + hdr.dx[limit_level][cx]/2,
                         hdr.geo_high[cx] - hdr.dx[limit_level][cx]/2,
                         hdr.grid_sizes[limit_level][cx])

    y_grid = np.linspace(hdr.geo_low[cy]  + hdr.dx[limit_level][cy]/2,
                         hdr.geo_high[cy] - hdr.dx[limit_level][cy]/2,
                         hdr.grid_sizes[limit_level][cy])


    # Array output
    if args.format == "array":
        # Interpolation bounds

        # Case when we slice all the values
        if variable == 'all':
            pass
            # To be added

        # Make a dict with output
        output = {'x': x_grid,
                  'y': y_grid,
                  variable: data}

        # Pickle into the jar
        pfile = open(outfile + ".pkl", 'wb')
        pickle.dump(output, pfile)
        if verbose > 0:
            print("Time to save uniform grid:", 
                  np.around(time.time() - output_start, 2), "s")

    # Image output
    else:
        # Pretty large figure
        fig = plt.figure(figsize=(8, 6))

        # Special case for grid_level 
        # (use the same number of colors as the number of levels)
        if args.variable == 'grid_level':
            raise ValueError('TODO')

        # Use log scale ?
        if args.log:
            norm = matplotlib.colors.LogNorm()
        else:
            norm = None

        if args.colormap in plt.colormaps():
            cmap = args.colormap
        else:
            cmap = 'jet'
        
        # Plot the slice
        plt.pcolormesh(x_grid,
                       y_grid,
                       data,
                       vmin=args.minimum,
                       vmax=args.maximum,
                       norm=norm,
                       cmap=cmap)

        # Add a colorbar  
        plt.colorbar(label=f'{args.variable}')
        ax = plt.gca()
        ax.set_aspect('equal')
        ax.set_xlabel(f"{coords_dict[cx]} [m]")
        ax.set_ylabel(f"{coords_dict[cy]} [m]")
        
        # save and close
        fig.savefig(outfile, dpi=500)
        plt.close(fig)

if __name__ == "__main__":
    main()
