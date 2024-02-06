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
from mandoline.utils import sanitize_field_name
from mandoline.data_reader import slice_box
from mandoline.slice_data import SliceData

matplotlib.use('Agg')

coords_dict = {0:'x',
               1:'y',
               2:'z'}

def main():
    """
    Main function running the mandoline Command Line Tool
    """

    # Argument parser
    parser = argparse.ArgumentParser(
            description="Fast  slices of (large) AMReX plotfiles")

    parser.add_argument(
            "plotfile", type=str,
            help="Path of the plotfile to slice")
    parser.add_argument(
            "--normal", "-n", type=int,
            help="Index of the coordinate normal to the slice x:0, y:1, z:2")
    parser.add_argument(
            "--position", "-p", type=float,
            help="position of the slice, defaults to domain center")
    parser.add_argument(
            "--variables", "-v", type=str, nargs='+',
            help=("variables names to slice, defaults to \"density\""
                  "\"all\" slices all the fields in the plotfile"
                  "\"grid_level\" outputs the AMR grid level data"))
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
            help="Output file name used to override the default value")
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
        verbose = 1
    else:
        verbose = args.verbose

    # Default output file
    if args.output is None:
        outfile_root = f"{args.plotfile}_{coords_dict[args.normal]}"
    else:
        outfile_root = args.output

    """
    Read the header files and define field and coords
    """
    # Header timer 
    header_start = time.time()

    # Define the data needed for the slice
    # This reads the plotfile and cell headers
    # And also stores and define the slice data
    # With default values to pass it to the 
    # Box reader multiprocessing function
    slc = SliceData(args.plotfile, 
                    fields=args.variables, 
                    normal=args.normal, 
                    pos=args.position,
                    limit_level=args.max_level)

    # Timing info
    if verbose > 0:
        print("Time to read header files:",
              np.around(time.time() - header_start, 2))

    """
    Reading the boxes intersecting with the plane
    """
    # Object to store the slices
    plane_data = {}
    # Multiprocessing pool input template
    # Passing the SliceData class copies it across
    # All subprocesses which is slow
    # For a given level
    for Lv in range(slc.limit_level + 1):
        # Box reading timer
        read_start = time.time()
        # Divide by level so we can stack later
        plane_data[Lv] = []
        # Multiprocessing inputs
        pool_inputs = []
        # For each box in that level
        for idx, box in enumerate(slc.boxes[f"Lv_{Lv}"]):
            # Check if the box intersect the slicing plane
            if (box[slc.cn][0] <= slc.pos and 
                box[slc.cn][1] >= slc.pos):
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
                p_in  = {'cx':slc.cx,
                         'cy':slc.cy,
                         'cn':slc.cn,
                         'dx':slc.dx,
                         'pos':slc.pos,
                         'limit_level':slc.limit_level,
                         'fidxs':slc.fidxs,
                         'Lv':Lv,
                         'idx':idx,
                         'indexes':slc.cells[f'Lv_{Lv}']['indexes'][idx],
                         'cfile':slc.cells[f'Lv_{Lv}']['files'][idx],
                         'offset':slc.cells[f'Lv_{Lv}']['offsets'][idx],
                         'box':box}
                pool_inputs.append(p_in) # Add to inputs

        # Read the data in parallel
        with multiprocessing.Pool() as pool:
            plane_data[Lv] = pool.map(slice_box, pool_inputs)
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
    left = {'data':[slc.limit_level_arr() for _ in range(slc.nfidxs)],
            'normal':slc.limit_level_arr()}
    # Array for the "right" or up side of the plane
    # The only convention is that "left" < slice_coordinate < "right"
    right = {'data':[slc.limit_level_arr() for _ in range(slc.nfidxs)],
            'normal':slc.limit_level_arr()}

    if slc.do_grid_level:
        grid_level = {'left':slc.limit_level_arr(),
                      'right':slc.limit_level_arr()}

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
                if slc.do_grid_level:
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
                if slc.do_grid_level:
                    grid_level['right'][xa:xo, ya:yo] = out['level']

    # Do the linear interpolation if normals are not the same
    # Empty arrays for the final data
    all_data = []
    # Boolean slicing array for when the points are the same
    bint = ~np.isclose(left['normal'], right['normal'])
    # Iterate with the number of fields
    for i in range(slc.nfidxs):
        data = slc.limit_level_arr()
        # Linear interpolation
        term1 = left['data'][i][bint] * (right['normal'][bint] - slc.pos) 
        term2 = right['data'][i][bint] * (slc.pos - left['normal'][bint])
        term3 = right['normal'][bint] - left['normal'][bint]
        data[bint] =  (term1 + term2) / term3
        # Could be either
        data[~bint] = right['data'][i][~bint]
        # For some reason
        all_data.append(data.T)

    if slc.do_grid_level:
        # Concatenate both sides
        all_levels = np.stack([grid_level['right'], grid_level['left']])
        # I guess the min is better for debuging
        # All things considered they should be pretty similar
        all_data.append(np.min(all_levels, axis=0).T)

    if verbose > 0:
        print("Interpolation time:", 
              np.around(time.time() - interp_start, 2), "s")

    """
    Output the slice
    """
    output_start = time.time()
    # Define the limit level coodinates vectors
    # TODO: move this to SliceData as a method
    x_grid = np.linspace(slc.geo_low[slc.cx]  + slc.dx[slc.limit_level][slc.cx]/2,
                         slc.geo_high[slc.cx] - slc.dx[slc.limit_level][slc.cx]/2,
                         slc.grid_sizes[slc.limit_level][slc.cx])

    y_grid = np.linspace(slc.geo_low[slc.cy]  + slc.dx[slc.limit_level][slc.cy]/2,
                         slc.geo_high[slc.cy] - slc.dx[slc.limit_level][slc.cy]/2,
                         slc.grid_sizes[slc.limit_level][slc.cy])

    # Case when we slice all the values
    if 'all' in args.variables:
        field_names = [name for name in slc.fields]
    else:
        # Treat grid_level differently
        field_names = [name for name in args.variables if name != 'grid_level']
        

    # Array output
    if args.format == "array":
        # Make a dict with output
        output = {'x': x_grid,
                  'y': y_grid,}
        # Store in output
        for i, name in enumerate(field_names):
            output[name] = all_data[i]

        # Works if only grid_level
        if slc.do_grid_level:
            output['grid_level'] = all_data[-1]

        # Add info to output if single field
        if len(args.variables) == 1:
            fname = sanitize_field_name(args.variables[0])
            outfile = '_'.join([outfile_root, fname])
        else:
            outfile = outfile_root
        # Pickle into the jar
        pfile = open(outfile + ".pkl", 'wb')
        pickle.dump(output, pfile)
        if verbose > 0:
            print("Time to save uniform grid:", 
                  np.around(time.time() - output_start, 2), "s")

    # Image output
    else:
        # Use log scale ?
        if args.log:
            norm = matplotlib.colors.LogNorm()
        else:
            norm = None

        # User supplied colormap
        if args.colormap in plt.colormaps():
            cmap = args.colormap
        else:
            cmap = 'jet'

        # Plots are the same with grid level
        if slc.do_grid_level:
            field_names.append('grid_level')
        # A figure per field
        for i, name in enumerate(field_names):

            if verbose > 1:
                plot_start = time.time()
                print(f"Plotting {name}...")
            # Pretty large figure
            fig = plt.figure(figsize=(8, 6))

            # Plot the slice
            plt.pcolormesh(x_grid,
                           y_grid,
                           all_data[i],
                           vmin=args.minimum,
                           vmax=args.maximum,
                           norm=norm,
                           cmap=cmap)

            # Add a colorbar  
            plt.colorbar(label=f'{name}')
            ax = plt.gca()
            ax.set_aspect('equal')
            ax.set_xlabel(f"{coords_dict[slc.cx]} [m]")
            ax.set_ylabel(f"{coords_dict[slc.cy]} [m]")

            if verbose > 1:
                print(f"Done! ({np.around(time.time() - plot_start, 2)} s)")

            if verbose > 1:
                save_start = time.time()
                print(f"Saving {name} plot...")
            
            # save and close
            outfile = '_'.join([outfile_root, sanitize_field_name(name)])
            fig.savefig(outfile, dpi=500)
            plt.close(fig)

            if verbose > 1:
                print(f"Done! ({np.around(time.time() - save_start, 2)} s)")
        

if __name__ == "__main__":
    main()
