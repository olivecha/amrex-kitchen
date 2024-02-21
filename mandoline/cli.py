import os
import sys
import time
import shutil
import argparse
import multiprocessing
import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
from mandoline import HeaderData
from mandoline.utils import sanitize_field_name, plotfile_ndims
from mandoline.data_reader import slice_box, plate_box
from mandoline.slice_data import SliceData


coords_dict = {0:'x',
               1:'y',
               2:'z',
               '2D':''}

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
                  ", \"all\" slices all the fields in the plotfile, "
                  "\"grid_level\" outputs the AMR grid level data"))
    parser.add_argument(
            "--max_level", "-L", type=int,
            help="Maximum AMR level loaded, defaults to finest level")
    parser.add_argument(
            "--serial", "-s", action='store_true',
            help="Flag to disable multiprocessing")
    parser.add_argument(
            "--format", "-f", type=str,
            help=("""Either "image", "array" or "plotfile".
                  image: creates and saves an image using matplotlib.
                  array: creates a numpy array with a uniform grid
                  with the resolution of --max_level and saves it
                  as a numpy file with separated fields.
                  plotfile: Keep the adaptive mesh refinement information
                  and save the slice and specified fields in a 2D amrex
                  plotfile"""))
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


    """
    Read the header files and define field and coords
    """
    # Header timer 
    header_start = time.time()

    # Class to handle slice parameters
    slc = SliceData(args.plotfile, 
                    fields=args.variables, 
                    normal=args.normal, 
                    pos=args.position,
                    limit_level=args.max_level)

    # Default output file
    if args.output is None:
        outbase = args.plotfile.split('/')[-1]
        if slc.ndims == 3:
            outfile_root = f"{coords_dict[slc.cn]}_{outbase}"
        else:
            outfile_root = f"2D_{outbase}"
    else:
        outfile_root = args.output

    # Define the workers pool if needed
    if args.serial:
        pass
    else:
        pool = multiprocessing.Pool()

    """
    Case for 2D plotfiles (no slicing required)
    """
    if slc.ndims == 2:
        #TODO: Use colander to remove fields from 2D plotfiles
        if args.format == "plotfile":
            raise NotImplementedError("TODO (you can use amrex-kitchen/"
                                      "colander")
        # The slice is just the header data
        # Object to store the slices
        plane_data = {}
        # For a given level
        for Lv in range(slc.limit_level + 1):
            # Box reading timer
            read_start = time.time()
            # Divide by level so we can stack later
            plane_data[Lv] = []
            # Multiprocessing inputs
            pool_inputs = []
            for indexes, cfile, offset, box in zip(slc.cells[f"Lv_{Lv}"]['indexes'],
                                                   slc.cells[f"Lv_{Lv}"]['files'],
                                                   slc.cells[f"Lv_{Lv}"]['offsets'],
                                                   slc.boxes[f"Lv_{Lv}"]):
                # Everything needed by the slice reader
                p_in  = {'cx':slc.cx,
                         'cy':slc.cy,
                         'dx':slc.dx,
                         'limit_level':slc.limit_level,
                         'fidxs':slc.fidxs,
                         'Lv':Lv,
                         'indexes':indexes,
                         'cfile':cfile,
                         'offset':offset,
                         'box':box}
                pool_inputs.append(p_in) # Add to inputs
            # Read the data in parallel
            if args.serial:
                plane_data[Lv] = list(map(plate_box, pool_inputs))
            else:
                plane_data[Lv] = pool.map(plate_box, pool_inputs)

            if verbose > 0:
                print(f"Time to read Lv {Lv}:", 
                      np.around(time.time() - read_start, 2))

        # Broadcast the data to empty arrays
        all_data = [slc.limit_level_arr() for _ in range(slc.nfidxs)]

        if slc.do_grid_level:
            grid_level = slc.limit_level_arr()
        else:
            grid_level = None

        # Parse the multiprocessing output
        # Do levels sequentially to update with finer data
        for Lv in plane_data:
            for out in plane_data[Lv]:
                # Add the slices if they are defined
                xa, xo = out['sx']  # x slice
                ya, yo = out['sy']  # y slice
                # add the field data
                for i in range(slc.nfidxs):
                    all_data[i][xa:xo, ya:yo] = out['data'][i]
                # add the normal coordinate
                # broadcast the grid level to the grid if needed
                if slc.do_grid_level:
                    grid_level[xa:xo, ya:yo] = out['level']

        if slc.do_grid_level:
            all_data.append(grid_level)

        for i, data in enumerate(all_data):
            all_data[i] = data.T
        

    # Case for 3D plotfiles (mandoline time)
    elif slc.ndims == 3:

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
                             'bidx':idx,
                             'indexes':slc.cells[f'Lv_{Lv}']['indexes'][idx],
                             'cfile':slc.cells[f'Lv_{Lv}']['files'][idx],
                             'offset':slc.cells[f'Lv_{Lv}']['offsets'][idx],
                             'box':box}
                    pool_inputs.append(p_in) # Add to inputs

            # Read the data in parallel (or serial)
            # The box reader returns the slice at both sides of the slicing
            # Plane if available (i.e. not a boundary, or between boxes). 
            # This allows handling the case when the slice is between boxes
            # This function interpolate between the two planes and returns
            if args.serial:
                plane_data[Lv] = list(map(slice_box, pool_inputs))
            else:
                plane_data[Lv] = pool.map(slice_box, pool_inputs)

            if verbose > 0:
                print(f"Time to read Lv {Lv}:", 
                      np.around(time.time() - read_start, 2))
        """
        Stacking and Interpolation
        """
        # Interpolation timer
        interp_start = time.time()
        # A list of the interpolated data arrays
        if (args.format == "array" or
            args.format == "image"):
            all_data = slc.reducemp_data_ortho(plane_data)

        elif args.format == "plotfile":
            all_data_bylevel, indexes, headers  = slc.interpolate_bylevel(plane_data)

        if verbose > 0:
            print("Interpolation time:", 
                  np.around(time.time() - interp_start, 2), "s")
    else:
        raise ValueError(f"Number of dimensions {ndims} not supported")

    """
    Output the slice
    """
    output_start = time.time()

    # Define the saved fields
    all_names = [name for name in slc.fields]
    field_names = [all_names[idx] for idx in slc.fidxs if idx is not None]
    
    if args.format == "plotfile":
        # Define the plotfile name
        if len(field_names) == 1:
            fname = sanitize_field_name(field_names[0])
            outfile = '_'.join([outfile_root, fname])
        else:
            outfile = outfile_root
        # Remove if exists ?
        if outfile in os.listdir():
            do_remove = input(f"{outfile} already exists remove ? (y/n) ")
            if do_remove == 'y':
                shutil.rmtree(outfile)
        # Create the plotfile dir
        os.mkdir(outfile)
        # Rewrite the header
        with open(os.path.join(outfile, "Header"), "w") as hfile:
            slc.write_2d_slice_global_header(hfile,
                                             field_names,
                                             indexes)
        # Write the level data
        

    # Array output
    elif args.format == "array":
        # Make a dict with output
        x_grid, y_grid = slc.slice_plane_coordinates()
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
            outfile = '_'.join([fname, outfile_root])
        else:
            outfile = outfile_root
        # Save with numpy so we can np.load
        np_file = outfile + ".npz"
        np.savez_compressed(np_file, **output)

        if verbose > 0:
            print("Time to save uniform grid:", 
                  np.around(time.time() - output_start, 2), "s")

    # Image output
    elif args.format == "image":
        # Compute grids
        x_grid, y_grid = slc.slice_plane_coordinates()
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
