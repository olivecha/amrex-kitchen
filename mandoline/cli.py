import os
import sys
import time
import argparse
import pickle
import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import matplotlib
from mandoline import HeaderData
matplotlib.use('Agg')

coords_dict = {0:'x',
               1:'y',
               2:'z'}
def main():
    """
    Main function running the mandoline Command Line Tool
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
    
    args = parser.parse_args()

    if args.output is None:
        filename = f"{args.plotfile}_{coords_dict[args.normal]}_{args.variable}"
    else:
        filename = args.output

    # Read the header
    hdr = HeaderData(args.plotfile, max_level=args.max_level)
    
    # Compute the box indexes in the slice
    slice_idx = hdr.slice_indexes(args.normal, args.position, args.max_level)

    points = hdr.points_from_indexes(slice_idx)
    
    # Read the data with indexes
    if args.variable is None:
        variable = 'density'
    else:
        variable = args.variable
    data = hdr.load_field_from_indexes(slice_idx, variable)
        
    # Find out what x and y are in the plot
    x_coord, y_coord = [i for i in range(3) if i != args.normal]
    
    # Array output
    if args.format == "array":
        # Interpolation bounds
        x_lo, x_hi = np.min(points[x_coord]), np.max(points[x_coord])
        y_lo, y_hi = np.min(points[y_coord]), np.max(points[y_coord])

        # Reshape points as tuples
        in_points = np.vstack([points[x_coord], points[y_coord]]).T

        # Output points
        x_out = np.linspace(x_lo, x_hi, 500)
        y_out = np.linspace(y_lo, y_hi, 500)
        out_points = np.meshgrid(x_out, y_out)

        # Interpolate on output points
        out_data = griddata(in_points, data, 
                            np.transpose([arr.flatten() for arr in out_points]),
                            fill_value=0.)

        # Make a dict with output
        output = {'x':x_out,
                  'y':y_out,
                  'var':out_data.reshape(500, 500)}

        # Pickle into the jar
        pfile = open(filename + ".pkl", 'wb')
        pickle.dump(output, pfile)
        print("Saved uniform grid!")

    # Image output
    else:
        fig = plt.figure(figsize=(8, 6))

        # Special case for grid_level 
        # (use the same number of colors as the number of levels)
        if args.variable == 'grid_level':
            if args.max_level is None:
                n_levels = np.arange(0, hdr.max_level + 2)
            else:
                n_levels = np.arange(0, args.max_level + 2)
        else:
            n_levels = 100  # Looks good

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
        plt.tricontourf(points[x_coord],
                        points[y_coord],
                        data,
                        n_levels,
                        vmin=args.minimum,
                        vmax=args.maximum,
                        norm=norm,
                        cmap=cmap)

        # Add a colorbar  
        plt.colorbar(label=f'{args.variable}')
        ax = plt.gca()
        ax.set_aspect('equal')
        ax.set_xlabel(f"{coords_dict[x_coord]} [m]")
        ax.set_ylabel(f"{coords_dict[y_coord]} [m]")
        
        # save and close
        fig.savefig(filename, dpi=300)
        plt.close(fig)

if __name__ == "__main__":
    main()
