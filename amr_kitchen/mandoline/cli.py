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
from amr_kitchen import HeaderData
from .utils import sanitize_field_name, plotfile_ndims
#from .blade import slice_box, plate_box
from .mandoline import Mandoline


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
            "--format", "-f", type=str, default="image",
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

    # Class to handle slice parameters
    # Define the data needed for the slice
    # This reads the plotfile and cell headers
    # And also stores and define the slice data
    # With default values
    mand = Mandoline(args.plotfile, 
                     fields=args.variables, 
                     limit_level=args.max_level,
                     serial=args.serial,
                     verbose=args.verbose)

    # Do one slice
    mand.slice(normal=args.normal, 
               pos=args.position, 
               outfile=args.output, 
               fformat=args.format,
               uselog=args.log,
               cmap=args.colormap,
               vmin=args.minimum,
               vmax=args.maximum)


if __name__ == "__main__":
    main()
