import os
import time
import numpy as np
import argparse
from .menu import Menu

def list_of_strings(arg):
    arg = [argument.strip() for argument in arg.split(",")]
    return arg

def main():
    # Argument parser
    parser = argparse.ArgumentParser(
            description="Displays fields and species of a plotfile")

    parser.add_argument(
            "plotfile", type=str,
            help="Path of the plotfile to read")
    
    parser.add_argument(
            "--has_var", "-hv", type=list_of_strings,
            help=("""Variables (between " ") to find in the plotfile"""))
    
    parser.add_argument(
            "--every", "-e", action='store_true',
            help="Flag to enable displaying every field in the database and if they are the plotfile or not")
    
    parser.add_argument(
            "--description", "-d", action='store_true',
            help="Flag to enable displaying the fields' descriptions")
    
    parser.add_argument(
            "--min_max", "-m", action='store_true',
            help="Flag to enable displaying the fields' absolute min and max throughout every level")
    
    parser.add_argument(
            "--finest_lv", "-f", action='store_true',
            help="(Equivalent to menu -m -f) Flag to enable displaying the fields' min and max at the finest level")


    args = parser.parse_args()
    """
    Input arguments sanity check
    """
    if args.plotfile is None:
        raise ValueError("Must specify a plotfile to filter")

    # Colander object
    Menu(plt_file=args.plotfile,
         has_var=args.has_var,
         every=args.every,
         description=args.description,
         min_max=args.min_max,
         finest_lv=args.finest_lv,
         )

if __name__ == "__main__":
    main()
