import os
import time
import numpy as np
import argparse
from .menu import Menu

def list_of_strings(arg):
    return arg.split(',')

def main():
    # Argument parser
    parser = argparse.ArgumentParser(
            description="Displays fields and species of a plotfile")

    parser.add_argument(
            "plotfile", type=str,
            help="Path of the plotfile to read")
    
    parser.add_argument(
            "--has_var", "-hv", type=list_of_strings,
            help=("Variable to find in the plotfile"))
    
    parser.add_argument(
            "--all", "-a", action='store_true',
            help="Flag to enable displaying every field in the database and if they are the plotfile or not")
    
    parser.add_argument(
            "--description", "-d", action='store_true',
            help="Flag to enable displaying the fields' descriptions")


    args = parser.parse_args()
    """
    Input arguments sanity check
    """
    if args.plotfile is None:
        raise ValueError("Must specify a plotfile to filter")

    # Colander object
    Menu(plt_file=args.plotfile,
         has_var=args.has_var,
         all=args.all,
         description=args.description)

if __name__ == "__main__":
    main()
