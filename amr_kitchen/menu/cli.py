import os
import time
import numpy as np
import argparse
from .menu import Menu


def main():
    # Argument parser
    parser = argparse.ArgumentParser(
            description="Reads fields and species of a plotfile")

    parser.add_argument(
            "plotfile", type=str,
            help="Path of the plotfile to read")

    parser.add_argument(
            "--fields", "-f", action='store_false',
            help="Flag to disable showing the fields in the plotfile")

    parser.add_argument(
            "--species", "-s", action='store_false',
            help="Flag to disable showing the species in the plotfile")
    
    parser.add_argument(
            "--has_var", "-hv", action='store_true',
            help="Flag to enable showing all the fields in the database and if they are the plotfile or not")


    args = parser.parse_args()
    """
    Input arguments sanity check
    """
    if args.plotfile is None:
        raise ValueError("Must specify a plotfile to filter")

    # Colander object
    Menu(plt_file=args.plotfile,
         variables=args.fields,
         species=args.species,
         has_var=args.has_var,)

if __name__ == "__main__":
    main()
