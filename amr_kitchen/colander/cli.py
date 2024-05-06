import os
import time
import numpy as np
import argparse
from .colander import Colander


def main():
    # Argument parser
    parser = argparse.ArgumentParser(
            description="Remove field and levels from AMReX plotfiles")

    parser.add_argument(
            "plotfile", type=str,
            help="Path of the plotfile to filter")

    parser.add_argument(
            "--variables", "-v", type=str, nargs='+',
            help=("Variable to keep in the filtered plotfile"))
    parser.add_argument(
            "--limit_level", "-l", type=int,
            help="Maximum AMR Level to keep in the filtered plotfile")
    parser.add_argument(
            "--serial", "-s", action='store_true',
            help="Flag to disable multiprocessing")
    parser.add_argument(
            "--output", "-o", type=str,
            help="Output path to store the filtered plotfile")

    args = parser.parse_args()
    """
    Input arguments sanity check
    """
    if args.plotfile is None:
        raise ValueError("Must specify a plotfile to filter")
    if len(args.variables) == 0:
        raise ValueError("Must keep at least one field")
    if args.output is None:
        raise ValueError("Must specify output path")

    # Colander object
    cld = Colander(plotfile=args.plotfile,
                   limit_level=args.limit_level,
                   output=args.output,
                   variables=args.variables)

    cld.strain()

if __name__ == "__main__":
    main()
