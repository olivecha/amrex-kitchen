import os
import time
import numpy as np
import argparse
from amr_kitchen.taste import Taster


def main():
    # Argument parser
    parser = argparse.ArgumentParser(
            description="Validate the sanity of AMReX plotfiles")

    parser.add_argument(
            "plotfile", type=str,
            help="Path of the plotfile to validate")
    parser.add_argument(
            "--limit_level", "-l", type=int,
            help="Limit level at which the validation will occur")
    parser.add_argument(
            "--boxes", "-b", action='store_false',
            help="Default to True. Call disables the validation of number and shape of boxes")
    parser.add_argument(
            "--maxmin", "-m", action='store_false',
            help="Default to True. Call disables the validation of maxs and mins")
    parser.add_argument(
            "--coordinates", "-c", action='store_false',
            help="Default to True. Call disables the validation of boxes' coordinates")
    parser.add_argument(
            "--nan", "-n", action='store_false',
            help="Default to True. Call disables the validation of NaNs")
    parser.add_argument(
            "--nofail", "-f", action='store_true',
            help="Default to False. Call bypasses all Errors fails and prints all Errors")
    
    args = parser.parse_args()
    """
    Input arguments sanity check
    """
    if args.plotfile is None:
        raise ValueError("Must specify a plotfile to validate")
    
    # Sanity validation 
    Taster(plt_file=args.plotfile,
           limit_level=args.limit_level,
           boxes=args.boxes,
           maxmin=args.maxmin,
           coordinates=args.coordinates,
           nan=args.nan,
           nofail=args.nofail,)

if __name__ == "__main__":
    main()
