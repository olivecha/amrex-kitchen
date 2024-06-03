import os
import time
import numpy as np
import argparse
from amr_kitchen.taste import tasting


def main():
    # Argument parser
    parser = argparse.ArgumentParser(
            description="Validate the sanity of AMReX plotfiles")

    parser.add_argument(
            "plotfile", type=str,
            help="Path of the plotfile to validate")
    
    args = parser.parse_args()
    """
    Input arguments sanity check
    """
    if args.plotfile is None:
        raise ValueError("Must specify a plotfile to validate")
    
    # Sanity validation 
    tasting(args.plotfile)

if __name__ == "__main__":
    main()
