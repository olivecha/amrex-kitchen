import os
import time
import numpy as np
import argparse
from .taste import Taster
from tqdm import tqdm


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
            "--no_bin_headers", "-nh", action='store_false',
            help="Do not validate the indexes in the binary headers")
    parser.add_argument(
            "--no_bin_shape", "-ns", action='store_false',
            help="Do not validate the shape of the binary data")
    parser.add_argument(
            "--bin_data", "-bd", action='store_true',
            help=("Validate the binary data against the max/mins in the"
                  " level headers and warn if NaNs are encountered"
                  " (This can take a while)"))
    parser.add_argument(
            "--box_coords", "-bc", action='store_true',
            help=("Validate that the boxes coordinates in the plotfile"
                  " header are consistent with the indices in the level"
                  " headers"))
    parser.add_argument(
            "--nofail", "-nf", action='store_true',
            help="Don't fail on the first error and perform all the checks")
    parser.add_argument(
            "--verbose", "-v", type=int,
            help="Verbosity level (Defaults to 1)")

    args = parser.parse_args()
    """
    Input arguments sanity check
    """
    if args.plotfile is None:
        raise ValueError("Must specify a plotfile to validate")

    Taster(plt_file=args.plotfile,
           limit_level=args.limit_level,
           binary_headers=args.no_bin_headers,
           binary_shape=args.no_bin_shape,
           binary_data=args.bin_data,
           boxes_coordinates=args.box_coords,
           nofail=args.nofail,
           verbose=args.verbose)

if __name__ == "__main__":
    main()
