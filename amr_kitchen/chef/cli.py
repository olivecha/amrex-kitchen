import os
import time
import numpy as np
import argparse
from .colander import Colander


def main():
    # Argument parser
    parser = argparse.ArgumentParser(
            description="Compute derived quantities from AMReX plotfiles")

    parser.add_argument(
            "plotfile", type=str,
            help="Path of the plotfile to cook")

    parser.add_argument(
            "--recipe", "-r", type=str, 
            help=("Recipe used to add new data, available recipes are:\n"
                  "- 'HRR': Heat release rate [W/m^3]\n"
                  "- 'ENT': Sensible enthalpy []\n"
                  "- 'SRi': Reaction rate of selected species"
                  "(must specify with --species argument)\n")
                  "- 'SDi': Diffusion coefficient of selected species"
                  "(must specify with --species argument)\n"))
    parser.add_argument(
            "--species", "-s", type=str, nargs='+',
            help="name of the species for which the recipe is prepared")

    parser.add_argument(
            "--mech", "-m", type=str,
            help=("Path to the Cantera mechanism to use for the derived"
                  " variables computations"))
    parser.add_argument(
            "--pressure", "-p", type=float,
            help="Pressure of the simulation data (atm)")

    args = parser.parse_args()
    """
    Input arguments sanity check
    """
    if args.plotfile is None:
        raise ValueError("Must specify a plotfile to cook")
    if args.mech is None:
        raise ValueError("Must specify a mechanism")
    if args.recipe not in ["HRR", "ENT", "SRi", "SDi"]:
        raise ValueError(f"Unknown recipe: {args.recipe}")
    if args.recipe in ['SRi', 'SDi']:
        if args.species is None:
            raise ValueError("Must specify species")

    # Colander object
    chef = Chef(plotfile=args.plotfile,
                recipe=args.recipe,
                species=args.species,
                mech=args.mech)

    cld.strain()

if __name__ == "__main__":
    main()
