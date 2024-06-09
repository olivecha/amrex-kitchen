import os
import time
import numpy as np
import argparse
from .chef import Chef
from tqdm import tqdm


def main():
    # Argument parser
    parser = argparse.ArgumentParser(
            description="Compute derived quantities from AMReX plotfiles")

    parser.add_argument(
            "plotfile", type=str,
            help="Path of the plotfile to cook")
    parser.add_argument(
            "--outdir", "-o", type=str,
            help=("Output path of the generated plotfile"
                 " (defaults to input_plotfile'_ck' "))
    parser.add_argument(
            "--recipe", "-r", type=str, 
            help=("Recipe used to add new data, available recipes are:\n"
                  "- 'HRR': Heat release rate [W/m^3]\n"
                  "- 'ENT': Sensible enthalpy []\n"
                  "- 'SRi': Reaction rate of selected species"
                  " (must specify with --species argument)\n"
                  "- 'SDi': Diffusion coefficient of selected species"
                  " (must specify with --species argument)\n"
                  "- 'RRi' : Rate of progress of selected reactions"
                  " (must specify with the --reactions argument)\n"
                  "- user defined: the path to a python file with"
                  " a function named 'recipe' in it\n"))
    parser.add_argument(
            "--species", "-s", type=str, nargs='+',
            help="name of the species for which the recipe is prepared")
    parser.add_argument(
            "--reactions", "-R", type=int, nargs='+',
            help="index of the reactions for which the recipe is prepared")
    parser.add_argument(
            "--mech", "-m", type=str,
            help=("Path to the Cantera mechanism to use for the derived"
                  " variables computations"))
    parser.add_argument(
            "--pressure", "-p", type=float,
            help="Pressure of the simulation data (atm)")

    args = parser.parse_args()

    # Chef object
    plot_chef = Chef(args.plotfile,
                     outfile=args.outdir,
                     recipe=args.recipe,
                     species=args.species,
                     reactions=args.reactions,
                     mech=args.mech,
                     serial=False,
                     pressure=args.pressure)
    plot_chef.cook()

if __name__ == "__main__":
    main()
