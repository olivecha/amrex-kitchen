import sys
import argparse
from .chk2plt import chk2plt


def main():
    # Argument parser
    parser = argparse.ArgumentParser(
            description="Convert PeleLmeX checkpoints to plotfiles")

    parser.add_argument(
            "--checkpoint", "-c", type=str,
            help="Path of the checkpoint to convert")

    parser.add_argument(
            "--plotfile_ref", "-p", type=str, 
            help="Path to a reference plotfile to read the species")
    parser.add_argument(
            "--species", "-s", type=int, nargs="+",
            help="Space separated species expected to be in the checkpoint")
    parser.add_argument(
            "--include_gradp", "-ip", action='store_false',
            help="Flag to include pressure gradient (defaults to true)")
    parser.add_argument(
            "--include_reactions", "-ir", action='store_true',
            help="Flag to include species reaction rates (defaults to false)")
    parser.add_argument(
            "--floor_massfracs", "-f", action='store_false',
            help="Flag to force species mass fractions to sum to 1.0")
    parser.add_argument(
            "--output", "-o", type=str,
            help="Output path to store the plotfile")
    args = parser.parse_args()
    """
    Input arguments sanity check
    """
    if args.checkpoint is None:
        raise ValueError("Must specify a checkpoint to convert")
    if (args.plotfile_ref is None and len(args.species) == 0):
        raise ValueError("Need reference plotfiles or species list to convert")

    # Colander object
    chk2plt(chkdir=args.checkpoint,
            target_plotfile=args.plotfile_ref,
            species=args.species,
            gradp=args.include_gradp,
            species_reactions=args.include_reactions,
            floor_massfracs=args.floor_massfracs,
            pltdir=args.output)

if __name__ == "__main__":
    main()

