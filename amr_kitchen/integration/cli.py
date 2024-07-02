import sys
from .integration import volume_integration
from amr_kitchen import PlotfileCooker
import argparse

def list_of_strings(arg):
    arg = [argument.strip() for argument in arg.split(",")]
    return arg

def main():
    # Argument parser
    parser = argparse.ArgumentParser(
            description="Prints the volume integral of the chosen field in a plotfile")

    parser.add_argument(
            "plotfile", type=str,
            help="Path of the plotfile to integrate")

    parser.add_argument(
            "--vars", "-v", type=list_of_strings,
            help=("""Variables (between " ") to integrate"""))

    args = parser.parse_args()
    """
    Input arguments sanity check
    """
    if args.plotfile is None:
        raise ValueError("Must specify a plotfile to integrate")
    

    # Creating a PlotfileCooker instance
    pck = PlotfileCooker(args.plotfile, ghost=True)
    # Integrating the chosen fields 
    for field in args.vars:
        integral = volume_integration(pck=pck,
                                      field=field)
        print(f"Volume integral of {field} in plotfile: {integral}")


if __name__ == "__main__":
    main()
