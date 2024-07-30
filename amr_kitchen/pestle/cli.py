import argparse
from pestle import volume_integral
from amr_kitchen import PlotfileCooker
from pestle import field_units


def main():
    # Argument parser
    parser = argparse.ArgumentParser(
            description="Prints the volume integral of the chosen field in a plotfile")

    parser.add_argument(
            "--variable", "-v", type=str,
            help="Variable to integrate")
    parser.add_argument(
            "--limit_level", "-l", type=int,
            help="Maximum AMR Level considered")
    parser.add_argument(
            "--volfrac", "-vf", action="store_true",
            help=("Use the volFrac field to obtain more"
                  " accurate integrals for plotfiles with embedded"
                  " boundaries. The contribution of a finite volume"
                  " to the integral is taken as value * dV * volFrac."))
    parser.add_argument(
            "plotfile", type=str,
            help="Path of the plotfile to integrate")

    args = parser.parse_args()
    """
    Input arguments sanity check
    """
    if args.plotfile is None:
        raise ValueError("Must specify a plotfile to integrate")

    # Creating a PlotfileCooker instance
    pck = PlotfileCooker(args.plotfile, ghost=True)
    # Integrating the chosen fields 
    integral = volume_integral(pck=pck,
                               field=args.variable,
                               limit_level=args.limit_level,
                               use_volfrac=args.volfrac )
    if args.variable in field_units:
        units = field_units[args.variable]
    elif "I_R" in args.variable:
        units = field_units["I_R"]
    else:
        units = "[-]"

    print(f"Volume integral of {args.variable} in plotfile: {integral:.15f} {units}")

if __name__ == "__main__":
    main()
