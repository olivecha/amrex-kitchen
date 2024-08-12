import sys
import argparse
from combine import combine
from amr_kitchen import PlotfileCooker

def main():
    # Argument parser
    parser = argparse.ArgumentParser(
            description="Combine two AMReX plotfiles")

    parser.add_argument(
            "--plotfile1", "-p1", type=str,
            help="Path of the first plotfile to combine")
    parser.add_argument(
            "--plotfile2", "-p2", type=str,
            help="Path of the second plotfile to combine")
    parser.add_argument(
            "--vars1", "-v1", type=list_of_strings,
            help=("Comma or space separated variables between quotes"
                  " to keep in the first plotfile"))
    parser.add_argument(
            "--vars2", "-v2", type=list_of_strings,
            help=("Comma or space separated variables between quotes"
                  " to keep in the second plotfile"))
    parser.add_argument(
            "--output", "-o", type=str,
            help="Output path to store the combined plotfile")
    parser.add_argument(
            "--serial", "-s", action='store_true',
            help="Flag to disable multiprocessing")
    
    args = parser.parse_args()
    plt1, plt2, pltout = sys.argv[1:]

    combine(PlotfileCooker(plt1), 
            PlotfileCooker(plt2), 
            pltout)

if __name__ == "__main__":
    main()
