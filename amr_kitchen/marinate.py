import sys
import pickle
from amr_kitchen import PlotfileCooker


def main():
    if sys.argv[1] in ['-h', '--help']:
        print(("Saves the PlotfileCooker object of a plotfile"
               " to a pickle file"))
    else:
        pck = PlotfileCooker(sys.argv[1],
                             maxmins=True,
                             ghost=True)
        with open(sys.argv[1] + '.pkl', 'wb') as pfile:
            pickle.dump(pck, pfile)

