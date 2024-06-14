import sys
from combine import combine
from amr_kitchen import PlotfileCooker

def main():
    plt1, plt2, pltout = sys.argv[1:]

    combine(PlotfileCooker(plt1), 
            PlotfileCooker(plt2), 
            pltout)

if __name__ == "__main__":
    main()
