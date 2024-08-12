import os
import sys

def main():
    if sys.argv[1] in ["-h", "--help"]:
        print("Display the current time step of a plotfile")
    else:
        pfdir = sys.argv[1]
        with open(os.path.join(pfdir, 'Header')) as hfile:
            hfile.readline()
            nfields = int(hfile.readline())
            for _ in range(nfields):
                hfile.readline()
            hfile.readline()
            print(f"Plotfile time =", float(hfile.readline()))


