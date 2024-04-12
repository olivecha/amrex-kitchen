import os
import sys
import argparse
import numpy as np
from amr_kitchen import HeaderData

def main():
    # Argument parser
    parser = argparse.ArgumentParser(
            description="A quick taste of the plotfile")

    parser.add_argument(
            "plotfile", type=str,
            help="Path of the plotfile to filter")

    parser.add_argument(
            "--variable", "-v", type=str, 
            help="Variable to taste")

    args = parser.parse_args()

    print(f"__ {os.path.split(args.plotfile)[1]} __")

    hdr = HeaderData(args.plotfile, header_only=True)

    fidx = hdr.fields[args.variable]

    for lv in range(hdr.limit_level + 1):
        cell_header = os.path.join(args.plotfile,
                                  f"Level_{lv}",
                                   "Cell_H")
        with open(cell_header, "r") as chead:
            for _ in range(4):
                chead.readline()
            nboxes = int(chead.readline().split(' ')[0].replace('(',''))
            for _ in range(3*nboxes + 7):
                chead.readline()
            vmaxes = []
            for _ in range(nboxes - 1):
                data = chead.readline().split(',')
                vmaxes.append(float(data[fidx]))
            lvmax = np.max(vmaxes)
            idxmax = np.argmax(vmaxes)
            ptmax = hdr.box_centers[lv][idxmax]
            boxmax = hdr.boxes[lv][idxmax]

            if lvmax > 1e4:
                print(f"Max {args.variable} at Level {lv}: {lvmax:.2e}")
                print(f"Max value position", ptmax)
            elif lvmax > 1e2:
                print(f"Max {args.variable} at Level {lv}: {lvmax:.1f}")
            elif lvmax > 1:
                print(f"Max {args.variable} at Level {lv}: {lvmax:.2f}")
            else:
                print(f"Max {args.variable} at Level {lv}: {lvmax}")


