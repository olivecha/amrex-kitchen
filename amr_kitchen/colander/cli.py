import os
from p_tqdm import p_map
from multiprocessing import Pool
from tqdm import tqdm
import time
import numpy as np
import argparse
from mandoline import HeaderData, parallel_cook, parallel_cook_byarray
from mandoline.colander import Colander


def main():
    # Argument parser
    parser = argparse.ArgumentParser(
            description="Remove field and levels from AMReX plotfiles")

    parser.add_argument(
            "plotfile", type=str,
            help="Path of the plotfile to filter")
    parser.add_argument(
            "--variables", "-r", type=str, nargs='+',
            help=("Variable to keep in the filtered plotfile"))
    parser.add_argument(
            "--limit_level", "-l", type=int,
            help="Maximum AMR Level to keep in the filtered plotfile")
    parser.add_argument(
            "--serial", "-s", action='store_true',
            help="Flag to disable multiprocessing")
    parser.add_argument(
            "--output", "-o", type=str,
            help="Output path to store the filtered plotfile")

    args = parser.parse_args()
    """
    Input arguments sanity check
    """
    if args.plotfile is None:
        raise ValueError("Must specify a plotfile to filter")
    if len(args.variables) == 0:
        raise ValueError("Must keep at least one field")
    if args.output is None:
        raise ValueError("Must specify output path")

    # Header data
    hdr = HeaderData(args.plotfile, limit_level=args.limit_level)
    # Create the output structure
    hdr.make_dir_tree(args.output)

    # Store new file offsets
    new_offsets = {}
    # START COOKING
    for Lv in range(hdr.limit_level + 1):
        level_files = np.array(hdr.cells[f"Lv_{Lv}"]["files"])
        cell_head_r = os.path.join(hdr.pfile,
                                   hdr.cell_paths[Lv] + '_H')
        # So we can use index
        call_args = []
        # For unique(cell_files)
        for bfile in np.unique(level_files):

            # Indexes of cells in the binary file
            cell_indexes = np.arange(level_files.shape[0])[level_files == bfile]
            # Path to the new binary file
            bfile_w = os.path.join(args.output,
                                   hdr.cell_paths[Lv].split('/')[0],
                                   bfile.split('/')[-1])
            mp_args = {"bfile":bfile,
                       "bfile_w":bfile_w,
                       "indexes":np.array(hdr.cells[f"Lv_{Lv}"]["indexes"])[cell_indexes],
                       "cell_indexes":cell_indexes,
                       "offsets_r":np.array(hdr.cells[f"Lv_{Lv}"]["offsets"])[cell_indexes],
                       "nvars":hdr.nvars,
                       "fields":hdr.fields,
                       "recipe":recipe,
                       "ncells":len(level_files)}
            call_args.append(mp_args)
        cook_start = time.time()
        print(call_args[0]["recipe"])
        print(f"Cooking Level {Lv}...")
        # One process per file
        new_offsets_stack = []
        pool = Pool()
        new_offsets_stack = pool.map(parallel_cook, call_args)
        #for arg in tqdm(call_args):
        #    new_offsets_stack.append(parallel_cook(arg))
        
        print(f"Done!", f"({np.around(time.time() - cook_start, 2)} s)")
        # Reduce arrays together
        new_offsets[Lv] = np.sum(new_offsets_stack, axis=0)

        # Update the new Cell header
        cell_head_w = os.path.join(args.output, hdr.cell_paths[Lv] + "_H")

        with open(cell_head_w, 'w') as ch_w, open(cell_head_r, 'r') as ch_r:
            # First two lines
            for i in range(2):
                l = ch_r.readline()
                ch_w.write(l)
            _ = ch_r.readline()
            ch_w.write("1\n")
            while True:
                l = ch_r.readline()
                if "FabOnDisk:" in l:
                    new_l = l.split()[:-1]
                    new_l.append(str(new_offsets[Lv][0]))
                    ch_w.write(' '.join(new_l) + "\n")
                    break
                else:
                    ch_w.write(l)
            for fst in new_offsets[Lv][1:]:
                l = ch_r.readline()
                new_l = l.split()[:-1]
                new_l.append(str(fst))
                ch_w.write(' '.join(new_l) + "\n")
            for l in ch_r:
                ch_w.write(l)

    # Update the plotfile header
    header_r = os.path.join(hdr.pfile, 'Header')
    header_w = os.path.join(args.output, 'Header')
    with open(header_r, 'r') as hr, open(header_w, 'w') as hw:
        # Version
        hw.write(hr.readline())
        _ = hr.readline()
        # Number of fields
        hw.write("1\n")
        for i in range(hdr.nvars):
            _ = hr.readline()
        # First word in recipe function
        hw.write(recipe.__doc__.split()[0] + '\n')
        # Rest of header (test if can be removed)
        for l in hr:
            if f"Level_{hdr.limit_level}/Cell" in l:
                hw.write(l)
                break
            else:
                hw.write(l)

if __name__ == "__main__":
    main()
