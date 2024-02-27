import os
import billiard as mp
from p_tqdm import p_map
from concurrent.futures import ProcessPoolExecutor as Pool
from tqdm import tqdm
import time
import numpy as np
import argparse


def main():
    # Argument parser
    parser = argparse.ArgumentParser(
            description="Post processing utilisty for AMReX plotfiles")

    parser.add_argument(
            "plotfile", type=str,
            help="Path of the plotfile to cook")
    parser.add_argument(
            "--recipe", "-r", type=str,
            help=("Path to a script containing a recipe function, or"
                  "key of a predetermine function (see --help)"))
    parser.add_argument(
            "--mech", "-m", type=str,
            help=("Path to a Cantera kinetics mechanism (if the recipe"
                  " uses Cantera)"))
    parser.add_argument(
            "--limit_level", "-l", type=int,
            help="Maximum AMR Level read in plotfile")
    parser.add_argument(
            "--num_cpus", "-np", type=int,
            help="Number of subprocesses to use")
    parser.add_argument(
            "--output", "-o", type=str,
            help="Output path to store the post processing")

    args = parser.parse_args()

    if args.mech is not None:
        os.environ["CANTERA_MECH"] = args.mech

    from mandoline import HeaderData, parallel_cook, parallel_cook_byarray
    from mandoline.colander import Colander

    # what are we cooking
    recipe = Colander(args.recipe)
    recipe.__doc__ = args.recipe

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
