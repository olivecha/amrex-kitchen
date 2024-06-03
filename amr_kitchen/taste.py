import os
import sys
import time
import multiprocessing
import pickle
import numpy as np
from tqdm import tqdm
from amr_kitchen import PlotfileCooker
from amr_kitchen.utils import TastesBadError

# This should live somewhere
def indexes_and_shape_from_header(header):
    """
    This takes the byte string of a box header in a plotfile binary
    file and infers the indexes of the box and the number of fields
    in the plotfile
    """
    h = header.decode('ascii')
    start, stop, _, nfields = h.split()[-4:]
    nfields = int(nfields)
    start = np.array(start.split('(')[-1].replace(')','').split(','), dtype=int)
    stop = np.array(stop.replace('(', '').replace(')', '').split(','), dtype=int)
    shape = stop - start + 1
    shape = [s for s in shape]
    shape.append(nfields)
    return [start, stop], tuple(shape)

# The prototype of this should also live somewhere as 
# it is faster than how binary files are read in colander.py
def mp_read_binfile_minmax(args):
    """
    Multiprocessing function to read binary file data
    with minimal input by using the ascii encoded
    box headers
    """
    bfilename = args[0]
    indexes = []
    min_vals = []
    max_vals = []
    with open(args, 'rb') as bfile:
        while True:
            try:
                h = bfile.readline()
                idx, shape = indexes_and_shape_from_header(h)
                arr = np.fromfile(bfile, 'float64', np.prod(shape))
                arr = arr.reshape(shape, order='F')
                indexes.append(idx)
                min_vals.append(np.min(arr, axis=len(shape)-1))
                max_vals.append(np.max(arr, axis=len(shape)-1))
            except Exception as e:
                break
    return indexes, min_vals, max_vals

def tasting(plt_file):

    # Create the PlotfileCooker class instance
    # With the min/max method it should be:
    # pck = PlotfileCooker(plt1, minmaxs=True)
    pck = PlotfileCooker(plt_file)

    # Check that the number of box coordinates matches the

    # **** My bad but we should always use the name "box" instead of cell ****
    # Faudrait combiner les attributs boxs et cells mais c'est un bon petit 
    # projet de refractoring genre pck.boxes devient pck.boxes['geo']
    # et pck.cells['files'] devient pck.boxes['files']

    # number of files in the cell header data
    for lv in range(pck.limit_level + 1):
        print(f"Checking level {lv}")
        # number of cells and fields according to the level header 
        nbr_cells_header = len(pck.cells[lv]['files'])
        # Same nbr of fields at each level 
        nbr_fields_header = pck.nfields

        # mins and maxs in the header data 
        mins_header, maxs_header = [], []
        all_mins_header, all_maxs_header = pck.min_max_header()
        for key in all_mins_header[lv]:
            # Ici tu calcule les valeurs min/maxs entre toutes les boxes
            # Pour le niveau courant
            mins_header.append(np.min(all_mins_header[lv][key]))
            maxs_header.append(np.max(all_maxs_header[lv][key]))
        # Ici tu calcule les valeurs max/mins parmis tout les fields
        # Je pense pas c'est ce que tu voulais faire...
        # Par contre c'est intéressant à avoir pour regarder s'il y a des NotANumber
        # Tu pourrais valider np.nanmax(mins_header) == np.max(mins_header)...
        # Tu pourrais profiler mais je pense c'est plus rapide que np.isnan(data).any()
        # Après je sais pas si le module qui sauvegarde les plotfiles filtre les nans 
        # Pour les max mins des boxes je vais regarder
        min_header, max_header = np.min(mins_header), np.max(maxs_header)

        nbr_box = 0 
        all_mins_cells, all_maxs_cells = [], []

        for cells in pck.bybinfile(lv):
            # Let's add the number of box in each binary file for the level 
            nbr_box += len(mp_read_binfile_minmax(cells[0])[0])

            # Let's collect the mins and maxs of every binary file for the level 
            mins_cell, maxs_cell = [], []
            idx, all_mins_box, all_maxs_box = mp_read_binfile_minmax(cells[0])
            # Ici all_mins_box c'est une liste la longueur du nombre de box dans la binary file
            # Et chaque élément c'est les valeurs mins/maxs pour chaque field
            # all_mins_box[100][3] c'est la valeur min dans la box 100 au niveau actuel du field avec 9
            # l'index 3
            for box in all_mins_box:
                mins_cell.append(np.min(box))
            all_mins_cells.append(np.min(mins_cell))

            for box in all_maxs_box:
                maxs_cell.append(np.max(box))
            all_maxs_cells.append(np.max(maxs_cell))

            # Let's read the binary files
            with open(cells[0], 'rb') as bfile:
                h = bfile.readline()
                # Index and shape of indiviual cells
                idx, shape = indexes_and_shape_from_header(h)
            
            # Shape according to the level header 
            shape_header = tuple(np.append(((idx[1]-idx[0])+1),nbr_fields_header))

            # Comparison between shapes
            try:
                raise TastesBadError((shape,shape_header))
            except TastesBadError as e:
                if e.value[0] != e.value[1]:
                    print(f"The shapes are not the same : {e.value[0]} != {e.value[1]} at cell {cells[0]}")
                    raise TastesBadError("Shape Error") 

        # Comparison between number of cells 
        try:
            raise TastesBadError((nbr_cells_header,nbr_box))
        except TastesBadError as e:
            if e.value[0] != e.value[1]:
                print(f"The number of cells is not the same : {e.value[0]} != {e.value[1]} at level {lv}")
                raise TastesBadError("Box Number Error") 
        #else:
            # Je suis pas sur encore mais j'avais plus pensé à une stratégie
            # no news is good news et le seul output c'est tqdm pour voir
            # qu'il se passe quelque chose quand la plotfile est longue
            # Aussi créer une Exception custom et raise la avec le plus d'information
            # possible dès qu'il y a un problème pour pas continuer à analyser pour rien
            # faudrait aussi que l'outil soit défini dans une fonction et qu'il retourne 
            # True/False si la plotfile est bonne ou pas avec un raise_error=True/False flag
            # Ça serait pour pouvoir loop sur toutes les plotfiles dans un directory dequoi demême
            #print(f"The number of cells in level {lv} is ok")

        # Comparison between the cells' and the level's mins and maxs :
        # np.isclose devrait fonctionner entre les deux 'mins'
        for m in all_mins_cells:
            try:
                raise TastesBadError((m,min_header))
            except TastesBadError as e:
                if e.value[0] < e.value[1]:
                    print(f"Level {lv}'s min : {e.value[1]} --> not absolute because {e.value[0]} < {e.value[1]} at cell {cells[0]}")
                    raise TastesBadError("Minimum Error") 
                
        for m in all_maxs_cells:
            try:
                raise TastesBadError((m,max_header))
            except TastesBadError as e:
                if e.value[0] > e.value[1]:
                    print(f"Level {lv}'s max : {e.value[1]} --> not absolute because {e.value[0]} > {e.value[1]} at cell {cells[0]}")
                    raise TastesBadError("Maximum Error") 
        print("Done!")
    

