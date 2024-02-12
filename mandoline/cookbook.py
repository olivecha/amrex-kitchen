import numpy as np
import cantera as ct
import multiprocessing as mp
import os

MECH = os.environ["CANTERA_MECH"]
GAS = ct.Solution(MECH)
SPECIES = {f"Y({sp.name})":sp.name for sp in GAS.species()}

def hrr_recipe(arrays):
    """
    HeatRelease
    """
    out = np.empty_like(arrays['temp'])
    for i in range(out.shape[0]):
        for j in range(out.shape[1]):
            for k in range(out.shape[2]):
                Y = {}
                for sp in SPECIES:
                    Y[SPECIES[sp]] = arrays[sp][i, j, k]
                # Handle Embedded boundaries
                if np.isclose(np.sum(list(Y.values())), 0.):
                    out[i, j, k] = 0. 
                else:
                    GAS.TPY = arrays['temp'][i, j, k], ct.one_atm, Y
                    out[i, j, k] = GAS.heat_release_rate
    return out

# For import
recipes = {'HeatRelease':hrr_recipe}
