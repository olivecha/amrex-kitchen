import sys
import pickle
from amr_kitchen import PlotfileCooker

pck = PlotfileCooker(sys.argv[1],
                     maxmins=True,
                     ghost=True)

with open(sys.argv[1] + '.pkl', 'wb') as pfile:
    pickle.dump(pck, pfile)

