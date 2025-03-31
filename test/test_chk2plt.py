import os
import shutil
import unittest
from amr_kitchen.taste import Taster
from amr_kitchen.chk2plt import CheckpointReader, chk2plt

class Testchk2plt(unittest.TestCase):
    chk3d = os.path.join("test_assets", "example_chk_3d")
    plt_ref = os.path.join("test_assets", "example_plt_3d")

    def test_checkpoint_reader(self):

        chk1 = CheckpointReader(self.chk3d)
        chk2 = CheckpointReader(self.chk3d, maxmins=True)

    def test_chk2plt(self):

        chk2plt(self.chk3d, self.plt_ref, species=[],
                gradp=True, species_reactions=False, 
                floor_massfracs=True, pltdir=os.path.join('test','plt_tmp'))

        Taster(os.path.join('test', 'plt_tmp'))

        shutil.rmtree(os.path.join("test","plt_tmp"))



