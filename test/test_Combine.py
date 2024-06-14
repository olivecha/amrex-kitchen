
import os
import shutil
import unittest
from amr_kitchen.combine import combine
from amr_kitchen import PlotfileCooker
from amr_kitchen.taste import Taster

class TestCombine(unittest.TestCase):
    pfile1 = os.path.join("test_assets", "plt1_Y")
    pfile2 = os.path.join("test_assets", "plt2_F")

    def test_combine(self): 
        plt1 = PlotfileCooker(self.pfile1)
        plt2 = PlotfileCooker(self.pfile2)
        pltout = os.path.join("test", "plt_combine_out")
        combine(plt1, plt2, pltout)
        self.assertTrue(Taster(pltout, nofail=True))
        shutil.rmtree(pltout)
