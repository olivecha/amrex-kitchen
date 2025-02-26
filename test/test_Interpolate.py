import os
import shutil
import unittest
from amr_kitchen import PlotfileCooker

class TestSlice(unittest.TestCase):
    pfile2d = os.path.join("test_assets", "example_plt_2d")
    pfile3d = os.path.join("test_assets", "example_plt_3d")

    @unittest.skip("Not implemented yet")
    def test_interpolation_call(self): 
        pck = PlotfileCooker(self.pfile3d)
        print(pck["temp"](0.012,0.012,0.012))
