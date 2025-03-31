import os
import shutil
import unittest
from amr_kitchen import PlotfileCooker

class TestSlice(unittest.TestCase):
    pfile2d = os.path.join("test_assets", "example_plt_2d")
    pfile3d = os.path.join("test_assets", "example_plt_3d")

    def test_interpolation_call(self): 
        pck = PlotfileCooker(self.pfile3d)
        # Test point inside box case
        print(pck["temp"](0.002,0.002,0.002))
        # Test point at box boundary case
        print(pck["temp"](0.012,0.012,0.012))
