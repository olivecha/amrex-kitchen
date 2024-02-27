import os
import unittest

from amr_kitchen import HeaderData

class TestSliceData(unittest.TestCase):
    pfile2d = "test_assets/example_plt_2d"
    pfile3d = "test_assets/example_plt_3d"

    def test_load2d(self):
        for Lv in [0, 1]:
            hdr = HeaderData(self.pfile2d, limit_level=Lv)
            self.assertIsInstance(hdr, HeaderData)
            self.assertEqual(hdr.ndims, 2)
            self.assertEqual(hdr.limit_level, Lv)

    def test_load3d(self):
        for Lv in [0, 1, 2]:
            hdr = HeaderData(self.pfile3d, limit_level=Lv)
            self.assertIsInstance(hdr, HeaderData)
            self.assertEqual(hdr.ndims, 3)
            self.assertEqual(hdr.limit_level, Lv)



