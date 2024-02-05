import os
import unittest

from mandoline.slice_data import SliceData

class TestSliceData(unittest.TestCase):

    def test_find_field(self):
        test_pfile = "test_assets/example_plt_2d"
        Lv = 1

        c = SliceData(test_pfile, 'temp', 0, 0., limit_level=Lv)
        self.assertListEqual(c.fidxs, [0])

        c = SliceData(test_pfile, ['temp', 'mag_vort'], 0, 0., limit_level=Lv)
        self.assertListEqual(c.fidxs, [0, 1])

        c = SliceData(test_pfile, 'grid_data', 0, 0., limit_level=Lv)
        self.assertListEqual(c.fidxs, [None])

        c = SliceData(test_pfile, 'all', 0, 0., limit_level=Lv)
        self.assertListEqual(c.fidxs, [0, 1, 2, 3, 4, 5, 6, 7, None])

        self.assertRaises(ValueError, SliceData,
                          test_pfile, 'wrong_field', 0, 0.,
                          limit_level=Lv)
            
