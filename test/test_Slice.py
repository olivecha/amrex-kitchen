import os
import unittest

from amr_kitchen.mandoline import Slice

class TestSlice(unittest.TestCase):
    pfile2d = "test_assets/example_plt_2d"
    pfile3d = "test_assets/example_plt_3d"

    def test_find_field(self):
        Lv = 1

        c = Slice(self.pfile2d, 'temp', 0, 0., 
                  limit_level=Lv)
        self.assertListEqual(c.fidxs, [0])

        c = Slice(self.pfile2d, ['temp', 'mag_vort'], 
                  0, 0., limit_level=Lv)
        self.assertListEqual(c.fidxs, [0, 1])

        c = Slice(self.pfile2d, 'grid_level', 
                  0, 0., limit_level=Lv)
        self.assertListEqual(c.fidxs, [None])

        c = Slice(self.pfile2d, 'all', 
                  0, 0., limit_level=Lv)
        self.assertListEqual(c.fidxs, 
                             [0, 1, 2, 3, 4, 
                              5, 6, 7, None])

        self.assertRaises(ValueError, Slice,
                          self.pfile2d, 
                          'wrong_field', 0, 0.,
                          limit_level=Lv)

    def test_limit_level_array_2d(self):
        c = Slice(self.pfile2d, "temp")
        arr = c.limit_level_arr()
        self.assertEqual(512, arr.shape[0])
        self.assertEqual(512, arr.shape[1])
            
    def test_limit_level_array_3d(self):
        c = Slice(self.pfile3d, "temp", 1)
        arr = c.limit_level_arr()
        self.assertEqual(32, arr.shape[0])
        self.assertEqual(32, arr.shape[1])

    def test_slice_plane_coordinates(self):
        c = Slice(self.pfile3d, "temp", 1)
        sx, sy = c.slice_plane_coordinates()
        self.assertTrue(sx[0] > c.geo_low[1])
        self.assertTrue(sx[-1] < c.geo_high[1])
        self.assertTrue(sy[0] > c.geo_low[2])
        self.assertTrue(sy[-1] < c.geo_high[2])
        
