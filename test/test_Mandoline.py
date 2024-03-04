import os
import unittest
import numpy as np

from amr_kitchen.mandoline import Mandoline

class TestMandoline(unittest.TestCase):
    pfile2d = "test_assets/example_plt_2d"
    pfile3d = "test_assets/example_plt_3d"
    ref_s3d = "test_assets/ref_slices_3d"

    def test_find_field(self):
        Lv = 1

        c = Mandoline(self.pfile2d, 'temp', limit_level=Lv)
        self.assertListEqual(c.fidxs, [0])

        c = Mandoline(self.pfile2d, ['temp', 'mag_vort'], limit_level=Lv)
        self.assertListEqual(c.fidxs, [0, 1])

        c = Mandoline(self.pfile2d, 'grid_level', limit_level=Lv)
        self.assertListEqual(c.fidxs, [None])

        c = Mandoline(self.pfile2d, 'all', limit_level=Lv)
        self.assertListEqual(c.fidxs, 
                             [0, 1, 2, 3, 4, 
                              5, 6, 7, None])

        self.assertRaises(ValueError, 
                          Mandoline,
                          self.pfile2d, 
                          'wrong_field',
                          limit_level=Lv)

    def test_limit_level_array_2d(self):
        c = Mandoline(self.pfile2d, "temp")
        arr = c.limit_level_arr()
        self.assertEqual(512, arr.shape[0])
        self.assertEqual(512, arr.shape[1])
            
    def test_limit_level_array_3d(self):
        c = Mandoline(self.pfile3d, "temp")
        arr = c.limit_level_arr()
        self.assertEqual(32, arr.shape[0])
        self.assertEqual(32, arr.shape[1])

    def test_slice_plane_coordinates(self):
        c = Mandoline(self.pfile3d, "temp", 1)
        sx, sy = c.slice_plane_coordinates()
        self.assertTrue(sx[0] > c.geo_low[1])
        self.assertTrue(sx[-1] < c.geo_high[1])
        self.assertTrue(sy[0] > c.geo_low[2])
        self.assertTrue(sy[-1] < c.geo_high[2])
        
    def test_slice_data_2d_serial(self):
        m = Mandoline(self.pfile2d, "temp", serial=True)
        out = m.slice(fformat="return")
        self.assertAlmostEqual(np.min(out["temp"]), 299.9403143172542)
        self.assertAlmostEqual(np.max(out["temp"]), 1517.690555398803) 

    def test_slice_data_2d_multiprocessing(self):
        m = Mandoline(self.pfile2d, "temp", serial=False)
        out = m.slice(fformat="return")
        self.assertAlmostEqual(np.min(out["temp"]), 299.9403143172542)
        self.assertAlmostEqual(np.max(out["temp"]), 1517.690555398803) 

    def test_slice_data_3d_serial(self):
        m = Mandoline(self.pfile3d, "temp", verbose=0, serial=True)

        # x normal slice
        out = m.slice(normal=0, pos=0, fformat="return")
        self.assertAlmostEqual(np.min(out["temp"]), 297.99999999999994)
        self.assertAlmostEqual(np.max(out["temp"]), 1579.8536855390914)
        out = m.slice(normal=0, pos=m.geo_high[0], fformat="return")
        self.assertAlmostEqual(np.max(out["temp"]),  1579.8536855390923)
        self.assertAlmostEqual(np.min(out["temp"]), 297.99999999999994)
        # y normal slice
        out = m.slice(normal=1, pos=0, fformat="return")
        self.assertAlmostEqual(np.min(out["temp"]), 298.0)
        self.assertAlmostEqual(np.max(out["temp"]), 1579.7695283039325)
        out = m.slice(normal=1, pos=m.geo_high[1], fformat="return")
        self.assertAlmostEqual(np.max(out["temp"]), 1579.8536855390937)
        self.assertAlmostEqual(np.min(out["temp"]), 298.0)
        # z normal slice
        out = m.slice(normal=2, pos=0, fformat="return")
        self.assertAlmostEqual(np.min(out["temp"]), 298.0040211808012)
        self.assertAlmostEqual(np.max(out["temp"]), 298.0044818553631)
        out = m.slice(normal=2, pos=m.geo_high[2], fformat="return")
        self.assertAlmostEqual(np.max(out["temp"]), 1579.8536855390937)
        self.assertAlmostEqual(np.min(out["temp"]), 1579.7695283039257)

    def test_slice_data_3d_multiprocessing(self):
        m = Mandoline(self.pfile3d, "temp", verbose=0, serial=False)
        # x normal slice
        out = m.slice(normal=0, pos=0, fformat="return")
        self.assertAlmostEqual(np.min(out["temp"]), 297.99999999999994)
        self.assertAlmostEqual(np.max(out["temp"]), 1579.8536855390914)
        out = m.slice(normal=0, pos=m.geo_high[0], fformat="return")
        self.assertAlmostEqual(np.max(out["temp"]),  1579.8536855390923)
        self.assertAlmostEqual(np.min(out["temp"]), 297.99999999999994)

    def test_slice_to_array_3d(self):
        m = Mandoline(self.pfile3d, "temp", verbose=0, serial=True)
        # X normal
        out = m.slice(normal=0, fformat="return")
        out_ref = np.load(os.path.join(self.ref_s3d, "Sxtempref.npz"))
        self.assertTrue(np.allclose(out_ref["temp"], out["temp"]))
        # Y normal
        out = m.slice(normal=1, fformat="return")
        out_ref = np.load(os.path.join(self.ref_s3d, "Sytempref.npz"))
        self.assertTrue(np.allclose(out_ref["temp"], out["temp"]))
        # Z normal
        out = m.slice(normal=2, fformat="return")
        out_ref = np.load(os.path.join(self.ref_s3d, "Sztempref.npz"))
        self.assertTrue(np.allclose(out_ref["temp"], out["temp"]))

