import os
import shutil
import unittest
from amr_kitchen import PlotfileCooker
from amr_kitchen.pestle import volume_integral

class Pestle(unittest.TestCase):
    example2d = os.path.join("test_assets", "example_plt_2d")
    example3d = os.path.join("test_assets", "example_plt_3d")
    exampleEB = os.path.join("test_assets", "plt_eb_3d")
    ref_vals_3d = {"density" : 2.459271208e-06,
                   "rhoh" : -0.2807197514,
                   "temp" : 0.00408613148,
                   "mag_vort" : 0.0001975304934,
                   "Y(H)" : 4.770825947e-12,
                   "Y(O2)" : 6.566789561e-07,
                   "Y(CO2)" : 1.446602546e-07,
                   "mixture_fraction" : 8.358649731e-08,}
    ref_regtest_volfract = 0.000072843417981
    ref_amrex_novolfrac = 7.315143997e-05


    def test_fail_on_2d(self):
        with self.assertRaises(ValueError):
            pck = PlotfileCooker(self.example2d, ghost=True)
            _ = volume_integral(pck, "temp")

    def test_with_ref_3d(self):
        pck = PlotfileCooker(self.example3d, ghost=True)
        for ky in self.ref_vals_3d:
            self.assertAlmostEqual(volume_integral(pck, ky),
                                   self.ref_vals_3d[ky])

    def test_reg_with_volfrac(self):
        pck = PlotfileCooker(self.exampleEB, ghost=True)
        self.assertAlmostEqual(volume_integral(pck, "density", use_volfrac=True),
                               self.ref_regtest_volfract)
        self.assertAlmostEqual(volume_integral(pck, "density", use_volfrac=False),
                               self.ref_amrex_novolfrac)

