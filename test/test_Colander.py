import os
import shutil
import unittest

from amr_kitchen.colander import Colander

class TestSlice(unittest.TestCase):
    pfile2d = os.path.join("test_assets", "example_plt_2d")
    pfile3d = os.path.join("test_assets", "example_plt_3d")

    def test_backward_compatible_2d(self): 
        """
        Keep all fields at all level and test
        that the resulting files are the same
        """
        try:
            all_variables=["temp",
                           "mag_vort",
                           "HeatRelease",
                           "Y(OH)",
                           "Y(CO)",
                           "Y(CH2O)",
                           "Y(C2H5)",
                           "Y(O)"]

            cld = Colander(self.pfile2d,
                           output="temp2d",
                           variables=all_variables)
            cld.strain()

            chdr_new = open(os.path.join("temp2d","Level_0","Cell_H")).read()
            chdr_ref = open(os.path.join(self.pfile2d, "Level_0", "Cell_H")).read()
            self.assertEqual(chdr_new, chdr_ref)
            shutil.rmtree("temp2d")
        except Exception as e:
            shutil.rmtree("temp2d")
            raise(e)



