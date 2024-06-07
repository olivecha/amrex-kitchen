import os
import shutil
import unittest
from amr_kitchen.taste import Taster
from amr_kitchen.utils import TastesBadError

class TestTaste(unittest.TestCase):
    goodplotfile_2d = os.path.join("test_assets", 
                                   "example_plt_2d")

    goodplotfile_3d = os.path.join("test_assets", 
                                   "example_plt_3d")

    badbinshape_2d = os.path.join("test_assets",
                                  "bad_plotfiles",
                                  "plt_2d_badbinaryshape")

    badindexes_2d = os.path.join("test_assets",
                                 "bad_plotfiles",
                                 "plt_2d_badindexes")

    missingbindata_2d = os.path.join("test_assets",
                                     "bad_plotfiles",
                                     "plt_2d_missingbindata")

    missingindexes_2d = os.path.join("test_assets",
                                     "bad_plotfiles",
                                     "plt_2d_missingindexes")

    badbinshape_3d = os.path.join("test_assets",
                                  "bad_plotfiles",
                                  "plt_3d_badbinaryshape")

    badindexes_3d = os.path.join("test_assets",
                                 "bad_plotfiles",
                                 "plt_3d_badindexes")

    missingbindata_3d = os.path.join("test_assets",
                                     "bad_plotfiles",
                                     "plt_3d_missingbindata")

    missingindexes_3d = os.path.join("test_assets",
                                     "bad_plotfiles",
                                     "plt_3d_missingindexes")

    def test_nofail_on_good_2d(self): 
        """
        Should not raise anything if the plotfile is
        good
        """
        # Assert the Taster evaluate as True and no
        # Errors are raised
        self.assertTrue(Taster(self.goodplotfile_2d))

    def test_nofail_on_good_3d(self): 
        """
        Should not raise anything if the plotfile is
        good
        """
        # Assert the Taster evaluate as True and no
        # Errors are raised
        self.assertTrue(Taster(self.goodplotfile_3d))

    def test_badindexes_in_level_header_2d(self):
        """
        Should raise a TastesBadError
        """
        with self.assertRaises(TastesBadError):
            Taster(self.badindexes_2d)
        # This should just return False
        self.assertFalse(Taster(self.badindexes_2d,
                                nofail=True))
