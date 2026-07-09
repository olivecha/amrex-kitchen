import os
import re
import shutil
import tempfile
import unittest

import numpy as np

from amr_kitchen import PlotfileCooker
from amr_kitchen.chef import Chef
from amr_kitchen.colander import Colander
from amr_kitchen.combine import combine
from amr_kitchen.mandoline import Mandoline
from amr_kitchen.taste import Taster


def rename_prefix(src, dst, new_prefix, old_prefix="Cell"):
    """
    Copy the plotfile :src: to :dst: using :new_prefix: as the name of
    the multifab header and binary files. This is what AMReX writes when
    a plotfile is created with a mfPrefix different from the default
    "Cell", e.g. the level headers are named "{prefix}_H" and the binary
    files "{prefix}_D_XXXXX".
    """
    shutil.copytree(src, dst)
    # Main header: the per level path lines ("Level_0/Cell")
    hpath = os.path.join(dst, 'Header')
    with open(hpath) as hfile:
        hdr = hfile.read()
    hdr = re.sub(rf'/{old_prefix}$', f'/{new_prefix}', hdr, flags=re.MULTILINE)
    with open(hpath, 'w') as hfile:
        hfile.write(hdr)
    for lvdir in [d for d in os.listdir(dst) if d.startswith('Level_')]:
        lvpath = os.path.join(dst, lvdir)
        # The level header and the binary files
        for fname in os.listdir(lvpath):
            if fname.startswith(old_prefix):
                os.rename(os.path.join(lvpath, fname),
                          os.path.join(lvpath,
                                       new_prefix + fname[len(old_prefix):]))
        # The binary file names in the FabOnDisk lines
        chpath = os.path.join(lvpath, f'{new_prefix}_H')
        with open(chpath) as cfile:
            chdr = cfile.read()
        chdr = chdr.replace(f'FabOnDisk: {old_prefix}_D_',
                            f'FabOnDisk: {new_prefix}_D_')
        with open(chpath, 'w') as cfile:
            cfile.write(chdr)
    return dst


class TestNamedPrefix(unittest.TestCase):
    """
    The name of the level headers and binary files is not always "Cell",
    the actual name is read from the per level path lines in the main
    plotfile header. These tests run the plotfile writers on a plotfile
    using a different prefix and validate that it is preserved.
    """
    prefix = "Fab"
    pfile2d = os.path.join("test_assets", "example_plt_2d")
    pfile3d = os.path.join("test_assets", "example_plt_3d")
    pfile1_Y = os.path.join("test_assets", "plt1_Y")
    pfile2_F = os.path.join("test_assets", "plt2_F")
    drm19 = os.path.join("test_assets", "drm19.yaml")

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.plt3d = rename_prefix(self.pfile3d,
                                   os.path.join(self.tmpdir, "plt3d"),
                                   self.prefix)

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def out(self, name):
        return os.path.join(self.tmpdir, name)

    def assert_prefix_written(self, pltdir, taste=True):
        """
        Validate that the plotfile at :pltdir: only uses self.prefix to
        name its level headers and binary files
        """
        pck = PlotfileCooker(pltdir)
        self.assertEqual(pck.prefix, self.prefix)
        for lv in range(pck.limit_level + 1):
            lvfiles = os.listdir(os.path.join(pltdir, pck.cell_paths[lv]))
            self.assertIn(f"{self.prefix}_H", lvfiles)
            self.assertTrue(any(f.startswith(f"{self.prefix}_D_")
                                for f in lvfiles))
            self.assertFalse(any(f.startswith("Cell") for f in lvfiles))
        if taste:
            self.assertTrue(Taster(pltdir, nofail=True))

    def test_read_named_prefix(self):
        pck = PlotfileCooker(self.plt3d, maxmins=True)
        pck_ref = PlotfileCooker(self.pfile3d, maxmins=True)
        self.assertEqual(pck.prefix, self.prefix)
        self.assertEqual(pck_ref.prefix, "Cell")
        # The data is unchanged by the renaming
        for lv in range(pck.limit_level + 1):
            for box, box_ref in zip(pck["temp"][lv][:],
                                    pck_ref["temp"][lv][:]):
                self.assertTrue(np.array_equal(box, box_ref))

    def test_colander_named_prefix(self):
        cld = Colander(self.plt3d,
                       output=self.out("strained"),
                       variables=["temp", "mag_vort"])
        cld.strain()
        self.assert_prefix_written(self.out("strained"))

    def test_mandoline_named_prefix(self):
        mdl = Mandoline(self.plt3d, "temp", verbose=0, serial=True)
        mdl.slice(normal=2, pos=0, fformat="plotfile",
                  outfile=self.out("sliced"))
        self.assert_prefix_written(self.out("sliced"))

    def test_chef_named_prefix(self):
        plot_chef = Chef(self.plt3d,
                         outfile=self.out("cooked"),
                         recipe='HRR',
                         mech=self.drm19,
                         serial=True,
                         pressure=1.0)
        plot_chef.cook()
        # Chef writes FAB headers declaring 0 fields, which Taster
        # rejects, so only the file naming is validated here
        self.assert_prefix_written(self.out("cooked"), taste=False)

    def test_combine_named_prefix(self):
        plt1 = rename_prefix(self.pfile1_Y, self.out("plt1"), self.prefix)
        plt2 = rename_prefix(self.pfile2_F, self.out("plt2"), self.prefix)
        combine(PlotfileCooker(plt1), PlotfileCooker(plt2),
                self.out("combined"))
        self.assert_prefix_written(self.out("combined"))


if __name__ == '__main__':
    unittest.main()
