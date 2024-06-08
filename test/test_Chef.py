import os
import shutil
import unittest
from amr_kitchen.chef import Chef

class TestSlice(unittest.TestCase):
    pfile2d = os.path.join("test_assets", "example_plt_2d")
    pfile3d = os.path.join("test_assets", "example_plt_3d")
    drm19 = os.path.join("test_assets", "drm19.yaml")

    def test_heat_release_recipe(self): 
        plot_chef = Chef(self.pfile3d,
                         outfile=os.path.join("test","chef_test_hrr"),
                         recipe='HRR',
                         mech=self.drm19,
                         serial=True,
                         pressure=1.0)
        plot_chef.cook()
        # TODO: validate HRR with existing plotfile
        with open(os.path.join("test",
                               "chef_test_hrr",
                               "Header")) as hfile:
            hfile.readline()
            hfile.readline()
            self.assertTrue("HeatRelease" in hfile.readline())

        shutil.rmtree(os.path.join("test","chef_test_hrr"))


    def test_enthalpy_recipe(self): 
        plot_chef = Chef(self.pfile3d,
                         outfile=os.path.join("test","chef_test_ent"),
                         recipe='ENT',
                         mech=self.drm19,
                         serial=True,
                         pressure=1.0)
        plot_chef.cook()
        # TODO: validate data with existing plotfile
        with open(os.path.join("test",
                               "chef_test_ent",
                               "Header")) as hfile:
            hfile.readline()
            hfile.readline()
            self.assertTrue("Enthalpy" in hfile.readline())

        shutil.rmtree(os.path.join("test","chef_test_ent"))


    def test_reaction_rate_recipe(self): 
        plot_chef = Chef(self.pfile3d,
                         outfile=os.path.join("test","chef_test_Ri"),
                         recipe='SRi',
                         species=['O2', 'H2'],
                         mech=self.drm19,
                         serial=True,
                         pressure=1.0)
        plot_chef.cook()
        # TODO: validate data with existing plotfile
        with open(os.path.join("test",
                               "chef_test_Ri",
                               "Header")) as hfile:
            hfile.readline()
            hfile.readline()
            self.assertTrue("IR" in hfile.readline())

        shutil.rmtree(os.path.join("test","chef_test_Ri"))

    def test_1sp_reaction_rate_recipe(self): 
        plot_chef = Chef(self.pfile3d,
                         outfile=os.path.join("test","chef_test_Ri"),
                         recipe='SRi',
                         species=['O2'],
                         mech=self.drm19,
                         serial=True,
                         pressure=1.0)
        plot_chef.cook()
        # TODO: validate data with existing plotfile
        with open(os.path.join("test",
                               "chef_test_Ri",
                               "Header")) as hfile:
            hfile.readline()
            hfile.readline()
            self.assertTrue("IR" in hfile.readline())

        shutil.rmtree(os.path.join("test","chef_test_Ri"))

    def test_diffusion_coefficients_recipe(self): 
        plot_chef = Chef(self.pfile3d,
                         outfile=os.path.join("test","chef_test_Di"),
                         recipe='SDi',
                         species=["H", "H2", "O2"],
                         mech=self.drm19,
                         serial=True,
                         pressure=1.0)
        plot_chef.cook()
        # TODO: validate data with existing plotfile
        with open(os.path.join("test",
                               "chef_test_Di",
                               "Header")) as hfile:
            hfile.readline()
            hfile.readline()
            self.assertTrue("DI" in hfile.readline())

        shutil.rmtree(os.path.join("test","chef_test_Di"))

    def test_user_sarray_recipe(self):
        plot_chef = Chef(self.pfile3d,
                         outfile=os.path.join("test","chef_test_user_sarr"),
                         recipe=os.path.join('test_assets',
                                             'user_recipes',
                                             'fuel_oxy_ratio.py'),
                         mech=self.drm19,
                         serial=True,
                         pressure=1.0)
        plot_chef.cook()
        # TODO: validate data with existing plotfile
        with open(os.path.join("test",
                               "chef_test_user_sarr",
                               "Header")) as hfile:
            hfile.readline()
            hfile.readline()
            self.assertTrue("omega_ratio_H2_O2" in hfile.readline())
        shutil.rmtree(os.path.join("test","chef_test_user_sarr"))

    def test_user_plotfile_recipe(self):
        plot_chef = Chef(self.pfile3d,
                         outfile=os.path.join("test","chef_test_user_pfile"),
                         recipe=os.path.join('test_assets',
                                             'user_recipes',
                                             'mass_frac_ratio.py'),
                         mech=self.drm19,
                         serial=True,
                         pressure=1.0)
        plot_chef.cook()
        # TODO: validate data with existing plotfile
        with open(os.path.join("test",
                               "chef_test_user_pfile",
                               "Header")) as hfile:
            hfile.readline()
            hfile.readline()
            self.assertTrue("Y_ratio_H2_CH4" in hfile.readline())
        shutil.rmtree(os.path.join("test","chef_test_user_pfile"))

    def test_byreaction_field(self):
        plot_chef = Chef(self.pfile3d,
                         outfile=os.path.join("test","chef_test_byreaction"),
                         recipe="RRi",
                         mech=self.drm19,
                         serial=True,
                         reactions=[0, 4, 5],
                         pressure=1.0)
        plot_chef.cook()
        # TODO: validate values inside the plotfile
        with open(os.path.join("test",
                               "chef_test_byreaction",
                               "Header")) as hfile:
            hfile.readline()
            hfile.readline()
            fields = []
            for ri in range(3):
                fields.append(hfile.readline())
            for ridx, field in zip([0, 4, 5],
                                   fields):
                self.assertTrue(field == f"R{ridx}\n")
        #shutil.rmtree(os.path.join("test", "chef_test_byreaction"))

