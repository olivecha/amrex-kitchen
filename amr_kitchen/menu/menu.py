import os
import sys
import time
import traceback
import multiprocessing
import pickle
import numpy as np
import re
from tqdm import tqdm

from amr_kitchen import PlotfileCooker

# TODO: not sure this needs to be a class
# could be a single script parsing the arguments 
class Menu(PlotfileCooker):
    """
    A class to read fields and species of a plotfile
    """

    def __init__(self, plt_file, variables=True, species=True,):
        """
        Constructor for the Menu
        """
        # TODO: The plotfile is already stored in the pfile attribute
        # when you call super().__init__ (which instantiate PlotfileCooker)
        self.plt_file = plt_file
        self.variables = variables
        self.species = species

        super().__init__(plt_file,)

        # Database in {Acronym: (Regular Expression, Definition)} format
        self.field_info = {"velocity": (r"^\w+_velocity$",
                                        "Velocity [units]"),
                            "density": (r"^density$",
                                        "Density [units]"),
                            "rhoh": (r"^rhoh$",
                                     "... [units]"),
                            "temp": (r"^temp$",
                                     "Temperature [units]"),
                            "RhoRT": (r"^RhoRT$",
                                      "... [units]"),
                            "divu": (r"^divu$",
                                     "... [units]"),
                            "gradp": (r"^gradp+\w$",
                                      "... [units]"),
                            "I_R": (r"^I_R\(.+\)$",
                                    "Species reaction rates [kg / s m^3]"),
                            "FunctCall": (r"^FunctCall$",
                                          "... [units]"),
                            "HeatRelease": (r"^HeatRelease$",
                                            "Heat release rate from chem. reactions [W / s m^3]"),
                            "avg_pressure": (r"^avg_pressure$",
                                             "Average Pressure [units]"),
                            "mag_vort": (r"^mag_vort$",
                                         "Vorticity magnitude [units]"),
                            "Y": (r"^Y\(.+\)$","... [units]"),
                            "mass_fractions": (r"^mass_fractions$",
                                               "Species mass fractions [units]"),
                            "mole_fraction": (r"^mole_fraction$",
                                              "Species mole fractions [units]"),
                            "diffcoeff": (r"^diffcoeff$",
                                          "Species mixture-averaged diffusion coefficients [units]"),
                            "lambda": (r"^lambda$",
                                       "Thermal diffusivity [units]"),
                            "viscosity": (r"^viscosity$",
                                          "Mixture viscosity [units]"),
                            "mixture_fraction": (r"^mixture_fraction$",
                                                 "Mixture fraction based on Bilger’s element formulation [units]"),
                            "progress_variable": (r"^progress_variable$",
                                                  "Progress variable based on a linear combination of Ys, T [units]"),
                            "avg_pressure": (r"^avg_pressure$",
                                             "Cell-averaged pressure [units]"),
                            "vorticity": (r"^vorticity$",
                                          "VortZ (2D) or VortX, VortY, VortZ (3D) [units]"),
                            "Qcrit": (r"^Qcrit$","Q-Criterion [units]"),
                            "kinetic_energy": (r"^kinetic_energy$",
                                               "Kinetic energy [units]"),
                            "enstrophy": (r"^enstrophy$","Enstrophy [units]"),
                            "rhominsumrhoY": (r"^rhominsumrhoY$",
                                              "Rho minus sum of rhoYs [units]"),
                            "DistributionMap": (r"^DistributionMap$",
                                                "The MPI-rank of each box [units]"),
                            }
        self.regexp_species = re.compile(self.field_info["Y"][0])

        self.menu()

    # TODO: I guess a class is okay if you do it like this
    # It'll prevent having a 200+ lines function which is
    # not ideal
    def menu(self):
        # TODO: the nice formatting I was talking about:
        # Check if there is a field defined for multiple species
        """if self.has_mass_fracs:
            # Get the length of the string of each species
            sp_lens = [len(sp) for sp in Y_species]
            # Amount of spaces padding so each string is the 
            # same length if max(sp_lens) = 3:
            # "CH4" stays the same "H" becomes "H  "
            pad = [l + (max(sp_lens) - l) for l in sp_lens]
            # Do the padding using format string in the format specifyer
            # (just found out you can do this)
            sp_padded = [f"{sp: <{p}}" for sp, p in zip(species, pad)]
            # Here this creates nice printable lines with 8 species
            # TODO: Up to a max of like 10 species per line (do some tests)
            # try to find a number of species per line that gives the fullest
            # block i.e. the last lines has close to the number for each line
            # Its gonna be some logic with % and // operators"""
			
			#Could be different be *creative*
			#$____________________ Plotfile Menu ____________________$
			#HeatRelease       : Heat Release rate [W / s m^3]
			#progress_variable : User defined reaction progress
			#I_R               : Species reaction rates [kg / s m^3
			

        if self.variables:
            self.variables_finder()
        if self.species:
            self.species_finder()
    


    def variables_finder(self):
        # Let's find all the fields in the plot file
        list_fields = []
        for field in list(self.fields):
            for key in self.field_info:
                regexp = re.compile(self.field_info[key][0])
                if regexp.search(field):
                    if key not in list_fields:
                        list_fields.append(key)
                        break
                    else:
                        break
            # A plot_file's field is not in the database
            else:
                if field not in list_fields:
                    list_fields.append(field)

        # Finding the maximum lenght of printed lines
        len_of_fields = [len(field) for field in list_fields]
        descriptions = []
        for field in list_fields:
            descriptions.append(self.field_info[field][1])
        len_of_description = [len(description) for description in descriptions]

        max_field = np.max(len_of_fields)  
        max_description = np.max(len_of_description)
        total_lenght = max_field + max_description + 2 

        title = "\n" + (total_lenght//3)*(" ")+"Fields found in file:"
        cap = "+"+("-"*total_lenght)+"+"
        print(title)
        print(cap)
        # Printing out the Fields on the menu
        for i in range(len(list_fields)):
            name = list_fields[i]
            spacing = max_field-len_of_fields[i]
            descriptions = self.field_info[name][1]
            print(name+(" ")*spacing+" : "+descriptions)
        print(cap+"\n")

    
    def species_finder(self):
        Y_species = []
        for field in self.fields:
            if self.regexp_species.search(field):
                 Y_species.append(field)
        Y_species = [re.sub(r"^Y\(", '', sp) for sp in Y_species]
        Y_species = [re.sub(r"\)$", '', sp) for sp in Y_species]
        Y_species.sort()

        # Get the length of the string of each species
        sp_lens = [len(sp) for sp in Y_species]
        # Amount of spaces padding so each string is the 
        # same length if max(sp_lens) = 3:
        # "CH4" stays the same "H" becomes "H  "
        pad = [l + (max(sp_lens) - l) for l in sp_lens]
        sp_padded = [f"{sp: <{p}}" for sp, p in zip(Y_species, pad)]
        sp_lines = [' '.join(sp_padded[i:i+8]) for i in range(0, len(sp_padded), 8)]
        line_lenght = len(sp_lines[0])
        cap = "+"+("-"*(line_lenght-2))+"+"
        title = "\n" + (line_lenght//3)*(" ")+"Species found in file:"
        # Printing out the Species on the menu
        print(title)
        print(cap)
        for l in sp_lines:
            print(l) 
        print(cap+"\n")





