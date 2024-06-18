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

class Menu(PlotfileCooker):
    """
    A class to read fields and species of a plotfile
    """
    # TODO: Mets la variable "field_info" ici
    # field_info = {...}
    # Elle va s'ajouter a chaque classe vu que c'est un attribut
    # qui change pas

    # TODO: Je pense pas c'est nécessaire de mettre l'option de pas show les
    # Variables ou les species ça rentre pas mal toujours dans une fenêtre de
    # terminal. Ça serait bien de pouvoir enlever les descriptions par contre
    def __init__(self, plt_file, variables=True, species=True, has_var=False):
        """
        Constructor for the Menu
        """
        # TODO: The plotfile is already stored in the pfile attribute
        # when you call super().__init__ (which instantiate PlotfileCooker)
        self.plt_file = plt_file
        self.variables = variables
        self.species = species
        self.has_var = has_var

        super().__init__(plt_file,)
        # TODO: Fais le __init__ header_only=True pour pas lire tout les
        # Cell_H juste pour avoir le nom des fields
        # Le best ça serait de pas créer de PlotfileCooker instance et
        # juste lire les lignes avec les fields dans plt00000/Header

        # Database in {Acronym: (Regular Expression, Definition)} format
        self.field_info = {"velocity": (r"^\w+_velocity$",
                                        "Velocity [m/s]"),
                            "density": (r"^density$",
                                        "Density [kg/m^3]"),
                            "rhoh": (r"^rhoh$",
                                     "Density * Specific Enthalpy [J / m^3]"),
                            "temp": (r"^temp$",
                                     "Temperature [K]"),
                            "RhoRT": (r"^RhoRT$",
                                      "Density * (R / M_bar) * Temperature [Pa]"),
                            "divu": (r"^divu$",
                                     "Divergence of the velocity field [1 / s]"),
                            "gradp": (r"^gradp+\w$",
                                      "Local pressure gradient [Pa / m]"),
                            "I_R": (r"^I_R\(.+\)$",
                                    "Species reaction rates [kg / s m^3]"),
                            "FunctCall": (r"^FunctCall$",
                                          "Number of function calls to chemistry [-]"),
                            "HeatRelease": (r"^HeatRelease$",
                                            "Heat release rate from chem. reactions [W / m^3]"),
                            "avg_pressure": (r"^avg_pressure$",
                                             "Cell-averaged pressure from the node-centers [Pa]"),
                            "mag_vort": (r"^mag_vort$",
                                         "Vorticity magnitude [1 / s]"),
                            "Y": (r"^Y\(.+\)$","Species mass fractions [-]"),
                            "mass_fractions": (r"^mass_fractions$",
                                               "Species mass fractions [-]"),
                            "X": (r"^X\(.+\)$","Species mole fractions [-]"),
                            "mole_fraction": (r"^mole_fraction$",
                                              "Species mole fractions [-]"),
                            "diffcoeff": (r"^diffcoeff$",
                                          "Mixture-averaged diffusion coefficients (mass) [m^2/s]"),
                            "lambda": (r"^lambda$",
                                       "Thermal diffusivity [W / m / K]"),
                            "viscosity": (r"^viscosity$",
                                          "Mixture viscosity [Pa-s]"),
                            "mixture_fraction": (r"^mixture_fraction$",
                                                 "Mixture fraction based on Bilger’s formulation [-]"),
                            "progress_variable": (r"^progress_variable$",
                                                  "User defined progress variable [-]"),
                            "vorticity": (r"^vorticity$",
                                          "Vortricity components [1 / s]"),
                            "Qcrit": (r"^Qcrit$","Q-Criterion [-]"),
                            "kinetic_energy": (r"^kinetic_energy$",
                                               "Kinetic energy as 0.5 * Rho * |U^2| [kg / m s^2]"),
                            "enstrophy": (r"^enstrophy$",
                                          "Enstrophy as 0.5 * Rho * |omega^2| [kg / m s^2]"),
                            #"rhominsumrhoY": (r"^rhominsumrhoY$",
                            #                  "Rho minus sum of rhoYs [units]"),
                            "volFrac": (r"^volFrac$",
                                         "Volume fraction at embedded boundaries [-]"),
                            "DistributionMap": (r"^DistributionMap$",
                                                "The MPI-rank of each box [units]"),
                            }
        self.regexp_species = re.compile(self.field_info["Y"][0])

        self.menu()

    def menu(self):
        # TODO: Remove stuff like this this is what git is for
        # Its still gonna be in the version history if you need it
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
			
            # TODO: Remove this
			#Could be different be *creative*
			#$____________________ Plotfile Menu ____________________$
			#HeatRelease       : Heat Release rate [W / s m^3]
			#progress_variable : User defined reaction progress
			#I_R               : Species reaction rates [kg / s m^3

        # TODO: Voir comment au dessus de __init__
        if self.variables:
            self.variables_finder()
        if self.species:
            self.species_finder()
    
    def show_variables(self, list_variables):
        # Let's print the Fields
        # Finding the maximum lenght of printed lines
        if self.has_var:
            len_of_fields = [len(key) for key in self.field_info]
            len_of_description = [len(self.field_info[key][1]) for key in self.field_info]
            yes_no = []
            for key in self.field_info:
                if key in list_variables:
                    yes_no.append("(Yes)")
                else:
                    yes_no.append("(No)")
        else:
            descriptions = []
            len_of_fields = [len(field) for field in list_variables]
            for field in list_variables:
                descriptions.append(self.field_info[field][1])
            len_of_description = [len(description) for description in descriptions]

        max_field = np.max(len_of_fields)  
        max_description = np.max(len_of_description)
        total_lenght = max_field + max_description + 2 
        if self.has_var:
            total_lenght += len("(yes)")
            title = "\n" + ((total_lenght//3)+2)*(" ")+"All known fields:"
        else:
            title = "\n" + (total_lenght//3)*(" ")+"Fields found in file:"
        cap = "+"+("-"*total_lenght)+"+"
        print(title)
        print(cap)
        # Printing out the Fields on the menu
        if self.has_var:
            print("Name"+" "*((max_field-len_of_fields[1])+3)+" "+"Present"+" "*2+"Description")
            print(cap)
            for i in range(len(list(self.field_info))):
                name = list(self.field_info)[i]
                spacing = (max_field-len_of_fields[i]) + 2 
                if yes_no[i] == "(yes)":
                    spacing -= 1
                description = self.field_info[name][1]
                print(name+(" ")*spacing+yes_no[i]+" : "+description)
        else:
            print("Name"+" "*((max_field-len_of_fields[1])+6)+"Description")
            print(cap)
            for i in range(len(list_variables)):
                name = list_variables[i]
                spacing = max_field-len_of_fields[i]
                description = self.field_info[name][1]
                print(name+(" ")*spacing+" : "+description)
        print(cap+"\n")

    def show_species(self, list_species):
        # Let's print the Species
        # Get the length of the string of each species
        sp_lens = [len(sp) for sp in list_species]
        # Amount of spaces padding so each string is the 
        # same length if max(sp_lens) = 3:
        # "CH4" stays the same "H" becomes "H  "
        pad = [l + (max(sp_lens) - l) for l in sp_lens]
        sp_padded = [f"{sp: <{p}}" for sp, p in zip(list_species, pad)]
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
        self.show_variables(list_fields)
    
    def species_finder(self):
        Y_species = []
        for field in self.fields:
            if self.regexp_species.search(field):
                 Y_species.append(field)
        Y_species = [re.sub(r"^Y\(", '', sp) for sp in Y_species]
        Y_species = [re.sub(r"\)$", '', sp) for sp in Y_species]
        Y_species.sort()
        self.show_species(Y_species)

        





