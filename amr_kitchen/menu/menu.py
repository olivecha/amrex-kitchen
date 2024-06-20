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
    # Database in {Acronym: (Regular Expression, Definition, Units)} format
    field_info = {"avg_pressure": (r"^avg_pressure$",
                                    "Cell-averaged pressure from the node-centers ", "[Pa]"),
                "density": (r"^density$",
                             "Density ", "[kg/m^3]"),
                "diffcoeff": (r"^D_+.",
                               "Mixture-averaged diffusion coefficients (mass) ", "[m^2/s]"),
                "DistributionMap": (r"^DistributionMap$",
                                     "The MPI-rank of each box ", "[-]"),
                "divu": (r"^divu$",
                          "Divergence of the velocity field ", "[1 / s]"),
                "enstrophy": (r"^enstrophy$",
                               "Enstrophy as 0.5 * Rho * |omega^2| ", "[kg / m s^2]"),
                "FunctCall": (r"^FunctCall$",
                               "Number of function calls to chemistry ", "[-]"),
                "gradp": (r"^gradp+\w$",
                           "Local pressure gradient ", "[Pa / m]"),
                "HeatRelease": (r"^HeatRelease$",
                                 "Heat release rate from chem. reactions ", "[W / m^3]"),
                "I_R": (r"^I_R\(.+\)$",
                         "Species reaction rates ", "[kg / s m^3]"),
                "kinetic_energy": (r"^kinetic_energy$",
                                    "Kinetic energy as 0.5 * Rho * |U^2| ", "[kg / m s^2]"),
                "lambda": (r"^lambda$",
                            "Thermal diffusivity ", "[W / m / K]"),
                "mag_vort": (r"^mag_vort$",
                              "Vorticity magnitude ", "[1 / s]"),
                "mass_fractions": (r"^mass_fractions$",
                                    "Species mass fractions ", "[-]"),
                "mixture_fraction": (r"^mixture_fraction$",
                                      "Mixture fraction based on Bilger’s formulation ", "[-]"),
                "mole_fraction": (r"^mole_fraction$",
                                   "Species mole fractions ", "[-]"),
                "progress_variable": (r"^progress_variable$",
                                       "User defined progress variable ", "[-]"),
                "Qcrit": (r"^Qcrit$",
                           "Q-Criterion ", "[-]"),
                "rhoh": (r"^rhoh$",
                          "Density * Specific Enthalpy ", "[J / m^3]"),
                "RhoRT": (r"^RhoRT$",
                           "Density * (R / M_bar) * Temperature ", "[Pa]"),
                "temp": (r"^temp$",
                          "Temperature ", "[K]"),
                "velocity": (r"^\w+_velocity$",
                              "Velocity ", "[m/s]"),
                "viscosity": (r"^viscosity$",
                               "Mixture viscosity ", "[Pa-s]"),
                "volFrac": (r"^volFrac$",
                             "Volume fraction at embedded boundaries ", "[-]"),
                "vorticity": (r"^vorticity$",
                               "Vortricity components ", "[1 / s]"),
                "X": (r"^X\(.+\)$",
                       "Species mole fractions ", "[-]"),
                "Y": (r"^Y\(.+\)$",
                       "Species mass fractions ", "[-]"),}

    def __init__(self, plt_file, has_var=None, every=False, description=False, min_max=False, finest_lv=False):
        """
        Constructor for the Menu\n
        +--------------------------------------------------------------------+\n
        plt_file --> Path of the plotfile to analyse\n
        has_var --> List of variables to find in the plotfile\n
        every --> Flag to enable displaying every field in the database\n
                        and if they are the plotfile or not\n
        description --> Flag to enable displaying the fields' descriptions\n
        min_max --> Flag to enable displaying the fields' absolute\n
                              min and max throughout every level\n
        finest_lv --> Flag to enable displaying the fields' min and max\n 
                            at the finest level\n
        +--------------------------------------------------------------------+\n
        """
        self.to_find = has_var
        self.every = every 
        self.description = description
        self.min_max = min_max
        self.finest_lv = finest_lv

        if self.min_max or self.finest_lv:
            super().__init__(plt_file,
                             maxmins=True,)
        else:
            self.plt_file = plt_file
            # Let's find the plotfile's fields
            filepath = os.path.join(self.plt_file, 'Header')
            with open(filepath) as hfile:
                self.version = hfile.readline()
                # field names
                self.nvars = int(hfile.readline())
                self.fields = {}
                for i in range(self.nvars):
                    self.fields[hfile.readline().replace('\n', '')] = i

        self.regexp_species = re.compile(self.field_info["Y"][0])
        self.menu()

    def menu(self):
        variables = self.variables_finder()
        if self.min_max or self.finest_lv:
            min_max_data = self.find_min_max() 
            self.show_min_max(min_max_data)
        if self.to_find:
            self.finding(variables) 
        if self.description or self.every:
            self.show_variables_description(variables)
        else:
            self.show_variables(variables)
        species = self.species_finder()
        self.show_species(species)
    
    def find_min_max(self):
        min_and_max = {}
        for field in self.fields:
            # Let's find the units 
            units = "[...]"
            for key in self.field_info:
                regexp = re.compile(self.field_info[key][0])
                if regexp.search(field):
                    units = self.field_info[key][2]
                    break
            if self.finest_lv:
                # Finest min and max 
                minimum = np.min(self.cells[self.limit_level]["mins"][field])
                maximum = np.max(self.cells[self.limit_level]["maxs"][field])
            else:
                # Absolute min and max 
                minimum = np.min([self.cells[lv]["mins"][field].min() for lv in range(self.limit_level + 1)]) 
                maximum = np.max([self.cells[lv]["maxs"][field].max() for lv in range(self.limit_level + 1)])
            min_and_max[field] = (str("{:.3e}".format(minimum)), str("{:.3e}".format(maximum)), units)
        return min_and_max
    
    def show_min_max(self, min_max_data):
        # To have two equal columns
        if not len(min_max_data)//2:
            min_max_data[""] = ("","","")
        middle = int(len(min_max_data)/2)
        # Every lenghts of both columns 
        all_lenght = [(len(field),
                       len(min_max_data[field][0]),
                       len(min_max_data[field][1]),
                       len(min_max_data[field][2])) for field in list(min_max_data)]
        # For both columns:
        max_name = np.max([all_lenght[i][0] for i in range(len(all_lenght))])
        max_min = np.max([all_lenght[i][1] for i in range(len(all_lenght))])
        max_max = np.max([all_lenght[i][2] for i in range(len(all_lenght))])
        max_unit = np.max([all_lenght[i][3] for i in range(len(all_lenght))])
        max_lenght = max_name + max_min + max_max + max_unit + 3 
        
        cap = "+"+"-"*max_lenght+"+\t+"+"-"*max_lenght+"+"
        title = " "*((len(cap)//2)-8)+"Fields' Mins and Maxs:"
        header1 = "Fields"+" "*(max_name-3)+"Min."+" "*(max_min-3)+"Max."+" "*(max_max-3)+"Units"
        header2 = " "*(max_unit-5)+"\tFields"+" "*(max_name-3)+"Min."+" "*(max_min-3)+"Max."+" "*(max_max-3)+"Units"
        print(title)
        print(cap)
        print(header1+header2)
        print(cap)
        # Printing the mins and maxs of the plotfile
        for i in range(middle):
            field1, field2 = list(min_max_data)[i], list(min_max_data)[i+middle]
            mini1, maxi1, units1 = min_max_data[field1]
            mini2, maxi2, units2 = min_max_data[field2]
            spacing_name1, spacing_name2 = " "*(max_name-all_lenght[i][0]), " "*(max_name-all_lenght[i+middle][0])
            spacing_min1, spacing_min2 = " "*(max_min-all_lenght[i][1]), " "*(max_min-all_lenght[i+middle][1])
            spacing_max1, spacing_max2 = " "*(max_max-all_lenght[i][2]), " "*(max_max-all_lenght[i+middle][2])
            spacing_units = " "*(max_unit-all_lenght[i][3])

            print((f"{field1}{spacing_name1} : {mini1}{spacing_min1} {maxi1}{spacing_max1} {units1}{spacing_units}\t"
                   f"{field2}{spacing_name2} : {mini2}{spacing_min2} {maxi2}{spacing_max2} {units2}"))
        print(cap)
        print("\n")

    def finding(self,list_variables):
        lines = []
        not_found = False
        for to_find in self.to_find:
            if to_find in list_variables:
                lines.append(f"{to_find} found!")
            else:
                lines.append(f"{to_find} not found!")
                not_found = True
        max_lenght = np.max([len(line) for line in lines])
        title = " "*((max_lenght//2)-7)+"Search results:"
        cap = "+"+("-"*(max_lenght-2)+"+")
        print(title)
        print(cap)
        for line in lines:
            print(line)
        if not_found:
            print(cap)
            print(("WARNING! Perharps fields were mispelled...\n"
                   "Please refer to the database by using flag --every -e\n"))
        else:
            print(cap+"\n")

    def show_variables_description(self, list_variables):
        # Let's print the Fields
        # Finding the maximum lenght of printed lines
        if self.every:
            len_of_fields = [len(key) for key in self.field_info]
            len_of_description = [len(self.field_info[key][1]+self.field_info[key][2]) for key in self.field_info]
            yes_no = []
            for key in self.field_info:
                if key in list_variables:
                    yes_no.append("Yes")
                else:
                    yes_no.append(" No")
        else:
            descriptions = []
            len_of_fields = [len(field) for field in list_variables]
            for field in list_variables:
                descriptions.append(self.field_info[field][1]+self.field_info[field][2])
            len_of_description = [len(description) for description in descriptions]

        max_field = np.max(len_of_fields)  
        max_description = np.max(len_of_description)
        total_lenght = max_field + max_description + 2 
        if self.every:
            total_lenght += 3
            title = ((total_lenght//3)+2)*(" ")+"All known fields:"
        else:
            title = (total_lenght//3)*(" ")+"Fields found in file:"
        cap = "+"+("-"*total_lenght)+"+"
        print(title)
        print(cap)
        # Printing out the Fields on the menu
        if self.every:
            print("Name"+" "*((max_field-len_of_fields[1])+2)+" "+"Present"+" Description")
            print(cap)
            for i in range(len(list(self.field_info))):
                name = list(self.field_info)[i]
                spacing = (max_field-len_of_fields[i]) + 2
                description = self.field_info[name][1]+self.field_info[name][2]
                # TODO: you could use the space padding formatter like for the species: {name:<{spacing}}
                print(name+(" ")*spacing+yes_no[i]+" : "+description)
        else:
            print("Name"+" "*(max_field-len("Name"))+" "*3+"Description")
            print(cap)
            for i in range(len(list_variables)):
                name = list_variables[i]
                spacing = max_field-len_of_fields[i]
                description = self.field_info[name][1]+self.field_info[name][2]
                print(name+(" ")*spacing+" : "+description)
        print(cap+"\n")

    def show_species(self, list_species):
        # Let's print the Species
        # Get the length of the string of each species
        sp_lens = [len(sp) for sp in list_species]
        # Amount of spaces padding so each string has the same lenght
        pad = [l + (max(sp_lens) - l) for l in sp_lens]
        sp_padded = [f"{sp: <{p}}" for sp, p in zip(list_species, pad)]
        sp_lines = [' '.join(sp_padded[i:i+8]) for i in range(0, len(sp_padded), 8)]
        line_lenght = len(sp_lines[0])
        cap = "+"+("-"*(line_lenght-2))+"+"
        title = (line_lenght//3)*(" ")+"Species found in file:"
        # Printing out the Species on the menu
        print(title)
        print(cap)
        for l in sp_lines:
            print(l) 
        print(cap+"\n")

    def show_variables(self, list_variables):
        # Let's print the Variables
        # Get the length of the string of each variable
        sp_lens = [len(sp) for sp in list_variables]
        # Amount of spaces padding so each string has the same lenght 
        pad = [l + (max(sp_lens) - l) for l in sp_lens]
        sp_padded = [f"{sp: <{p}}" for sp, p in zip(list_variables, pad)]
        field_per_line = 4
        sp_lines = [' '.join(sp_padded[i:i+field_per_line]) for i in range(0, len(sp_padded), field_per_line)]
        line_lenght = len(sp_lines[0]) 
        cap = "+"+("-"*(line_lenght-2))+"+"
        title = ((line_lenght//2)-10)*(" ")+"Fields found in file:"
        # Printing out the variables on the menu
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
                    self.field_info[field] = (field, "Units unknown", "[...]")
        list_fields.sort(key=str.lower)
        return list_fields
    
    def species_finder(self):
        Y_species = []
        for field in self.fields:
            if self.regexp_species.search(field):
                 Y_species.append(field) 
        # Et s'il n'y a aucun Field Y(...) dans la plotfile ?? 
        Y_species = [re.sub(r"^Y\(", '', sp) for sp in Y_species]
        Y_species = [re.sub(r"\)$", '', sp) for sp in Y_species]
        Y_species.sort()
        return Y_species
