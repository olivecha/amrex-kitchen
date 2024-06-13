import os
import sys
import time
import traceback
import multiprocessing
import pickle
import numpy as np
from tqdm import tqdm
from amr_kitchen import PlotfileCooker

class Menu(PlotfileCooker):
    """
    A class to read and present fields of a plotfile
    """

    def __init__(self, plt_file, variables=True, species=True,):
        """
        Constructor for the Menu
        """
        self.plt_file = plt_file
        self.variables = variables
        self.species = species

        super().__init__(plt_file,)

        self.menu()

    def menu(self):
        if self.variables:
            self.variables_finder()
        if self.species:
            self.species_finder()
    
    def variables_finder(self):
        # Possible fields in the form "acronym" : "name"
        possible_fields = ["velocity",
                           "density",
                           "rhoh",
                           "temp",
                           "RhoRT",
                           "divu",
                           "gradp",
                           "I_R",
                           "FunctCall",
                           "HeatRelease",
                           "avg_pressure",
                           "mag_vort",
                           "Y",               
                        ]

        list_fields = []
        for field in list(self.fields):
            for name in possible_fields:
                if name in field:
                    if name not in list_fields:
                        list_fields.append(name)
                        next
                    else:
                        pass
                # N'arrive encore pas à avoir la bonne logique 
                """elif name not in field:
                    if field not in list_fields:
                        print(f"'{field}' is not in the field database")
                        new_field = input(f"Do you want '{field}' added to the list of field ? (y/n): ").lower()
                        if new_field == "y":
                            list_fields.append(field)
                            print("Added!")
                        else:
                            pass
                    else:
                        pass"""
                
        # Printing out the fields on the menu
        separator = "-"
        reponse = "yes"
        num_separator = 20
        width = num_separator + len(reponse) + 2
        cap = separator*width
        print(f"+{cap}+")
        #print(f"|{space*round((width/2)-1)}MENU")
        for field in list_fields:
            print(f"|{field} {(20-len(field))*separator} {reponse}|")
        print(f"+{cap}+")

    
    def species_finder(self):
        ...

#plt_file = os.path.join(os.getcwd(),'plt00000')
#Menu(plt_file)
    
