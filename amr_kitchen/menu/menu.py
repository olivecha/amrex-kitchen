import os
import sys
import time
import traceback
import multiprocessing
import pickle
import numpy as np
from tqdm import tqdm
from amr_kitchen import PlotfileCooker

# TODO: not sure this needs to be a class
# could be a single script parsing the arguments 
class Menu(PlotfileCooker):
    """
    A class to read and present fields of a plotfile
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

        self.menu()

    # TODO: I guess a class is okay if you do it like this
    # It'll prevent having a 200+ lines function which is
    # not ideal
    def menu(self):
        # TODO: you could use regular expressions to find
        # fields with a shared pattern (google it its like
        # simple but complicated at the same time)
        regexp_species = r"^Y\(.+\)$"
        # This matches:
        # - the starting "Y(" : ^Y\(   you have to escape (
        # - any characters like N2 : ".+"
        # - the tailling ) : \)$  with escaping ) also
        # Find where the regex matches
        matches = [re.search(regexp_species, f) for f in pck.fields]
        # Remove non matched fields and access the matched string
        mass_fracs = [rm.string for rm in matches if rm is not None]
        if len(mass_fracs) > 0:
            # This should be defined by a method from __init__
            self.has_mass_fracs = True
            # reuse the regex to remove the "Y("
            Y_species = [re.sub(r"^Y\(", '', sp) for sp in mass_fracs]
            # Same for the SP")"
            Y_species = [re.sub(r"\)$", '', sp) for sp in Y_species]
            # This is nice becomes sometimes species have the form:
            # Y(CHN(S)) so if you just remove all the '()' you
            # break the species string

        # TODO: the nice formatting I was talking about:
        # Check if there is a field defined for multiple species
        if self.has_mass_fracs:
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
            # Its gonna be some logic with % and // operators
            sp_lines = [' '.join(sp_padded[i:i+8]) for i in range(0, len(sp_padded), 8)]
            print("Species found in file:")
            for l in sp_lines:
                print(l)
			# This looks like:
			""" # You could sort by size of elements to make it nicer
			N2     H2     H      O2     O      H2O    OH     H2O2  
			HO2    CO     CO2    CH4    CH3    CH3O2H CH3O2  CH3O  
			CH2OH  CH2O   HCO    C2H6   C2H5   C2H4   C2H3   CH2CHO
			"""
			# TODO: Tu feras un truc similaire pour les fields infos
			# pour que ça fasse dequoi comme:
			"""
			Could be different be *creative*
			$____________________ Plotfile Menu ____________________$
			HeatRelease       : Heat Release rate [W / s m^3]
			progress_variable : User defined reaction progress
			I_R               : Species reaction rates [kg / s m^3
			"""

        if self.variables:
            self.variables_finder()
        if self.species:
            self.species_finder()
    
    def variables_finder(self):
        # Possible fields in the form "acronym" : "name"
        # TODO: you can put this just under the class definition
        # And its gonna be accessible as an attribute
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
		# TODO: also do a dict using the PeleLMeX documentation like:
		field_info = {"mag_vort": "Vortricity magnitude [units]"
		 			  "I_R": "Species reaction rates [kg / s m^3]"}

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
				# #TODO: Je parserait les fields comme ça:
				# 1. Diviser en single/multiple fields
				# 	- single tu fais juste print l'item dans le dict
				# 	  plus haut genre print(field, ':', field_info[field]) 
				# 	- multiple tu set un flag genre self.has_vel = True
				#     Tu pourrait search les regex de chaque species_field
				# 	  sur chaque field et set les flags des species fields:
				self.species_fields = set()	
				for field in self.fields: # Ça itère déjà sur les keys
					for sp_rgx sp_field_name in zip([r"^Y\(.+\)$", r"^I_R\(.+\)$", ...], # Définie sous class(...)
									  				["Y", "I_R", "D_I"]):
						if re.search(sp_rgx) is not None:
							# The set doesnt allow duplicates
							self.species_fields.add(sp_field_name)
				# self.species_fields = {'Y', 'I_R'} (genre)
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
		# TODO: Voir le comment sur le pretty formatting
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

    
    # TODO: not sure why you need to separate that
    def species_finder(self):
        pass
        #...

#plt_file = os.path.join(os.getcwd(),'plt00000')
#Menu(plt_file)
    
