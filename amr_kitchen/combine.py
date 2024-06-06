import os
import sys
import time
import multiprocessing
import pickle
import numpy as np
from tqdm import tqdm
from amr_kitchen import PlotfileCooker
from amr_kitchen.colander.colander import *

def combine(plotfile1, plotfile2, fields1=None, fields2=None, output=None):
		"""
		plotfile1: path to the first plotfile
		plotfile2: path to the second plotfile
		fields1: fields to keep in plotfile1 (defaults to all)
		fields2: fields to keep in plotfile2 (defaults to all)
		output: path of the new plotfile (defaults to both names jointed 
										  and in the same directory as plotfile1)
		"""
		# Handle the case when both plotfiles share some fields (don't keep both)
		# First step
		pck1, pck2 = PlotfileCooker(plotfile1), PlotfileCooker(plotfile2)
		assert pck1 == pck2

		# Output's file name 
		if output == None:
			path1,tail1 = os.path.split(plotfile1)
			tail2 = os.path.basename(plotfile2)

			output_path = os.path.join(path1,tail1+"-"+tail2)
		else:
			output_path = output

		# List of all combine fields in the two plot file without duplicates 
		output_fields = []

		# List of fields in the first plot file
		if fields1 == None:
			plt1_fields = list(pck1.fields)
		else:
			plt1_fields = fields1

		# List of fields in the second plot file
		if fields2 == None:
			plt2_fields = list(pck2.fields)
		else:
			plt2_fields = fields2

		# Creation of output_fields
		output_fields = plt1_fields.copy()
		for field in plt2_fields:
			if field in output_fields:
				pass 
			else:
				output_fields.append(field)

		"""
		Limit level, boxes, cell indexes are the same for 
		plotfile1 and plotfile2 (within a tolerance)
		"""                  
		# Limit level of output
		output_limit_level = pck1.limit_level

		output_boxes = []
		output_idx = []
		for lv in range(output_limit_level + 1):
			# Boxes of output
			level_boxes = []
			for box in pck1.boxes[lv]:
				level_boxes.append(box)
			for box in pck2.boxes[lv]:
				level_boxes.append(box)
			output_boxes.append(level_boxes)

			# Cell indexes of output
			level_idx = []
			for idx in pck1.cells[lv]['indexes']:
				level_idx.append(idx)
			for idx in pck2.cells[lv]['indexes']:
				level_idx.append(idx)
			output_idx.append(level_idx)
			



# J'ai ça (PlotfileCooker(pck2d).cells[0]['indexes']) pour avoir tous les indexes par niveau. 
#                                                     Comment savoir ils appartiennent à quel field ? 
# CE N'EST SUREMENT PAS GRÂVE 

# QUESTION IMPORTANTE :
# Est-ce que je fais en sorte que tout aille dans une nouvelle plofile ou est-ce que 
# je fais juste de nouveaux arrays avec les informations des 2 plotfiles combinées ?? 
		
		                                                         
"""pck2d = os.path.join(os.getcwd(),"test_assets","example_plt_2d")
pck3d = os.path.join(os.getcwd(),"test_assets","example_plt_3d")
print(combine(pck2d,pck2d))

"""

def main():
	cld = Colander(plotfile=os.path.join(os.getcwd(),"test_assets","example_plt_2d"),
                   limit_level=0,
                   output=os.path.join("C:\\","Users","franc","Desktop","Stage-POLY-2024","Data","test"),
                   variables=["temp"])

	cld.strain()

	'''cld = Colander(plotfile=os.path.join(os.getcwd(),"test_assets","example_plt_2d"),
                   limit_level=0,
                   output=os.path.join("C:\\","Users","franc","Desktop","Stage-POLY-2024","Data","test"),
                   variables=["mag_vort"])

	cld.strain()'''

if __name__ == "__main__":
    main()


