# Welcome to the AMReX Kitchen

<img width="450" alt="image" src="https://github.com/olivecha/amrex-kitchen/assets/78630053/f2115bef-887a-4408-8d20-6be15bb0a4a4">

This is a collection of tools to perform post processing and visualization tasks 
on the AMReX plotfile format. For the moment all tools are written in python
and make use of the `multiprocessing` module when advantageous to do so. 

A kitchen themed naming convention is employed as it helps remembering the
command names and conveys a positive sentiment. It also makes sense of the
colloquial phrase *"Let him cook"* and the Breaking Bad reference *"we need to
cook"* which is well suited to the analysis of reacting flows.

## Overview of the tools

- **chef:** Compute (cook) derived thermochemical quantities (recipes) using Cantera
            `SolutionArrays`.

- **colander:** Strains out fileds or levels from a plotfile to make the file-size more manageable.

- **combine:** Combine fields from two plotfile with the same geometric mesh data.

- **mandoline:** Create (fast) slices of 3D plotfiles or covered grids of 2D plotfiles

- **menu:** Print out field information and min/max values by only reading the plotfile headers

- **pestle:** Grind and combine, fast and low memory volume integrals for multi-level plotfiles.

- **taste:** Parse the plotfile data to ensure it is not corrupted and that no files are missing.

- **whip:** Increased volume, 3D uniform grids from adaptive mesh refinement data.
            Smaller data types are supported (`float32`, ...) to reduce memory usage

- **minuterie:** Its all about timing, prints out the time of the plotfile,
                 much faster than amrex/Tools/Plotfile/ftime.cpp. (work in progress)

- **PlotfileCooker:** The base python class used by most tools to process the
                  plotfile data. Also makes available methods which provide
                  iterators over the data in the plotfile by each AMR level.
  
- **marinate:** Save the `PlotfileCooker` instance to binary using the python `pickle` library.
                This allows exploring the adaptive mesh refinement without needing the whole plotfile. 
  

## Installation

The python pakage will eventually be distributed with pip, but in the mean time
it can be installed manually by cloning the repository:

Clone the repository
```
git clone https://github.com/olivecha/amrex-kitchen.git
```
Go into it
```
cd amrex-kitchen
```

(Optional) Activate your python virtual or conda environment.

Install the requirements:
```
pip install -r requirements.txt
```

Then install the amrex-kitchen module:
```
pip install -e .
```
The `-e` flag makes the source files editable so reinstalling after each `git pull` is not necessary.

## Performance

All tools are tested on large (> 1TB) plotfiles with around 40 fields and 5 AMR Levels.
The computationaly intensive tools were optimized to work best on a single node
of a CPU cluster (Niagara) with a large number of CPUs (40), and a large amount of RAM (> 100 GB).

Compared to the tools proposed in the AMReX repository, what is proposed here is arguably a bit slower when a large amount of computing
power is available, as the `multiprocessing` module is limited to a single compute node. However, the tools proposed here work well with 
arbitrarly large plotfiles, and scale linearly with plotfile size (with some exceptions).

# Documentation

## Chef

The **chef** has entered the Amrex kitchen: with this tool you can apply **recipes** to plotfiles.

<img width="300" alt="image" src="https://github.com/user-attachments/assets/d0f20c96-a34a-4118-80d8-a75dacd468a5">

### Usage:

There are predefined recipes like Enthalpy and HeatRelease:
```
chef --outdir plthrr --recipe HRR --mech mechanism.yaml --pressure 1 plt00000
```
You can apply recipes on selected species, like the reaction rate of H2 and O2:
```
chef --outdir pltomega_sp --recipe SRi --species H2 O2 --mech mechanism.yaml --pressure 1 plt00000
```
Or selected individual reactions like the net rate of progress:
```
chef --outdir pltnetrxrate --recipe RRi --reactions 0 1 2 3 --mech mechanism.yaml --pressure 1 plt00000
```
But you can also be the chef and define your own recipes in an arbitrary `.py` file and pass it as an argument. 
The function must take two or three arguments. if it takes three they are a dictionnary containing the indexes
of the fields in the plotfile, a box array with shape (nx, ny, nz, nfields) and a Cantera SolutionArray with shape
(nx, ny, nz). In the case of two arguments the SolutionArray is omitted which is faster. The docstring of the 
function is taken as the field name in the plotfile.

For example using the SolutionArray to access the reactions rates and compute the ratio between the production
rates of two species:
```python
# my_recipe.py

def recipe(field_indexes, box_array, sol_array):
    """
    fuel_oxy_omega_ratio
    """
    id_H2 = sol_array.species_index('H2')
    id_O2 = sol_array.species_index('O2')
    omega_H2 = sol_array.net_production_rates[:, :, :, id_H2]
    omega_O2 = sol_array.net_production_rates[:, :, :, id_O2]
    return omega_H2/omega_O2
```
Then using the following command:
```
$ chef --output pltrxratio --recipe my_recipe.py --mech mechanism.yaml --pressure 1 plt00000
```

Or computing the difference between the temperature and an arbitrary value wich is function of 
the progress variable:
```
# my_recipe2.py

T_vs_C = a_user_defined_function...

def recipe(field_indexes, box_array):
    """
    T_minus_T_C
    """
    temp = box_array[:, :, :, field_indexes["temp"]]
    C = box_array[:, :, :, field_indexes["progress_variable"]]
    return temp - T_vs_C(C)
```
Which is computed using a similar command (without needing mechanism or pressure inputs):
```
$ chef --output pltTdiff --recipe my_recipe2.py plt00000
```

All available options and flags are as follows:
```
usage: chef [-h] [--outdir OUTDIR] [--recipe RECIPE] [--species SPECIES [SPECIES ...]] [--reactions REACTIONS [REACTIONS ...]]
            [--mech MECH] [--pressure PRESSURE]
            plotfile

positional arguments:
  plotfile              Path of the plotfile to cook

options:
  -h, --help            show this help message and exit
  --outdput, -o OUTDIR   Output path of the generated plotfile (defaults to input_plotfile'_ck'
  --recipe, -r RECIPE   Recipe used to add new data, available recipes are: - 'HRR': Heat release rate [W/m^3] - 'ENT': Sensible
                        enthalpy [] - 'SRi': Reaction rate of selected species (must specify with --species argument) - 'SDi':
                        Diffusion coefficient of selected species (must specify with --species argument) - 'RRi' : Rate of
                        progress of selected reactions (must specify with the --reactions argument) - user defined: the path to a
                        python file with a function named 'recipe' in it
  --species, -s SPECIES [SPECIES ...]
                        name of the species for which the recipe is prepared
  --reactions, -R REACTIONS [REACTIONS ...]
                        index of the reactions for which the recipe is applied
  --mech, -m MECH       Path to the Cantera mechanism to use for the derived variables computations
  --pressure, -p PRESSURE  Pressure of the simulation data (atm)
```

## Colander

<img width="300" alt="colander" src="https://github.com/olivecha/amrex-kitchen/assets/78630053/aec452e7-520e-4f2c-bbb2-44c79dd0c6ca">

Strain out variables or levels from plotfiles:

### Usage:

```
usage: colander [-h] [--variables VARIABLES [VARIABLES ...]] [--limit_level LIMIT_LEVEL] [--serial] [--output OUTPUT] plotfile

Filter out field and remove levels from AMReX plotfiles

positional arguments:
  plotfile              Path of the plotfile to filter

options:
  -h, --help            show this help message and exit
  --variables, -v VARIABLES [VARIABLES ...]   Variable to keep in the filtered plotfile
                        
  --limit_level, -l LIMIT_LEVEL   Maximum AMR Level to keep in the filtered plotfile
  --serial, -s          Flag to disable multiprocessing
  --output, -o OUTPUT   Output path to store the filtered plotfile
```

The output field is mandatory, keeping only the first adaptive mesh refinement level would yield the command:
```
colander --output plt_level0 --variables all --limit_level 0 plt50000
```

Creating a plotfile of the temperature field with the finest mesh refinement level:
```
colander --variables temp --output plt_temp_fine  plt50000
```

## Combine

Combine fields from two plotfiles into a new one. This can be usefull to add back a field produced by `Chef` into the original plotfile.

### Usage:
```
usage: combine [-h] [--plotfile1 PLOTFILE1] [--plotfile2 PLOTFILE2] [--vars1 VARS1] [--vars2 VARS2] [--output OUTPUT] [--serial]

Combine two AMReX plotfiles

options:
  -h, --help            show this help message and exit
  --plotfile1, -p1 PLOTFILE1   Path of the first plotfile to combine
  --plotfile2, -p2 PLOTFILE2   Path of the second plotfile to combine
  --vars1, -v1 VARS1    Comma or space separated variables between quotes to keep in the first plotfile
  --vars2, -v2 VARS2    Comma or space separated variables between quotes to keep in the second plotfile
  --output, -o OUTPUT   Output path to store the combined plotfile
  --serial, -s          Flag to disable multiprocessing
```

## Mandoline

<img width="350" alt="image" src="https://github.com/olivecha/mandoline/assets/78630053/857636f4-e49d-41d2-b428-6e12b6874157">

Fast slices of (large) AMReX plotfiles. This is equivalent to
`amrex/Tools/Plotfile/fsnapshot.cpp`, and retains the command line arguments names when available, 
but it only loads the data needed for the slice, so large plotfiles can be
visualized rapidly. This tools also supports multiple output formats (`.png`
image, `.npy` array and 2D AMReX plotfiles). Slicing a 2D plotfile will simply
create a uniform grid for the data, and put it in the output format. If the
`plotfile` format is used with a 2D plotfile the binary files will be rewritten
with the requested fields in the output file.


### Usage

```
usage: mandoline [-h] [--normal NORMAL] [--position POSITION] [--variables VARIABLES [VARIABLES ...]] [--max_level MAX_LEVEL]
                 [--serial] [--format FORMAT] [--output OUTPUT] [--colormap COLORMAP] [--minimum MINIMUM] [--maximum MAXIMUM]
                 [--log] [--verbose VERBOSE]
                 plotfile

Fast slices of (large) AMReX plotfiles

positional arguments:
  plotfile              Path of the plotfile to slice

options:
  -h, --help            show this help message and exit

  --normal, -n NORMAL   Index of the coordinate normal to the slice x:0, y:1, z:2

  --position, -p POSITION   position of the slice in mesh coordinates, defaults to domain center.

  --variables, -v VARIABLES [VARIABLES ...]
                        variables names to slice, defaults to "density", "all" slices all the fields in the plotfile, "grid_level"
                        outputs the AMR grid level data.

  --max_level, -L MAX_LEVEL   Maximum AMR level loaded, defaults to finest level.

  --serial, -s          Flag to disable multiprocessing, this is usefull when processing many small plotfiles in a bash loop.

  --format, -f FORMAT   Either "image", "array" or "plotfile". image: creates and saves an image using matplotlib. array: creates
                        a numpy array with a uniform grid at the resolution of --max_level and saves it as a numpy file with
                        separated fields. plotfile: Keep the adaptive mesh refinement information and save the slice and specified
                        fields in a 2D amrex plotfile

  --output, -o OUTPUT   Output file name used to override the default value.

  --verbose, -V VERBOSE   Verbosity level, defaults to 1

 # Image output flags
  --colormap, -c COLORMAP   A named matplotlib colormap, defaults to jet.
  --minimum, -m MINIMUM     Minimum value used in the colormap
  --maximum, -M MAXIMUM.    Maximum value used in the colormap
  --log, -l                 Flag to use log scale in the plotted image
```
### More information

AMReX blocks intersecting the slicing plane are read from the plotfile. Arrays
with the closests points at either side of the plane are constructed for each
level and a linear interpolation is performed to compute the data in the slice. 

Support for non-orthogonal slicing planes, and downsampling of higher level
data to increase speed and reduce output size is planned.

**Known issue:** when the slice position is between the last point of a higher
level box and the fist point of a lower level box it is possible that the region
will only be covered with lower level data. Moving the slice by the grid
resolution of the highest level is a temporary fix.

## Menu

Display information about what fields are in a plotfile and their units, descriptions and min/max values. Units are taken from the PeleLMeX code.

```
usage: menu [-h] [--has_var HAS_VAR] [--every] [--description] [--min_max] [--finest_lv] plotfile

Displays fields and species of a plotfile

positional arguments:
  plotfile              Path of the plotfile to read

options:
  -h, --help            show this help message and exit
  --every, -e           Flag to enable displaying every field in the database and if they are the plotfile or not
  --description, -d     Flag to enable displaying the fields' descriptions
  --min_max, -m         Flag to enable displaying the fields' absolute min and max throughout every level
  --finest_lv, -f       Flag to enable displaying the fields' min and max only at the finest level
```

## Pestle

Volume integrals of multilevel plotfile. This tool was profiled againts the `fvolumesum.cpp` script in AMReX tools and performed better with the same number of cores. Additionnaly, the `volFrac` field defined for embedded boundaries in PeleLmeX can be taken into acount to produce more accurate integrals.

```
usage: pestle [-h] [--variable VARIABLE] [--limit_level LIMIT_LEVEL] [--volfrac] plotfile

Prints the volume integral of the chosen field in a plotfile

positional arguments:
  plotfile              Path of the plotfile to integrate

options:
  -h, --help            show this help message and exit
  --variable, -v VARIABLE Variable to integrate
  --limit_level, -l LIMIT_LEVEL Maximum AMR Level considered
  --volfrac, -vf        Use the volFrac field to obtain more accurate integrals for plotfiles with embedded boundaries. The
                        contribution of a finite volume to the integral is taken as: value * dV * volFrac.
```

## Taste

Validate the sanity of AMReX plotfiles. This tool was created to help validate the other tools as they were developped. Also, the tool first check if the file hierarchy of the plotfile is complete, which is usefull to check if a plotfile was purged by computing clusters temporary disk space policies.

```
usage: taste [-h] [--limit_level LIMIT_LEVEL] [--no_bin_headers] [--no_bin_shape] [--bin_data] [--box_coords] [--nofail]
             [--verbose VERBOSE]
             plotfile

positional arguments:
  plotfile              Path of the plotfile to validate

options:
  -h, --help            show this help message and exit
  --limit_level, -l LIMIT_LEVEL
                        Limit level at which the validation will occur
  --no_bin_headers, -nh
                        Do not validate the indexes in the binary headers
  --no_bin_shape, -ns   Do not validate the shape of the binary data
  --bin_data, -bd       Validate the binary data against the max/mins in the level headers and warn if NaNs are encountered (This
                        can take a while)
  --box_coords, -bc     Validate that the boxes coordinates in the plotfile header are consistent with the indices in the level
                        headers
  --nofail, -nf         Don't fail on the first error and perform all the checks
  --verbose, -v VERBOSE
                        Verbosity level (Defaults to 1)
```

## Whip

Create a uniform 3D array of a plotfile field and save it to the compressed numpy format. Lower level data is broadcasted to the finest AMR Level. A prompt requires the user to validate if the expected size of the uniform grid can fit into memory. 

```
usage: whip [-h] [--variable VARIABLE] [--limit_level LIMIT_LEVEL] [--outfile OUTFILE] [--dtype DTYPE] [--nochecks] plotfile

positional arguments:
  plotfile              Path of the plotfile used to make the uniform grid

options:
  -h, --help                      show this help message and exit
  --variable, -v VARIABLE         Variable used to create the uniform grid
  --limit_level, -l LIMIT_LEVEL   Maximum AMR Level considered
  --outfile, -o OUTFILE           Output file to override the default
  --dtype, -d DTYPE               Data type used (defaults to float64), but float32 or integer types can be used to save space (this is
                                  slower as the data has to be recasted by numpy)
  --nochecks, -y                  Do not prompt if the required memory is acceptable
```

## Minuterie 

This tool simply displays the time of a plotfile. This is useful to create timeseries plot on the fly without reading the whole plotfile by piping the output from multiple plotfiles into a text file. A minimal amount of data is read from the header to obtain the plotfile time.

```
$ minuterie test_assets/example_plt_2d
Plotfile time = 0.4994722514455662
```

## Marinate

Saves the PlotfileCooker object of a plotfile to a pickle file. This allows working with the `PlotfileCooker` object of very large plotfiles (> 1 Tb) on a desktop computer, as the object should take less than 1 Gb. The pickle file has the name of the plotfile directory with the extension `.pkl`.

