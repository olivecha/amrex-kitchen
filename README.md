# Welcome to the AMReX Kitchen

<img width="450" alt="image" src="https://github.com/olivecha/amrex-kitchen/assets/78630053/f2115bef-887a-4408-8d20-6be15bb0a4a4">

This is a collection of tools to perform post processing and visualization tasks 
on the AMReX plotfile format. For the moment all tools are written in python
and make use of the `multiprocessing` module when advantageous to do so. 
All tools are tested on large (> 1TB) plotfiles with multiple fields and AMR Levels.
The computationaly intensive tools were optimized to work best on a single node
of a CPU cluster with a large number of CPUs (40), and a decent amount of RAM (> 100 GB). 

A kitchen themed naming convention is employed as it helps remembering the
command names and conveys a positive sentiment. It also makes sense of the
colloquial phrase *"Let him cook"* and the Breaking Bad reference *"we need to
cook"* which is well suited to the analysis of reacting flows.

## Overview of the tools

- **mandoline:** Create (fast) slices of 3D plotfiles or covered grids of 2D plotfiles

- **colander:** Strains out fileds or levels from a plotfile to make the file-size more manageable.

- **spoon:** A quick taste, print out basic plotfile information. (work in progress)

- **pestle:** Grind and combine, volume integrals using a multilevel covering grid. (work in progress)

- **whip:** Increased volume, 3D covering grids and on-disk arrays for large
            plotfiles. (work in progress)

- **minuterie:** Its all about timing, prints out the time of the plotfile,
                 much faster than amrex/Tools/Plotfile/ftime.cpp. (work in progress)

- **pantry:** Lists the fields in the plotfile with human readable output (planned work)

- **chef:** Compute (cook) derived thermo chemical quantities using Cantera
            `SolutionArrays`. (work in progress)

- **HeaderData:** The base python class used by most tools to process the
                  plotfile data. Also makes available methods which provide
                  iterators over the data in the plotfile by each AMR level.
                  (will probably be renamed with a kitchen inspired name)

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

Install the requirements (most likely you already have them):
```
pip install -r requirements.txt
```

Then install the amrex-kitchen module:
```
pip install -e .
```
The `-e` flag makes the source files editable so reinstalling after each `git pull` is not necessary.

# Documentation

## mandoline

<img width="350" alt="image" src="https://github.com/olivecha/mandoline/assets/78630053/857636f4-e49d-41d2-b428-6e12b6874157">

Fast slices of (large) AMReX plotfiles. This is equivalent to
`amrex/Tools/Plotfile/fsnapshot.cpp`, and retains the command line arguments names, 
but it only loads the data needed for the slice, so large plotfiles can be
visualized rapidly. This tools also supports multiple output formats (`.png`
image, `.npy` array and 2D AMReX plotfiles). Slicing a 2D plotfile will simply
create a uniform grid for the data, and put it in the output format. If the
`plotfile` format is used with a 2D plotfile the binary files will be rewritten
with the requested fields in the output file.


### Usage

```
usage: mandoline [-h] [--normal NORMAL] [--position POSITION]
                 [--variable VARIABLE] [--max_level MAX_LEVEL]
                 [--format FORMAT] [--output OUTPUT] [--colormap COLORMAP]
                 [--minimum MINIMUM] [--maximum MAXIMUM] [--log]
                 [--verbose VERBOSE]
                 plotfile

Fast approximate slices of AMReX plotfiles

positional arguments:
  plotfile              Path of the plotfile to slice

optional arguments:
  -h, --help            show this help message and exit
  --normal NORMAL, -n NORMAL
                        Coordinate normal to the slice x:0, y:1, z:2, defaults
                        to x.
  --position POSITION, -p POSITION
                        position of the slice, defaults to mid plane.
  --variable VARIABLE, -v VARIABLE
                        variable names, defaults to "density", "all" keeps all
                        the variables in the plotfile
  --max_level MAX_LEVEL, -L MAX_LEVEL
                        Maximum AMR level loaded, defaults to finest level.
  --format FORMAT, -f FORMAT
                        Either image, array or plotfile 
                        image: creates and saves an image using matplotlib
                        array: creates a numpy array
                        with am uniform grid and saves it with the numpy
                        format (can be loaded with `np.load`).
                        plotfile: creates a 2D plotfile at the slicing plane
                        with the requested fields. This is useful to process a
                        single data slice with AMReX compatible tools
  --output OUTPUT, -o OUTPUT
                        File name used to override the default value
  --verbose VERBOSE, -V VERBOSE
                        Verbosity level, defaults to 1

  #### Arguments for image format output ####

  --colormap COLORMAP, -c COLORMAP
                        A named matplotlib colormap, defaults to jet
  --minimum MINIMUM, -m MINIMUM
                        Minimum value used in the colormap
  --maximum MAXIMUM, -M MAXIMUM
                        Maximum value used in the colormap
  --log, -l             Flag to use log scale in the plotted image
```
### How it works

AMReX blocks intersecting the slicing plane are read from the plotfile. Arrays
with the closests points at either side of the plane are constructed for each
level and a linear interpolation is performed to compute the data in the slice. 

Support for non-orthogonal slicing planes, and downsampling of higher level
data to increase speed and reduce output size is planned.

**Known issue:** when the slice position is between the last point of a higher
level box and the fist point of a lower level box it is possible that the region
will only be covered with lower level data. Moving the slice by the grid
resolution of the highest level is a temporary fix. 


## 2. Colander

<img width="300" alt="colander" src="https://github.com/olivecha/amrex-kitchen/assets/78630053/aec452e7-520e-4f2c-bbb2-44c79dd0c6ca">

Strain out variables or levels from plotfiles, see:

```
$ colander --help
```



