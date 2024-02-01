# mandoline
A command line tool to create fast slices of (large) AMReX plotfiles.

<img width="350" alt="image" src="https://github.com/olivecha/mandoline/assets/78630053/857636f4-e49d-41d2-b428-6e12b6874157">

## Instalation

First install the requirements (most likely you already have them):
```
pip install -r requirements.txt
```
Then install the mandoline executable:
```
pip install -e .
```
The `-e` flag makes the source files editable so reinstalling after each `git pull` is not necessary.

## Usage

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
                        Coordinate normal to the slice x:0, y:1, z:2
  --position POSITION, -p POSITION
                        position of the slice, defaults to mid plane
  --variable VARIABLE, -v VARIABLE
                        variable name, defaults to "density"
  --max_level MAX_LEVEL, -L MAX_LEVEL
                        Maximum AMR level loaded, defaults to finest level
  --format FORMAT, -f FORMAT
                        Either image or array image: creates and saves an
                        image using matplotlib array: creates a numpy array
                        with a 500x500 uniform grid and saves it as a pickle
                        file of a python dict
  --output OUTPUT, -o OUTPUT
                        File name used to override the default value
  --colormap COLORMAP, -c COLORMAP
                        A named matplotlib colormap, defaults to jet
  --minimum MINIMUM, -m MINIMUM
                        Minimum value used in the colormap
  --maximum MAXIMUM, -M MAXIMUM
                        Maximum value used in the colormap
  --log, -l             Flag to use log scale in the plotted image
  --verbose VERBOSE, -V VERBOSE
                        Verbosity level, defaults to 1
```
## How it works

1. The block information is read from the plotfile header
2. For each level up to the specified max level boxes intersecting with the plane are read, but only the points closest to either side of the slicing plane are kept, along with their distance from the plane.
3. Empty arrays at the maximum level are created for both side of the plane.
4. The array are populated with the slice data for each grid level, overwriting lower level data with higher level data.
5. A linear interpolation is performed between the arrays representing both sides of the slicing plane.
6. The final data is either plotted or saved.


