# mandoline

A command line tool to create fast slices of (large) AMReX plotfiles.

## Usage

```
mandoline [args] plotfile
          -n --normal_dir {Normal coordinate to the slice x:0, y:1, z:2}
          -c --coordinate {Coordinate of the slice, defaults to mid plane}
          -v --variable {variable name, defaults to "density"}
          -L --max_level {Maximum AMR level used, defaults to finest level}
          -o --output {Either "image" or "array"
                       plot: creates a image saved as {plotfile}_{dir}_{variable}.png
                       array: creates a numpy array saved as {plotfile}_{dir}_{variable}.npy}
          -f --file {File name used to override the default value}
          -- Image only options --
          -c --colormap {A named matplotlib colormap, defaults to viridis}
          -m --minimum {Minimum value used in the colormap}
          -M --maximum {Maximum value used in the colormap}
          -l --log {Flag to use log scale in the plot}
```

## Caveats

`mandoline` infers the mesh geometry from the grid box data in the plotfile header file. 
This means that the number of point in a given slice is dependent on the `amr.blocking_factor` parameter.
On a given level the number of points in the slice will equal or lower than the number of cells divided by `amr.blocking_factor`.
This means that for a problem with 256 cells at level 0 and the default `amr.blocking_factor` of 8 a slice at level 0 would have 32x32 points. 
With `amr.ref_ratio = 2`, slice resolutions at levels 1, 2, 3 and 4 would be 64x64, 128x128, 256x256 and 512x512, a not so bad resolution for a 4" image. 

