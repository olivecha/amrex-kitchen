# Mandoline

A command line tool to create fast slices of (large) AMReX plotfiles.

## Usage

```
mandoline [args] plotfile
          -v --variable {variable name}
          -c --colormap {A named matplotlib colormap, defaults to viridis}
          -L --max_level {Maximum AMR level used default to finest level}
          -n --normal_dir {Normal coordinate to the slice x:0, y:1, z:2}
          -c --coordinate {Coordinate of the slice, default to mid plane}
```

