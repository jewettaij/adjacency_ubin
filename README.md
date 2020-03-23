[![GitHub repo size](https://img.shields.io/github/repo-size/jewettaij/adjacency_ubin)]()
[![License](https://img.shields.io/badge/License-GPL2-green.svg)]()

adjacency_ubin
===========
Generate a (sparse) adjacency matrix for a point cloud.  This program
generates a list of all pairs of points that lie within a user-specified
cutoff distance. This program uses a simple uniform binning algorithm.

## Typical Usage:

```
adjacency_ubuin -r rmax [optional arguments..] < coordinate_file > matrix_file
```

This simple program reads a file containing X Y Z coordinates of a list of
particles and prints out a list indicating which pairs of particles who are
nearby (ie. within a distance of "rmax").

Uniform binning (the simplest possible algorithm) is used to reduce
the running time from *O(n^2) to *O(n+nx\*ny\*nz)*
(where *n* is the number of points,
and *nx*, *ny*, *nz* are the number of bins in the x,y,z directions).
(This can optionally be controlled using the "-bins" argument.)
This crude program is not well suited for highly non-uniform distributions
of points, or point clouds in high dimensions.  For these applications,
use a more sophisticated spatial subdivision algorithm.

## Input format:
The file containing the coordinates of the points in the point cloud
should be a simple 3-column text file.
```
x1 y1 z1
x2 y2 z2
x3 y3 z3
x4 y4 z4
:  :  :
xn yn zn
```

## Output format:
All pairs of points which are in contact are reported using
the following format.
```
i1 j1 1
i2 j2 1
i3 j3 1
:  :  :
```
Here "i" and "j" refer to the ID numbers of two different particles
that are nearby each other (separated by a distance of less than "rmax").
*(The 3rd column ("1") is only useful when computing contact frequency
  from multiple different point clouds.  See below.)

### Multiple point clouds (trajectories"):

This list of coordinates can be optionally divided into multiple groups or 
"snapshots", which can be thought of as frames from a trajectory.
This file can read in RAW (simple 3-colum ASCII, shown above) or XYZ format.
This program generates a sparse matrix indicating which pairs of particles
are within a certain distance (rmax) relative to each other.
When multiple snapshots are present, it will keep track of how
frequently those particles were nearby eachother.


### Input format (.RAW file format):

The coordinate file can be divided into multiple groups or snapshots
These are separated/delimited by blank lines (or lines with only comments).

```
# 1st snapshot:
x11 y11 z11
x12 y12 z12
x13 y13 z13
 :  :  :
x1n y1n z1n

# 2nd snapshot:
x21 y21 z21
x22 y22 z22
x23 y23 z23
 :  :  :
x2n y2n z2n

# 3rd snapshot:
x31 y31 z31
x32 y32 z32
x33 y33 z33
 :  :  :
x3n y3n z3n

 :  :   :
```
Each non-blank line of this file contains 3 numbers which are the coordinates
of a particle. That particle is assigned an ID number equal to the number of 
lines read since the beginning of the current snapshot (starting at 1). 

### Output format for trajectories:
```
i1 j1 count1
i2 j2 count2
i3 j3 count3
 :  :   :
```
Here "i" and "j" refer to the ID numbers of two different particles
and "count" refers to the number of times (over all the snapshots) that 
those two particles were found to be nearby each other
(ie, found to be separated by a distance of less than "rmax").
("count" is 1 if the coordinate file contains only a single snapshot.
 Pairs "i" and "j" which are never in contact (count=0) are not reported.
 Pairs where i > j are redundant and are not reported.)

## Optional arguments:
```
  -r rmax      <- specify the threshold contact distance (1.5 by default).

  -raw         <-The default input file format: 3-column space-
                 delimited text file with blank lines separating frames.
                 (RAW format is assumed by default.)
  -xyz         <-Instead, assume the input coordinate file uses ".xyz" format.
                 (This feature is not well tested.  Hopefully it works.)
  -bins N      <-Specify the number of bins
                 Alternatley, you can specify 3 numbers (nx ny nz)
                 if you want it to vary in different directions.
                 Binning is only used to speed up the computation, so
                 these numbers should not effect the result, as long as 
                       (xmax-xmin)/nx >= rmax
                       (ymax-ymin)/ny >= rmax
                       (zmax-zmin)/nz >= rmax
                 where xmax, xmin, ymax, ymin, zmax, zmin are the maximum and
                 minimum coordinates in the input file.
                 However using a small number of bins may reduce memory usage.
                 (NOTE: The "-bins" option has not been tested and may not work)
```


## Changing the number of dimensions:

*Note:* This program can be compiled to work with coordinate files for particles
in dimensions other than 3 by changing the compiler flags at compile time
from *"-DG_DIM=3"* to *"-DG_DIM=4"*. (The number of dimensions is 3 by default.)
*(One way to do this editing the *"setup_gcc.sh"* file.  See below.)*

This program works well for 2D and 3D point clouds.
However this simple program is *particularly ill-suited* for high
dimensional point clouds.  Uniform binning requires exponentially
more space and time as the number of dimensions increases.


## Compilation

### Linux and Apple macOS:

```
cd src
source setup_gcc.sh
make
```

*(Note:  If you are not using the bash shell,
enter "bash" into the terminal beforehand.)*

*(Note: Apple users will need to install the gcc compiler
and other build tools using Xcode or brew.)*

*(Note: If you receive an error regarding "omp" or "OpenMP", then use
"setup_gcc_serial.sh" instead.  This may be necessary for apple users.)*

### Windows 10:

Install HyperV (with linux), or the Windows Subsystem for Linux (WSL) and run

```
sudo apt-get install build-essential
```
and then follow the instructions above.
(Older windows users can install Cygwin or MinGW, or linux via virtualbox.)


## Development Status: *alpha*

This program was quickly written has only tested on a few point clouds.
Some of the arguments documented above may not yet work.


## License

*adjacency_ubin* is available under the terms of the [GPL-2.0 license](LICENSE.txt).
