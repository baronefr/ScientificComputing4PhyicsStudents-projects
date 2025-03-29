# task 03 handout

handout from Francesco Pio Barone

## Submission comments

### 1: generate arrays

Source file: `generator.c`

### 2: daxpy

Source file: `daxpy_from_config.c`

The settings are specified in `daxpy.conf`.

### 3: makefile

The makefile is able to compile the source file for both (1) and (2). Run `make run` to build and execute the code. Eventually, run `make clean` to get rid of the executables, then `make cleandata` to remove the data files too.

### 4: HDF5

The directory `hdf5/` contains the previous source files, but edited to write and read file in HDF5 format.
Please check that HDF5 lib is installed in your system and edit the tags `-L/usr/lib64 -lhdf5_hl -lhdf5` in the makefile according to your installation.

### 5: GSL

I wrote the code to perform $y \leftarrow ax+y$ with [GSL](https://www.gnu.org/software/gsl/) in the directory `gsl/`. This is only a quick demo to show how to allocate the vectors and perform the operation with native GSL routines.

## External sources

`parser.h` is taken from [github.com/welljsjs/Config-Parser-C](https://github.com/welljsjs/Config-Parser-C), as suggested from the task sheet.
