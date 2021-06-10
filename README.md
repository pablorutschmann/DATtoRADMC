# JUPYTER .dat to RADMC3D input files conversion class

`v1.0.2`

Rutschmann Pablo, June 2021

rupablo@student.ethz.ch

This is a class that can convert the hydrodynamical simulation output files of the JUPYTER (file format .dat) to input files for RADMC3D. It can also extend the coordinates and data files by a user set number of cells. 

Based on a project by Evert Nasedkin @nenasedk , 2018: [JUPITER_VTKFileConversion](https://github.com/nenasedk/JUPITER_VTKFileConversion)

## Requirements
* Python 3.8
* numpy
* scipy
* astropy


## Installation

Download the DATtoRADMC.py python file and import the module into another python file, where you initialise the class and set up the conversion.

```python
import DATtoRADMC as *

foo = DATtoRADMC()
```


## Usage of the class `DATtoRADMC()`

To convert the JUPYTER files, you first have to give several user inputs.

### Provide information about the simulation

* Setting the simulation number, where `sim_num` is an integer.

```python
foo.SetOutNumber(sim_num)
```

* Setting the number of layers including the base layer, where `grid_levels` is an integer.

```python
foo.SetLevel(grid_levels)
```

* Setting the number of refinement levels, where `ref_levels` is an integer.

```python
foo.SetRefinement(ref_levels)
```

* Whether to mirror the vertical axis, where `mir` is a boolean, default is `True`.
```python
foo.SetMirror(mir)
```

* Setting the number of cells to extend by in the vertical direction in each direction(so the length of the vertical axis will increase by `2 * n_ext`), where `n_ext` is an integer, default is `30`. `0` for no extension. 

```python
foo.SetN_ext(n_ext)
```

* Setting the orbital radius of the planet in Au, where `radius` is a float.
```python
foo.SetRadius(radius)
```

* Setting the mass in solar masses, where `mass` is a float.
```python
foo.SetMass(mass)
```

* Setting the Filepath Information of the Jupyter Output files folder, where `filepath` is a string.
```python
foo.SetInDir(filepath)
```

* Setting the list of hydrodynamical fields to convert, where `field_list` is a list of strings, default is `'all'`, which fetches all the features from files it can find in the data directory.

```python
foo.SetFeatures(field_list)
```

* Run the conversion. Once all inputs were given, the conversion can be started by calling the `Wrapper()` function and the files will be converted.

```python
foo.Wrapper()
```

## Usage with `script_convert.py` in the command line

With the python script `script_convert.py` you can easily run the conversion directly from the command line. The order of the parameters is the same as above. The following lists the required arguments. Note that the field list argument MUST be the final argument in the command.
```
python3 script_convert.py [sim_num] [grid_levels] [ref_level] [radius] [mass]
```

or with full options:
```
python3 script_convert.py [sim_number] [grid_level] [ref_level] -s [mir] -e [n_ext] [radius] [mass] -d [directory] -f [field_list]
```

If no list of fields is given it defaults to `'all'`, which fetches all the files in the data directory.

The directory should contain a folder labelled `outputXXXXX`, where `XXXXX` is the simulation output number padded to 5 digits.

An example for output 235, with 5 grid levels, 4 refinement levels, radius of 50 Au and 1 solar mass, converting the gas- and dustdensity:
```
python3 script_convert.py 235 5 4 50 1 -f gasdensity dustdensity
```


## Process

The `Wrapper()` function includes the whole conversion process. 
It starts by setting up all the necessary directories and calculates the numerical constants.

Then it reads in the descriptor file `Descriptor.dat` of the simulation and builds the grid. In this step in removes the innermost radial cell if there are an odd number of radial cells. If specified, it also extend the vertical axis. Once the grid layers are built, the grid parent information is calculated, including the location in the parent layer where the current layer and the size of the layer as measured in units of the parent layer starts.

After this the data file conversion starts. Here the `Wrapper()` function iterates over the list of features and converts one file after the other. The data file, containing all data points in one line of values for each layer, is read in. Each line is reshaped into a 3D array and, if specified, extended(See below for details). Then the iteration order is adjusted and the reordered lines are cached for later.
Next, the data file for the current feature is written in the appropriate structure. The dust data files are generated from the gas data files. The density is divided by 100 and all other are just copied.

Once all files have been successfully converted, you can find the files in a folder `RADMC3DXXXXX`, where `XXXXX` is the zero-padded simulation number.

### Extension

The extension of the data files depends on the hydrodynamical field. 

* gas/dust density
The gas and dust density data is extended by fitting a gaussian on the initial data and then evaluated on the extended grid.

* all other fields
For all the other fields the extension works by copying the lowest and highest the specified number of times and appending them.

 

## Roadmap

* make a user friendly python script

* create pip install functionality


## License