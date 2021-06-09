# JUPYTER .dat to RADMC3D input files conversion class

'v1.0'

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

```
import DATtoRADMC as *

foo = DATtoRADMC()
```


## Usage

To convert the JUPYTER files, you first have to give several user inputs.

### Provide information about the simulation

* Setting the simulation number, where i_sim is an integer.

```
foo.SetOutNumber(i_sim)
```

* Setting the number of layers including the base layer, where n_levels is an integer.

```
foo.SetLevel(n_levels)
```

* Setting the number of refinement levels, where n_ref is an integer.

```
foo.SetRefinement(n_ref)
```

* Whether to mirror the vertical axis, default is True.
```
foo.SetMirror(True)
```

* Whether to expand the vertical axis, default is True.
```
foo.SetExtend(True)
```

* Setting the number of cells to extend by in the vertical direction, where n_ext is an integer, default is 30.

```
foo.SetN_ext(n_ext)
```

* Setting the radius of the orbit in Au, where rad is a float.
```
foo.SetRadius(rad)
```

* Setting the mass in jupiter masses, where m is a float.
```
foo.SetMass(m)
```

* Setting the Filepath Information of the Jupyter Output files folder, where filepath is a string.
```
foo.SetInDir(filepath)
```

* Setting the list of hydrodynamical fields to convert, where fields_list is a list of strings.

```
foo.SetFeatures(fields_list)
```

* Run the conversion. Once all inputs were given, the conversion can be started and the files will be converted.

```
foo.Wrapper()
```

## Process

The Wrapper() function includes the whole conversion process. 
It starts by setting up all the necessary directories and calculates the numerical constants.
Then it reads in the descriptor file Descriptor.dat of the simulation and builds the grid. In this step in removes the innermost radial cell if there are an odd number of radial cells. If specified, it also extend the vertical axis. Once the grid layers are built, the grid parent information is calculated, including the location in the parent layer where the current layer and the size of the layer as measured in units of the parent layer starts.
After this the data file conversion starts. Here the Wrapper() function iterated over the list of features and converts one file after the other. The data file, containing all data points in one line of values for each layer, is read in. Each line is reshaped into a 3D array and, if specified, extended(See below for details). Then the iteration order is adjusted and the reordered lines are cached for later.
Next, the data file for the current feature is written in the appropriate structure. The dust data files are generated from the gas data files. The density is divided by 100 and all other are just copied.
Once all files have been successfully converted, you can find the files in a folder 'RADMC3DXXXXX', where 'XXXXX' is the zero-padded simulation number.

### Extension

The extension of the data files depends on the hydrodynamical field. 

* gas/dust density
The gas and dust density data is extended by fitting a gaussian on the initial data and then evaluated on the extended grid.

* all other fields
For all the other fields the extension works by copying the lowest and highest the specified number of times and appending them.

 




## License
[MIT](https://choosealicense.com/licenses/mit/)