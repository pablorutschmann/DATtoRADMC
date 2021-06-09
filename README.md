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

#Setting the mass in jupiter masses, where m is a float.
```
foo.SetMass(m)
```


## License
[MIT](https://choosealicense.com/licenses/mit/)