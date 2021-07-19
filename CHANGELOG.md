#v1.6.1
* fixed total number of cells for velocity files


#v1.6
* truncated the vertices to 5 spaces after the comma for the phi direction. Last vertice is now exactly 2 pi of RADMC3D.
* hardcoded the half pi in the theta axis with the exact value according to RADMC3D
* New feature: Writing data files in ascii.

#v1.5
* changed the initial parameters of the gaussian fit to estimators: max , mean and standart deviation.
* added dust evaporation feature

#v1.4.1
* changed the input of the path information to just the basepath of the parent directory.

#v1.4
* actually implemented a working `force` feature with the help of a new `fetch_features` function.

#v1.3
* added `force` feature. If `True`, it generates the dust files from the gas file in any case. If `False` it only generates it if it the dust data file does not exist.

#v1.2.1
*  actually extend all fields with repeating outermost layers, except for density, which is still extended with a gaussian.

#v1.2
* added feature, which allows to fetch all the features/fields from the data directory.

#v1.0.2
* fixed a bug where if you extend but not mirror it would crash.

#v1.0.1
* changed the way you set the extension. Now only need to input the number of cells to extend by. `0` corresponds to no extension. No longer a need for boolean input.

#v1.0
