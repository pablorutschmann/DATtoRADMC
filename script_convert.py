# This is a script to automate the conversion of dat to RADMC3D files
#
# default usage:
# python script_convert.py [sim_number] [grid_level] [ref_level] [radius] [mass]
#
# with optional arguments:
# python script_convert.py [sim_number] [grid_level] [ref_level] -s [mir] -e [n_ext] [radius] [mass] -d [directory] -f [force] -l [field_list]
#
# The list of fields MUST always be the final argument.
#
# If converting multiple files, it is recommended to write a bash script to
# automate calling this module multiple times.

import DATtoRADMC
import argparse
import sys
import string

parser = argparse.ArgumentParser(description='dat to RADMC3D file conversion')

# Output numbers
parser.add_argument('o',
                    action = 'store',
                    nargs = 1,
                    type = int,
                    metavar = 'sim_num',
                    help= 'Simulation number')

#Number of Grid Levels including the base layer
parser.add_argument('n',
                    action = 'store',
                    nargs = 1,type = int,
                    metavar = 'grid_level',
                    help= 'Number of grid levels including base layer')

#Number of mesh refinement levels
parser.add_argument('g',
                    action = 'store',
                    nargs = 1,type = int,
                    metavar = 'ref_level',
                    help= 'Number of grid refinement levels')

#Whether to mirror the vertical axis
parser.add_argument('-s',
                    action = 'store',
                    nargs = 1,
                    type = bool,
                    metavar = 'mir_bool',
                    help= 'Boolean for mirroring, default is True')

#Setting the number of cells to extend
parser.add_argument('-e',
                    action = 'store',
                    nargs = 1,
                    type = int,
                    metavar = 'n_ext',
                    help= 'Number of cells to extend by on both sides, default is 30')

# Science
# Radius
parser.add_argument('r',
                    action = 'store',
                    nargs = 1,
                    type = float,
                    metavar = 'radius',
                    help= 'Planetary orbital radius in AU')
#Stellar mass
parser.add_argument('m',
                    action = 'store',
                    nargs = 1,
                    type = float,
                    metavar = 'mass',
                    help= 'Stellar mass in solar mass units.')

parser.add_argument('-d',
                    "--directory",
                    action = 'store',
                    type = str,
                    nargs = 1,
		            required = False,
                    help= "Directory to output folders, should contain a folder labelled outputXXXXX, where XXXXX is the simulation output number padded to 5 digits.")

parser.add_argument('-f',
                    '--force',
                    action = 'store',
                    type = bool,
                    nargs = 1,
		            required = False,
                    help= 'If True it ignores the found dust files and generates them form the gas files. If False it generates only the dust files that have no data file.')


parser.add_argument('-l',
                    '--listoffields',
                    action = 'append',
                    type = str,
                    nargs = argparse.REMAINDER,
                    help= 'List of hydrodynamic fields to convert, default is all. which fetches the features from the data directory',
	   	            required = False)


args = parser.parse_args()

print("Converting outputs from " + str(args.o[0]).zfill(5))



dv = DATtoRADMC.DATtoRADMC()
#Setting the simulation number
dv.SetSimNumber(args.o[0])

#Setting the number of layers including the base layer
dv.SetLevel(args.n[0])


#Setting the number of refinement levels
dv.SetRefinement(args.g[0])

#Whether to mirror the vertical axis
if args.s is not None:
    dv.SetMirror(args.s[0])

#Setting the number of cells to extend by
if args.e is not None:
    dv.SetExtend(args.e[0])

#Setting the radius
dv.SetRadius(args.r[0])

#Setting the mass
dv.SetMass(args.m[0])

#Setting the parent directory of the folder with the JUPITER output folder in it.
if args.force is not None:
    dv.SetForce(args.force[0])

if args.directory is not None:
    dv.SetBasePath(args.d[0])

#Setting the features
if args.listoffields is not None:
    dv.SetFeatures(args.listoffields[0])
    if 'all' in args.listoffields[0]:
        print('Converting all files in the data directory')
    else:
        print("Converting fields " + str(args.listoffields[0]))

else:
    print('Converting all files in the data directory')

dv.Wrapper()

