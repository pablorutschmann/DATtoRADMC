import os, sys
import numpy as np
from scipy.optimize import curve_fit
import astropy.units as u

class DATtoRADMC:
    """
    Convert JUPITER .DAT files to binary Files for RADMC3D

    Reads in a Descriptor.dat file and a hydrodynamic field file (.dat)
    Outputs a
    """

    def __init__(self):
        # Simulation Information
        self.outNumber = -1  # Simulation Number
        self.nLevels = -1  # Number of Mesh Levels (with Base Grid)
        self.nRefinements = -1 # Refinement levels
        self.mirror = False # Mirror theta axis
        self.n_extend = 30 # Number of cells to extend by. Default: 30. 0 for no extension
        self.features = []  # All to be converted Hydrofields
        self.__feature = 'notafeat' # Currently converted __feature, should not be accessed by user.
        self.__cur00 = 0.0
        self.nLevelCoords = []  # Number of Vertices along each axis at each refinement level [[len(phi[0]),len(r[0]),len(th[0])],[...],...]
        self.oddR = False # Oddness of Base Grid Radial Cell Number
        self.precis = 8 # Floating point precision of the output files
        self.nrspec = 1 # Number of fluids

        self.featlist = ['gasdensity', 'gastemperature', 'gasenergy',
                         'gaserad', 'gasopacity', 'gaspotential',
                         'gasstheat', 'gastau', 'gastaucell',
                         'gasvelocity', 'dustdensity', 'dusttemperature',
                         'dustenergy', 'dusterad', 'dustopacity',
                         'dustpotential', 'duststheat', 'dusttau',
                         'dusttaucell', 'dustvelocity']
        self.out_featlist = ['gas_density', 'gas_temperature', 'gas_energy',
                         'gas_erad', 'gas_opacity', 'gas_potential',
                         'gas_stheat', 'gas_tau', 'gas_taucell',
                         'gas_velocity', 'dust_density', 'dust_temperature',
                         'dust_energy', 'dust_erad', 'dust_opacity',
                         'dust_potential', 'dust_stheat', 'dust_tau',
                         'dust_taucell', 'dust_velocity']

        self.completed = []

        # Filepath information
        self.dataDir = 'notadir'  # Jupyter output folder
        self.dataOutPath = 'notapath'  # Destination Folder
        self.inFilename = 'notaname'  # Filename constructed from __feature, outNumber and nLevels (Changes for every level)
        self.outFilename = 'notaname'
        self.descriptorName = 'notadesc' # Name of the descriptor file
        self.BASEPATH = os.getcwd() + '/' # Current Directory

        # Science information
        self.radius = -1.
        self.rcgs = -1. * u.g
        self.mass = -1.
        self.mcgs = -1. * u.g
        # Science Constants
        self.TEMP = -1.
        self.DENS = -1.
        self.PERIOD = -1.
        self.VEL = -1.
        self.units = "CGS"

        # Grid Information

        # Stores coordinates of grid points at each refinement level. Example: LevelCoords['r'][0] is the base grid r array
        self.LevelCoords = {'r': [],
                            'th': [],
                            'phi': []}

        self.ncells = []  # Number of Grid Cells in each layer: [[n_r_0, n_th_0, n_phi_0],[n_r_1,...]...]
        self.ncells_filt = {}  # Number of points in unfiltered in each mesh refinement level
        self.cellist = []  # Which cells were filtered out
        # Caches the Parent Information in a Dictionary
        self.pi_r = {}
        self.pi_th = {}
        self.pi_phi = {}
        # Caches the converted data of each layer in a list for the current feature. Is reset before converting next __feature.
        self.converted = []

        self.mins = []  # Min boundaries of a each mesh level
        self.maxs = []  # Max boundaries of a each mesh level

# ------------------------------------------------------------------------------------------------

    # Basic user Set functions
    def SetOutNumber(self, n):
        self.outNumber = n

    def SetLevel(self,level):
        self.nLevels = level

    def SetRefinement(self,ref):
        self.nRefinements = ref

    def SetMirror(self,bool):
        self.mirror = bool

    def SetExtend(self,n_ext):
        self.n_extend = n_ext

    def SetFeatures(self, feats):
        for feat in feats:
            if feat in self.featlist:
                self.features.append(feat)
            else:
                print(feat + " is not a recognised input, may not be implemented. Continuing...")
                self.features.remove(feat)
        self.__feature = self.features[0]

    def SetBasePath(self,path):
        self.BASEPATH = path
    def SetInDir(self,path):
        self.dataDir = path
    def SetOutDir(self,path):
        self.dataOutPath = path
    def SetInFilename(self,fname):
        self.inFilename = fname

    def SetUnits(self, units):
        ulist = ["CGS","cgs","AU","Au","au"]
        if units in ulist:
            self.units = units
        else:
            print("Please select a valid unit system.")
        return
    # Science Set Functions
    def SetRadius(self, rad):
        """
        Orbital radius of companion
        """
        self.radius = rad*u.AU
        self.rcgs = self.radius.to(u.cm)
        if(self.mass>0. and self.TEMP <0):
            self.SetConstants()

    def SetMass(self,mass):
        """
        Mass of star in solar masses
        """
        self.mass = mass*u.M_sun
        self.mcgs = self.mass.to(u.g)
        if(self.radius > 0. and self.TEMP < 0):
            self.SetConstants()

    def SetConstants(self):
        # 6.67e-8 - cgs G
        # 8.314e7 - cgs gas constant
        if self.radius <= 0:
            raise Exception('Radius not set! Please set a radius.')
        if self.mass <= 0:
            raise Exception('Mass not set! Please set a Mass.')
        if self.radius > 0 and self.mass > 0:
            self.TEMP = ((self.rcgs.value)/(np.sqrt((self.rcgs.value)**3 / 6.67259e-8 / (self.mcgs.value))))**2 / 8.314e7
            self.DENS = self.mcgs.value/(self.rcgs.value)**3
            self.PERIOD = 2*np.pi*np.sqrt((self.rcgs.value)**3 / (6.67259e-8 * self.mcgs.value))
            self.VEL = self.rcgs.value/(self.PERIOD/2*np.pi)


    # ----------------------------------------------------------------------------------------

    # Directory and file Setup
    def SetupDirs(self):
        # Error checking
        if self.outNumber < 0:
            print("Please set the simulation number")
            return

        # Create directory paths
        self.dataDir = self.BASEPATH + "output" + str(self.outNumber).zfill(5) + "/"
        self.dataOutPath = self.BASEPATH + "RADMC" + str(self.outNumber).zfill(5) + "/"
        self.descriptorName = "Descriptor" + str(self.outNumber) + ".dat"
        self.outFilename = self.__feature + str(self.outNumber).zfill(3) + "_" + str(self.nLevels) + ".dat"

        # Check existance
        print("Data directory is: " + self.dataDir)
        print("If this is incorrect, please set a base path to the data directory.")
        if not os.path.isdir(self.dataDir):
            print("ERROR: data directory " + self.dataDir + " does not exist!")
            return

        # Make out dir
        if not os.path.isdir(self.dataOutPath):
            os.mkdir(self.dataOutPath)

    def SetupNames(self, level):
        if self.outNumber < 0:
            print("Please set the simulation number")
            return
        if self.__feature == 'notafeat':
            print("Please input a __feature (e.g. gasvelocity)")
            return
        if self.nLevels < 0:
            print("Please input the number of mesh refinement levels")
            return

        self.inFilename = self.__feature + str(self.outNumber) + "_" + str(level) + "_" + str(level) + ".dat"

        try:
            assert (os.path.isfile(self.dataDir + self.inFilename))
        except AssertionError:
            print(self.dataDir + self.inFilename + " Does not exist. Please enter a valid filename or filepath.")
            sys.exit(1)

    # ---------------------------------------------------------------------------------------------
    # Important part starts here
    # ---------------------------------------------------------------------------------------------

    # ---------------------------------------------------------------------------------------------
    # Wrapper
    #
    # Wraps the whole conversion process under one function. Once all user inputs were given this function can be called and it will:
    # 1. Setup the Directories
    # 2. Build the Base Grid and calculates/stores the Parent Information
    # 3. Goes through all the features and performs the conversion process:
    #     a. Prepares the converted data by reordering it and possibly mirrors/extends it
    #     b. Writes the converted data into the RADMC3D input write_data_file
    # 4. Checks the dimensions of the data and grid file
    # 5. Finally writes the Grid File
    # ---------------------------------------------------------------------------------------------

    def Wrapper(self):
        self.SetupDirs()
        self.SetConstants()
        self.GetCoordinates()

        self.completed = []

        for feat in self.features:
            if 'dust' in feat:
                pass
            else:
                print(feat)
                self.__feature = feat
                self.ncells_filt[self.__feature] = []
                self.ConvertFiles()
                self.write_data_file()
                self.completed.append(feat)

        self.write_grid_file()
        print('Conversion Completed, converted ' + str(len(self.completed)) + '/' + str(len(self.features)))

    # ---------------------------------------------------------------------------------------------
    # Convert Files
    #
    # Converts all Levels of Jupyter Output Data to 1D lists in the right order.
    #   1. reads in the data files for the current __feature and reshapes it into 3D-Array
    #   2. mirrors/extend the 3D-Array
    #   3. reorders each level
    #   3. appends the converted 1D lists to the Cache List 'converted'.
    # ---------------------------------------------------------------------------------------------
    def ConvertFiles(self):
        converted_dat = []
        self.converted = []
        for i in range(self.nLevels):
            self.SetupNames(i)
            in_dat = np.fromfile(self.dataDir + self.inFilename, dtype='double')

            if "velocity" in self.__feature:
                data2 = np.array([], dtype=np.float64)
                data2 = np.append(data2, in_dat.astype(np.float64))
                data2 = np.reshape(data2, (3, -1))
                # 3 lines with entire data: theta, radial, phi
                data3 = np.column_stack((data2[0], data2[1], data2[2]))

                reordered_th = self.reorder_one_line(data3[:,0],i)
                reordered_r = self.reorder_one_line(data3[:,1], i)
                reordered_phi = self.reorder_one_line(data3[:,2],i)

                reordered_dat = []

                for i in range(len(reordered_th)):
                    reordered_dat.append(reordered_r[i])
                    reordered_dat.append(reordered_th[i])
                    reordered_dat.append(reordered_phi[i])

                print(len(reordered_dat))
                reordered_dat = [x * self.VEL for x in reordered_dat]

            else:
                reordered_dat = self.reorder_one_line(in_dat,i)
                if ("density" in self.__feature):
                    reordered_dat = [x * self.DENS for x in reordered_dat]
                if ("temperature" in self.__feature):
                    reordered_dat = [x * self.TEMP for x in reordered_dat]

            self.converted.append(reordered_dat)


    # ---------------------------------------------------------------------------------------------
    # GetCoordinates

    # This function generates all the Information for the amr_grid.inp file:
    #   1. Reads in the Descriptor file
    #   2. Stores the grid points in the LevelCoords Dictionary
    #   3. Stores the Number of Coordinates of each axis at each level in nLevelCoords
    #   4. Calls get_parent_info to store the parent information for each axis in their dictionaruies
    #
    # ---------------------------------------------------------------------------------------------
    def GetCoordinates(self):
        phi = []
        r = []
        th = []
        self.ncells = []

        if self.nLevels < 0:
            raise Exception('Please set the number of mesh levels!')

        if self.nRefinements < 0:
            raise Exception('Please set the number of refinement levels!')

        dsc = open(self.dataDir + self.descriptorName)

        for i in range(self.nLevels):
            self.SetupNames(i)
            dsc = open(self.dataDir + self.descriptorName)
            for j, line in enumerate(dsc):
                if j == 6 + (i*11):
                    n_phi, n_r, n_th = [int(x) for x in line.split()]
                    if i == 0 and n_r % 2 == 1:
                        self.oddR = True
                    self.ncells.append([n_r, n_th, n_phi])
                # Invert theta, phi so that coordinate system is right handed
                # (needed for cell orientation)
                if j == 8 + (i*11):
                    cur = [np.float64(x) for x in line.split()]  # Import one line into list of values
                    cur.pop(0)  # First and last two points are 'ghost points'
                    cur.pop(0)
                    cur.pop()
                    cur.pop()
                    if i == 0:
                        self.__cur00 = cur[0]
                        #print(round(cur[0] + np.pi,8))
                    # shift by pi
                    # print('Cur: ')
                    # print(cur)
                    # print('Cur + pi: ')
                    # print([x - self.__cur00 for x in cur])
                    phi.append([x - self.__cur00 for x in cur])
                    cur = []
                if j == 9 + (i*11):
                    cur = [np.float64(x) for x in line.split()]
                    cur.pop(0)  # First and last two points are 'ghost points'
                    cur.pop(0)
                    #pop first radial vertice for odd number of cells
                    if self.oddR == True and i == 0:
                        cur.pop(0)
                    cur.pop()
                    cur.pop()
                    #ToDo Scale Radius to right Units
                    cur_scaled = [np.float64(x * (self.rcgs.value))  for x in cur]
                    r.append(cur_scaled)
                    cur = []
                if j == 10 + (i*11):
                    cur = [np.float64(x) for x in line.split()]
                    cur.pop(0)  # First and last two points are 'ghost points'
                    cur.pop(0)
                    cur.pop()
                    cur.pop()
                    th.append((self.extend_th_coords(cur,i)))
                    cur = []
            self.nLevelCoords.append([len(phi[i]),len(r[i]),len(th[i])])


            dsc.seek(0)

        #Put phi,r,th array into LevelCoords
        self.LevelCoords['phi'] = phi
        self.LevelCoords['r'] = r
        self.LevelCoords['th'] = th
        self.get_parent_info()
        # ToDo Check if nLevelCoords of base grid matches self.ncells
        # ToDo Print Length of Base Grid Arrays if they match
        if self.ncells[0][0] == self.nLevelCoords[0][0]:
            raise Exception('Number of cells and number of vertices in the Base-Grid do not match!')
        else:
            print(str(self.nLevelCoords[0][0]) + ' x ' + str(self.nLevelCoords[0][1]) + ' x ' + str(self.nLevelCoords[0][2]) + " Vertices in filtered Base-Grid.")

    # ---------------------------------------------------------------------------------------------
    # get_parent_information
    # used in get_coordinates
    # Calculates the amr_grid.inp Grid File Parent Information and stores it in the Dictionaries.
    # Actually just a wrapper for the pi_one_axis function.
    # ---------------------------------------------------------------------------------------------
    def get_parent_info(self):
        # -----------------------------------------------------------------------------------------
        # pi_one_axis
        # Input: list of coords for each layer for one axis
        # For a list of Coords for each layer it:
        #   1. Find the Minima and Maxima
        #   2. Searches for the index closest to the minima of the parent coords from the front
        #   3. Searches for the index closest to the maxima of the parent coords from the back
        #   4. Stores all the information in a dictionary
        # Output: Dictionary with the mins, maxs, first and last shared and relative sizes.
        # -----------------------------------------------------------------------------------------
        def pi_one_axis(n_coords):
            dict = {}
            dict['mins'] = []
            dict['maxs'] = []
            dict['first shared'] = []
            dict['last shared'] = []
            dict['sizes'] = []
            #Starts at level 1
            for i in range(1, len(n_coords)):
                dict['mins'].append(min(n_coords[i]))
                dict['maxs'].append(max(n_coords[i]))
                #FIXED: iterate from front for th centered around 0
                smallest_diff = np.abs(n_coords[i - 1][0] - n_coords[i][0])
                closest_index = 0
                for j,item in enumerate(n_coords[i-1][1:], 1):
                    if abs(item - n_coords[i][0]) < smallest_diff:
                        closest_index = j
                        smallest_diff = np.abs(n_coords[i - 1][j] - n_coords[i][0])
                # +1 because radmc3d counts index from 1
                first_shared = closest_index + 1
                #first_shared = np.argmin(np.abs(n_coords[i - 1] - n_coords[i][0]))+1
                #print(n_coords)
                #FIXED: Iterate from Back for th centered around 0
                smallest_diff = np.abs(n_coords[i - 1][-1] - n_coords[i][-1])
                closest_index = len(n_coords[i - 1]) - 1
                for j, item in reversed(list(enumerate(n_coords[i - 1][0:-1]))):
                    if abs(item - n_coords[i][-1]) < smallest_diff:
                        closest_index = j
                        smallest_diff = np.abs(n_coords[i - 1][j] - n_coords[i][-1])
                # +1 because radmc3d counts index from 1
                last_shared = closest_index + 1
                #last_shared = np.argmin(np.abs(n_coords[i - 1] - n_coords[i][-1]))+1
                dict['first shared'].append(first_shared)
                dict['last shared'].append(last_shared)
                dict['sizes'].append(last_shared - first_shared)
            return dict
        self.pi_r = pi_one_axis(self.LevelCoords['r'])
        self.pi_th= pi_one_axis(self.LevelCoords['th'])
        self.pi_phi = pi_one_axis(self.LevelCoords['phi'])

    # ---------------------------------------------------------------------------------------------
    # extend_th_coords
    # used in get_coordinates
    # Input: Theta Array of Grid Points for one layer
    #   1. Extends, if wanted, the theta axis by n_extend cells
    #   2. Mirrors, if wanted, the theta axis along the midplane
    # Output: extended and mirrored theta array
    # ---------------------------------------------------------------------------------------------
    def extend_th_coords(self,th_array, index):
        len_th = len(th_array)
        th_diff = th_array[1] - th_array[0]
        #if ext == True, extends the theta array by 30 before mirroring
        if self.n_extend > 0 and index == 0:
            for i in range(0, self.n_extend):
                th_array.insert(0, th_array[0] - th_diff)
            len_th = len(th_array)

        flipped_th = []
        #Creates the mirror array
        if self.mirror == True:
            for i in range(1, len_th):
                flipped_th.append(np.float64(th_array[-1] + (i * th_diff)))
            return th_array + flipped_th
        else:
            return th_array

    # ---------------------------------------------------------------------------------------------
    # reorder_one_line
    # This function is used in ConvertFiles
    # It takes the 1D-list from one layer of the data file and reorders the coordinates from (phi,r,th) r,th,phi) fastest to slowest iterated.
    #   1. Reshapes the data array into a 3D-Array with Axis [th, r, phi]
    #   2. Mirrors, if wanted, the theta direction
    #   3. Extends, if wanted, the theta direction by n_extend cells
    #   4. Reorders 3D_array into a 1D List with (r,th,phi) fastest to slowest iterated.
    # Output: Reordered and extended/mirrored 1D list for that refinement layer
    # ---------------------------------------------------------------------------------------------
    def reorder_one_line(self,array,index):

        num_r, num_th, num_phi = self.ncells[index]

        reshaped_dat = np.reshape(array, [num_th, num_r, num_phi])

        if self.oddR == True and index == 0:
            # print(np.shape(reordered_dat[0,:,:]))
            reshaped_dat = reshaped_dat[:, 1:, :]
            num_r -= 1

        if self.mirror == True:
            flipped = np.flip(reshaped_dat,axis=0)
            reshaped_dat = np.concatenate((reshaped_dat,flipped),axis=0)
            num_th = num_th * 2
            #self.ncells[index][1] = num_th

        if self.n_extend > 0 and index == 0 and 'density' in self.__feature:

            def g(x, amp, mean, stddev):
                z = (x - mean) / stddev
                y = amp * np.exp(-z ** 2 / 2)
                return y


            ynew = self.LevelCoords['th'][0]
            yold = ynew[self.n_extend:-self.n_extend]

            #Calculating the effective coordinates of the data points
            def get_data_coords(vertices):
                data_coords = []
                for i in range(np.size(vertices, axis=0)-1):
                    data_coords.append(np.mean([vertices[i],vertices[i+1]]))
                return data_coords

            data_coords_new = get_data_coords(ynew)
            data_coords_old = data_coords_new[self.n_extend:-self.n_extend]

            extended_dat = np.zeros((num_th + 2 * self.n_extend,num_r,num_phi))

            # fixing the boundaries to be continous
            sigma = np.ones(len(data_coords_old))
            sigma[0] = 0.01
            sigma[-1] = 0.01

            # do fitting and extrapolation
            for r in range(num_r):
                for phi in range(num_phi):
                    th_old = reshaped_dat[:,r,phi]
                    p_init = [th_old.max(), 1.57, 0.01]
                    coeff, cov = curve_fit(g, data_coords_old, th_old, p0=p_init, sigma=sigma, maxfev=3000)
                    th_new_bottom = g(data_coords_new[:self.n_extend], coeff[0], coeff[1], coeff[2])
                    th_new_top = g(data_coords_new[-self.n_extend:], coeff[0], coeff[1], coeff[2])

                    th_new = np.concatenate([th_new_bottom,th_old,th_new_top])
                    extended_dat[:, r, phi] = th_new


            reshaped_dat = extended_dat
            num_th += 2 * self.n_extend
                #self.ncells[index][1] = num_th

        if self.n_extend and index == 0 and 'temperature' in self.__feature:
            temp_top = reshaped_dat[-1,:,:]
            temp_bot = reshaped_dat[0,:,:]
            temp_top = np.tile(temp_top, (self.n_extend,1,1))
            temp_bot = np.tile(temp_bot, (self.n_extend,1,1))
            reshaped_dat = np.concatenate([temp_bot, reshaped_dat, temp_top],axis=0)
            num_th += 2 * self.n_extend

        #update ncells with the new numbers
        n_th,n_r,n_phi = np.shape(reshaped_dat)
        self.ncells_filt[self.__feature].append([n_r, n_th, n_phi])
        if n_th == num_th and n_r == num_r and n_phi == num_phi:
            print('Layer ' + str(index) + ' check Passed: Cell Numbers match!')
        else:
            print('Check Failed: Cell Numbers dont match!')

        #self.ncells[index] = [n_r,n_th,n_phi]

        reordered = []
        for k in range(num_phi):
            for i in range(num_th):
                for j in range(num_r):
                    reordered.append(reshaped_dat[i, j, k])

        return reordered

    # ---------------------------------------------------------------------------------------------
    # write_grid_file
    # This function is used in Wrapper
    # writes out the amr_grid.imp file
    # ---------------------------------------------------------------------------------------------

    def write_grid_file(self):
        outfile = open(self.dataOutPath + 'amr_grid.inp', 'w')
        try:
            #ToDo Add iputs for first few lines
            outfile.write('1' + '\n')
            if self.nLevels == 0:
                outfile.write('0' + '\n')
            else:
                outfile.write('10' + '\n')
            outfile.write('100' + '\n')
            outfile.write('0 \n')
            outfile.write('1 1 1' + '\n')
            outfile.write(" ".join(str(item) for item in self.ncells_filt[self.features[0]][0]))
            outfile.write('\n')
            if self.nLevels > 1:
                outfile.write(f'{self.nLevels - 1} {self.nRefinements}\n')
            outfile.write(" ".join(str(item) for item in self.LevelCoords['r'][0]))
            outfile.write('\n')
            outfile.write(" ".join(str(item) for item in self.LevelCoords['th'][0]))
            outfile.write('\n')
            outfile.write(" ".join(str(item) for item in self.LevelCoords['phi'][0]))
            outfile.write('\n')
            if self.nLevels > 1:
                for i in range(len(self.pi_r['first shared'])):
                    outfile.write(f'{str(i)} {str(self.pi_r["first shared"][i])} {str(self.pi_th["first shared"][i])} {str(self.pi_phi["first shared"][i])} {str(self.pi_r["sizes"][i])} {str(self.pi_th["sizes"][i])} {str(self.pi_phi["sizes"][i])}\n')

        finally:
            outfile.close()

    # ---------------------------------------------------------------------------------------------
    # write_data_file
    # This function is used in Wrapper
    # writes out the data file for the current __feature
    # ---------------------------------------------------------------------------------------------

    def write_data_file(self):
        n_tot_filt = 0
        for i in range(len(self.ncells_filt[self.__feature])):
            n_tot_filt += np.prod(self.ncells_filt[self.__feature][i])


        array_dbl = np.concatenate(self.converted)
        n_tot = len(array_dbl)
        if n_tot == n_tot_filt:
            print('Total Number of Cells match!')
            print('n_tot: ' + str(n_tot))
            print('n_tot_filt: ' + str(n_tot_filt))
        else:
            print('Error: Total Number of Cells dont match!')
            print('n_tot: ' + str(n_tot))
            print('n_tot_filt: ' + str(n_tot_filt))

        array_int = [1,self.precis,n_tot,self.nrspec]
        #long integerf
        array_int = np.array(array_int,dtype='int64')

        if 'gas' in self.__feature:
            feat_dust = 'dust' + self.__feature[3:]
            if feat_dust in self.features:
                if 'density' in feat_dust:
                    array_dust_dbl = array_dbl * 0.01
                else:
                    array_dust_dbl = array_dbl
                feature_dust = self.out_featlist[self.featlist.index(feat_dust)]
                if 'temperature' in feature_dust:
                    outfile_dust = open(self.dataOutPath + feature_dust + '.bdat', 'wb')
                else:
                    outfile_dust = open(self.dataOutPath + feature_dust + '.binp', 'wb')

                try:
                    array_int.tofile(outfile_dust)
                    array_dust_dbl.tofile(outfile_dust)
                    print(feature_dust + ' file was generated from the gas' + self.__feature[3:] + '.dat file')
                    self.completed.append(feat_dust)
                finally:
                    outfile_dust.close()


        array_dbl = np.array(array_dbl, dtype = 'float64')
        feature =  self.out_featlist[self.featlist.index(self.__feature)]
        if 'temperature' in self.__feature:
            outfile = open(self.dataOutPath + feature + '.bdat', 'wb')
        else:
            outfile = open(self.dataOutPath + feature +'.binp', 'wb')

        try:
            array_int.tofile(outfile)
            array_dbl.tofile(outfile)

        finally:
            outfile.close()