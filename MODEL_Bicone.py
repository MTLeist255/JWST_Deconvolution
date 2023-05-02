# 2022 July 14: cleaned up code to generate bicone models
# 1) Create bicone model -> add to square
# 2) Create point source -> add to square + bicone
# 3) Create several point sources -> add to square + bicone
# 4) Create Gaussian background -> add to square + bicone + point source
# 5) Create galaxy emission -> add to square + bicone + point source
# 6) save FITS output

# set the input/output path
from scipy.ndimage import interpolation

path = 'Images/JWST/1_Input_Models/Models/'

from Convenience_Functions import *

# create a 1024x1024 array of zeros
naxis = (1024,1024)
square = np.zeros(naxis)

# set the simulated wavelength
wave = 10
# set name index
name_idx = 3
# create single pixel point source
point_source = 7e6
# create the scaling factor for the bicone
scale1 = 2.49

# set the names for the output files
name = ['Model_'+str(wave)+'um_single_pixel',
        'Model_'+str(wave)+'um_single_pixel_disk',
        'Model_'+str(wave)+'um_single_pixel_disk_bicone',
        'Model_'+str(wave)+'um_AGN_complicated',
        'Model_'+str(wave)+'um_galaxy_test',
        'Model_Background'
        ]

#************************************* 1) Create 2D- bicone model *****************************************************#

# create 2d array and radial component
# LARGE bicone: side = 80
# SMALL bicone: side = 40
side = 32
rad = int(side/2)
array2d = np.zeros(shape = (rad,rad))

# #ORIGINAL -> DO NOT EDIT ONLY UN-COMMENT
# #create upper left quadrant of image
for x in np.arange(rad, 0, -1):
    # divide (rad -x/n) by some number to change the opening angle of the cone
    zeros = np.zeros(int(rad-x/1))
    # create decreasing arrays
    # divide (rad -x/n) by some number to change the opening angle of the cone
    array = np.arange(rad, rad-x/1, -1)

    # fill the front end of the decreasing arrays with 0's to get a square
    # 2darray
    full = np.append(zeros,array)
    # flux ratio goals: point source / TOTAL flux in bicone
    array2d[rad-x,:] = full*scale1

# mirror upper left quadrant to create upper cone
array2d = np.concatenate((array2d,np.flip(array2d, axis = 1)), axis = 1)

# mirror entire upper cone to create bottom cone
# array2d = squareD
array2d = np.concatenate((array2d,np.flip(array2d, axis=0)), axis = 0)
coneshape = array2d

# add radial dependence utilizing code from last week
array1d = np.arange(-side/2, side/2)
array2d = np.tile(array1d, (side, 1))
for x in array1d:
    dist = np.sqrt(array1d**2+ x**2)
    array2d[:,int(x+side/2)] = dist

radial = array2d

#combine
final_product = coneshape**3 / np.sqrt(np.sqrt(radial+1))

# SMALL bicone: size = 40
final_product2 = final_product

# Set sharp edges of cone manually
# bottom cone
final_product2[1:0, 38:2] = 0
final_product2[2:1, 31:9] = 0
final_product2[3:2, 27:13] = 0
final_product2[4:3, 23:17] = 0
final_product2[5:4, 19:21] = 0

# top cone
final_product2[2:38,0:1] = 0
final_product2[9:31,1:2] = 0
final_product2[13:27,2:3] = 0
final_product2[17:23,3:4] = 0
final_product2[21:19,4:5] = 0

#********************************** 2) Create 2D-bicone by hand **************************************************#

# generate the bicone by hand?
arrayI = np.zeros(shape=(11, 11))
# generate the quadrant I of the bicone by hand (annoying)
# row 1
count = 0
value = [2.96e5, 1e5, 9e3, 4.5e3, 4e3, 1,1,1,1,1,1]
for i in range (0,int(len(value))):
    print('(', count, ',', count + 1, ') | value = ', value[i])
    arrayI[0:1, count:count+1] += value[i]
    count += 1

# row 2
arrayI[1:2, 0:1] = 1e3
count = 1
value2 = [5.87e5, 2.93e5, 1.46e5, 7.3e4, 3.6e4, 1.8e4, 10,10,10,10]
for i in range (0,int(len(value2))):
    print('(', count, ',', count + 1, ') | value = ', value2[i])
    arrayI[1:2, count:count+1] += value2[i]
    count += 1

# row 3
arrayI[2:3, 0:1] = 100
arrayI[2:3, 1:2] = 1e3
count = 2
value3 = [8.78e5, 4.5e5, 2.2e5, 9e4, 7e4, 5e4, 3e4, 100,100]
for i in range (0,int(len(value3))):
    print('(', count, ',', count + 1, ') | value = ', value3[i])
    arrayI[2:3, count:count+1] += value3[i]
    count += 1


# row 4
count = 3
trim = 0
arrayI[3:4, 0:1] += 0
arrayI[3:4, 1:2] += 0
arrayI[3:4, 2:3] += 0
for i in range (0,7):
    value = 1.17e6-trim
    arrayI[3:4, count:count+1] += value
    trim += 2.2e5
    count += 1

arrayI[3:4, 9:10] = 2e4

# row 5
count = 4
trim = 0
arrayI[4:5, 0:1] += 0
arrayI[4:5, 1:2] += 0
arrayI[4:5, 2:3] += 0
arrayI[4:5, 3:4] += 0
for i in range (0,7):
    value = 1.46e6-trim
    arrayI[4:5, count:count+1] += value
    trim += 2.345e5
    count += 1

# row 6
count = 5
trim = 0
arrayI[5:6, 0:1] += 0
arrayI[5:6, 1:2] += 0
arrayI[5:6, 2:3] += 0
arrayI[5:6, 3:4] += 0
arrayI[5:6, 4:5] += 0
for i in range (0,6):
    value = 1.75e6-trim
    arrayI[5:6, count:count+1] += value
    trim += 3.3e5
    count += 1

# SWITCH -> fill the exterior columns in a for loop and the interior columns manually
# row 7
count = 0
arrayI[6:7, 6:7] += 2.04e6
arrayI[6:7, 7:8] += 1.53e6
arrayI[6:7, 8:9] += 1e6
arrayI[6:7, 9:10] += 8e5
for i in range (0,6):
    value = 0
    arrayI[6:7, count:count+1] += value
    count +=1

# row 8
count = 0
arrayI[7:8, 7:8] += 2.33e6
arrayI[7:8, 8:9] += 1.8e6
arrayI[7:8, 9:10] += 1.3e6
for i in range (0,7):
    value = 0
    arrayI[7:8, count:count + 1] += value
    count += 1

# row 9
count = 0
arrayI[8:9, 8:9] += 2.9e6
arrayI[8:9, 9:10] += 2.7e6
for i in range (0,8):
    value = 0
    arrayI[8:9, count:count + 1] += value
    count += 1

# row 10
count = 0
arrayI[9:10, 9:10] += 3e6
for i in range (0,9):
    value = 0
    arrayI[9:10, count:count + 1] += value
    count += 1

arrayII = np.concatenate((arrayI, np.flip(arrayI, axis=0)), axis=0)
arrayIIb = np.delete(arrayII, 11, 0)
arrayIIc = np.delete(arrayIIb, 10, 1)

# rotate along x-axis?
arrayIII = np.concatenate((arrayIIc, np.flip(arrayIIc, axis=1)), axis=1)

# add arrayIII to larger array?
arrayV = np.zeros(shape=(21, 21))
arrayV[0:21, 0:20]+=arrayIII
arrayX = np.insert(arrayV, 10, 0, 1)
arrayY = np.delete(arrayX, 21, 1)

# fill in the column values
final_arr = arrayY
final_arr[0:1, 10:11] += 1
final_arr[1:2, 10:11] += 10
final_arr[2:3, 10:11] += 100
final_arr[3:4, 10:11] += 1e4
final_arr[4:5, 10:11] += 7e4
final_arr[5:6, 10:11] += 2.5e5
final_arr[6:7, 10:11] += 5e5
final_arr[7:8, 10:11] += 1e6
final_arr[8:9, 10:11] += 2e6
final_arr[9:10, 10:11] += 2.5e6
final_arr[10:11, 10:11] += point_source
final_arr[11:12, 10:11] += 2.5e6
final_arr[12:13, 10:11] += 2e6
final_arr[13:14, 10:11] += 1e6
final_arr[14:15, 10:11] += 5e5
final_arr[15:16, 10:11] += 2.5e5
final_arr[16:17, 10:11] += 7e4
final_arr[17:18, 10:11] += 1e4
final_arr[18:19, 10:11] += 100
final_arr[19:20, 10:11] += 10
final_arr[20:21, 10:11] += 1

#******************************* 3) Create Gaussian background component *********************************************#

# add extended emission
# gaussiad2D(peak, x_center, y_center, x_sigma, y_sigma)
# NOTE: X-, y-sigma = fwhm / 2*np.sqrt(2ln(2))
disk = Gaussian2D(3e4, 10, 10, 0.5,1)
#disk2 = Gaussian2DKernel(0.5)
dummy_background = Gaussian2D(30, 495, 488, 200, 200)

# total background component
a = Gaussian2D(175, 128, 128, 100, 100, theta = 180. * np.pi / 180)
a1 = Gaussian2D(155, 128, 128, 25, 65, theta = 75. * np.pi / 180)

# upper and lower wing component
b = Gaussian2D(450, 68, 178, 30, 7, theta = 30. * np.pi / 180)
c = Gaussian2D(300, 158, 88, 30, 9, theta = 175. * np.pi / 180)
c1 = Gaussian2D(200, 185, 95, 15, 7, theta = 75. * np.pi / 180)

# far-lower right extension
d = Gaussian2D(400, 128, 30, 10, 5, theta = 180. * np.pi / 180)
e = Gaussian2D(250, 153, 15, 35, 7, theta = 170. * np.pi / 180)
f = Gaussian2D(250, 175, 30, 15, 7, theta = 75. * np.pi / 180)

# Random extended components
g = Gaussian2D(220, 135, 180, 7, 5, theta = 65. * np.pi / 180)
h = Gaussian2D(220, 128, 220, 25, 6, theta = 165. * np.pi / 180)
i = Gaussian2D(220, 25, 50, 10, 20, theta = 90. * np.pi / 180)
j = Gaussian2D(190, 220, 200, 10, 20, theta = 50. * np.pi / 180)
k = Gaussian2D(210, 35, 80, 5,6, theta = 50. * np.pi / 180)
l = Gaussian2D(250, 35, 110, 25,6, theta = 65. * np.pi / 180)
m = Gaussian2D(200, 67, 74, 5,4, theta = 65. * np.pi / 180)
n = Gaussian2D(175, 224, 134, 15,4, theta = 85. * np.pi / 180)
n1 = Gaussian2D(205, 70, 230, 17,9, theta = 85. * np.pi / 180)

# components to remove
o = Gaussian2D(175, 180, 125, 40,10, theta = 65. * np.pi / 180)
p = Gaussian2D(155, 69, 122, 75,10, theta = 60. * np.pi / 180)
q = Gaussian2D(120, 78, 29, 15,10, theta = 60. * np.pi / 180)
r = Gaussian2D(125, 205, 135, 10,10, theta = 65. * np.pi / 180)

# build the Gaussian image
ny = nx = 20
y1, x1 = np.mgrid[0:ny, 0:nx]

# function to add background component
g_ext = disk(x1,y1)

# function to add extended components
g_rem = o(x1, y1) + p(x1, y1) + q(x1, y1) + r(x1, y1)

# Generate a dummy background to subtract from galaxy background
ny1 = nx1 = 256
y2, x2 = np.mgrid[0:982, 0:1017]

# function to add background componenet
g_ext1 = dummy_background(x2,y2)

#******************************** 4) Create galaxy emission component ************************************************#

# Import the background galaxy model
model = 'ngc5728_HLA.fits'

# read in data
A = get_pkg_data_filename(path + 'Galaxy_Models/'+model)
A_list = fits.open(A)
galaxy = A_list[1].data
print(galaxy.shape)

# subtract dummy background?
galaxy = galaxy - g_ext1
galaxy[galaxy < 0] = 1

rot = ndimage.rotate(galaxy, 90, reshape=True)
galaxy[galaxy<0] = 0

# scale
galaxy = rot[395:651, 361:617]

# rescale galaxy: HST -> JWST
galaxy2 = resize_psf(galaxy, input_pixel_scale = 0.13, output_pixel_scale = 0.111)
galaxy2[galaxy2<0] = 0

# resize rescaled image
galaxy = galaxy2[20:276, 20:276]
galaxy[galaxy<0] = 0

# 2b) scale?
norm_galaxy = galaxy/galaxy.max()
galaxy = norm_galaxy*1000

# mask negative pixels
galaxy[galaxy<0] = 0

#*************************************** 5) save FITS output ********************************************************#
squareE = square

if name_idx == 0:
    # add bicone + point source + galaxy emission only
    # add point source to galaxy image
    squareE[512:513, 511:512] += point_source

elif name_idx == 1:
    # add bicone + point source + galaxy emission only
    # add point source to galaxy image
    squareE[511:512, 512:513] += point_source
    # add disk
    squareE[501:521, 502:522] += g_ext

elif name_idx == 2:
    # Add bicone #1 to galaxy image
    final_product2 = final_product2[5:27, 5:27]
    # shift the cone by 1/2 a pixel down and to the left?
    final_product2 = shift(final_product2, (0.5, -0.5), mode='nearest')
    squareE[501:523, 501:523] += final_product2

    # add bicone #2 -> comment above if user-selected here
    #squareE[501:522, 501:522] += final_arr

    # add disk
    squareE[502:522, 501:521] += g_ext
    # add point source
    squareE[511:512, 511:512] += point_source

elif name_idx == 3:
    # Add bicone #1 to galaxy image -> subpixel shift?
    final_product2 = final_product2[5:27, 5:27]
    # shift the cone by 1/2 a pixel down and to the left?
    final_product2 = shift(final_product2, (0.5, -0.5), mode='nearest')
    squareE[501:523, 501:523] += final_product2

    # add bicone #2 -> comment above if user-selected here
    #squareE[501:522, 501:522] += final_arr
    squareE[501:521, 501:521] += g_ext

    # add galaxy background
    squareE[511:512, 511:512] += point_source
    squareE[384:640, 384:640] += galaxy

save_FITS(squareE, filename = None, name = str(name[name_idx]), output_dir = path)


