# 1 Feb 2022: Source code to format the input model into the MIRISim requirements
# This code gives the minimal example to create an input .fits file. Populate the DATA array with the desired data
# Original Creator: Christophe COSSOU
# source: https://wiki.miricle.org//pub/Internal/Software/ChristopheCossou/MIRISim_in_Python_userguide.pdf

#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Make a minimal SkyCube that works
'''
from Convenience_Functions import *

# set the input path
path = 'Images/JWST/1_Input_Models/Models/'
wave = 10
name = ['Model_'+str(wave)+'um_single_pixel',
        'Model_'+str(wave)+'um_single_pixel_disk',
        'Model_'+str(wave)+'um_single_pixel_disk_bicone',
        'Model_'+str(wave)+'um_AGN_complicated',
        'Model_'+str(wave)+'um_galaxy_test',
        'Model_Background'
        ]

# set name index
name_idx = 3

# Read in FITS input model -> bicone
model = get_pkg_data_filename(path + name[name_idx] + '.fits')
model_list = fits.open(model)
model_data = model_list[0].data
model_list.close()
#
# # Reshape the 2D array to 3D
data3d = model_data.reshape((model_data.shape[0], model_data.shape[1],-1))

# Establish array size
# (z, y, x) = (y_lim, wavelength coverage, x_lim)
# wavelength coverage -> the number of different wavelengths for MIRISim to simulate. This number is iterated using
# dwave as a reference, i.e., dwave = 0.5, z = 40 -> wavelength coverage between 5-25 um, incrementing every 0.25 um
naxis = (1024,1024,40)

# Populate data array -> ADD TO THIS LINE WITH THE DESIRED DATA
data = np.zeros(naxis)
# read-in data
data[0:1024, 0:1024, 0:40] = data3d

# size of a pixel in arcsec
dpix = 0.1110

# Micron, must be constant throughout spectrum
dwave = 0.5

# Index of the reference pixel -> starts at 1
crpix = (0.5, 0.5, 0.5)

# value of the reference pixel
# Axis 1 & 2 ar in RA & Dec
# axis 3 is in micron
# FOV: (RA, Dec, Micron)
crval = (-57., -57., 5)

# Pixel size: RA, Dec, Wavelength
cdelt = (dpix, dpix, dwave)

# Prepare data to correct dimension order
data = np.moveaxis(data, -1, 0)
# flips the data top -> bottom
data = np.flip(data, axis =1)

# Write data to FITS
hdu = fits.PrimaryHDU(data)

# create necessary headers -> RA
# position of the reference pixel in the skycube
# coordinate system
hdu.header['CRPIX1'] = crpix[0]
# position of the reference pixel on the sky
hdu.header['CRVAL1'] = crval[0]
# size of the pixel
hdu.header['CDELT1'] = cdelt[0]
# unit for each dimension
hdu.header['CUNIT1'] = u.arcsec.to_string()
hdu.header['CTYPE1'] = 'RA---TAN'

# create necessary headers -> DEC
hdu.header['CRPIX2'] = crpix[1]
hdu.header['CRVAL2'] = crval[1]
hdu.header['CDELT2'] = cdelt[1]
hdu.header['CUNIT2'] = u.arcsec.to_string()
hdu.header['CTYPE2'] = 'DEC--TAN'

# create necessary headers -> WAVE
hdu.header['CRPIX3'] = crpix[2]
hdu.header['CRVAL3'] = crval[2]
hdu.header['CDELT3'] = cdelt[2]
hdu.header['CUNIT3'] = u.um.to_string()
hdu.header['CTYPE3'] = 'WAVE'

hdu.header['UNITS'] = (u.microjansky/u.arcsec**2).to_string()

# poppy compliant
hdu.header['PIXELSCL'] = 0.1110

# Title of output FITS file
outFitsName = name[name_idx] + '.fits'
hdu1 = fits.HDUList([hdu])

# define output directory
output_dir = '../../../mirisim/FITS/'
#output_dir = 'Images/JWST/Inputs/Models/Final_Models'
outFits = os.path.join(output_dir, outFitsName)
# create the header file
# create the FITS file and send to mirisim directory
hdu1.writeto(outFits, overwrite = True)