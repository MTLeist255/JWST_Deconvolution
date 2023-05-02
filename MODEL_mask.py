# 16 Nov 2022: simple code to read in model images, mask bad pixels, and save output FITS file

from Convenience_Functions import *

# set the model type
model_type = 'Ideal_Model'
# set the model name
model = 'Model_AGN_complicated'
# set the simulated filter
filter = 'F1800W'
# set the simulated wavelength
obs_wave = 18.0

# set the import path for the model image FITS file
model_path = 'Images/JWST/2_MIRISim_outputs/Deconvolution_testing/Residual_outputs/FITS/' + str(model_type) + '/' + str(model) + '/'
# read in the image data
image = get_pkg_data_filename(model_path + str(filter) + '_' + str(model) + '.fits')
# get the necessary image data
image_list = fits.open(image)
# copy the header information: PRIMARY + SCI
prime_hdr = image_list[0].header.copy()
sci_hdr = image_list[1].header.copy()
# get the image data
image_data = image_list[1].data

# mask bad pixels
x_arr = [0]
y_arr = [0]
for i in range(0, int(len(x_arr))):
    # change bad pixels to 0's
    image_data[y_arr[i]:y_arr[i]+1, x_arr[i]:x_arr[i]+1] = 0

# save the FITS header information for each filter
# save outputs to single FITS file with appropriate FITS header information
base_filename = 'MASKED_' + str(filter) + '_' + str(model) + '.fits'
outfile = os.path.join(model_path, base_filename)
hdu = fits.HDUList()

sci_hdr['BUNIT'] = 'DN'
sci_hdr[''] = ' Modified for Kraken MFBD deconvolution'
sci_hdr['PIXELSCL'] = (0.1110, 'pixel scale')
sci_hdr['WAVELEN'] = (obs_wave, 'Observed wavelength')
sci_hdr['FILTER'] = (str(filter), 'Observation filter')
sci_hdr['INSTRUM'] = ('MIRI', 'Observation instrument')
sci_hdr['METHOD'] = ('Richardson-Lucy', 'Deconvolution method')
sci_hdr['NAXIS1'] = (256, 'X-axis length')
sci_hdr['NAXIS2'] = (256, 'Y-axis length')

# append the simulated image
# append new header information -> single image
hdu.append(fits.PrimaryHDU(image_data, header=prime_hdr))
hdu.append(fits.ImageHDU(header=sci_hdr))
hdu.writeto(outfile, overwrite=True)
image_list.flush()

