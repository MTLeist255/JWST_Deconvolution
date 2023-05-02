# 8 mar 2023: format the JWST/MIRI observations to the Kraken input format

from Convenience_Functions import *

# set the input and output path
im_path = 'Images/JWST/5_ERS_Data/4_STAGE3_outputs/NGC5728_obs2/'
out_path = im_path + 'Kraken_Formatted/'

# set the observation filter
filter = 'F2100W'

# read in the MIRI data
image = get_pkg_data_filename(im_path + 'NGC-5728-2_'+str(filter)+'_final.fits')
image_list = fits.open(image)
image_data = image_list[1].data

# save the FITS header info
prime_hdr = image_list[0].header.copy()
sci_hdr = image_list[1].header.copy()

# crop the image array to 256x256
print('Image shape before trim: ', image_data.shape)
# reshpae image array: (y,x)
trim_data = image_data[50:306, 42:298]
print('Image shape after trim: ', trim_data.shape)

# pad negative pixel values = 1
trim_data[trim_data <0] += 0
pad_data = trim_data

# save the formatted image to a FITS file
base_filename = 'PAD_NGC5728_' + str(filter) + '_obs2.fits'
outfile = os.path.join(out_path, base_filename)
hdu = fits.HDUList()

sci_hdr['NAXIS'] = (256, 'New image array x-axis length')
sci_hdr['NAXIS'] = (256, 'New image array y-axis length')

# append the simulated image
# append new header information -> single image
# append new header information -> single image
hdu.append(fits.PrimaryHDU(pad_data, header=prime_hdr))
hdu.append(fits.ImageHDU(header=sci_hdr))
hdu.writeto(outfile, overwrite=True)
image_list.flush()