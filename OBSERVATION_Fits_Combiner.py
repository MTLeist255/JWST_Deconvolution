# 20 Nov 2022: code to combine the individual Kraken deconvolved FITS files for the JWST/MIRI observations

from Convenience_Functions import *

# set the iterations number
iter = 50
# set the number to start iterating through (i.e., skip how many at the start?)
skip = 0
# Set the observation filter used
filter = 'F2100W'
# set the observation title
obs = 'PAD_NGC5728_'+str(filter)+'_obs2'
# set the import and output paths
obs_path = 'Images/JWST/5_ERS_Data/4_Stage3_Outputs/NGC5728_obs2/'
fits_path = 'Images/JWST/4_Deconvolution/5_Kraken/NGC5728/obs2/FITS/RAW_FITS/deconv_tests/psfs/' + str(filter) + '/kraken_'+str(filter)+'_complex_vars_deconv/ch1_hvar2/'
out_path = 'Images/JWST/4_Deconvolution/5_Kraken/NGC5728/obs2/FITS/'

# Read in the observed image
image = get_pkg_data_filename(obs_path + '/kraken_Formatted/' + obs + '.fits')
# read in image data for deconvolution and saving
image_list = fits.open(image)
image_data = image_list[0].data
# copy the primary header information
prime_hdr = image_list[0].header.copy()
# copy the SCI extension header information
sci_hdr = image_list[1].header.copy()

# Iterate through each FITS file in the deconvolution directory, appending the image data to a single array
decon_im = []

for i in range(skip, iter):
    # read in the image data
    decon = get_pkg_data_filename(fits_path + 'object_iter' + str(i) + '.fits')
    # get the necessary image data
    decon_list = fits.open(decon)
    # get the image data
    decon_data = decon_list[0].data

    # resize images -> ONLY FOR KRAKEN: 1024x1024 -> 256x256
    decon_trim = decon_data[385:641, 385:641]

    # append the image data
    decon_im.append(decon_trim)

# Append the MIRISim image to the start of this array
img_array = [image_data] + decon_im

# Save all images in the array to a single FITS file
# save outputs to single FITS file with appropriate FITS header information
base_filename = 'Kraken_update_NGC5728_'+str(filter)+'_obs2_DET_SAMP.fits'
outfile = os.path.join(out_path, base_filename)
hdu = fits.HDUList()

count = 0
for k in img_array:

    if count == 0:
        # append the primary header
        hdu.append(fits.PrimaryHDU(header=prime_hdr))

    elif count == 1:
        # append the original image and science header
        sci_hdr = image_list[1].header.copy()

        # Append the primary header info to extension 1
        sci_hdr['BUNIT'] = 'mJy/sr'
        sci_hdr[''] = ' and Kraken MFBD Deconvolution Parameters'
        sci_hdr['SAMP'] = ('DET_SAMP', 'PSF extension used for deconvolution')
        sci_hdr['METHOD'] = ('Kraken MFBD', 'Deconvolution method')
        sci_hdr['ITERATS'] = (str(iter), 'Number of iterations')
        hdu.append(fits.ImageHDU(img_array[count-1], header=sci_hdr))
        fits.writeto(out_path + 'test.fits', img_array[count], overwrite=True)
    else:
        # append each deconvolved image
        hdu.append(fits.ImageHDU(img_array[count], header=sci_hdr))

    count += 1

hdu.writeto(outfile, overwrite=True)





