# Nov 7 2022: a simple code to combine individual FITS files iteratively

from Convenience_Functions import *

# set the number to start iterating through (i.e., skip how many at the start?)
skip = 0
# set the iterations number
iter = 150

# set the directory
Dir = '5_Kraken'

# set the model type
model_type = 'Ideal_Model'

# set the model
model = 'Model_AGN_complicated'

# set the observation filter
filter = 'F2100W'

# Set the FITS import path for the MIRISim image
mirisim_path = 'Images/JWST/2_MIRISim_outputs/Deconvolution_testing/Residual_outputs/FITS/' +str(model_type) + '/' + str(model)

# Set the FITS import path for the deconvolution data
fits_path = 'Images/JWST/4_Deconvolution/'+str(Dir)+'/FITS/'+str(model_type) + '/RAW_FITS/deconv_tests/psfs/21new_' + str(filter) + '/kraken_' +str(filter)+ '_complex_vars_deconv/ch1_hvar2/'

# set the output path for the combined FITS files
out_path = 'Images/JWST/4_Deconvolution/'+str(Dir)+'/FITS/'+str(model_type)+'/' +str(model) + '/21new/'

# read in the MIRISim image, saving the primary and science FITS header
image = get_pkg_data_filename(mirisim_path + '/SCALED_' + str(filter) + '_' + str(model) + '.fits')
# get the necessary image data
image_list = fits.open(image)
# get the image data
image_data = image_list[0].data

# set the simulated wavelength arbitrarily
if filter == 'F560W':
    obs_wave = 5.6
elif filter == 'F1000W':
    obs_wave = 10.0
elif filter == 'F1500W':
    obs_wave = 15.0
elif filter == 'F1800W':
    obs_wave = 18.0
else:
    obs_wave = 21.0

# set the deconvolution output name
if Dir == '4_AIDA':
    method = 'AIDA'
elif Dir == '5_Kraken':
    method = 'Kraken MFBD'
elif Dir == '6_IWFT':
    method = 'Iterative Wiener Filtering and Thresholding (IWFT)'

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
    if Dir == '5_Kraken':
        # note update this based on object centroid
        decon = decon_data[385:641, 385:641]
    else:
        decon = decon_data

    # append the image data
    decon_im.append(decon)


# Append the MIRISim image to the start of this array
img_array = [image_data] + decon_im

# Save all images in the array to a single FITS file
# save outputs to single FITS file with appropriate FITS header information
base_filename = 'SCALED_' + str(filter) + '_' + str(model) + '_DET_SAMP.fits'
outfile = os.path.join(out_path, base_filename)
hdu = fits.HDUList()

count = 0
for k in img_array:
    # copy the primary header information
    prime_hdr = image_list[1].header.copy()
    sci_hdr = image_list[0].header.copy()

    if count == 0:
        hdu.append(fits.PrimaryHDU(header=prime_hdr))

    elif count == 1:
        # copy the SCI extension header information
        sci_hdr['BUNIT'] = 'DN'
        sci_hdr[''] = ' and ' + str(method) + ' Deconvolution Parameters'
        sci_hdr['PIXELSCL'] = (0.1110, 'pixel scale')
        sci_hdr['WAVELEN'] = (str(obs_wave), 'Observed wavelength')
        sci_hdr['FILTER'] = (str(filter), 'Observation filter')
        sci_hdr['INSTRUM'] = ('MIRI', 'Observation instrument')
        sci_hdr['SAMP'] = ('DET_SAMP', 'PSF extension used for deconvolution')
        sci_hdr['METHOD'] = (str(method), 'Deconvolution method')
        sci_hdr['ITERATS'] = (str(int(iter - skip)), 'Number of iterations')
        hdu.append(fits.ImageHDU(img_array[count-1], header=sci_hdr))

    else:
        # append each deconvolved image
        hdu.append(fits.ImageHDU(img_array[count-1], header=sci_hdr))

    count += 1

hdu.writeto(outfile, overwrite=True)
image_list.flush()





