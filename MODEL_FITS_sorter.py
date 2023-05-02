# 12 Feb 2023: short file to read in  a FITS cube and exclude a user-defined number of FITS entries, then resave
# the data into a new array

from Convenience_Functions import *

# set the iterations number
iter = 125
# set the number to start iterating through (i.e., skip how many at the start?)
skip = 0

# set the observation filter
filter = 'F2100W'

# set the directory
Dir = '7_Sparse_Condat_Vu'

# set the model type
model_type = 'Ideal_Model'

# set the model
model = 'Model_AGN_complicated'

# Set the FITS import path for the MIRISim image
mirisim_path = 'Images/JWST/2_MIRISim_outputs/Deconvolution_testing/Residual_outputs/FITS/'+str(model_type)+'/' + str(model) + '/'

# Set the FITS import path for the deconvolution data
fits_path = 'Images/JWST/4_Deconvolution/'+str(Dir)+'/FITS/'+str(model_type) + '/' + str(model) + '/Original'

# set the output path for the combined FITS files
out_path = 'Images/JWST/4_Deconvolution/'+str(Dir)+'/FITS/'+str(model_type)+'/' +str(model)

# read in the MIRISim image, saving the primary and science FITS header
# read in the image data
image = get_pkg_data_filename(mirisim_path + 'SCALED_' + str(filter) + '_' + str(model) + '.fits')
# get the necessary image data
image_list = fits.open(image)
# copy the header information: PRIMARY + SCI
prime_hdr = image_list[1].header.copy()
sci_hdr = image_list[0].header.copy()
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
elif Dir == '7_Sparse_Condat_Vu':
    method = 'Sparse Condat-Vu'

# Iterate through each FITS file in the deconvolution directory, appending the image data to a single array
decon_im = []

# read in the deconvolution data array
decon = get_pkg_data_filename(fits_path + '/Scaled_' + str(filter) + '_' + str(model) + '_DET_SAMP.fits')
# get the necessary image data
decon_list = fits.open(decon)

for i in range(skip, iter):
    # get the image data
    decon_data = decon_list[i].data
    decon = decon_data


    if i == 1+skip:
        print('This iteration skipped: #', i)
    else:
        # append the image data
        decon_im.append(decon)

# Append the MIRISim image to the start of this array
img_array = [image_data] + decon_im

# Save all images in the array to a single FITS file
# save outputs to single FITS file with appropriate FITS header information
base_filename = 'SCALED_'+ str(filter) + '_' + str(model) + '_DET_SAMP.fits'
outfile = os.path.join(out_path, base_filename)
hdu = fits.HDUList()

count = 0
for k in img_array:

    if count == 0:
        # append the simulated image
        hdu.append(fits.PrimaryHDU(header=prime_hdr))
    elif count == 1:
        # Append the primary header info to extension 0
        sci_hdr['BUNIT'] = 'DN'
        sci_hdr[''] = ' and '+str(method)+' Deconvolution Parameters'
        sci_hdr['PIXELSCL'] = (0.1110, 'pixel scale')
        sci_hdr['WAVELEN'] = (str(obs_wave), 'Observed wavelength')
        sci_hdr['FILTER'] = (str(filter), 'Observation filter')
        sci_hdr['INSTRUM'] = ('MIRI', 'Observation instrument')
        sci_hdr['SAMP'] = ('DET_SAMP', 'PSF extension used for deconvolution')
        sci_hdr['METHOD'] = (str(method), 'Deconvolution method')
        sci_hdr['ITERATS'] = (str(int(iter-skip)), 'Number of iterations')

        # Append the science header info and MIRISim image data to the extension 1
        hdu.append(fits.ImageHDU(img_array[count-1], header=sci_hdr))

    else:
        # append each deconvolved image
        hdu.append(fits.ImageHDU(img_array[count-1], header=sci_hdr))

    count += 1

hdu.writeto(outfile, overwrite=True)
image_list.flush()



