# 20 July 2023: simple code to covert the individual Kraken deconvolved FITS files form 1024x1024 images arrays -> 256x256

from Convenience_Functions import *

# set the model type
model_type = 'Ideal_Model'
# set the model name(s)
model = ['TEST3_NOISEx2', 'TEST4_NOISEx3', 'TEST5_NOISEx4', 'TEST6_NOISEx5', 'TEST7_NOISEx6']
# set the simulated filter(s)
filter = ['F560W', 'F1000W', 'F1500W', 'F1800W', 'F2100W']

# 2. read in data by model type
for i in range(0, int(len(model))):
    # 2a. read in data by filter
    for j in range(0, int(len(filter))):
        # 2. generate a 1024x1024 array of zeros
        naxis = (256, 256)
        square = np.zeros(naxis)
        # set the simulated wavelength
        if j == 0:
            obs_wave = 5.6
        elif j == 1:
            obs_wave = 10.0
        elif j == 2:
            obs_wave = 15.0
        elif j == 3:
            obs_wave = 18.0
        else:
            obs_wave = 21.0

        # set the import path for the model image FITS file
        model_path = 'Images/JWST/2_MIRISim_outputs/Deconvolution_testing/Residual_outputs/FITS/' + str(model_type) + '/Model_AGN_complicated/20July2023_TESTS/' + str(model[i]) + '/'
        # set the output path to save the modifed FITS file
        out_path = 'Images/JWST/2_MIRISim_outputs/Deconvolution_testing/Residual_outputs/FITS/' + str(model_type) + '/Model_AGN_complicated/20July2023_TESTS/' + str(model[i]) + '/'
        # read in the image data
        image = get_pkg_data_filename(model_path + 'Kraken_' + str(filter[j]) + '_DECONV.fits')
        # get the necessary image data
        image_list = fits.open(image)
        # copy the header information: PRIMARY + SCI
        #prime_hdr = image_list[0].header.copy()
        sci_hdr = image_list[0].header.copy()
        # get the image data
        image_data = image_list[0].data
        print('#' + str(i) + ') ' + ' model, filter: ' + str(model[i]) +' | ' + str(filter[j]))

        # 3. resize the old image and add to a new array
        # NOTE: this is based on the objects centroid position (rounded up to the nearest whole pixel value)
        trim_image_data = image_data[384:640, 384:640]
        square[0:256, 0:256] += trim_image_data

        # 4. save the FITS header information for each filter
        # save outputs to single FITS file with appropriate FITS header information
        base_filename = str(filter[j]) + '_DECONV.fits'
        outfile = os.path.join(out_path, base_filename)
        hdu = fits.HDUList()

        sci_hdr['BUNIT'] = 'DN'
        sci_hdr[''] = ' Modified for Kraken MFBD deconvolution'
        sci_hdr['PIXELSCL'] = (0.1110, 'pixel scale')
        sci_hdr['WAVELEN'] = (obs_wave, 'Observed wavelength')
        sci_hdr['FILTER'] = (str(filter[j]), 'Observation filter')
        sci_hdr['INSTRUM'] = ('MIRI', 'Observation instrument')
        sci_hdr['NAXIS'] = (256, 'New image array x-axis length')
        sci_hdr['NAXIS'] = (256, 'New image array y-axis length')
        sci_hdr['TEST'] = (1, str(model[i]))

        # append the simulated image
        # append new header information -> single image
        hdu.append(fits.PrimaryHDU(square, header=sci_hdr))
        #hdu.append(fits.ImageHDU(header=sci_hdr))
        hdu.writeto(outfile, overwrite=True)
        image_list.flush()