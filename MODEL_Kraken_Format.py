# Nov 7 2022: simple code to convert 256x256 -> 1024x1024 padded arrays to work with the Kraken MFBD MATLAB code

from Convenience_Functions import *

# set the model type
model_type = 'Ideal_Model'
# set the model name(s)
model = ['Model_single_pixel', 'Model_single_pixel_disk', 'Model_single_pixel_disk_bicone', 'Model_AGN_complicated']
# set the simulated filter(s)
filter = ['F560W', 'F1000W', 'F1500W', 'F1800W', 'F2100W']

# 2. read in data by model type
for i in range(0, int(len(model))):
    # 2a. read in data by filter
    for j in range(0, int(len(filter))):
        # 2. generate a 1024x1024 array of zeros
        naxis = (1024, 1024)
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
        model_path = 'Images/JWST/2_MIRISim_outputs/Deconvolution_testing/Residual_outputs/FITS/' + str(model_type) + '/' + str(model[i]) + '/'
        # set the output path to save the modifed FITS file
        out_path = 'Images/JWST/2_MIRISim_outputs/Deconvolution_testing/Residual_outputs/FITS/' + str(model_type) + '/Kraken_Modified/' + str(model[i]) + '/'
        # read in the image data
        image = get_pkg_data_filename(model_path + 'SCALED_' + str(filter[j]) + '_' + str(model[i]) + '.fits')
        # get the necessary image data
        image_list = fits.open(image)
        # copy the header information: PRIMARY + SCI
        prime_hdr = image_list[0].header.copy()
        sci_hdr = image_list[1].header.copy()
        # get the image data
        image_data = image_list[1].data
        print('model, filter: ', model[i], filter[j])

        # 3. add the model data to the new array
        # NOTE: this is based on the objects centroid position (rounded up to the nearest whole pixel value)
        square[384:640, 384:640] += image_data

        # 4. save the FITS header information for each filter
        # save outputs to single FITS file with appropriate FITS header information
        base_filename = 'Kraken_' + str(filter[j]) + '_' + str(model[i]) + '.fits'
        outfile = os.path.join(out_path, base_filename)
        hdu = fits.HDUList()

        sci_hdr['BUNIT'] = 'DN'
        sci_hdr[''] = ' Modified for Kraken MFBD deconvolution'
        sci_hdr['PIXELSCL'] = (0.1110, 'pixel scale')
        sci_hdr['WAVELEN'] = (obs_wave, 'Observed wavelength')
        sci_hdr['FILTER'] = (str(filter[j]), 'Observation filter')
        sci_hdr['INSTRUM'] = ('MIRI', 'Observation instrument')
        sci_hdr['NAXIS'] = (1024, 'New image array x-axis length')
        sci_hdr['NAXIS'] = (1024, 'New image array y-axis length')

        # append the simulated image
        # append new header information -> single image
        hdu.append(fits.PrimaryHDU(square, header=prime_hdr))
        hdu.append(fits.ImageHDU(header=sci_hdr))
        hdu.writeto(outfile, overwrite=True)
        image_list.flush()