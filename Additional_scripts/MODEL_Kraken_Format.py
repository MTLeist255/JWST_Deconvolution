# Nov 7 2022: simple code to convert 256x256 -> 1024x1024 padded arrays to work with the Kraken MFBD MATLAB code

from Convenience_Functions import *

# set the model type
model_type = 'Ideal_Model'
# set the model name(s)
sample = ['standard']
# set the simulated filter(s)
filter = ['F560W', 'F1000W', 'F1500W', 'F1800W', 'F2100W']

# 2. read in data by model type
for i in range(0, int(len(sample))):
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
        mirisim_path = 'Images/JWST/2_MIRISim_outputs/Deconvolution_testing/Residual_outputs/FITS/' + str(model_type) + '/Model_AGN_complicated/current_12Aug2023/' + str(sample[i]) + '/'
        out_path = 'Images/JWST/5_ERS_Data/5_Deconvolution_results/NGC5728/Deconvolution/Padded_Images/toy_AGN_model/'
        # read in the image data
        image = get_pkg_data_filename(mirisim_path + str(filter[j]) + '_' + str(sample[i]) + '_model.fits')
        # get the necessary image data
        image_list = fits.open(image)
        # copy the header information: PRIMARY + SCI
        #prime_hdr = image_list[0].header.copy()
        sci_hdr = image_list[0].header.copy()
        # get the image data
        image_data = image_list[0].data
        print('#' + str(i) + ') ' + ' model, filter: ' + str(sample[i]) +' | ' + str(filter[j]))

        # 3. add the model data to the new array
        # NOTE: this is based on the objects centroid position (rounded up to the nearest whole pixel value)
        square[384:640, 384:640] += image_data

        # 4. save the FITS header information for each filter
        # save outputs to single FITS file with appropriate FITS header information
        base_filename = str(filter[j]) + '_' + str(sample[i]) + '_model_1024.fits'
        outfile = os.path.join(out_path, base_filename)
        hdu = fits.HDUList()

        sci_hdr[''] = ' Modified for Kraken MFBD deconvolution'
        sci_hdr['NAXIS'] = (1024, 'New image array x-axis length')
        sci_hdr['NAXIS'] = (1024, 'New image array y-axis length')

        # append the simulated image
        # append new header information -> single image
        hdu.append(fits.PrimaryHDU(square, header=sci_hdr))
        hdu.writeto(outfile, overwrite=True)
        image_list.flush()