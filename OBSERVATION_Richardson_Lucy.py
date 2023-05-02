# 2022 Nov 19: code to deconvolved JWST/MIRI images using non-circulant Richardson-Lucy

from Convenience_Functions import *

# display the deconvolved image?
display = False

# set the number of iterations
iter = 251

# set the observation title
obs = ['F2100W']

# set the import and output paths
im_path = 'Images/JWST/5_ERS_Data/4_STAGE3_outputs/NGC5728_obs1/RL_Formatted/'
out_path = 'Images/JWST/4_Deconvolution/1_Richardson_Lucy/NGC5728/obs1/FITS/'

# read in the image data
# copy the FITS header info
for i in range(0, int(len(obs))):
    # Read in the observed image
    image = get_pkg_data_filename(im_path + 'RL_NGC5728_' + str(obs[i]) + '_obs1.fits')

    # read in image data for deconvolution and saving
    image_list = fits.open(image)
    image_data = image_list[1].data

    # copy the primary header information
    prime_hdr = image_list[0].header.copy()

    # copy the SCI extension header information
    sci_hdr = image_list[1].header.copy()

    # read the filter used and set the reference PSF accordingly
    filter = image_list[0].header['FILTER']
    if filter == 'F560W':
        obs_wave = 5
    elif filter == 'F1000W':
        obs_wave = 10
    elif filter == 'F1500W':
        obs_wave = 15
    elif filter == 'F1800W':
        obs_wave = 18
    else:
        obs_wave = 21

    print('Check reference PSF is correct: ', obs_wave)

    # Read in the reference PSF
    ref_psf = get_pkg_data_filename(im_path + 'RL_reference_PSF_' + str(obs[i]) + '_DET_DIST_FULL_original2.fits')

    # read in image data for deconvolution
    psf_list = fits.open(ref_psf)
    # ext 0- SAMP_PSF: PSF binned to detector sampling
    # ext 1- DET_SAMP: PSF binned to detector sampling
    # ext 2- OVERDIST: PSF modified based on detector geometric distortions
    # ext 3- DET_DIST: PSF modified based on detector geometric distortions binned to detector sampling <- USE ME
    psf = psf_list[3].data

    # For me -> keep track of progress
    plt.subplot(121)
    plt.imshow(image_data, origin='lower', norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=image_data.max()),cmap='turbo')
    plt.title('NGC5728 ' + str(obs[i]) + ' image')
    plt.subplot(122)
    plt.imshow(psf, origin='lower', norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=psf.max()),cmap='turbo')
    plt.title(str(obs[i]) + ' reference PSF')
    plt.show()

    # Deconvolve the observed image using Richardson-Lucy, saving the deconvolved images
    decon_im = []

    # set the number of deconvolution iterations: range(starting value, ending value, incrementing value)
    for j in range(1, iter, 1):
        # richardson_lucy_np(observation, PSF, regularization parameter)
        deconvolution = richardson_lucy_np(image_data, psf, j)
        print(obs[i] + ' iter #: ' + str(j))

        # append the deconvolved image
        decon_im.append(deconvolution)

    # save the observed image and deconvolved images to a single FITS file
    final_arr = [image_data] + decon_im

    # save outputs to single FITS file with appropriate FITS header information
    base_filename = 'RL_NGC5728_' + str(obs[i]) + '_obs1_DET_SAMP_'+str(coord)+'.fits'
    outfile = os.path.join(out_path, base_filename)
    hdu = fits.HDUList()

    count = 0
    for k in image:

        if count == 0:

            # copy the primary header information
            prime_hdr = image_list[0].header.copy()

            # copy the SCI extension header information
            sci_hdr = image_list[1].header.copy()
            sci_hdr['NAXIS1'] = (130, 'x-axis limits')
            sci_hdr['NAXIS2'] = (130, 'y-axis limits')


            sci_hdr[''] = ' and Deconvolution Parameters'
            sci_hdr['METHOD'] = ('Richardson-Lucy', 'Deconvolution method')
            sci_hdr['ITERATS'] = (str(iter), 'Number of iterations')
            # append the observed image
            # append new header information -> single image
            hdu.append(fits.PrimaryHDU(header=prime_hdr))
            hdu.append(fits.ImageHDU(final_arr[count], header=sci_hdr))
            hdu.writeto(outfile, overwrite=True)
            image_list.flush()

        else:
            # append each deconvolved image
            hdu.append(fits.ImageHDU(image[count]))

        count += 1

    hdu.writeto(outfile, overwrite=True)

    if display:
        for l in range(0, int(len(final_arr)), 10):
            image1 = image[l]
            title = str('iteration: ' + str(l))
            im = plt.imshow(image1, origin='lower',norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=image1.max()),cmap='viridis')
            plt.title(title)
            plt.xlim(0, 130)
            plt.ylim(0,130)
            add_colorbar(im, label='Flux density (mJy/sr)')
            plt.savefig(filter + '_iter_' + str(l))
            plt.show()


    # close image imports
    image_list.close()
    psf_list.close()