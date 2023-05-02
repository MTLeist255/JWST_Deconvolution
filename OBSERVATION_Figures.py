# 16 Mar 2023: code to plot various figures of JWST/MIRI data

from Convenience_Functions import *

# Plot the merit function grid
merit_grid = False
# normalize the aperture flux
normalize = True
# plot the aperture flux as a delta function
# NOTE: only one of these can be toggled at a time
delta = False
# toggle the diffraction limit display?
diff_lim = True

# Plot the SED grid
SED = False

# plot the observed vs deconvolved image grid
image_grid = False

# lot the residual image grid
residual_grid = False

# plot RGN images
RGB = False

# set the deconvolution iterations globally
merit_iter = [6,18,19,20,18]

# set the x-lim for the merit_grid plot
xlim = 30

# arbitraily set the filters used
filter = ['F560W', 'F1000W', 'F1500W', 'F1800W', 'F2100W']

# set the observed wave
obs_wave = [5.6, 10.0, 15.0, 18.0, 21.0]

# set the deconvolution type directory
dir = '5_Kraken'
# set the output path
out_path = 'Images/JWST/4_Deconvolution/'+str(dir)+'/NGC5728/obs2/Images/'

# loop through image data, rotate by proper WCS coordinates, trim, and save for later plotting
obs_imgs = []
dec_imgs = []
cenx_obs = []
ceny_obs = []
cenx_dec = []
ceny_dec = []
for i in range(0, int(len(filter))):
    # read in the image data
    # set the import path
    im_path = 'Images/JWST/4_Deconvolution/'+str(dir)+'/NGC5728/obs2/FITS/28Mar2023_Kraken_old/'
    image = get_pkg_data_filename(im_path + 'kraken_NGC5728_'+str(filter[i])+'_obs2_DET_SAMP_XY.fits')
    img_list = fits.open(image)
    img_data = img_list[1].data

    # copy the primary header information
    prime_hdr = img_list[0].header.copy()
    # copy the SCI extension header information
    sci_hdr = img_list[1].header.copy()

    # set the deconvolved image data
    dec_data = img_list[merit_iter[i]+1].data

    # save each image to a dummy FITS file with the necessary FITS info for WCS rotation
    fits.writeto(im_path + 'Dummy_fits_obs.fits', img_data, sci_hdr, overwrite=True)
    fits.writeto(im_path + 'Dummy_fits_dec.fits', dec_data, sci_hdr, overwrite=True)

    # read in reference ALMA image for WCS coordinate rotation
    rot_im = get_pkg_data_filename('ngc5728_alma_full_pyspeckit_fit_single_gaussian_28-08-2018_physical_params.fits')
    rot_list = fits.open(rot_im)

    # read in the dummy FITS file for both images
    dummy_im = get_pkg_data_filename(im_path + 'Dummy_fits_obs.fits')
    dummy_list = fits.open(dummy_im)
    dummy_dec = get_pkg_data_filename(im_path + 'Dummy_fits_dec.fits')
    dummy_dec_list = fits.open(dummy_dec)

    # re-project both images to match ALMA WCS
    obs_img = reproject_interp(dummy_list[0], rot_list[0].header)
    rot_image_obs = obs_img[0]
    dec_img = reproject_interp(dummy_dec_list[0], rot_list[0].header)
    rot_image_dec = dec_img[0]

    # trim each rotated image (observed/deconvolved) from 490x490 array to 256x256, centering image at 128x128
    if i == 0:
        # step 1) determine the centroid for both images
        # step 2) trim 100 pixels off each side from both images based on centroid location

        # shift the deconvolved image?
        dec = shift(rot_image_dec, (-0.5,-0.5), mode='nearest')
        center_obs = centroid_quadratic(rot_image_obs,xpeak=214,ypeak=258,fit_boxsize=20)
        center_dec = centroid_quadratic(rot_image_dec, xpeak=216, ypeak=257, fit_boxsize=20)

        cenx_obs.append(center_obs[0])
        ceny_obs.append(center_obs[1])
        cenx_dec.append(center_obs[0])
        ceny_dec.append(center_obs[1])

        print('F560W obs centroid: ', round(center_obs[0]), round(center_obs[1]))
        print('F560W dec centroid: ', round(center_dec[0]), round(center_dec[1]))

        # save manual inputs as backup -> Richardson-Lucy: obs2
        # obs_trim = rot_image_obs[159:359, 115:315]
        # dec_trim = rot_image_dec[158:358, 116:316]

        # save manual inputs as backup -> Kraken (old): obs2
        obs_trim = rot_image_obs[158:358, 114:314]
        dec_trim = rot_image_dec[157:357, 115:315]
    elif i == 1:
        center_obs = centroid_quadratic(rot_image_obs,xpeak=215,ypeak=259,fit_boxsize=20)
        center_dec = centroid_quadratic(rot_image_dec, xpeak=217, ypeak=258, fit_boxsize=20)

        cenx_obs.append(center_obs[0])
        ceny_obs.append(center_obs[1])
        cenx_dec.append(center_obs[0])
        ceny_dec.append(center_obs[1])

        print('F1000W obs centroid: ', round(center_obs[0]), round(center_obs[1]))
        print('F1000W dec centroid: ', round(center_dec[0]), round(center_dec[1]))

        # # save manual inputs as backup -> Richardson-Lucy: obs2
        # obs_trim = rot_image_obs[159:359, 114:314]
        # dec_trim = rot_image_dec[159:359, 116:316]

        # save manual inputs as backup -> Kraken (old): obs2
        obs_trim = rot_image_obs[158:358, 113:313]
        #obs_trim = rot_image_obs[130:386, 85:341]
        dec_trim = rot_image_dec[160:360, 115:315]

    elif i == 2:
        center_obs = centroid_quadratic(rot_image_obs,xpeak=172,ypeak=274,fit_boxsize=20)
        center_dec = centroid_quadratic(rot_image_dec, xpeak=173, ypeak=274, fit_boxsize=20)

        cenx_obs.append(center_obs[0])
        ceny_obs.append(center_obs[1])
        cenx_dec.append(center_obs[0])
        ceny_dec.append(center_obs[1])

        print('F1500W obs centroid: ', round(center_obs[0]), round(center_obs[1]))
        print('F1500W dec centroid: ', round(center_dec[0]), round(center_dec[1]))

        # # save manual inputs as backup -> Richardson-Lucy: obs2
        # obs_trim = rot_image_obs[174:374, 73:273]
        # dec_trim = rot_image_dec[174:374, 74:274]

        # save manual inputs as backup -> Kraken (old): obs2
        obs_trim = rot_image_obs[173:373, 72:272]
        dec_trim = rot_image_dec[175:375, 74:274]
    elif i == 3:
        center_obs = centroid_quadratic(rot_image_obs,xpeak=173,ypeak=273,fit_boxsize=20)
        center_dec = centroid_quadratic(rot_image_dec, xpeak=173, ypeak=273, fit_boxsize=20)

        cenx_obs.append(center_obs[0])
        ceny_obs.append(center_obs[1])
        cenx_dec.append(center_obs[0])
        ceny_dec.append(center_obs[1])

        print('F1800W obs centroid: ', round(center_obs[0]), round(center_obs[1]))
        print('F1800W dec centroid: ', round(center_dec[0]), round(center_dec[1]))

        # # save manual inputs as backup -> Richardson-Lucy: obs2
        # obs_trim = rot_image_obs[174:374, 73:273]
        # dec_trim = rot_image_dec[174:374, 74:274]

        # save manual inputs as backup -> Kraken (old): obs2
        obs_trim = rot_image_obs[173:373, 72:272]
        dec_trim = rot_image_dec[175:375, 74:274]
    else:
        # arbitrarily set the image max for plotting later
        max21 = dec_data.max()
        center_obs = centroid_quadratic(rot_image_obs,xpeak=173,ypeak=275,fit_boxsize=20)
        center_dec = centroid_quadratic(rot_image_dec, xpeak=173, ypeak=274, fit_boxsize=20)

        cenx_obs.append(center_obs[0])
        ceny_obs.append(center_obs[1])
        cenx_dec.append(center_obs[0])
        ceny_dec.append(center_obs[1])

        print('F2100W obs centroid: ', round(center_obs[0]), round(center_obs[1]))
        print('F2100W dec centroid: ', round(center_dec[0]), round(center_dec[1]))

        # # save manual inputs as backup -> Richardson-Lucy: obs2
        # obs_trim = rot_image_obs[175:375, 73:273]
        # dec_trim = rot_image_dec[175:375, 74:274]

        # save manual inputs as backup -> Kraken (old): obs2
        obs_trim = rot_image_obs[174:374, 72:272]
        dec_trim = rot_image_dec[174:374, 73:273]

    # set the centroids for each image
    # obs_x, obs_y = int(center_obs[0]), int(center_obs[1])
    # dec_x, dec_y = int(center_dec[0]), int(center_dec[1])
    # print(str(filter[i]) + ' obs. centroid: ' + str(obs_x) + ', ' + str(obs_y))
    # print(str(filter[i]) + ' dec. centroid: ' + str(dec_x) + ', ' + str(dec_y))
    #
    # # trim the image array down to a 200x200 array
    # #obs_trim = rot_image_obs[obs_y-100:obs_y+100, obs_x-100:obs_x+100]
    # #dec_trim = rot_image_dec[dec_y-100:dec_y+100, dec_x-100:dec_x+100]
    #
    # # for testing only -> save images to individual FITS files
    fits.writeto(im_path + 'Dummy_after_rot_obs_'+str(filter[i])+'.fits', obs_trim, sci_hdr, overwrite=True)
    fits.writeto(im_path + 'Dummy_after_rot_dec_'+str(filter[i])+'.fits', dec_trim, sci_hdr, overwrite=True)
    # #print('trimmed rotated obs img shape: ', obs_trim.shape)
    # #print('trimmed rotated dec img shape: ', dec_trim.shape)

    # append each image
    #obs_trim = rot_image_obs
    #dec_trim = rot_image_dec
    obs_imgs.append(obs_trim)
    dec_imgs.append(dec_trim)

    # for testing only
    # fig, ax = plt.subplots(nrows=2, ncols=3, figsize=(12,6), tight_layout=True)
    # ax[0][0].imshow(img_data, origin='lower', norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=img_data.max()), cmap='turbo')
    # ax[0][0].set_title('Original image')
    # ax[0][1].imshow(rot_image_obs, origin='lower', norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=img_data.max()),cmap='turbo')
    # ax[0][1].set_title('Reprojected image')
    # ax[0][2].imshow(obs_trim, origin='lower', norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=img_data.max()),cmap='turbo')
    # ax[0][2].set_title('Reprojected and\n cropped image')
    # ax[1][0].imshow(dec_data, origin='lower', norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=dec_data.max()), cmap='turbo')
    # ax[1][0].set_title('Original deconvolved image')
    # ax[1][1].imshow(rot_image_dec, origin='lower', norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=dec_data.max()),cmap='turbo')
    # ax[1][1].set_title('Reprojected deconvolved image')
    # ax[1][2].imshow(dec_trim, origin='lower', norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=dec_data.max()),cmap='turbo')
    # ax[1][2].set_title('Reprojected and cropped\n deconvolved image')
    # plt.show()


if merit_grid:
    # plot the merit functions in a 2x1 grid
    # set the output path for the merit function measurements
    mes_path = 'Images/JWST/4_Deconvolution/5_Kraken/NGC5728/obs2/Measurements/'
    # set the save file name
    save_merit = out_path + '/Merit_Function_NGC5728_obs2'

    # iterate through each filter measurement, extracting the necessary plotting data
    nint_full = []
    flux_full = []
    fwhm_full = []
    for i in range(0, int(len(filter))):
        # read in the merit function data for each waveband
        nint, flux, fwhm = np.loadtxt(mes_path +str(filter[i])+'/RL_Image_Plot_Kraken_NGC5728_'+str(filter[i])+'_obs2_DET_SAMP_XY_21um_aper.txt',unpack=True)

        # dummy index
        # append flux values to dummy index
        dummyA = []

        # normalize the flux counts?
        if normalize:
            A = flux

            flux_norm = []
            for x in range(0, xlim):
                #print('before normalizing: ', A[x])
                normA = A[x] / A[0]
                #print('after normalizing: ', normA)

                flux_norm.append(normA)

            flux = flux_norm

        nintB = []
        flux_B = []
        fwhm_B = []
        for i in range(0, xlim):
            # convert nint from float -> int
            nint1 = nint[i]
            nint1 = int(nint1)

            # read in FWHM data
            fwhmA = fwhm[i]

            # trim float values
            Afwhm = '{:.3f}'.format(fwhmA)

            # read in flux data
            fluxA = flux[i]

            # toggle between absolute flux values or delta flux values
            if delta:
                if i == 0:
                    # append current flux values to dummy index
                    dummyA.append(fluxA)
                    fluxA = 0
                else:
                    # append current flux values to dummy index
                    dummyA.append(fluxA)

                    # calculate delta values
                    ratA = np.abs(fluxA / dummyA[i - 1])
                    #print('delta ratio: ', ratA)

                    fluxA = ratA

            # trim float value for plotting merit functions
            Aflux = '{:.3f}'.format(fluxA)

            # # append the deconvolution iteration, flux value, and FWHM value for this waveband
            nintB.append(nint1)
            flux_B.append(fluxA)
            fwhm_B.append(fwhmA)

        # append the deconvolution iteration, flux value, and FWHM value for all wavebands
        nint_full.append(nintB)
        flux_full.append(flux_B)
        fwhm_full.append(fwhm_B)

    # Create the 1x2 grid
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(11,6), tight_layout=True)

    # set the flux conservation y-lims based on user-input
    if normalize == True:
        ax[1].set_ylabel('Normalized Aperture Flux', fontsize=15)
        ax[1].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        # RL obs2 y-lim
        # ax[1].set_ylim(1, 1.145)
        save_merit = save_merit + '_norm'
    elif delta == True:
        ax[1].set_ylabel('ΔAperture Flux', fontsize=15)
        ax[1].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        save_merit = save_merit + '_delta'
    else:
        ax[1].set_ylabel('Aperture Flux (mJy/sr)', fontsize=15)
        ax[1].yaxis.set_major_formatter(FormatStrFormatter('%.2e'))
        save_merit = save_merit

    # plot the FWHM and flux values
    for i in range(0, int(len(filter))):
        color = ['b', 'g', 'y', 'darkorange', 'r']
        # plot the FWHM values
        ax[0].plot(nint_full[i], fwhm_full[i], color=color[i], label= str(filter[i]) + ' (' + str(merit_iter[i]) + ' iterations)', linewidth=1)
        ax[0].plot(merit_iter[i], fwhm_full[i][merit_iter[i]], color=color[i], marker='X', markersize=15)

        # plot the flux values
        ax[1].plot(nint_full[i], flux_full[i], color=color[i], label= str(filter[i]) + ' (' + str(merit_iter[i]) + ' iterations)', linewidth=1)
        #ax[1].plot(merit_iter[i], flux_full[i][merit_iter[i]], color=color[i], marker='X', markersize=15)

    if diff_lim:
        # set the calculated diffraction limits
        diff_lim = [0.207, 0.328, 0.488, 0.591, 0.674]
        ax[0].hlines(0.197, 0, 3, color='k')
        ax[0].text(3.5, 0.21, 'F560W PSF FWHM', color='b', weight='bold')
        ax[0].hlines(0.317, 3, 6, color='k')
        ax[0].text(6.5, 0.317, 'F1000W PSF FWHM', color='g', weight='bold')
        ax[0].hlines(0.451, 1, 4, color='k')
        ax[0].text(4.5, 0.451, 'F1500W PSF FWHM', color='y', weight='bold')
        ax[0].hlines(0.516, 1, 4, color='k')
        ax[0].text(4.5, 0.516, 'F1800W PSF FWHM', color='darkorange', weight='bold')
        ax[0].hlines(0.571, 1, 4, color='k')
        ax[0].text(4.5, 0.571, 'F2100W PSF FWHM', color='r', weight='bold')

    # set the global figure items
    ax[0].text(2, 0.75, 'Kraken', fontsize=20)
    ax[0].set_ylabel('FWHM (")', fontsize = 15)
    ax[0].set_xlim(0, xlim)
    ax[0].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax[0].set_title('FWHM Limit', fontsize = 20)
    ax[0].legend(loc='upper right')
    ax[1].set_title('Flux Conservation', fontsize=20)
    ax[1].set_xlim(0, xlim)
    fig.supxlabel('Iterations', fontsize=15)
    #plt.suptitle('Model', fontsize=30, y=0.99)
    plt.suptitle('NGC 5728 Merit Functions', fontsize=25, y=0.99)
    plt.savefig(save_merit + '.png')
    plt.show()

if SED:
    # set the save file name
    save_merit = out_path + '/SED_NGC5728_obs2'

    # plot the deconvolved individual FWHM/flux SEDs
    obs_fwhm = [0.287, 0.559, 0.632, 0.746, 0.826]
    obs_flux = [0.025, 0.029, 0.161, 0.204, 0.266]

    # create a line for each value
    obs_fwhm_line = [obs_fwhm[0], obs_fwhm[1], obs_fwhm[2], obs_fwhm[3], obs_fwhm[4]]
    obs_flux_line = [obs_flux[0], obs_flux[1], obs_flux[2], obs_flux[3], obs_flux[4]]

    # plot the deconvolved individual FWHM/flux SEDs
    dec_fwhm = [0.185, 0.265, 0.295, 0.341, 0.369]
    dec_flux = [0.029, 0.032, 0.182, 0.237, 0.315]

    # create a line for each value
    dec_fwhm_line = [dec_fwhm[0], dec_fwhm[1], dec_fwhm[2], dec_fwhm[3], dec_fwhm[4]]
    dec_flux_line = [dec_flux[0], dec_flux[1], dec_flux[2], dec_flux[3], dec_flux[4]]

    # plot the deconvolved individual FWHM/flux SEDs
    rl_fwhm = [0.166, 0.316, 0.323, 0.353, 0.393]
    rl_flux = [0.029, 0.032, 0.179, 0.231, 0.304]

    # create a line for each value
    rl_fwhm_line = [rl_fwhm[0], rl_fwhm[1], rl_fwhm[2], rl_fwhm[3], rl_fwhm[4]]
    rl_flux_line = [rl_flux[0], rl_flux[1], rl_flux[2], rl_flux[3], rl_flux[4]]

    # plot the SED
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(12,6), tight_layout=True)

    # A) plot the FWHM values
    # plot the single line FWHM values
    ax[0].plot(obs_wave, obs_fwhm_line, color='k', linestyle='--', label='JWST/MIRI')
    ax[0].plot(obs_wave, dec_fwhm_line, color='r', label='Kraken')
    #ax[0].plot(obs_wave, rl_fwhm_line, color='g', label='Richardson-Lucy')

    # plot the single line flux values
    ax[1].semilogy(obs_wave, obs_flux_line, color='k', linestyle='--', label='JWST/MIRI')
    ax[1].semilogy(obs_wave, dec_flux_line, color='r', label='Kraken')
    #ax[1].semilogy(obs_wave, rl_flux_line, color='g', label='Richardson-Lucy')

    for i in range(0, int(len(obs_wave))):
        # plot the FWHM values
        # plot the miri fwhm values
        ax[0].plot(obs_wave[i], obs_fwhm[i], color='k', marker='X')
        # plot the kraken fwhm values
        ax[0].plot(obs_wave[i], dec_fwhm[i], color='r', marker='o')
        # plot the RL fwhm values
        #ax[0].plot(obs_wave[i], rl_fwhm[i], color='g', marker='o')

        # plot the flux values
        # plot the miri flux values
        ax[1].plot(obs_wave[i], obs_flux[i], color='k', marker='X')
        # plot the kraken flux values
        ax[1].plot(obs_wave[i], dec_flux[i], color='r', marker='o')
        # plot the RL flux values
        #ax[1].plot(obs_wave[i], rl_flux[i], color='g', marker='o')


    # add the diffraction limits?
    # set the calculated diffraction limits
    diff_lim = [0.207, 0.328, 0.488, 0.591, 0.674]
    ax[0].hlines(diff_lim[0], 4.88, 6.42, color='k', linestyle = '--')
    ax[0].text(3.1, diff_lim[0] + 0.007, 'F560W PSF FWHM', color='b', weight='bold')

    ax[0].hlines(diff_lim[1], 8.76, 11.1, color='k', linestyle='--')
    ax[0].text(7, diff_lim[1] + 0.007, 'F1000W PSF FWHM', color='g', weight='bold')

    ax[0].hlines(diff_lim[2], 13.1, 17.1, color='k', linestyle='--')
    ax[0].text(11.2, diff_lim[2] + 0.007, 'F1500W PSF FWHM', color='y', weight='bold')

    ax[0].hlines(diff_lim[3], 16.0, 20.3, color='k', linestyle='--')
    ax[0].text(14.3, diff_lim[3] +0.007, 'F1800W PSF FWHM', color='darkorange', weight='bold')

    ax[0].hlines(diff_lim[4], 17.8, 24.4, color='k', linestyle='--')
    ax[0].text(17.5, diff_lim[4] +0.007, 'F2100W PSF FWHM', color='r', weight='bold')

    # set the plot specific params
    # set the plot specific params
    ax[0].set_title('FWHM Limit', fontsize=15)
    ax[0].set_xlim(3, 25)
    ax[0].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax[0].xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax[0].set_ylabel('FWHM (")', fontsize=15)
    ax[0].legend(loc='upper left')
    ax[1].set_yscale('log')

    # # set the y-tcik labels manually?
    # labels = ['8.0e+04', '1.0e+05', '5.0e+5', '1.0e+06']
    # positions = [8e4, 1e5, 5e5, 1e6]
    # ax[1].yaxis.set_major_locator(ticker.FixedLocator(positions))
    # ax[1].yaxis.set_major_formatter(ticker.FixedFormatter(labels))
    # ax[1].set_ylim(8e4, 1.1e6)
    # ax[1].set_yticklabels([])

    ax[1].set_title('Flux Conservation', fontsize=15)
    ax[1].set_xlim(3, 25)
    # set the RA and Dec labels
    labels = ['$3x10^{-1}$', '$1x10^{-1}$', '$2x10^{-2}$']
    positions = [0.3, 0.1, 0.025]
    ax[1].yaxis.set_major_locator(ticker.FixedLocator(positions))
    ax[1].yaxis.set_major_formatter(ticker.FixedFormatter(labels))
    ax[1].set_ylim(0.02, 0.35)
    #ax[1].yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
    ax[1].xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax[1].set_ylabel('Aperture Flux (Jy)', fontsize=15)
    ax[1].legend(loc='upper left')

    fig.suptitle('NGC 5728', fontsize=20)
    fig.supxlabel('Observed Wavelength (μm)', fontsize=15)
    plt.savefig(save_merit + '.png')
    plt.show()

if image_grid:
    # set the save file name
    save_merit = out_path + '/Image_grid_NGC5728_obs2'

    # set the RA and Dec labels
    labels = ['-10', '-5', '0', '5', '10']
    positions = [10, 50, 100, 150, 190]

    # read in WCS rotated/cropped images form above
    fig, ax = plt.subplots(nrows=2, ncols=5, figsize=(12, 5), constrained_layout = True)
    for i in range(0, int(len(obs_wave))):
        im1 = ax[0][i].imshow(obs_imgs[i], origin = 'lower',norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=max21), cmap='RdYlBu_r' )
        ax[0][i].set_title(str(filter[i]), fontsize = 15)
        im2 = ax[1][i].imshow(dec_imgs[i], origin='lower',norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=max21), cmap='RdYlBu_r')
        ax[1][i].text(43, 175, str(merit_iter[i]) + ' Iterations', color='w', fontsize=12)

        if i == 0:
            ax[0][i].tick_params(direction='in', color='w', length=6,  left=True, top=False, right=False,bottom=True, labelleft=True, labeltop=False, labelright=False, labelbottom=False)
            ax[1][i].tick_params(direction='in', color='w', length=6)
        else:
            ax[0][i].tick_params(direction='in', color='w', length=6, axis='both', left=True, top=False, right=False,bottom=True, labelleft=False, labeltop=False, labelright=False, labelbottom=False)
            ax[1][i].tick_params(direction='in', color='w', length=6, axis='both', left=True, top=False, right=False,bottom=True, labelleft=False, labeltop=False, labelright=False, labelbottom=True)

        # set the RA/DEC labels-x axis
        ax[0][i].xaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[0][i].xaxis.set_major_formatter(ticker.FixedFormatter(labels))
        ax[1][i].xaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[1][i].xaxis.set_major_formatter(ticker.FixedFormatter(labels))


    # Set orientation compass: x,y
    # create the height of the arrow
    # size = [txt1-x, txt1-y, txt2-x, txt2-y,
    # arrow1-x, arrow1-y, arrow1-xlen, arrow1-ylen,
    # arrow2-x, arrow2-y, arrow2-xlen, arrow2-ylen]
    size = [175, 57, 127, 3,
            183, 10, 0, 35,
            183, 10, -35, 0
            ]
    # add text
    ax[1][4].text(size[0], size[1], 'N', color='w', fontsize=12)
    ax[1][4].text(size[2], size[3], 'E', color='w', fontsize = 12)
    arrow1 = plt.arrow(size[4], size[5], size[6], size[7], width=1,color='w')
    # create the length of the arrow
    arrow2 = plt.arrow(size[8], size[9], size[10], size[11], width=1,color='w')

    # add arrows
    ax[1][4].add_patch(arrow1)
    ax[1][4].add_patch(arrow2)

    # set the plot specific values
    ax[0][0].set_ylabel('JWST/MIRI', fontsize = 12)
    ax[1][0].set_ylabel('Kraken', fontsize=12)

    # set the RA/DEC labels-y axis
    ax[0][0].yaxis.set_major_locator(ticker.FixedLocator(positions))
    ax[0][0].yaxis.set_major_formatter(ticker.FixedFormatter(labels))
    ax[1][0].yaxis.set_major_locator(ticker.FixedLocator(positions))
    ax[1][0].yaxis.set_major_formatter(ticker.FixedFormatter(labels))

    fig.supxlabel('RA-offset (")', fontsize=15)
    fig.supylabel('Dec-offset (")', fontsize=15)
    fig.suptitle('NGC 5728', fontsize=20)

    fig.colorbar(im2, ax=ax.flat, label='Observed Flux (mJy/sr)', format=FormatStrFormatter('%.1e'), pad=0)
    fig.subplots_adjust(hspace=0, wspace=0)
    plt.savefig(save_merit + '.png')
    plt.show()

if residual_grid:
    # set the save file name
    save_merit = out_path + '/Residual_grid_NGC5728_obs2'

    # set the RA and Dec labels
    #labels = ['-10', '-5', '0', '5', '10']
    #positions = [10, 50, 100, 150, 190]

    labels = ['-5', '-2.5', '0', '2.5', '5']
    positions = [2, 25, 50, 75, 98]

    # read in WCS rotated/cropped images form above
    fig, ax = plt.subplots(nrows=3, ncols=5, figsize=(12,8), constrained_layout=True)


    for i in range(0, int(len(obs_wave))):
        # create residual images: (deconv. - observ.)
        residual = dec_imgs[i] - obs_imgs[i]

        # plot each image: obs - dec - res
        im1 = ax[0][i].imshow(obs_imgs[i][51:151, 51:151], origin='lower', interpolation = 'nearest', norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=dec_imgs[i].max()), cmap='turbo')
        ax[0][i].text(28,85,str(filter[i]),color = 'w', fontsize=15, weight='bold')
        im2 = ax[1][i].imshow(dec_imgs[i][51:151, 51:151], origin='lower', interpolation = 'nearest', norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=dec_imgs[i].max()), cmap='turbo')
        ax[1][i].text(25, 85, str(merit_iter[i]) + ' Iterations', color = 'w', fontsize = 12)
        im3 = ax[2][i].imshow(residual[51:151, 51:151], origin='lower', interpolation = 'nearest',norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=dec_imgs[i].max()), cmap='turbo')

        ax[0][0].yaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[0][0].yaxis.set_major_formatter(ticker.FixedFormatter(labels))
        ax[1][0].yaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[1][0].yaxis.set_major_formatter(ticker.FixedFormatter(labels))
        ax[2][0].yaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[2][0].yaxis.set_major_formatter(ticker.FixedFormatter(labels))
        ax[2][i].xaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[2][i].xaxis.set_major_formatter(ticker.FixedFormatter(labels))

        ax[0][0].tick_params(direction='in', color='w', length=8, left=True, top=False, right=False, bottom=True,labelleft=True, labeltop=False, labelright=False, labelbottom=False)
        ax[1][0].tick_params(direction='in', color='w', length=8, left=True, top=False, right=False, bottom=True,labelleft=True, labeltop=False, labelright=False, labelbottom=False)
        ax[2][0].tick_params(direction='in', color='w', length=8, left=True, top=False, right=False, bottom=True,labelleft=True, labeltop=False, labelright=False, labelbottom=True)
        ax[0][i].tick_params(direction='in', color='w', length=8, left=True, top=False, right=False, bottom=True,labelleft=False, labeltop=False, labelright=False, labelbottom=False)
        ax[1][i].tick_params(direction='in', color='w', length=8, left=True, top=False, right=False, bottom=True,labelleft=False, labeltop=False, labelright=False, labelbottom=False)
        ax[2][i].tick_params(direction='in', color='w', length=8, left=True, top=False, right=False, bottom=True,labelleft=False, labeltop=False, labelright=False, labelbottom=True)


        #fig.colorbar(im3, label='Observed Flux (mJy/sr)', format=FormatStrFormatter('%.1e'), pad=0.01)

    # Set orientation compass: x,y
    # create the height of the arrow
    # size = [txt1-x, txt1-y, txt2-x, txt2-y,
    # arrow1-x, arrow1-y, arrow1-xlen, arrow1-ylen,
    # arrow2-x, arrow2-y, arrow2-xlen, arrow2-ylen]

    # for a 200x200 image grid
    # size = [175, 57, 127, 3,
    #         183, 10, 0, 35,
    #         183, 10, -35, 0
    #         ]

    # for a 100x100 image grid
    size = [91, 31, 68, 7,
            95, 10, 0, 15,
            95, 10, -15, 0
            ]

    # add text
    ax[2][4].text(size[0], size[1], 'N', color='w', fontsize=12)
    ax[2][4].text(size[2], size[3], 'E', color='w', fontsize = 12)
    arrow1 = plt.arrow(size[4], size[5], size[6], size[7], width=1,color='w')
    # create the length of the arrow
    arrow2 = plt.arrow(size[8], size[9], size[10], size[11], width=1,color='w')

    # add arrows
    ax[2][4].add_patch(arrow1)
    ax[2][4].add_patch(arrow2)

    # set the plot specific values
    ax[0][0].set_ylabel('JWST/MIRI', fontsize=15)
    ax[1][0].set_ylabel('Deconvolved', fontsize=15)
    ax[2][0].set_ylabel('Residual', fontsize=15)
    fig.supxlabel('RA-offset (")', fontsize=15)
    fig.supylabel('Dec-offset (")', fontsize=15)
    fig.suptitle('NGC 5728', fontsize=20)
    fig.subplots_adjust(hspace=0, wspace=0)
    plt.savefig(save_merit + '.png')
    plt.show()

if RGB:
    # plot 3-color RGB image grid (2x3) comparing the observed vs deconvolved images
    # column 1: F560W, F1000W, F1500W
    # column 2: F1000W, F1500W, F1800W
    # column 3: F1500W, F1800W, F2100W
    # set the save file name
    save_merit = out_path + '/RGB_grid_NGC5728_obs2'

    # set the RA and Dec labels
    labels = ['-10', '-5', '0', '5', '10']
    positions = [10, 50, 100, 150, 190]

    # read in WCS rotated/cropped images from above
    fig, ax = plt.subplots(nrows=2, ncols=3, figsize=(9,6), constrained_layout=True)

    for i in range(0, 3):

        # set the observation filter
        blue = filter[i]
        green = filter[int(i+1)]
        red = filter[int(i+2)]

        # generate the observ. RGB image
        obs_img_rgb = np.zeros((obs_imgs[i].shape[0], obs_imgs[i].shape[1], 3), dtype=float)
        obs_img_rgb[:, :, 0] = img_scale.log(obs_imgs[i+2], scale_min=0, scale_max=obs_imgs[i+2].max())
        obs_img_rgb[:, :, 1] = img_scale.log(obs_imgs[i+1], scale_min=0, scale_max=obs_imgs[i+1].max())
        obs_img_rgb[:, :, 2] = img_scale.log(obs_imgs[i], scale_min=0, scale_max=obs_imgs[i].max())

        # generate the deconv. RGB image
        dec_img_rgb = np.zeros((dec_imgs[i].shape[0], dec_imgs[i].shape[1], 3), dtype=float)
        dec_img_rgb[:, :, 0] = img_scale.log(dec_imgs[i+2], scale_min=0, scale_max=dec_imgs[i+2].max())
        dec_img_rgb[:, :, 1] = img_scale.log(dec_imgs[i+1], scale_min=0, scale_max=dec_imgs[i+1].max())
        dec_img_rgb[:, :, 2] = img_scale.log(dec_imgs[i], scale_min=0, scale_max=dec_imgs[i].max())

        # plot the RGB images
        ax[0][i].imshow(obs_img_rgb, origin='lower')
        ax[1][i].imshow(dec_img_rgb, origin='lower')

        mono = {'family': 'monospace'}
        # observation tags
        ax[0][i].text(5, 10, str(blue), color='royalblue', fontsize=10)
        ax[0][i].text(5, 26, str(green), color='chartreuse', fontsize=10)
        ax[0][i].text(5, 41, str(red), color='r', fontsize=10)

        # deconvolved tags
        ax[1][i].text(5, 10, str(blue), color='royalblue', fontsize=10)
        ax[1][i].text(5, 26, str(green), color='chartreuse', fontsize=10)
        ax[1][i].text(5, 41, str(red), color='r', fontsize=10)


    # set the RA/DEC labels-y axis
    ax[0][0].yaxis.set_major_locator(ticker.FixedLocator(positions))
    ax[0][0].yaxis.set_major_formatter(ticker.FixedFormatter(labels))
    ax[1][0].yaxis.set_major_locator(ticker.FixedLocator(positions))
    ax[1][0].yaxis.set_major_formatter(ticker.FixedFormatter(labels))

    #  set the RA/DEC labels-x axis
    ax[1][0].xaxis.set_major_locator(ticker.FixedLocator(positions))
    ax[1][0].xaxis.set_major_formatter(ticker.FixedFormatter(labels))
    ax[1][1].xaxis.set_major_locator(ticker.FixedLocator(positions))
    ax[1][1].xaxis.set_major_formatter(ticker.FixedFormatter(labels))
    ax[1][2].xaxis.set_major_locator(ticker.FixedLocator(positions))
    ax[1][2].xaxis.set_major_formatter(ticker.FixedFormatter(labels))

    # set the image params
    ax[0][0].tick_params(direction='in', color='w', length=6, left=True, top=False, right=False, bottom=True, labelleft=True, labeltop=False, labelright=False, labelbottom=False)
    ax[0][1].tick_params(direction='in', color='w', length=6, left=True, top=False, right=False, bottom=True,labelleft=False, labeltop=False, labelright=False, labelbottom=False)
    ax[0][2].tick_params(direction='in', color='w', length=6, left=True, top=False, right=False, bottom=True,labelleft=False, labeltop=False, labelright=False, labelbottom=False)
    ax[1][0].tick_params(direction='in', color='w', length=6, left=True, top=False, right=False, bottom=True,labelleft=True, labeltop=False, labelright=False, labelbottom=True)
    ax[1][1].tick_params(direction='in', color='w', length=6, left=True, top=False, right=False, bottom=True,labelleft=False, labeltop=False, labelright=False, labelbottom=True)
    ax[1][2].tick_params(direction='in', color='w', length=6, left=True, top=False, right=False, bottom=True,labelleft=False, labeltop=False, labelright=False, labelbottom=True)

    # Set orientation compass: x,y
    # create the height of the arrow
    # size = [txt1-x, txt1-y, txt2-x, txt2-y,
    # arrow1-x, arrow1-y, arrow1-xlen, arrow1-ylen,
    # arrow2-x, arrow2-y, arrow2-xlen, arrow2-ylen]
    size = [177, 50, 137, 5,
            183, 10, 0, 30,
            183, 10, -30, 0
            ]
    # add text
    ax[1][2].text(size[0], size[1], 'N', color='w', fontsize=12)
    ax[1][2].text(size[2], size[3], 'E', color='w', fontsize = 12)
    arrow1 = plt.arrow(size[4], size[5], size[6], size[7], width=1,color='w')
    # create the length of the arrow
    arrow2 = plt.arrow(size[8], size[9], size[10], size[11], width=1,color='w')

    # add arrows
    ax[1][2].add_patch(arrow1)
    ax[1][2].add_patch(arrow2)

    fig.supxlabel('RA-offset (")', fontsize=15)
    fig.supylabel('Dec-offset (")', fontsize=15)
    fig.suptitle('NGC 5728', fontsize=20)
    ax[0][0].set_ylabel('JWST/MIRI', fontsize=12)
    ax[1][0].set_ylabel('Deconvolved', fontsize=12)
    fig.subplots_adjust(hspace=0.1, wspace=0.1)
    plt.savefig(save_merit + '.png')
    plt.show()





