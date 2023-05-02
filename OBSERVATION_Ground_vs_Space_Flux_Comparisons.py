# 2022 Oct 19: Quick code to compare the measured flux in a similiar aperture size between observations of NGC 5728
# made by JWST and those made by VLT and Gemini South

from Convenience_Functions import *

# set which images you want displayed
jwst = True
RL_decon = False
kraken = True
ground_img = True

# set the deconvolved index values
RL_idx = [8, 14, 27, 23, 22]
K_idx = [14, 12, 16, 17, 19]

aper_size = 3.783

# Set the JWST images path
space_path = 'Images/JWST/5_ERS_Data/1_Raw_data/ngc5728_MIRI_Imaging/'
# set the ground-based images import path
ground_path = 'NGC5728_Deconvolution_Rosario/'
# set the Richardson-Lucy deconvolved image path
RL_path = 'Images/JWST/4_Deconvolution/1_Richardson_Lucy/NGC5728/FITS/'
# set the Richardson-Lucy deconvolved image path
kraken_path = 'Images/JWST/4_Deconvolution/5_Kraken/NGC5728/FITS/'

# set the JWST/MIRIimage filter name
filter = ['F560W', 'F1000W', 'F1500W', 'F1800W', 'F2100W']

# set the ground based image names
ground = ['8_74um_NGC5728_Si2_2005-07-09T00-00',
          '11-25um_NGC5728_PAH2_2010-03-12T07-11',
          '12_27um_NGC5728_NEII_1_2010-03-12T07-19',
          '13-04um_NGC5728_NEII_2_2008-03-21T05-36',
          '18_3um_NGC5728_Qa_2005-07-09T00-05'
          ]

# set the title for the ground observations
ground_title = ['Gemini South/T-ReCS 8.7 μm',
                'VLT/VISIR 11.3 μm',
                'VLT/VISIR 12.3 μm',
                'VLT/VISIR 13.0 μm',
                'Gemini South/T-ReCS 18.3 μm',
                ]

# set the titles for the JWST/MIRI observations
space_title = ['JWST/MIRI 5.6 μm', 'JWST/MIRI 10.0 μm', 'JWST/MIRI 15.0 μm', 'JWST/MIRI 18.0 μm', 'JWST/MIRI 21.0 μm']

# annulus data
space_anul = []
RL_anul = []
K_anul = []
ground_anul = []

# aperture data
space_aper = []
RL_aper = []
K_aper = []
ground_aper = []

# ground-based error bars
aperLOW = []
aperHIGH = []
anulLOW = []
anulHIGH = []
count = 0
# 1) read in the image data
for i in range(0, int(len(filter))):
    if i == 0:
        # F560W
        Space_xpeak = 151
        Space_ypeak = 138
        wave = 5.6

        # deconvolved image
        dec_x = 66
        dec_y = 66

        # 8.74 um: 8_74um_NGC5728_Si2_2005-07-09T00-00
        Ground_xpeak = 28
        Ground_ypeak = 26
        wave1 = 8.7
        px_scl = 1

        box = 10
        threshold = 0.5

    elif i == 1:
        # F1000W
        Space_xpeak = 151
        Space_ypeak = 139
        wave = 10.0

        # deconvolved image
        dec_x = 64.9
        dec_y = 64.9

        # 11.25 um: 11-25um_NGC5728_PAH2_2010-03-12T07-11
        Ground_xpeak = 27
        Ground_ypeak = 27
        wave1 = 11.3
        px_scl = 1

        box = 10
        threshold = 0.5

    elif i == 2:
        # F1500W
        Space_xpeak = 174
        Space_ypeak = 116
        wave = 15.0

        # deconvolved image
        dec_x = 66
        dec_y = 64.9

        # 12.27 um: 12_27um_NGC5728_NEII_1_2010-03-12T07-19
        Ground_xpeak = 26
        Ground_ypeak = 26
        wave1 = 12.3
        px_scl = 1

        box = 10
        threshold = 0.5

    elif i == 3:
        # F1800W
        Space_xpeak = 197
        Space_ypeak = 86
        wave = 18.0

        # deconvolved image
        dec_x = 66
        dec_y = 66

        # 13.04 um: 13-04um_NGC5728_NEII_2_2008-03-21T05-36
        Ground_xpeak = 27
        Ground_ypeak = 27
        wave1 = 13.0
        px_scl = 1

        box = 10
        threshold = 0.5

    elif i == 4:
        # F2100W
        Space_xpeak = 221
        Space_ypeak = 55
        wave = 21.0

        # deconvolved image
        dec_x = 80.9
        dec_y = 54.9

        # 18.3 um:
        Ground_xpeak = 28
        Ground_ypeak = 26
        wave1 = 18.3
        px_scl = 1

        box = 10
        threshold = 0.5

    # read in the JWST/MIRI images
    obs = get_pkg_data_filename(space_path + 'NGC5728_' + filter[i] + '_final.fits')
    obs_list = fits.open(obs)
    obs_im2 = obs_list[1].data
    obs_im2 = obs_im2[Space_ypeak - 25:Space_ypeak+25, Space_xpeak-25:Space_xpeak+25]

    # convert from mJy/sr -> mJy
    factor = ((0.1110 ** 2) * 2.350443e-5) * 1000
    obs_im2 = obs_im2 * factor

    # mask negative pixel values
    obs_im2[obs_im2<0] += 0

    # save the image to a place-holder fits file
    save_FITS(obs_im2, filename='residual',
              name='RL_Merit_function_holder',
              output_dir=ground_path, pixel_scale=0.1110,
              filter=filter[i], wavelength=wave, instrument='MIRI')

    # Read in the formatted JWST/MIRI images FITS file
    obs_ims = get_pkg_data_filename(ground_path + 'residual_RL_Merit_function_holder.fits')
    obs_list2 = fits.open(obs_ims)
    obs_im = obs_list2[0].data

    # read in the RL deconvolved images
    RL = get_pkg_data_filename(RL_path + 'NGC5728_' + filter[i] + '_final_DET_SAMP.fits')
    RL_list = fits.open(RL)
    RL_im2 = RL_list[RL_idx[i]].data
    RL_im2 = RL_im2[int(round(dec_y)) - 25:int(round(dec_y))+25, int(round(dec_x))-25:int(round(dec_x))+25]

    # convert from mJy/sr -> mJy
    RL_im2 = RL_im2 * factor

    # mask negative pixel values
    RL_im2[RL_im2<0] += 0

    # save the image to a place-holder fits file
    save_FITS(RL_im2, filename='residual',
              name='RL_Merit_function_holder1',
              output_dir=ground_path, pixel_scale=0.1110,
              filter=filter[i], wavelength=wave, instrument='MIRI')

    # Read in the formatted RL images FITS file
    RL_ims = get_pkg_data_filename(ground_path + 'residual_RL_Merit_function_holder1.fits')
    RL_list2 = fits.open(RL_ims)
    RL_im = RL_list2[0].data

    # read in the kraken deconvolved images
    K = get_pkg_data_filename(kraken_path + 'NGC5728_' + filter[i] + '_final_DET_SAMP.fits')
    K_list = fits.open(K)
    K_im2 = RL_list[K_idx[i]].data
    K_im2 = K_im2[int(round(dec_y)) - 25:int(round(dec_y))+25, int(round(dec_x))-25:int(round(dec_x))+25]

    # convert from mJy/sr -> mJy
    K_im2 = K_im2 * factor

    # mask negative pixel values
    K_im2[K_im2<0] += 0

    # save the image to a place-holder fits file
    save_FITS(K_im2, filename='residual',
              name='RL_Merit_function_holder2',
              output_dir=ground_path, pixel_scale=0.1110,
              filter=filter[i], wavelength=wave, instrument='MIRI')

    # Read in the formatted RL images FITS file
    K_ims = get_pkg_data_filename(ground_path + 'residual_RL_Merit_function_holder2.fits')
    K_list2 = fits.open(K_ims)
    K_im = K_list2[0].data


    # Read in the ground-based images
    gro = get_pkg_data_filename(ground_path + str(ground[i]) + '.fits')
    gro_list = fits.open(gro)
    gro_im = gro_list[0].data

    # mask negative pixel values
    gro_im[gro_im < 0] += 0

    # save the image to a place-holder fits file
    save_FITS(gro_im, filename='residual',
              name='RL_Merit_function_holder3',
              output_dir=ground_path, pixel_scale=px_scl,
              filter=filter[i], wavelength=wave1, instrument='MIRI')

    # Read in the deconvolved FITS file
    gro_ims = get_pkg_data_filename(ground_path + 'residual_RL_Merit_function_holder3.fits')
    gro_list2 = fits.open(gro_ims)
    ground_im = gro_list2[0].data

    # 3) set the object centroids
    # # measure the centroid: JWST/MIRI
    center_quad = centroid_quadratic(obs_im, xpeak=25, ypeak=25, fit_boxsize=box)
    center_space = (center_quad[0], center_quad[1])

    # # measure the centroid: RL
    center_quad = centroid_quadratic(RL_im, xpeak=25, ypeak=25, fit_boxsize=box)
    center_RL = (center_quad[0], center_quad[1])

    # # measure the centroid: kraken
    center_quad = centroid_quadratic(K_im, xpeak=25, ypeak=25, fit_boxsize=box)
    center_K = (center_quad[0], center_quad[1])

    # # measure the centroid: ground
    center_quad = centroid_quadratic(ground_im, xpeak=Ground_xpeak, ypeak=Ground_ypeak, fit_boxsize=box)
    center_ground = (center_quad[0], center_quad[1])

    # measure the flux of the point source based on the aperture size
    position_space = [(center_space[0] - 0.2), (center_space[1] - 0.2)]
    position_RL = [(center_RL[0] - 0.2), (center_RL[1] - 0.2)]
    position_K = [(center_K[0] - 0.2), (center_K[1] - 0.2)]
    if i == 0:
        position_ground = [(center_ground[0]+1), (center_ground[1]+1)]
    elif i == 2:
        position_ground = [(center_ground[0] - 1), (center_ground[1] - 1)]
    else:
        position_ground = [(center_ground[0]), (center_ground[1])]


    # measure the JWST/MIRI images-aperture
    aperture_jwst = photutils_CircularAperture(position_space, r=aper_size)
    phot_table = aperture_photometry(obs_im, aperture_jwst, method='exact')
    phot_table['aperture_sum'].info.format = '%.8g'
    aperture_sum = phot_table['aperture_sum'][0]
    aper_jwst = float(aperture_sum)

    # measure the RL images-aperture
    aperture_RL = photutils_CircularAperture(position_RL, r=aper_size)
    phot_table2 = aperture_photometry(RL_im, aperture_RL, method='exact')
    phot_table2['aperture_sum'].info.format = '%.8g'
    aperture_sum2 = phot_table2['aperture_sum'][0]
    aper_RL = float(aperture_sum2)

    # measure the Kraken images-aperture
    aperture_K = photutils_CircularAperture(position_K, r=aper_size)
    phot_table3 = aperture_photometry(K_im, aperture_K, method='exact')
    phot_table3['aperture_sum'].info.format = '%.8g'
    aperture_sum3 = phot_table3['aperture_sum'][0]
    aper_K = float(aperture_sum3)

    # measure the ground images-aperture
    aperture_gro = photutils_CircularAperture(position_ground, r=aper_size)
    phot_table4 = aperture_photometry(gro_im, aperture_gro, method='exact')
    phot_table4['aperture_sum'].info.format = '%.8g'
    aperture_sum4 = phot_table4['aperture_sum'][0]
    aper_gro = float(aperture_sum4)

    # measure the extended emission for each image
    # Perform annulus photometry and remove the background from the total flux measured above
    # define the inner and outer annulus radius
    outer_rad = 12.153

    # measure the JWST/MIRI images-anullus
    anullus_jwst = photutils_CircularAperture(position_space, r=outer_rad)
    phot_table5 = aperture_photometry(obs_im, anullus_jwst, method='exact')
    phot_table5['aperture_sum'].info.format = '%.8g'
    aperture_sum5 = phot_table5['aperture_sum'][0]
    anul_jwst = float(aperture_sum5)

    # measure the RL images-anullus
    anullus_RL = photutils_CircularAperture(position_RL, r=outer_rad)
    phot_table6 = aperture_photometry(RL_im, anullus_RL, method='exact')
    phot_table6['aperture_sum'].info.format = '%.8g'
    aperture_sum6 = phot_table6['aperture_sum'][0]
    anul_RL = float(aperture_sum6)

    # measure the Kraken images-anullus
    anullus_K = photutils_CircularAperture(position_K, r=outer_rad)
    phot_table7 = aperture_photometry(K_im, anullus_K, method='exact')
    phot_table7['aperture_sum'].info.format = '%.8g'
    aperture_sum7 = phot_table7['aperture_sum'][0]
    anul_K = float(aperture_sum7)

    # measure the ground images-anullus
    anullus_gro = photutils_CircularAperture(position_ground, r=outer_rad)
    phot_table8 = aperture_photometry(obs_im, anullus_gro, method='exact')
    phot_table8['aperture_sum'].info.format = '%.8g'
    aperture_sum8 = phot_table8['aperture_sum'][0]
    anul_gro = float(aperture_sum8)

    # subtract
    subtract_jwst = anul_jwst - aper_jwst
    subtract_RL = anul_RL - aper_RL
    subtract_K = anul_K - aper_K
    subtract_gro = anul_gro - aper_gro

    if i ==4:
        # the longest waveband gets a 20% error
        # calculate ground_img errors- aperture
        aper = aper_gro * 0.2
        ap_low = aper_gro - aper
        ap_high = aper_gro + aper
        # calculate ground_img errors- annulus
        anul = subtract_gro * 0.2
        an_low = subtract_gro - anul
        an_high = subtract_gro + anul
    else:
        # all else get a 10% error band
        aper = aper_gro * 0.1
        ap_low = aper_gro - aper
        ap_high = aper_gro + aper
        # calculate ground_img errors- annulus
        anul = subtract_gro * 0.1
        an_low = subtract_gro - anul
        an_high = subtract_gro + anul


    # 5) append photometry measurements (for comparison plot later)
    space_anul.append(subtract_jwst)
    space_aper.append(aper_jwst)
    RL_anul.append(subtract_RL)
    RL_aper.append(aper_RL)
    K_anul.append(subtract_K)
    K_aper.append(aper_K)
    ground_anul.append(subtract_gro)
    ground_aper.append(aper_gro)

    # append the error bars
    aperLOW.append(ap_low)
    aperHIGH.append(ap_high)
    anulLOW.append(an_low)
    anulHIGH.append(an_high)
    count += 2

    print('| Space-based Flux (mJy) (aperture | annulus): ', aper_jwst, ' | ', subtract_jwst,
          #'| Kraken (mJy) (aperture | annulus): ', aper_K, ' | ', subtract_K,
          '| RL (mJy) (aperture | annulus): ', aper_RL, subtract_RL,
          '| Ground-based Flux (mJy) (aperture | annulus): ', aper_gro, ' | ', subtract_gro
          )

    # display and save the JWST images?
    if jwst:
        # display the apertures to make sure the're fitting correctly
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 7), tight_layout=True)

        # trim decimals
        total = '{:.2f}'.format(aper_jwst)
        sub = '{:.2f}'.format(subtract_jwst)

        # set the image positions
        x_positions = [2.5, 16, 25, 34, 47.5]
        y_positions = [2.5, 16, 25, 34, 47.5]
        labels = ['-2.5', '-1', '0', '1', '2.5']

        ap_patches = aperture_jwst.plot(color='r', lw=2)
        ap_patches2 = anullus_jwst.plot(color='w', lw=2, linestyle='--')
        im00 = plt.imshow(obs_im, origin='lower', cmap='viridis',
                          norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=obs_im.max()))
        plt.title(space_title[i], fontsize=20)
        ax.xaxis.set_major_locator(ticker.FixedLocator(x_positions))
        ax.xaxis.set_major_formatter(ticker.FixedFormatter(labels))
        ax.yaxis.set_major_locator(ticker.FixedLocator(y_positions))
        ax.yaxis.set_major_formatter(ticker.FixedFormatter(labels))
        ax.set_ylabel('Y (")', fontsize=15)
        ax.set_xlabel('X (")', fontsize=15)
        add_colorbar(im00, label='Log-scaled Flux Density (mJy)', format=FormatStrFormatter('%.1e'))
        #plt.savefig(ground_path + 'Space_' + filter[i] + '_OBS_21um.png')

        plt.grid(False)
        #plt.show()

    # display/save the RL decon images images?
    if RL_decon:
        # display the apertures to make sure the're fitting correctly
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 7), tight_layout=True)

        # trim decimals
        total1 = '{:.2f}'.format(aper_RL)
        sub1 = '{:.2f}'.format(subtract_RL)

        # set the image positions
        x_positions = [2.5, 16, 25, 34, 47.5]
        y_positions = [2.5, 16, 25, 34, 47.5]
        labels = ['-2.5', '-1', '0', '1', '2.5']

        ap_patches3 = aperture_RL.plot(color='r', lw=2)
        ap_patches4 = anullus_RL.plot(color='w', lw=2, linestyle='--')
        im00 = plt.imshow(RL_im, origin='lower', cmap='viridis',
                          norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=RL_im.max()))
        plt.title(str(filter[i]) + ' Richardson-Lucy', fontsize=20)
        ax.xaxis.set_major_locator(ticker.FixedLocator(x_positions))
        ax.xaxis.set_major_formatter(ticker.FixedFormatter(labels))
        ax.yaxis.set_major_locator(ticker.FixedLocator(y_positions))
        ax.yaxis.set_major_formatter(ticker.FixedFormatter(labels))
        ax.set_ylabel('Y (")', fontsize=15)
        ax.set_xlabel('X (")', fontsize=15)
        add_colorbar(im00, label='Log-scaled Flux Density (mJy)', format=FormatStrFormatter('%.1e'))
        #plt.savefig(ground_path + 'RL_' + str(filter[i]) + '_OBS_21um.png')

        plt.grid(False)
        #plt.show()

    # display/save the kraken decon images images?
    if kraken:
        # display the apertures to make sure the're fitting correctly
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 7), tight_layout=True)

        # trim decimals
        total1 = '{:.2f}'.format(aper_K)
        sub1 = '{:.2f}'.format(subtract_K)

        # set the image positions
        x_positions = [2.5, 16, 25, 34, 47.5]
        y_positions = [2.5, 16, 25, 34, 47.5]
        labels = ['-2.5', '-1', '0', '1', '2.5']

        ap_patches3 = aperture_K.plot(color='r', lw=2)
        ap_patches4 = anullus_K.plot(color='w', lw=2, linestyle='--')
        im00 = plt.imshow(K_im, origin='lower', cmap='viridis',
                          norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=K_im.max()))
        plt.title(str(filter[i]) + ' Kraken', fontsize=20)
        ax.xaxis.set_major_locator(ticker.FixedLocator(x_positions))
        ax.xaxis.set_major_formatter(ticker.FixedFormatter(labels))
        ax.yaxis.set_major_locator(ticker.FixedLocator(y_positions))
        ax.yaxis.set_major_formatter(ticker.FixedFormatter(labels))
        ax.set_ylabel('Y (")', fontsize=15)
        ax.set_xlabel('X (")', fontsize=15)
        add_colorbar(im00, label='Log-scaled Flux Density (mJy)', format=FormatStrFormatter('%.1e'))
        #plt.savefig(ground_path + 'Kraken_' + str(filter[i]) + '_OBS_21um.png')

        plt.grid(False)
        #plt.show()

    # display/save the ground based images?
    if ground_img:
        # display the apertures to make sure the're fitting correctly
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 7), tight_layout=True)

        # trim decimals
        total1 = '{:.2f}'.format(aper_gro)
        sub1 = '{:.2f}'.format(subtract_gro)

        # set the image positions
        x_positions = [2.5, 16, 25, 34, 47.5]
        y_positions = [2.5, 16, 25, 34, 47.5]
        labels = ['-2.5', '-1', '0', '1', '2.5']

        ap_patches3 = aperture_gro.plot(color='r', lw=2)
        ap_patches4 = anullus_gro.plot(color='w', lw=2, linestyle='--')
        im00 = plt.imshow(ground_im, origin='lower', cmap='viridis',
                          norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=ground_im.max()))
        plt.title(ground_title[i], fontsize=20)
        ax.xaxis.set_major_locator(ticker.FixedLocator(x_positions))
        ax.xaxis.set_major_formatter(ticker.FixedFormatter(labels))
        ax.yaxis.set_major_locator(ticker.FixedLocator(y_positions))
        ax.yaxis.set_major_formatter(ticker.FixedFormatter(labels))
        ax.set_ylabel('Y (")', fontsize=15)
        ax.set_xlabel('X (")', fontsize=15)
        add_colorbar(im00, label='Log-scaled Flux Density (mJy)', format=FormatStrFormatter('%.1e'))
        #plt.savefig(ground_path + 'Ground_' + ground[i] + '_OBS_21um.png')

        plt.grid(False)
        #plt.show()

# plot flux comparison plot (mJy vs wavelength)
space_wave = [5.6, 10, 15, 18, 21]
ground_wave = [8.74, 11.25, 12.27, 13.04, 18.3]

# uncertainties
# aper_high = [aperLOW[0], aperLOW[1], aperLOW[0], aperLOW[0], aperLOW[0], ]
# aper_low = [12.78, 18.74, 47.59, 83.79, 98.83]
# error_high = [7.59, 2.19, 102.45, 98.62, 114.36]
# error_low = [6.21, 1.79, 83.83, 80.53, 76.24]

# convert aperture and annulus size from pixel -> arcsec and trim to 2 decimals
ap = ((aper_size)*2)*0.1110
anul = 2.7-((aper_size*2)*0.1110)

# trim decimals
apert = '{:.2f}'.format(ap)
anul = '{:.2f}'.format(anul)

fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(16, 7), tight_layout=True)
ax[0].set_title('JWST/MIRI', fontsize=15)
#ax[0].plot(5,5, color = 'w', label = 'aperture size = '+str(apert)+' (")')
ax[0].plot(space_wave, space_aper, color = 'k', label = 'Unresolved core')
ax[0].plot(space_wave, space_anul, color = 'r', label = 'Resolved extended')
ax[0].legend(loc = 'upper left')
ax[0].set_ylim(0, 210)
ax[2].set_title('Ground-Based', fontsize=15)
ax[2].plot(ground_wave, ground_aper, color = 'k', label = 'Aperture')
ax[2].plot(ground_wave, ground_anul, color = 'r', label = 'Annulus')
ax[2].set_xlim(6, 21)

if kraken:
    ax[1].plot(space_wave, K_aper, color='k')
    ax[1].plot(space_wave, K_anul, color='r')
    ax[1].set_title('Kraken Deconvolved')
    #ax[1].legend(loc='upper left')
    ax[1].set_ylim(0, 210)
    # plot the error bars individually
    for i in range(len(space_wave)):
        ax[1].plot(space_wave[i], K_anul[i], color='r', marker='o')
        ax[1].plot(space_wave[i], K_aper[i], color='k', marker='o')

# plot the error bars individually
for i in range(len(space_wave)):
    ax[0].plot(space_wave[i], space_anul[i], color='r', marker = 'o')
    ax[0].plot(space_wave[i], space_aper[i], color='k', marker='o')
    ax[2].plot(ground_wave[i], ground_anul[i], color='r', marker = 'X', markersize = 5)
    ax[2].plot(ground_wave[i], ground_aper[i], color='k', marker='X', markersize=5)

    # aperture error bars
    ax[2].vlines(x=ground_wave[i], ymin=aperLOW[i], ymax=aperHIGH[i], color = 'b')
    ax[2].hlines(y=aperLOW[i], xmin=ground_wave[i]-0.1, xmax=ground_wave[i]+0.1, colors='b')
    ax[2].hlines(y=aperHIGH[i], xmin=ground_wave[i] - 0.1, xmax=ground_wave[i] + 0.1, colors='b')

    # annulus error bars
    ax[2].vlines(x=ground_wave[i], ymin=anulLOW[i], ymax=anulHIGH[i], color = 'b')
    ax[2].hlines(y=anulLOW[i], xmin=ground_wave[i]-0.1, xmax=ground_wave[i]+0.1, colors='b')
    ax[2].hlines(y=anulHIGH[i], xmin=ground_wave[i] - 0.1, xmax=ground_wave[i] + 0.1, colors='b')


fig.suptitle('Photometric Comparisons', fontsize=20)
fig.supylabel('Flux Density (mJy)')
fig.supxlabel('Observed Wavelength (μm)')
if RL_decon == True:
    plt.savefig(ground_path + 'RL_Ground_Space_6Dec2022_uncertanties.png')
elif kraken == True:
    plt.savefig(ground_path + 'Kraken_Ground_Space_21um_13dec2022_uncertanties.png')
else:
    plt.savefig(ground_path + 'Ground_Space_8Dec2022_uncertanties_RL_21um_swapped.png')
plt.show()




