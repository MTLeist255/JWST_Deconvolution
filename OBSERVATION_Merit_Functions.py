# 2022 Nov 19: code to measure the merit fuctions of the JWST/MIRI observations

from Convenience_Functions import *

# Plot the merit function results?
plot = True
# Perform aperture photometry on other regions in the image array?
photom_test = True
# if photom_test: fit an aperture across the background galaxy
flux_galaxy = True
# if photom_test: fit an aperture based on the F2100W PSF Airy disk
aper_21 = True

# set the deconvolution path
Dir = '5_Kraken'
# set image specific global params:
# image filter
filter = 'F2100W'
# set the observation title
obs = 'Kraken_update_NGC5728_'+str(filter)+'_obs2_DET_SAMP'

# set global measurement parameters
# centroid guess
guess = (127,127)
# centroid box size
box = 10
# photometry aperture radius
aper_rad = 17.8
# image pad level
pad = 0.5
# FWHM threshold
threshold = 0
# galaxy aperture
galax_aper = 60

# FOR ADDITIONAL TESTING ONLY
if filter == 'F2100W':
    # bicone semi-major axis
    semi_major = 17
    # bicone semi-minor axis
    semi_minor = 15
    # point source aperture size
    PS_aper = 9
elif filter == 'F1800W':
    semi_major = 16
    semi_minor = 14
    PS_aper = 7.5
else:
    semi_major = 15
    semi_minor = 13
    if filter == 'F1500W':
        PS_aper = 6
    elif filter == 'F1000W':
        PS_aper = 5
    else:
        PS_aper = 3.5



# set import paths
path = 'Images/JWST/4_Deconvolution/'+str(Dir)+'/NGC5728/obs2/FITS/'
out_path = 'Images/JWST/4_Deconvolution/'+str(Dir)+'/NGC5728/obs2/Measurements/Update_25Apr2023/' + str(filter)


# read in deconvolved data and set the pixel scale and observation wavelength
decon_ims = get_pkg_data_filename(path + obs +  '.fits')
decon_list = fits.open(decon_ims)

# set the observations pixel scale
pixel_scale = float(0.1110)

# set the observations wavelength
if filter == 'F560W':
    obs_wave = 5.6
elif filter == 'F1000W':
    obs_wave = 10.0
elif filter == 'F1500W':
    obs_wave = 15.0
elif filter == 'F1800W':
    obs_wave = 18.0
elif filter == 'F2100W':
    obs_wave = 21.0

# set the data arrays
iter = []
fwhm_px = []
fwhm_arc = []
flux = []
flux_sum = []
cenX = []
cenY = []
change = []
change_fwhm = []
count = 0

# iterate through deconvolved data
# NOTE: index 0 only contains header info, no data. Start at index 1
for i in range(1, int(len(decon_list))):
    # # read in image and save to place-holder fits file
    decon = decon_list[i].data

    # save the image to a place-holder fits file
    save_FITS(decon, filename='residual',
              name='RL_Merit_function_holder',
              output_dir=path, pixel_scale=pixel_scale,
              filter=filter, wavelength=obs_wave, instrument='MIRI')

    # Read in the deconvolved FITS file
    decon_ims = get_pkg_data_filename(path + '/residual_RL_Merit_function_holder.fits')
    decon_list2 = fits.open(decon_ims)
    decon2 = decon_list2[0].data

    # measure centroid using a quadratic fit
    center_quad = centroid_quadratic(decon2,
                                     xpeak=guess[0],
                                     ypeak=guess[1],
                                     fit_boxsize=box)
    center = (center_quad[0], center_quad[1])


    # measure FWHM
    fwhm = MAE_measure_FWHM_gaussian(HDUlist=decon_ims,
                                     ext=0,
                                     threshold=threshold,
                                     plot=False,
                                     print_params=False,
                                     centroid=center)

    # measure photometry
    decon_flux = MAE_Aperture_Photometry(HDUlist=decon_ims,
                                   ext=0,
                                   radius=aper_rad,
                                   normalize=False,
                                   centroid=center,
                                   pad=pad)

    # 9. append measurements ALL
    # append the iteration
    iter.append(count)
    # append the FWHM (pixel)
    fwhm_px.append(fwhm[1])
    # append the FWHM (arcsec)
    fwhm_arc.append(fwhm[0])
    # append the flux measurements (ALL)
    flux.append(decon_flux)
    # append the total image flux
    flux_sum.append(decon.sum())
    # append the measured centroid
    cenX.append(center[0])
    cenY.append(center[1])

    # check progress on measurements
    if i==1:
        fluxP = decon_flux[0]
        change.append(fluxP)
        change_fwhm.append(fwhm[0])
        new_fwhm = fwhm[0]
        print('Iteration / centroid (x,y) / fwhm (") / flux (mJy/sr): ', i, center, new_fwhm, fluxP)
    else:
        fluxP = decon_flux[0]
        percent = np.abs(((fluxP - change[0]) / change[0]) * 100)
        change_fwhm.append(fwhm[0])
        percent_fwhm = np.abs(((change_fwhm[i - 1] - change_fwhm[i - 2]) / change_fwhm[i - 2]) * 100)
        new_fwhm = fwhm[0]
        print('Iteration / centroid (x,y) / fwhm (") / % change / flux (mJy/sr) / % change: ', i, center, new_fwhm, percent_fwhm, fluxP, percent)

    # increase the count
    count += 1

# write everything to output files:
lineA = []
lineB = []

# append the merit function measures
merit_nint = []
merit_flux = []
merit_fwhm = []
for i in range(0, len(iter)):
    # save the total selection criteria data: nint, aperture flux, fwhm (pix)
    line_select = [str(iter[i]) + '\t'
                   + str(flux[i][0]) + '\t\t'
                   + str(fwhm_arc[i]) + '\t\t'
                   ]

    # save the total data: nint, total flux, aperture flux, fwhm ("), fwhm (pixel), cenX, cenY
    line_total = [str(iter[i]) + '\t'
                  + str(flux_sum[i]) + '\t\t'
                  + str(flux[i][0]) + '\t'
                  + str(fwhm_arc[i]) + '\t'
                  + str(fwhm_px[i]) + '\t'
                  + str(cenX[i]) + '\t'
                  + str(cenY[i])
                  ]

    lineA.append(line_select)
    lineB.append(line_total)

    # save the merit function data: nint, flux, fwhm (")
    merit_nint.append(iter[i])
    merit_flux.append(flux[i][0])
    merit_fwhm.append(fwhm_arc[i])

# save selection merits to a single text file
with open(out_path + '/' + 'RL_SELECT_' + str(obs) + '_21um_aper.txt','w') as f:
     f.write('\n'.join(['\n'.join(lines) for lines in lineA]) + '\n')

# save total measurements to a single text file
with open(out_path + '/'  + 'RL_TOTAL_' + str(obs) + '_21um_aper.txt','w') as f:
     f.write('\n'.join(['\n'.join(lines) for lines in lineB]) + '\n')

#save the merit functions
data2 = ''
with open (out_path + '/' + 'RL_Image_Plot_' + str(obs) + '_21um_aper.txt', 'w') as fh:
   for a, b, c in zip(merit_nint, merit_flux, merit_fwhm):
       print('%.0f  %.5f  %.5f ' % (a, b, c), file = fh)

with open(out_path + '/' + 'RL_Image_Plot_' + str(obs) + '_21um_aper.txt', 'r') as fp:
   data2 = fp.read()

with open(out_path + '/' + 'RL_Image_Plot_' + str(obs) + '_21um_aper.txt', 'w') as fp:
   fp.write(data2)


# plot simple comparison
if plot:
    # for testing purposes
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(15, 10), tight_layout=True)
    ax1.set_title('FWHM')
    ax1.plot(iter, merit_fwhm, marker='o', color='b')
    ax1.plot(iter, merit_fwhm, color='r')
    ax1.set_ylabel('FWHM (")')
    ax2.set_title('Flux')
    ax2.plot(iter, merit_flux, marker='o', color='b')
    ax2.plot(iter, merit_flux, color='r')
    ax2.set_ylabel('Flux Density (mJy/sr)')
    fig.supxlabel('Iterations')
    fig.suptitle(str(obs))
    plt.show()

# Additional photometric testing
if photom_test:
    # set the data arrays
    iter = []
    flux_PS = []
    flux_bicone = []
    flux_galaxy = []
    aper21 = []
    count = 0

    # iterate through deconvolved data
    # NOTE: index 0 only contains header info, no data. Start at index 1
    #for i in range(1, int(len(decon_list))):
    for i in range(2, 251):
        # read in image and save to place-holder fits file
        photom = decon_list[i].data

        # save the image to a place-holder fits file
        save_FITS(photom, filename='residual',
                  name='RL_Merit_function_holder_PhotomTest',
                  output_dir=path, pixel_scale=pixel_scale,
                  filter=filter, wavelength=obs_wave, instrument='MIRI')

        # Read in the deconvolved FITS file
        photom_ims = get_pkg_data_filename(path + '/residual_RL_Merit_function_holder_PhotomTest.fits')
        photom_list = fits.open(photom_ims)
        photom_im = photom_list[0].data

        # measure centroid using a quadratic fit
        center_quad = centroid_quadratic(photom_im,
                                         xpeak=guess[0],
                                         ypeak=guess[1],
                                         fit_boxsize=box)
        center = (center_quad[0], center_quad[1])

        # measure photometry -> central point source
        PS_flux = MAE_Aperture_Photometry(HDUlist=photom_ims,
                                             ext=0,
                                             radius=aper_rad,
                                             normalize=False,
                                             centroid=center,
                                             pad=pad)

        bicone_flux = MAE_Aperture_Photometry(HDUlist=photom_ims,
                                       ext=0,
                                       radius=aper_rad,
                                       normalize=False,
                                       centroid=center,
                                       aperture='Ellipse',
                                       semi_major=semi_major,
                                       semi_minor=semi_minor,
                                       theta=55.0,
                                       pad=pad
                                       )

        # append measurements ALL
        # append the iteration
        iter.append(count)
        # append the flux measurements (point source)
        flux_PS.append(PS_flux)
        # append the flux measurements (bicone)
        flux_bicone.append(PS_flux)
        # increase the count
        count += 1

        # toggle whether to include galaxy measurements or not
        if flux_galaxy:
            # measure photometry -> galaxy
            galax_flux = MAE_Aperture_Photometry(HDUlist=photom_ims,
                                                 ext=0,
                                                 radius=galax_aper,
                                                 normalize=False,
                                                 centroid=center,
                                                 pad=pad)

            # append the flux measurements (bicone)
            flux_galaxy.append(galax_flux[0])
        else:
            flux_galaxy.append(0)

        # toggle whether to include the 21 um aperture size measurements
        if aper_21:
            # measure photometry using a fixed aperture -> set by 21 um image
            aper_flux21 = MAE_Aperture_Photometry(HDUlist=photom_ims,
                                                  ext=0,
                                                  radius=17.8,
                                                  normalize=False,
                                                  centroid=center,
                                                  pad=pad)

            # append the flux measurements (bicone)
            aper21.append(aper_flux21[0])
        else:
            aper21.append(0)

    # write everything to an output file
    lineC = []
    lineD = []

    # append the merit function measures
    merit_nint2 = []
    merit_PS_flux = []
    merit_bicone_flux = []
    merit_galaxy_flux = []
    merit_aper21_flux = []
    for i in range(0, len(iter)):
        # save the merit function data: nint, PS_flux, bicone_flux, galaxy_flux
        merit_nint2.append(iter[i])
        merit_PS_flux.append(flux_PS[i][0])
        merit_bicone_flux.append(flux_bicone[i][0])
        merit_galaxy_flux.append(flux_galaxy[i])
        merit_aper21_flux.append(aper21[i])


    # save the merit functions
    data2 = ''
    with open(out_path + '/' + '/Photometry_Tests/RL_PHOTOM_Testing_' +str(obs)+ '.txt', 'w') as fh:
        for a, b, c, d, e in zip(merit_nint2, merit_PS_flux, merit_bicone_flux, merit_galaxy_flux, merit_aper21_flux):
            print('%.0f  %.5f  %.5f  %.5f %.5f' % (a, b, c, d, e), file=fh)

    with open(out_path + '/' + '/Photometry_Tests/RL_PHOTOM_Testing_' +str(obs)+ '.txt', 'r') as fp:
        data2 = fp.read()

    with open(out_path + '/' + '/Photometry_Tests/RL_PHOTOM_Testing_' +str(obs)+ '.txt', 'w') as fp:
        fp.write(data2)

# close the FITS importers
decon_list.close()
decon_list2.close()
photom_list.close()