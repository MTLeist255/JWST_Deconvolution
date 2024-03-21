# 2022 Oct 27: Updated merit functions code to reflect current paper progress

from Convenience_Functions import *

# Plot the merit function results?
plot = True
# Perform aperture photometry on other regions in the image array?
photom_test = False
# if photom_test: fit an aperture across the background galaxy
flux_galaxy = False

# set image specific global params:
# set the deconvolution path
Dir = '5_Kraken'
# set the deconvolved surname
# sample type
sample = 'standard'
# image filter
filter = 'F1500W'
# model title
model = 'Model_AGN_complicated'
# model type
model_type = 'Ideal_Model'

# set measurement specific params
# centroid guess
guess = (127,127)
# centroid box size
box = 15
# FWHM threshold
threshold = 0
# photometry aperture radius
aper_rad = 17.8
# image pad level
pad = 0.5
# galaxy aperture
galax_aper = 75

# set wavelength specific params
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
# set the import path for the deconvolved image
path = 'Images/JWST/4_Deconvolution/'+str(Dir)+'/FITS/' + str(model_type) + '/' + str(model) + '/current_12Aug2023/'+ str(sample) + '/'
# set the output path for the merit function measurements
out_path = 'Images/JWST/4_Deconvolution/'+str(Dir)+'/Measurements/' + str(model_type) + '/' + str(model) + '/current_12Aug2023/' + str(sample) + '/'

# read in deconvolved data and set the pixel scale and observation wavelength
decon_ims = get_pkg_data_filename(path +str(filter)+'_'+str(sample)+'_model.fits')
#decon_ims = get_pkg_data_filename(path + 'objects_2100W_deconv_mrl_ovar2.fits')
decon_list = fits.open(decon_ims)
print('Deconvolution array list length: ', len(decon_list))
#dec_data = decon_list[0].data
#print(len(dec_data))

# set the observations pixel scale
pixel_scale = decon_list[1].header['PIXELSCL']
pixel_scale = float(pixel_scale)

# set the observations filter
filter = decon_list[1].header['FILTER']

# set the observations wavelength
obs_wave = decon_list[1].header['WAVELEN']
obs_wave = float(obs_wave)

# 4. set the data arrays
iter = []
fwhm_px = []
fwhm_arc = []
flux = []
flux_sum = []
cenX = []
cenY = []
count = 0
change_fwhm = []
change = []
fwhm_per = []
flux_per = []
# iterate through deconvolved data
# NOTE: index 0 only contains header info, no data. Start at index 1
for i in range(1,int(len(decon_list))):
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

    # 6. measure centroid using a quadratic fit
    center_quad = centroid_quadratic(decon2,
                                     xpeak=guess[0],
                                     ypeak=guess[1],
                                     fit_boxsize=box)
    center = (center_quad[0], center_quad[1])


    # 7. measure FWHM
    fwhm = MAE_measure_FWHM_gaussian(HDUlist=decon_ims,
                                     ext=0,
                                     threshold=threshold,
                                     plot=False,
                                     print_params=False,
                                     centroid=center)

    # 8. measure photometry
    decon_flux = MAE_Aperture_Photometry(HDUlist=decon_ims,
                                   ext=0,
                                   radius=aper_rad,
                                   normalize=False,
                                   centroid=center,
                                   pad=pad)

    # Print merit function progress
    if i==1:
        fluxP = decon_flux[0]
        change.append(fluxP)
        change_fwhm.append(fwhm[1])
        new_fwhm = fwhm[1]
        print('Iteration / centroid (x,y) / fwhm (pixel) / flux: ', i, center, new_fwhm, fluxP)
        fwhm_per.append(0)
        flux_per.append(0)
    else:
        fluxP = decon_flux[0]
        percent = np.abs(((fluxP - change[0]) / change[0]) * 100)
        change_fwhm.append(fwhm[1])
        percent_fwhm = np.abs(((change_fwhm[i - 1] - change_fwhm[i - 2]) / change_fwhm[i - 2]) * 100)
        new_fwhm = fwhm[1]
        print('Iteration / centroid (x,y) / fwhm (pixel) / % change / flux / % change: ', i, center, new_fwhm, percent_fwhm, fluxP, percent)
        # append the percent changes
        fwhm_per.append(percent_fwhm)
        flux_per.append(percent)

    # increase the count
    count += 1

    # append measurements ALL
    # append the iteration
    iter.append(count)
    # append the flux measurements (ALL)
    flux.append(decon_flux)
    # append the total image flux
    flux_sum.append(decon2.sum())
    # append the measured centroid
    cenX.append(center[0])
    cenY.append(center[1])
    # append the FWHM (pixel)
    fwhm_px.append(fwhm[1])
    # append the FWHM (arcsec)
    fwhm_arc.append(fwhm[0])


# 12. write everything to output files:
lineA = []
lineB = []

# append the merit function measures
merit_nint = []
merit_flux = []
merit_fwhm = []
for i in range(0, len(iter)):
    # save the total selection criteria data: nint, aperture flux, fwhm (pix)
    line_select = [str(iter[i]) + '\t'
                   + str(flux[i][0]) + '\t'
                   + str(flux_per[i]) + '\t'
                   + str(fwhm_px[i]) + '\t'
                   + str(fwhm_per[i])
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

    # append the measurements for writing
    lineA.append(line_select)
    lineB.append(line_total)

    # save the merit function data: nint, flux, fwhm (")
    merit_nint.append(iter[i])
    merit_flux.append(flux[i][0])
    merit_fwhm.append(fwhm_px[i])

# # save selection merits to a single text file
with open(out_path + str(filter) + '/SELECT_' + str(sample) + '.txt','w') as f:
     f.write('\n'.join(['\n'.join(lines) for lines in lineA]) + '\n')

# save total measurements to a single text file
with open(out_path +  str(filter) + '/TOTAL_' + str(sample) + '.txt','w') as f:
     f.write('\n'.join(['\n'.join(lines) for lines in lineB]) + '\n')

#save the merit functions
data2 = ''
with open (out_path + str(filter) + '/Image_Plot_' + str(sample) + '.txt', 'w') as fh:
   for a, b, c in zip(merit_nint, merit_flux, merit_fwhm):
       print('%.0f  %.5f  %.5f ' % (a, b, c), file = fh)

with open(out_path + str(filter) + '/Image_Plot_' + str(sample) + '.txt', 'r') as fp:
   data2 = fp.read()

with open(out_path + str(filter) + '/Image_Plot_'+str(sample) + '.txt', 'w') as fp:
   fp.write(data2)


# 11. plot simple comparison
if plot:
    # for testing purposes
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(15, 10), tight_layout=True)
    ax1.set_title('FWHM Limit')
    ax1.plot(iter, merit_fwhm, marker='o', color='b')
    ax1.plot(iter, merit_fwhm, color='r')
    ax1.set_ylabel('FWHM (pixel)')
    ax2.set_title('Flux Conservation')
    ax2.plot(iter, merit_flux, marker='o', color='b')
    ax2.plot(iter, merit_flux, color='r')
    ax2.set_ylabel('Aperture counts')
    fig.supxlabel('Iterations')
    fig.suptitle(str(filter) + ' ' + str(model) + '\nIters: ' + str(int(len(decon_list))-1), fontsize = 20)
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

    # 5. iterate through deconvolved data
    # NOTE: index 0 only contains header info, no data. Start at index 1
    for i in range(1, int(len(decon_list))):
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

        # measure the flux in the bicone
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

    # write everything to an output file
    lineC = []
    lineD = []

    # append the merit function measures
    merit_nint2 = []
    merit_PS_flux = []
    merit_bicone_flux = []
    merit_galaxy_flux = []
    for i in range(0, len(iter)):
        # save the merit function data: nint, PS_flux, bicone_flux, galaxy_flux
        merit_nint2.append(iter[i])
        merit_PS_flux.append(flux_PS[i][0])
        merit_bicone_flux.append(flux_bicone[i][0])
        merit_galaxy_flux.append(flux_galaxy[i])


    # save the merit functions
    data2 = ''
    with open(out_path + str(filter) + '/Photometry_tests/SCALED_PHOTOM_Testing_' +str(sample) + '.txt', 'w') as fh:
        for a, b, c, d in zip(merit_nint2, merit_PS_flux, merit_bicone_flux, merit_galaxy_flux):
            print('%.0f  %.5f  %.5f  %.5f' % (a, b, c, d), file=fh)

    with open(out_path + str(filter) + '/Photometry_tests/SCALED_PHOTOM_Testing_' + str(sample) + '.txt', 'r') as fp:
        data2 = fp.read()

    with open(out_path + str(filter) + '/Photometry_tests/SCALED_PHOTOM_Testing_' + str(sample) + '.txt', 'w') as fp:
        fp.write(data2)

# close the FITS importers
decon_list.close()
decon_list2.close()
photom_list.close()