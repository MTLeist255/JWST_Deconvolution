# 18 Jan 2020: Source code to read in the Kraken output files, measure the image statistics, determine the converged
# image criteria, and save all stats

from Convenience_Functions import *

# save the merit functions to a .txt file and corresponding plot?
save_merits = False
# plot the final merit functions?
plot_merits = False
# save the final FITS in x,y coordinates
save_final_FITS = False

# set the import filter
filter = ['F2100W']
# Set the user-defined radius equal to the filter PSF diffraction limit
# MIRI PSFS diffraction limits (in pixels): F560W (1.882), F1000W (2.982), F1500W (4.436), F1800W (5.373), F2100W (6.127)
# source: https://jwst-docs.stsci.edu/jwst-mid-infrared-instrument/miri-performance/miri-point-spread-functions
r = [6.127]

# set the number of times deconvolved
dec_iter = 102
skip = 2

# set the observation object name
object = 'NGC4388'
# set the deconvolution method used
method = 'Kraken_LSE'
# set the reference PSF used
psf = 'standard_psf'

# set the input/output paths
global_path = 'Images/JWST/5_ERS_Data/5_Deconvolution_results/'+object+'/'
im_path = global_path + 'Deconvolution/Padded_images/'
out_path = global_path + 'Deconvolution/Results/'+method+'/'+psf+'/'
final_dec_path = out_path + 'final_images/'
metric_path = out_path + 'merits/'

# orange-peel testing paths
#final_dec_path = global_path + 'Results/orange_peel_tests/'+method+'L/toy_AGN_model/'
#metric_path = final_dec_path + 'merits/'

# Iterate through each FITS file in the deconvolution directory, appending the image data to a single array
for i in range(0, int(len(filter))):
    # check-in
    print('\n######################### Check-in: Object: ', object,
          ' | Method: ', method,
          ' | Reference PSF: ', psf,
          ' | Observed filter: ', filter[i],
          '#########################')


    # NOTE: the below images have their merit functions measured
    # read in the padded reference MIRI data
    obj_name = object + '_' + filter[i] + '_1024.fits'
    image = get_pkg_data_filename(im_path + obj_name)
    image_list = fits.open(image)
    image_data = image_list[0].data
    # save the header info
    sci_hdr = image_list[0].header

    # read in the deconvolved data
    dec_name = object + '_' + filter[i] + '_DECONV.fits'
    hdul = fits.open(out_path + dec_name)

    # build arrays to save everything
    # image results
    image_arr = []
    final_fits = []

    # full deconvolution results
    iter = []
    ffwhm = []
    delta_B_fwhm = []
    delta_T_fwhm = []
    fflux = []
    delta_B_flux = []
    delta_T_flux = []
    s_param = []
    delta_B_S = []
    delta_T_S = []
    std_dev = []

    # delta index values
    idx_fwhm = []
    idx_flux = []
    idx_S = []

    # set the stopping criteria
    stop = False
    # NOTE: for Kraken LSE and AMRL the 1st iteration is a noise estimate image, the 2nd image is the original observation,
    # start iterations at index 2
    for j in range(skip, dec_iter):
        if j == skip:
            # NOTE: again we are skipping the 1st-2 iterations as these are not part of the deconvolved results
            # the very first measurement should be of the original image
            dec_data = image_data

            # the following params are for the SMoothness parameter (set once and hold fixed)
            # set the radius for the boxcar filter size (smoothness parameter)
            smooth_radius = r[i]

            # Compute the image centroid
            # NOTE: assumed to be the brightest pixel in the image array
            imy1, imx1 = np.unravel_index(np.argmax(dec_data, axis=None), dec_data.shape)

            # find the encircled energy of 99% of the flux of the original image
            # NOTE: sets the outer radius for the smoothness parameter
            cog = COG(dec_data, (imy1, imx1), max_level=99, params=False)
            boxsize = cog[0]

        else:
            data = hdul[0].data
            dec_data = data[j]

        # guess the centroid position based on the image peak
        im_y, im_x = np.unravel_index(np.argmax(dec_data, axis=None), dec_data.shape)
        # determine the centroid from a quadratic fit using the image peak as the initial guess
        xycent = centroid_quadratic(dec_data, xpeak=im_y, ypeak=im_x, fit_boxsize=15)

        # measure merit functions
        fwhm, flux, S, centroid = measure_merits(dec_data, 17.8, boxsize, smooth_radius, xycent, convert=True,sanity_check=False, plot=False)


        # determine stopping criteria
        if j == skip:
            image_data = dec_data

            # save the initial fwhm, flux, and S value
            initial_fwhm = fwhm
            initial_flux = flux
            initial_S = S

            # save the initial fwhm, flux, and S values to an array
            idx_fwhm.append(fwhm)
            idx_flux.append(flux)
            idx_S.append(S)

            # NOTE: to keep the .txt outputs consistent for index = 0 append a blank space in the delta measurements
            delta_fwhm_total = 0.0
            delta_fwhm_between = 0.0
            delta_flux_total = 0.0
            delta_flux_between = 0.0
            delta_S_total = 0.0
            delta_S_between = 0.0
            stddev = 0

            # print the initial results
            print(j - skip, ' iterations: ',
                  'centroid (y,x): ', centroid,
                  ' | FWHM ("): ', '%.5f' % fwhm,
                  ' | Δ_total (%): ', '%.5f' % delta_fwhm_total,
                  ' | Δ_between (%): ', '%.5f' % delta_fwhm_between,
                  ' | Aperture flux (Jy): ', '%.5f' % flux,
                  ' | Δ_total (%): ', '%.5f' % delta_flux_total,
                  ' | Δ_between (%): ', '%.5f' % delta_flux_between,
                  ' | S-index: ', '%.5f' % S,
                  ' | stddev (σ): ''%.5f' % stddev,
                  ' | Δ_total (%): ', '%.5f' % delta_S_total,
                  ' | Δ_between (%): ', '%.5f' % delta_S_between,
                  )
        else:
            # measure Δ
            # FWHM
            delta_fwhm_total = calculate_difference(initial_fwhm, fwhm)
            delta_fwhm_between = calculate_difference(idx_fwhm[j - skip], fwhm)
            # flux
            delta_flux_total = calculate_difference(initial_flux, flux)
            delta_flux_between = calculate_difference(idx_flux[j - skip], flux)
            # S
            delta_S_total = calculate_difference(initial_S, S)
            delta_S_between = calculate_difference(idx_S[j - skip], S)
            stddev = np.std(s_param[:j])

            # FWHM LIMIT MET
            # Save the FITS file and the final merit data if the flux criteria is met before the FWHM criteria
            if delta_fwhm_between < 0.1:
                if stop != True:
                    stop = True
                    print('\n ##################################### FWHM limit reached- iteration: ', j-skip,'#####################################')
                    print('\n', j - skip, ' iterations: ',
                          'centroid (y,x): ', centroid,
                          ' | FWHM ("): ', '%.5f' % fwhm,
                          ' | Δ_total (%): ', '%.5f' % delta_fwhm_total,
                          ' | Δ_between (%): ', '%.5f' % delta_fwhm_between,
                          ' | Aperture flux (Jy): ', '%.5f' % flux,
                          ' | Δ_total (%): ', '%.5f' % delta_flux_total,
                          ' | Δ_between (%): ', '%.5f' % delta_flux_between,
                          ' | S-index: ', '%.5f' % S,
                          ' | stddev (σ): ''%.5f' % stddev,
                          ' | Δ_total (%): ', '%.5f' % delta_S_total,
                          ' | Δ_between (%): ', '%.5f' % delta_S_between,
                          )
                    print('\n######################################################################################################################\n')

                    # append all the information and update the FITS header
                    final_fits.append(dec_data)
                    final_iter = j-skip
                    final_fwhm = '%.5f' % fwhm
                    final_flux = '%.5f' % flux
                    final_S = '%.5f' % S
                    sci_hdr['CONVCRT'] = ('FWHM', 'Criteria Criteria image converged to')

                else:
                    print(j - skip, ' iterations: ',
                          'centroid (y,x): ', centroid,
                          ' | FWHM ("): ', '%.5f' % fwhm,
                          ' | Δ_total (%): ', '%.5f' % delta_fwhm_total,
                          ' | Δ_between (%): ', '%.5f' % delta_fwhm_between,
                          ' | Aperture flux (Jy): ', '%.5f' % flux,
                          ' | Δ_total (%): ', '%.5f' % delta_flux_total,
                          ' | Δ_between (%): ', '%.5f' % delta_flux_between,
                          ' | S-index: ', '%.5f' % S,
                          ' | stddev (σ): ''%.5f' % stddev,
                          ' | Δ_total (%): ', '%.5f' % delta_S_total,
                          ' | Δ_between (%): ', '%.5f' % delta_S_between,
                          )

            # FLUX LIMIT MET
            # Save the FITS file and the final merit data if the flux criteria is met before the FWHM criteria
            elif delta_flux_total > 30:
                # check if criteria is met
                if stop != True:
                    stop = True
                    print('\n ##################################### Flux limit reached- iteration: ', j - skip, '#####################################')
                    print('\n', j - skip, ' iterations: ',
                          'centroid (y,x): ', centroid,
                          ' | FWHM ("): ', '%.5f' % fwhm,
                          ' | Δ_total (%): ', '%.5f' % delta_fwhm_total,
                          ' | Δ_between (%): ', '%.5f' % delta_fwhm_between,
                          ' | Aperture flux (Jy): ', '%.5f' % flux,
                          ' | Δ_total (%): ', '%.5f' % delta_flux_total,
                          ' | Δ_between (%): ', '%.5f' % delta_flux_between,
                          ' | S-index: ', '%.5f' % S,
                          ' | stddev (σ): ''%.5f' % stddev,
                          ' | Δ_total (%): ', '%.5f' % delta_S_total,
                          ' | Δ_between (%): ', '%.5f' % delta_S_between,
                          )
                    print(
                        '\n######################################################################################################################\n')

                    # append all the information and update the FITS header
                    final_fits.append(dec_data)
                    final_iter = j - skip
                    final_fwhm = '%.5f' % fwhm
                    final_flux = '%.5f' % flux
                    final_S = '%.5f' % S
                    sci_hdr['CONVCRT'] = ('FWHM', 'Criteria Criteria image converged to')

                else:
                    print(j - skip, ' iterations: ',
                          'centroid (y,x): ', centroid,
                          ' | FWHM ("): ', '%.5f' % fwhm,
                          ' | Δ_total (%): ', '%.5f' % delta_fwhm_total,
                          ' | Δ_between (%): ', '%.5f' % delta_fwhm_between,
                          ' | Aperture flux (Jy): ', '%.5f' % flux,
                          ' | Δ_total (%): ', '%.5f' % delta_flux_total,
                          ' | Δ_between (%): ', '%.5f' % delta_flux_between,
                          ' | S-index: ', '%.5f' % S,
                          ' | stddev (σ): ''%.5f' % stddev,
                          ' | Δ_total (%): ', '%.5f' % delta_S_total,
                          ' | Δ_between (%): ', '%.5f' % delta_S_between,
                          )

            # Save the FITS file and the final merit data if the flux criteria is met before the FWHM criteria
            elif stddev > 0.02:
                if stop != True:
                    stop = True
                    print('\n ##################################### S-index standard deviation > 0.02: ', j-skip,'#####################################')
                    print('\n', j - skip, ' iterations: ',
                          'centroid (y,x): ', centroid,
                          ' | FWHM ("): ', '%.5f' % fwhm,
                          ' | Δ_total (%): ', '%.5f' % delta_fwhm_total,
                          ' | Δ_between (%): ', '%.5f' % delta_fwhm_between,
                          ' | Aperture flux (Jy): ', '%.5f' % flux,
                          ' | Δ_total (%): ', '%.5f' % delta_flux_total,
                          ' | Δ_between (%): ', '%.5f' % delta_flux_between,
                          ' | S-index: ', '%.5f' % S,
                          ' | stddev (σ): ''%.5f' % stddev,
                          ' | Δ_total (%): ', '%.5f' % delta_S_total,
                          ' | Δ_between (%): ', '%.5f' % delta_S_between,
                          )
                    print('\n######################################################################################################################\n')

                    # append all the information and update the FITS header
                    final_fits.append(dec_data)
                    final_iter = j-skip
                    final_fwhm = '%.5f' % fwhm
                    final_flux = '%.5f' % flux
                    final_S = '%.5f' % S
                    sci_hdr['CONVCRT'] = ('S-deviation', 'Criteria Criteria image converged to')

                else:
                    print(j - skip, ' iterations: ',
                          'centroid (y,x): ', centroid,
                          ' | FWHM ("): ', '%.5f' % fwhm,
                          ' | Δ_total (%): ', '%.5f' % delta_fwhm_total,
                          ' | Δ_between (%): ', '%.5f' % delta_fwhm_between,
                          ' | Aperture flux (Jy): ', '%.5f' % flux,
                          ' | Δ_total (%): ', '%.5f' % delta_flux_total,
                          ' | Δ_between (%): ', '%.5f' % delta_flux_between,
                          ' | S-index: ', '%.5f' % S,
                          ' | stddev (σ): ''%.5f' % stddev,
                          ' | Δ_total (%): ', '%.5f' % delta_S_total,
                          ' | Δ_between (%): ', '%.5f' % delta_S_between,
                          )

            else:
                print(j - skip, ' iterations: ',
                      'centroid (y,x): ', centroid,
                      ' | FWHM ("): ', '%.5f' % fwhm,
                      ' | Δ_total (%): ', '%.5f' % delta_fwhm_total,
                      ' | Δ_between (%): ', '%.5f' % delta_fwhm_between,
                      ' | Aperture flux (Jy): ', '%.5f' % flux,
                      ' | Δ_total (%): ', '%.5f' % delta_flux_total,
                      ' | Δ_between (%): ', '%.5f' % delta_flux_between,
                      ' | S-index: ', '%.5f' % S,
                      ' | stddev (σ): ''%.5f' % stddev,
                      ' | Δ_total (%): ', '%.5f' % delta_S_total,
                      ' | Δ_between (%): ', '%.5f' % delta_S_between,
                      )

        # append all of the results
        # imaging
        image_arr.append(data)

        # full deconvolution results
        iter.append(j-skip)
        ffwhm.append(fwhm)
        delta_B_fwhm.append(delta_fwhm_total)
        delta_T_fwhm.append(delta_fwhm_between)
        fflux.append(flux)
        delta_B_flux.append(delta_flux_total)
        delta_T_flux.append(delta_flux_between)
        s_param.append(S)
        delta_B_S.append(delta_S_total)
        delta_T_S.append(delta_S_between)
        std_dev.append(stddev)

        # index for Δ calculations
        idx_fwhm.append(fwhm)
        idx_flux.append(flux)
        idx_S.append(S)


    if save_merits:
        # save the final merit functions to a text file
        # write the metrics to a .txt file for additional sanity checking
        lineA = []
        lineB = []

        # append the merit function measures
        merit_nint = []
        merit_flux = []
        merit_fwhm = []
        merit_smooth = []
        for k in range(0, len(iter)):
            # save the total selection criteria data: nint, aperture flux, fwhm (pix) + the delta (total and between_ values
            line_select = [str(iter[k]) + '\t'
                           + str(ffwhm[k]) + '\t'
                           + str(delta_B_fwhm[k]) + '\t'
                           + str(delta_T_fwhm[k]) + '\t'
                           + str(fflux[k]) + '\t'
                           + str(delta_B_flux[k]) + '\t'
                           + str(delta_T_flux[k]) + '\t'
                           + str(s_param[k]) + '\t'
                           + str(std_dev[k]) + '\t'
                           + str(delta_B_S[k]) + '\t'
                           + str(delta_T_S[k])
                           ]
            lineA.append(line_select)

            # save the merit function data: nint, flux, fwhm (")
            merit_nint.append(iter[k])
            merit_fwhm.append(ffwhm[k])
            merit_flux.append(fflux[k])
            merit_smooth.append(s_param[k])

        # save selection merits to a single text file
        with open(metric_path + '/' + str(filter[i]) + '_SELECT_metrics.txt', 'w') as f:
            f.write('\n'.join(['\n'.join(lines) for lines in lineA]) + '\n')

        # save the merit functions for later plotting
        data2 = ''
        with open(metric_path + '/' + str(filter[i]) + '_PLOT_metrics.txt', 'w') as fh:
            for a, b, c, d in zip(merit_nint, merit_flux, merit_fwhm, merit_smooth):
                print('%.0f  %.5f  %.5f %.5f' % (a, b, c, d), file=fh)

        with open(metric_path + '/' + str(filter[i]) + '_PLOT_metrics.txt', 'r') as fp:
            data2 = fp.read()

        with open(metric_path + '/' + str(filter[i]) + '_PLOT_metrics.txt', 'w') as fp:
            fp.write(data2)

    if plot_merits:
        # plot the final merit function results
        fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, figsize=(20, 6), tight_layout=True)
        ax1.set_title('FWHM Limit')
        ax1.plot(merit_nint, merit_fwhm, marker='o', color='b')
        ax1.plot(merit_nint, merit_fwhm, color='r')
        ax1.set_ylabel('FWHM (")')

        ax2.plot(merit_nint, merit_flux, marker='o', color='b')
        ax2.plot(merit_nint, merit_flux, color='r')
        ax2.set_title('Flux Conservation')
        ax2.set_ylabel('Aperture flux (Jy)')

        ax3.plot(merit_nint, merit_smooth, marker='o', color='b')
        ax3.plot(merit_nint, merit_smooth, color='r')
        ax3.set_title('S-index')
        ax3.set_ylabel('S')

        ax1.set_xlim(0, dec_iter)
        ax2.set_xlim(0, dec_iter)
        ax3.set_xlim(0, dec_iter)

        fig.supxlabel('Iterations')
        fig.suptitle(object + ' ' + filter[i] + ' ' + method, fontsize=20)
        plt.savefig(metric_path + filter[i] + '_metrics_plot.png')
        plt.show()

    if save_final_FITS:
        # Save the final image data
        dec_img_data = final_fits[0]

        # save the intial emrits to 5 decimal places
        obs_fwhm = '%.5f' % initial_fwhm
        obs_flux = '%.5f' % initial_flux
        obs_S = '%.5f' % initial_S

        # save the final input image
        obs_name = object + '_' + filter[i] + '_FINAL.fits'
        sci_hdr['FILTER'] = (filter[i], 'Observation filter')
        sci_hdr['FWHM'] = (float(obs_fwhm), 'Observed image FWHM (")')
        sci_hdr['FLUX'] = (float(obs_flux), 'Observed image aperture flux(Jy)')
        sci_hdr['S-INDEX'] = (float(obs_S), 'Observed image smoothness index')
        fits.writeto(final_dec_path + obs_name, image_data, sci_hdr, overwrite=True)

        # save the final deconvolved image
        sci_hdr['FILTER'] = (filter[i], 'Observation filter')
        sci_hdr['METHOD'] = (method, 'Deconvolution method used')
        sci_hdr['CONITER'] = (int(final_iter), 'Iteration number image converged at')
        sci_hdr['FWHM'] = (float(final_fwhm), 'Converged image FWHM (")')
        sci_hdr['FLUX'] = (float(final_flux), 'Converged image aperture flux(Jy)')
        sci_hdr['S-INDEX'] = (float(final_S), 'Converged image smoothness index')
        dec_name = object + '_' + filter[i] + '_FINAL_DECONV.fits'
        fits.writeto(final_dec_path + dec_name, dec_img_data, sci_hdr, overwrite=True)
