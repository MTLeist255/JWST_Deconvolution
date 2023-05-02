# 2022 Nov 3: code to generate (1) individual model SEDs and (2) 2x2 grid of SEDs for all four models

from Convenience_Functions import *

# plot individual SEDs
individual_SED = True

# plot flux measurements as a ratio instead of an SED
ratio = True

# set the input params
# model type
model_type = 'Ideal_Model'
# set the model index:
idx = 3
# model title: (1) single pixel, (2) single pixel + disk, (3) single pixel + disk + ionization cone, (4) complicated AGN
model = ['Model_single_pixel', 'Model_single_pixel_disk', 'Model_single_pixel_disk_bicone', 'Model_AGN_complicated']
# merit function time
title = ['Toy Unresolved Point Source Model', 'Toy Polar Dust and Unresolved Point Source Model',
         'Toy Extended Dust, Polar Dust, and Unresolved Point Source Model', 'Toy AGN Model']

# set the simulated wavelength
wave = [5.6, 10.0, 15.0, 18.0, 21.0]

sample = '_DET_SAMP_'

# Manually input the measurements: FWHM and aperture flux
if idx == 0:
    # model: single pixel
    # set the fwhm plotting values
    miri_fwhm = [1.76, 2.99, 4.48, 5.35, 6.2]
    RL_fwhm = [0.63, 0.65, 0.78, 1.44, 1.69]
    WH_fwhm = [0.9, 2.09, 3.47, 3.85, 4.29]
    ME_fwhm = [1.45, 2.31, 3.76, 4.37, 5.26]
    AIDA_fwhm = [0.9, 1.33, 2.98, 1.98, 2.64]
    Kraken_fwhm = [0.62, 0.65, 0.7, 0.72, 0.88]

    # set the flux plotting values
    miri_flux = [6.67e4, 8.65e5, 6.6e6, 6.9e6, 7.39e6]
    RL_flux = [7.51e4, 9.38e5, 7.26e6, 7.68e6, 8.29e6]
    WH_flux = [8.68e4, 1.13e6, 8.59e6, 8.84e6, 9.61e6]
    ME_flux = [8.67e4, 1.05e6, 8.12e6, 8.26e6, 9.29e6]
    AIDA_flux = [5.76e4, 7.63e5, 6.43e6, 6.53e6, 6.96e6]
    Kraken_flux = [7.77e4, 9.71e5, 7.51e6, 8.01e6, 8.65e6]
elif idx == 1:
    # model: single pixel and disk
    # set the fwhm plotting values
    miri_fwhm = [1.96, 3.18, 4.64, 5.44, 6.28]
    RL_fwhm = [0.69, 0.82, 1.66, 1.71, 2.11]
    WH_fwhm = [1.77, 2.48, 3.36, 3.88, 4.53]
    ME_fwhm = [1.66, 2.37, 3.41, 4.47, 5.04]
    AIDA_fwhm = [1.3, 1.75, 2.57, 3.18, 3.23]
    Kraken_fwhm = [0.68, 0.83, 1.91, 1.57, 1.98]

    # set the flux plotting values
    miri_flux = [7.56e4, 9.75e5, 7.48e6, 7.76e6, 8.35e6]
    RL_flux = [8.52e4, 1.05e6, 8.23e6, 8.27e6, 9.4e6]
    WH_flux = [9.1e4, 1.27e6, 9.73e6, 9.87e6, 1.07e7]
    ME_flux = [9.83e4, 1.16e6, 8.91e6, 9.56e6, 1.01e7]
    AIDA_flux = [7.13e4, 9.27e5, 7.19e6, 7.25e6, 8.11e6]
    Kraken_flux = [8.73e4, 1.08e6, 8.47e6, 8.93e6, 9.83e6]
elif idx == 2:
    # model: single pixel, disk and bicone
    # set the fwhm plotting values
    miri_fwhm = [2.58, 7.35, 11.37, 12.34, 13.44]
    RL_fwhm = [1.27, 2.34, 3.38, 3.78, 4.76]
    WH_fwhm = [1.65, 3.37, 6.59, 8.77, 10.59]
    ME_fwhm = [1.47, 2.46, 4.04, 4.90, 5.94]
    AIDA_fwhm = [1.54, 2, 4.05, 4.81, 6.29]
    Kraken_fwhm = [1, 1.83, 3.18, 3.49, 4.48]

    # set the flux plotting values
    miri_flux = [2.17e5, 2.8e6, 2.14e7, 2.29e7, 2.38e7]
    RL_flux = [2.49e5, 3.09e6, 2.41e7, 2.62e7, 2.77e7]
    WH_flux = [2.32e5, 3.33e6, 2.54e7, 2.72e7, 2.79e7]
    ME_flux = [2.8e5, 3.23e6, 2.5e7, 2.73e7, 2.87e7]
    AIDA_flux = [2.22e5, 2.85e6, 2.24e7, 2.44e7, 2.58e7]
    Kraken_flux = [2.54e5, 3.11e6, 2.48e7, 2.7e7, 2.86e7]
else:
    # model: toy AGN
    # set the fwhm plotting values
    miri_fwhm = [3.08, 4.87, 7.29, 8.34, 9.70]
    RL_fwhm = [1.21, 1.74, 2.95, 3.37, 4.37]
    WH_fwhm = [1.49, 2.43, 3.71, 4.28, 5.31]
    ME_fwhm = [1.43, 2.53, 3.93, 4.70, 5.68]
    AIDA_fwhm = [1.91, 3.02, 3.99, 4.87, 5.55]
    Kraken_fwhm = [1.19, 1.67, 2.92, 3.25, 4.22]

    # set the flux plotting values
    miri_flux = [2.81e5, 2.82e6, 2.75e7, 2.88e7, 3.06e7]
    RL_flux = [3.21e5, 3.09e6, 3.04e7, 3.26e7, 3.49e7]
    WH_flux = [3.66e5, 3.40e6, 3.36e7, 3.61e7, 3.88e7]
    ME_flux = [3.46e5, 3.18e6, 3.13e7, 3.36e7, 3.56e7]
    AIDA_flux = [2.84e5, 2.81e6, 2.81e7, 3.01e7, 3.25e7]
    Kraken_flux = [3.25e5, 3.12e6, 3.12e7, 3.33e7, 3.54e7]

# append the FWHM values to create a line
obs_line_fwhm = [miri_fwhm[0], miri_fwhm[1], miri_fwhm[2], miri_fwhm[3], miri_fwhm[4]]
RL_line_fwhm = [RL_fwhm[0], RL_fwhm[1], RL_fwhm[2], RL_fwhm[3], RL_fwhm[4]]
WH_line_fwhm = [WH_fwhm[0], WH_fwhm[1], WH_fwhm[2], WH_fwhm[3], WH_fwhm[4]]
MEM_line_fwhm = [ME_fwhm[0], ME_fwhm[1], ME_fwhm[2], ME_fwhm[3], ME_fwhm[4]]
AIDA_line_fwhm = [AIDA_fwhm[0], AIDA_fwhm[1], AIDA_fwhm[2], AIDA_fwhm[3], AIDA_fwhm[4]]
kraken_line_fwhm = [Kraken_fwhm[0], Kraken_fwhm[1], Kraken_fwhm[2], Kraken_fwhm[3], Kraken_fwhm[4]]

# append the flux values to create a line
obs_line_flux = [miri_flux[0], miri_flux[1], miri_flux[2], miri_flux[3], miri_flux[4]]
RL_line_flux = [RL_flux[0], RL_flux[1], RL_flux[2], RL_flux[3], RL_flux[4]]
WH_line_flux = [WH_flux[0], WH_flux[1], WH_flux[2], WH_flux[3], WH_flux[4]]
MEM_line_flux = [ME_flux[0], ME_flux[1], ME_flux[2], ME_flux[3], ME_flux[4]]
AIDA_line_flux = [AIDA_flux[0], AIDA_flux[1], AIDA_flux[2], AIDA_flux[3],AIDA_flux[4]]
kraken_line_flux = [Kraken_flux[0], Kraken_flux[1], Kraken_flux[2], Kraken_flux[3], Kraken_flux[4]]

# plot the individual SEDS
if individual_SED:
    # set the output path for the merit function measurements
    out_path = 'Images/JWST/4_Deconvolution/Total_Comparisons/'+str(model_type)+'/SEDs/Single/'

    # set the figure save title
    save_merit = out_path + '/SCALED_' + str(model[idx]) + '_DET_SAMP_All_diff_lim'

    # plot the SED
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(14, 7), tight_layout=True)

    # A) plot the FWHM values
    # plot the single line FWHM values
    ax[0].plot(wave, obs_line_fwhm, color='k', linestyle='--', label='MIRISim Toy Model')
    ax[0].plot(wave, RL_line_fwhm, color='r', label='Richardson-Lucy')
    ax[0].plot(wave, WH_line_fwhm, color='g', label='IWFT')
    ax[0].plot(wave, MEM_line_fwhm, color='b', label='Sparse')
    ax[0].plot(wave, AIDA_line_fwhm, color='y', label='AIDA')
    ax[0].plot(wave, kraken_line_fwhm, color='darkorange', label='Kraken')

    # plot the single line flux values
    ax[1].plot(wave, obs_line_flux, color='k', linestyle='--', label='MIRISim Toy Model')
    ax[1].plot(wave, RL_line_flux, color='r', label='Richardson-Lucy')
    ax[1].plot(wave, WH_line_flux, color='g', label='IWFT')
    ax[1].plot(wave, MEM_line_flux, color='b', label='Sparse')
    ax[1].plot(wave, AIDA_line_flux, color='y', label='AIDA')
    ax[1].plot(wave, kraken_line_flux, color='darkorange', label='Kraken')

    for i in range(0, int(len(wave))):
        # plot the FWHM values
        # plot the miri fwhm values
        ax[0].plot(wave[i], miri_fwhm[i], color='k', marker='X')
        # plot the RL fwhm values
        ax[0].plot(wave[i], RL_fwhm[i], color='r', marker='o')
        # plot the WH fwhm values
        ax[0].plot(wave[i], WH_fwhm[i], color='g', marker='o')
        # plot the MEM fwhm values
        ax[0].plot(wave[i], ME_fwhm[i], color='b', marker='o')
        # plot the AIDA fwhm values
        ax[0].plot(wave[i], AIDA_fwhm[i], color='y', marker='o')
        # plot the kraken FWHM values
        ax[0].plot(wave[i], Kraken_fwhm[i], color='darkorange', marker='o')

        # plot the flux values
        # plot the miri flux values
        ax[1].plot(wave[i], miri_flux[i], color='k', marker='X')
        # plot the RL flux values
        ax[1].plot(wave[i], RL_flux[i], color='r', marker='o')
        # plot the WH flux values
        ax[1].plot(wave[i], WH_flux[i], color='g', marker='o')
        # plot the MEM flux values
        ax[1].plot(wave[i], ME_flux[i], color='b', marker='o')
        # plot the AIDA flux values
        ax[1].plot(wave[i], AIDA_flux[i], color='y', marker='o')
        # plot the kraken flux values
        ax[1].plot(wave[i], Kraken_flux[i], color='darkorange', marker='o')

    # add the diffraction limits?
    # set the calculated diffraction limits
    diff_lim = [1.882, 2.982, 4.436, 5.373, 6.127]
    ax[0].hlines(diff_lim[0], 4.88, 6.42, color='k', linestyle='--')
    ax[0].text(6.9, diff_lim[0] - 0.055, 'F560W PSF FWHM', color='b', weight='bold')
    ax[0].hlines(diff_lim[1], 8.76, 11.1, color='k', linestyle='--')
    ax[0].text(3.1, diff_lim[1] - 0.075, 'F1000W PSF FWHM', color='g', weight='bold')
    ax[0].hlines(diff_lim[2], 13.1, 17.1, color='k', linestyle='--')
    ax[0].text(7, diff_lim[2] - 0.045, 'F1500W PSF FWHM', color='y', weight='bold')
    ax[0].hlines(diff_lim[3], 16.0, 20.3, color='k', linestyle='--')
    ax[0].text(10, diff_lim[3] - 0.06, 'F1800W PSF FWHM', color='darkorange', weight='bold')
    ax[0].hlines(diff_lim[4], 17.8, 24.4, color='k', linestyle='--')
    ax[0].text(11.8, diff_lim[4] - 0.045, 'F2100W PSF FWHM', color='r', weight='bold')

    # set the plot specific params
    # set the plot specific params
    ax[0].set_title('FWHM Limit', fontsize=15)
    ax[0].set_xlim(3, 25)
    ax[0].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax[0].xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax[0].set_ylabel('FWHM (pixel)', fontsize=15)
    ax[0].legend(loc='upper left')
    ax[1].set_title('Counts Conservation', fontsize=15)
    ax[1].set_xlim(3, 25)
    ax[1].yaxis.set_major_formatter(FormatStrFormatter('%.1e'))
    ax[1].xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax[1].set_ylabel('Aperture counts', fontsize=15)
    ax[1].legend(loc='upper left')

    fig.suptitle(title[idx], fontsize=20)
    fig.supxlabel('Simulated Wavelength (μm)', fontsize=15)
    plt.savefig(save_merit + '.png')
    plt.show()

# plot the flux values as a ratio instead of an SED
if ratio:
    # set the output path for the merit function measurements
    out_path = 'Images/JWST/4_Deconvolution/Total_Comparisons/'+str(model_type)

    # set the figure save title
    save_merit = out_path + '/SCALED_' + str(model[idx]) + '_DET_SAMP_All_diff_lim'

    # calculate the fwhm and flux ratios
    # calculate the flux/FWHM ratios
    RL_fwhm_rat = []
    WH_fwhm_rat = []
    ME_fwhm_rat = []
    AIDA_fwhm_rat = []
    Kraken_fwhm_rat = []
    fwhm_ylim_min = []
    fwhm_ylim_max = []

    RL_flux_rat = []
    WH_flux_rat = []
    ME_flux_rat = []
    AIDA_flux_rat = []
    Kraken_flux_rat = []
    flux_ylim_min = []
    flux_ylim_max = []

    # iterate through each wavelength
    for i in range(0,int(len(wave))):
        # fwhm ratio = technique/mirisim
        RL_fwhm_ratio = miri_fwhm[i]/RL_fwhm[i]
        WH_fwhm_ratio = miri_fwhm[i]/ WH_fwhm[i]
        ME_fwhm_ratio = miri_fwhm[i]/ ME_fwhm[i]
        AIDA_fwhm_ratio = miri_fwhm[i]/ AIDA_fwhm[i]
        Kraken_fwhm_ratio = miri_fwhm[i]/ Kraken_fwhm[i]

        # flux ratio = technique/mirisim
        RL_flux_ratio = miri_flux[i]/ RL_flux[i]
        WH_flux_ratio = miri_flux[i]/ WH_flux[i]
        ME_flux_ratio = miri_flux[i]/ME_flux[i]
        AIDA_flux_ratio = miri_flux[i]/ AIDA_flux[i]
        Kraken_flux_ratio = miri_flux[i]/ Kraken_flux[i]

        # append the fwhm ratios
        RL_fwhm_rat.append(RL_fwhm_ratio)
        WH_fwhm_rat.append(WH_fwhm_ratio)
        ME_fwhm_rat.append(ME_fwhm_ratio)
        AIDA_fwhm_rat.append(AIDA_fwhm_ratio)
        Kraken_fwhm_rat.append(Kraken_fwhm_ratio)

        # append the flux ratios
        RL_flux_rat.append(RL_flux_ratio)
        WH_flux_rat.append(WH_flux_ratio)
        ME_flux_rat.append(ME_flux_ratio)
        AIDA_flux_rat.append(AIDA_flux_ratio)
        Kraken_flux_rat.append(Kraken_flux_ratio)

        # set the y-lims for each plot based on the individual results
        total_fwhm = [RL_fwhm_ratio, WH_fwhm_ratio, ME_fwhm_ratio, AIDA_fwhm_ratio,Kraken_fwhm_ratio]
        total_flux = [RL_flux_ratio, WH_flux_ratio, ME_flux_ratio, AIDA_flux_ratio,Kraken_flux_ratio]

        # determine the minimum value
        y_min_fwhm = min(total_fwhm)-0.15
        y_min_flux = min(total_flux)-0.01

        # determine the maximum value
        y_max_fwhm = max(total_fwhm)+0.1
        y_max_flux = max(total_flux)+0.1

        # append the min/max values
        fwhm_ylim_min.append(y_min_fwhm)
        fwhm_ylim_max.append(y_max_fwhm)
        flux_ylim_min.append(y_min_flux)
        flux_ylim_max.append(y_max_flux)

    # set the min/max y-lims for each plot based on the overall results
    FWHM_min_y = min(fwhm_ylim_min)
    FWHM_max_y = max(fwhm_ylim_max)
    FLUX_min_y = min(flux_ylim_min)
    FLUX_max_y = max(flux_ylim_max)

    # plot the SED
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(14, 7), tight_layout=True)

    # iterate through each ratio value and plot
    for i in range(0, int(len(wave))):
        if i == 0:
            # plot the FWHM values
            # plot the RL fwhm ratios
            ax[0].plot(wave[i], RL_fwhm_rat[i], color='r', marker='o', label='Richardson-lucy', markersize=10,linestyle='None')
            # plot the WH fwhm ratios
            ax[0].plot(wave[i], WH_fwhm_rat[i], color='g', marker='D', label='IWFT', markersize=10, linestyle='None')
            # plot the MEM fwhm ratios
            ax[0].plot(wave[i], ME_fwhm_rat[i], color='b', marker='v', label='Sparse', markersize=10,
                       linestyle='None')
            # plot the AIDA fwhm ratios
            ax[0].plot(wave[i], AIDA_fwhm_rat[i], color='y', marker='s', label='AIDA', markersize=10, linestyle='None')
            # plot the Kraken fwhm ratios
            ax[0].plot(wave[i], Kraken_fwhm_rat[i], color='darkorange', marker='x', label = 'Kraken', markersize=10, linestyle='None')

            # plot the flux values
            # plot the RL flux ratios
            ax[1].plot(wave[i], RL_flux_rat[i], color='r', marker='o', label='Richardson-lucy', markersize=10,linestyle='None')
            # plot the WH flux ratios
            ax[1].plot(wave[i], WH_flux_rat[i], color='g', marker='D', label='IWFT', markersize=10, linestyle='None')
            # plot the MEM flux ratios
            ax[1].plot(wave[i], ME_flux_rat[i], color='b', marker='v', label='Sparse', markersize=10,
                       linestyle='None')
            # plot the AIDA flux ratios
            ax[1].plot(wave[i], AIDA_flux_rat[i], color='y', marker='s', label='AIDA', markersize=10, linestyle='None')
            # plot the Kraken flux ratios
            ax[1].plot(wave[i], Kraken_flux_rat[i], color='darkorange', marker='x', label = 'Kraken', markersize=10, linestyle='None')
        else:
            # plot the FWHM values
            # plot the RL fwhm ratios
            ax[0].plot(wave[i], RL_fwhm_rat[i], color='r', marker='o', markersize=10,linestyle='None')
            # plot the WH fwhm ratios
            ax[0].plot(wave[i], WH_fwhm_rat[i], color='g', marker='D', markersize=10, linestyle='None')
            # plot the MEM fwhm ratios
            ax[0].plot(wave[i], ME_fwhm_rat[i], color='b', marker='v', markersize=10,linestyle='None')
            # plot the AIDA fwhm ratios
            ax[0].plot(wave[i], AIDA_fwhm_rat[i], color='y', marker='s', markersize=10, linestyle='None')
            # plot the Kraken fwhm ratios
            ax[0].plot(wave[i], Kraken_fwhm_rat[i], color='darkorange', marker='x', markersize=10, linestyle='None')

            # plot the flux values
            # plot the RL flux ratios
            ax[1].plot(wave[i], RL_flux_rat[i], color='r', marker='o', markersize=10,linestyle='None')
            # plot the WH flux ratios
            ax[1].plot(wave[i], WH_flux_rat[i], color='g', marker='D', markersize=10, linestyle='None')
            # plot the MEM flux ratios
            ax[1].plot(wave[i], ME_flux_rat[i], color='b', marker='v', markersize=10, linestyle='None')
            # plot the AIDA flux ratios
            ax[1].plot(wave[i], AIDA_flux_rat[i], color='y', marker='s', markersize=10, linestyle='None')
            # plot the Kraken flux ratios
            ax[1].plot(wave[i], Kraken_flux_rat[i], color='darkorange', marker='x', markersize=10, linestyle='None')


    # set the plot specific params
    ax[0].set_title('FWHM Limit', fontsize=15)
    ax[0].set_xlim(3, 25)
    ax[0].set_ylim([FWHM_min_y, FWHM_max_y])
    ax[0].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax[0].xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax[0].set_ylabel('MIRISim / Deconvolved', fontsize = 15)
    ax[0].legend(loc='upper left')

    # set the plot specific params
    ax[1].set_title('Counts Conservation', fontsize = 15)
    ax[1].set_ylim([FLUX_min_y, FLUX_max_y])
    ax[1].set_xlim(3, 25)
    ax[1].set_ylabel('MIRISim / Deconvolved', fontsize = 15)
    ax[1].legend(loc='upper left')

    fig.suptitle(title[idx] + ' Ratios', fontsize = 20)
    fig.supxlabel('Simulated Wavelength (μm)', fontsize = 15)
    plt.savefig(save_merit + '.png')
    plt.show()
