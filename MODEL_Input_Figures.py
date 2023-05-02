# 17 Feb 2023: Simple code to plot a 2x2 figure of input model components and a 2x3 figure of reference PSFs

from Convenience_Functions import *

# display the model components
model = False

# Display the PSF grid
PSF = False

# set the output path
out_path = 'Images/JWST/4_Deconvolution/Total_Comparisons/'

# read in the input model data
mod_path = 'Images/JWST/1_Input_Models/Models/'

if model:
    # read in the model component data
    # single pixel model
    model1 = get_pkg_data_filename(mod_path + 'sp_component.fits')
    model1_list = fits.open(model1)
    model1_data = model1_list[0].data
    print('Point source total flux: ', model1_data.sum())

    # single pixel + disk model
    model2 = get_pkg_data_filename(mod_path + 'polar_component.fits')
    model2_list = fits.open(model2)
    model2_data = model2_list[0].data
    print('Polar dust total flux: ', model2_data.max())

    # single pixel + disk + bicone model
    model3 = get_pkg_data_filename(mod_path + 'bicone_component.fits')
    model3_list = fits.open(model3)
    model3_data = model3_list[0].data
    print('Bicone total flux: ', model3_data.sum())

    # toy AGN model
    model4 = get_pkg_data_filename(mod_path + 'galaxy_component2.fits')
    model4_list = fits.open(model4)
    model4_data = model4_list[0].data
    print('Galaxy total flux: ', model4_data.sum())

    # # plot the input model data
    # # plot the 2x2 grid image data

    fig, ax = plt.subplots(nrows=1, ncols=4, figsize=(9, 4), constrained_layout=True)
    # set the colorbar label
    label = 'Log-scaled Image Peak Counts'

    positions = [0.05, 2, 4.2]
    labels = ['-2', '-0', '2']
    positions2 = [0.5, 5, 9.5]
    labels2 = ['-5', '0', '5']
    positions3 = [1, 10, 20, 30, 38]
    labels3 = ['-20', '-10', '0', '10', '20']
    positions4 = [15, 65, 128, 193, 241]
    labels4 = ['-128', '-64', '0', '64', '128']

    # Complicated AGN
    im00 = ax[0].imshow(model4_data, origin='lower', norm=ImageNormalize(stretch=LogStretch(), vmin=5, vmax=4000),
                        cmap='afmhot')
    ax[0].set_xlabel('Galactic nuclei region\n ~5 kpc', fontsize=12)
    ax[0].hlines(108, 108, 148, color='chartreuse')
    ax[0].hlines(148, 108, 148, color='chartreuse')
    ax[0].vlines(108, 108, 148, color='chartreuse')
    ax[0].vlines(148, 108, 148, color='chartreuse')
    ax[0].xaxis.set_major_locator(ticker.FixedLocator(positions4))
    ax[0].xaxis.set_major_formatter(ticker.FixedFormatter(labels4))
    ax[0].yaxis.set_major_locator(ticker.FixedLocator(positions4))
    ax[0].yaxis.set_major_formatter(ticker.FixedFormatter(labels4))

    # Bicone
    im01 = ax[1].imshow(model3_data, origin='lower',
                        norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=model3_data.max()), cmap='afmhot')
    ax[1].set_xlabel('Ionization bicone\n ~400 pc', fontsize=12)
    ax[1].hlines(15, 15, 25, color='chartreuse')
    ax[1].hlines(25, 15, 25, color='chartreuse')
    ax[1].vlines(15, 15, 25, color='chartreuse')
    ax[1].vlines(25, 15, 25, color='chartreuse')
    ax[1].xaxis.set_major_locator(ticker.FixedLocator(positions3))
    ax[1].xaxis.set_major_formatter(ticker.FixedFormatter(labels3))
    ax[1].yaxis.set_major_locator(ticker.FixedLocator(positions3))
    ax[1].yaxis.set_major_formatter(ticker.FixedFormatter(labels3))

    # 2D-Gaussian disk
    im10 = ax[2].imshow(model2_data, origin='lower', cmap='afmhot',
                        norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=2e5), )
    ax[2].set_xlabel('Polar dust halo\n ~100 pc', fontsize=12)
    ax[2].hlines(2.5, 2.5, 7.5, color='chartreuse')
    ax[2].hlines(7.5, 2.5, 7.5, color='chartreuse')
    ax[2].vlines(2.5, 2.5, 7.5, color='chartreuse')
    ax[2].vlines(7.5, 2.5, 7.5, color='chartreuse')
    ax[2].xaxis.set_major_locator(ticker.FixedLocator(positions2))
    ax[2].xaxis.set_major_formatter(ticker.FixedFormatter(labels2))
    ax[2].yaxis.set_major_locator(ticker.FixedLocator(positions2))
    ax[2].yaxis.set_major_formatter(ticker.FixedFormatter(labels2))

    # Single pixel
    im11 = ax[3].imshow(model1_data, origin='lower', cmap='afmhot')
    # ax[0].text(7, 28, 'Dusty torus', color='w', fontsize=10)
    ax[3].set_xlabel('Central region\n ~17.5 pc', fontsize=12)
    ax[3].xaxis.set_major_locator(ticker.FixedLocator(positions))
    ax[3].xaxis.set_major_formatter(ticker.FixedFormatter(labels))
    ax[3].yaxis.set_major_locator(ticker.FixedLocator(positions))
    ax[3].yaxis.set_major_formatter(ticker.FixedFormatter(labels))

    ax[3].tick_params(direction='in', color='w', length=6, axis='both', left=True, top=False, right=False, bottom=True,
                      labelleft=True, labeltop=False, labelright=False, labelbottom=True, pad=-17, labelcolor='w',
                      labelsize='small')
    ax[2].tick_params(direction='in', color='w', length=6, axis='both', left=True, top=False, right=False, bottom=True,
                      labelleft=True, labeltop=False, labelright=False, labelbottom=True, pad=-17, labelcolor='w',
                      labelsize='small')
    ax[1].tick_params(direction='in', color='w', length=6, axis='both', left=True, top=False, right=False, bottom=True,
                      labelleft=True, labeltop=False, labelright=False, labelbottom=True, pad=-17, labelcolor='w',
                      labelsize='small')
    ax[0].tick_params(direction='in', color='w', length=6, axis='both', left=True, top=False, right=False, bottom=True,
                      labelleft=True, labeltop=False, labelright=False, labelbottom=True, pad=-20, labelcolor='w',
                      labelsize='small')

    fig.subplots_adjust(hspace=0, wspace=0)

    # fig.suptitle('Toy AGN Model Components', fontsize=15, y=0.95)
    fig.supxlabel('X (pixels)', fontsize=15)
    fig.supylabel('Y (pixels)', fontsize=15)
    plt.savefig(out_path + 'Model_Inputs_Grid_scaled4.png')
    plt.show()

if PSF:
    # read in the reference PSFs
    # 5.6 um
    psf_path = 'Images/JWST/1_Input_Models/PSFs/Reference_PSFs/OVERSAMP_4/'
    psf1 = get_pkg_data_filename(psf_path + 'WebbPSF_reference_PSF_F560W.fits')
    psf1_list = fits.open(psf1)
    psf1_data = psf1_list[3].data

    # 10 um
    psf2 = get_pkg_data_filename(psf_path + 'WebbPSF_reference_PSF_F1000W.fits')
    psf2_list = fits.open(psf2)
    psf2_data = psf2_list[3].data

    # 15 um
    psf3 = get_pkg_data_filename(psf_path + 'WebbPSF_reference_PSF_F1500W.fits')
    psf3_list = fits.open(psf3)
    psf3_data = psf3_list[3].data

    # 18 um
    psf4 = get_pkg_data_filename(psf_path + 'WebbPSF_reference_PSF_F1800W.fits')
    psf4_list = fits.open(psf4)
    psf4_data = psf4_list[3].data

    # 21 um
    psf5 = get_pkg_data_filename(psf_path + 'WebbPSF_reference_PSF_F2100W.fits')
    psf5_list = fits.open(psf5)
    psf5_data = psf5_list[3].data

    # plot the reference PSFs
    fig, ax = plt.subplots(nrows=1, ncols=5, figsize=(15, 4), tight_layout=True)

    # set the colorbar label
    label = 'Log-scaled Image Peak (counts/pixel)'

    positions4 = [5, 25, 50, 75, 100, 125, 145]
    labels4 = ['-75', '-50', '-25', '0', '25', '50', '75']

    # single pixel
    im00 = ax[0].imshow(psf1_data[53:203, 52:202], origin='lower',
                        norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=psf1_data.max()), cmap='afmhot')
    ax[0].text(54, 125, 'F560W', color='w', fontsize=15, weight='bold')
    ax[0].xaxis.set_major_locator(ticker.FixedLocator(positions4))
    ax[0].xaxis.set_major_formatter(ticker.FixedFormatter(labels4))
    ax[0].yaxis.set_major_locator(ticker.FixedLocator(positions4))
    ax[0].yaxis.set_major_formatter(ticker.FixedFormatter(labels4))
    # add_colorbar(im00, label=label, format=FormatStrFormatter('%.1e'))

    # single pixel + disk
    im01 = ax[1].imshow(psf2_data[53:203, 52:202], origin='lower',
                        norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=psf2_data.max()), cmap='afmhot')
    ax[1].text(50, 125, 'F1000W', color='w', fontsize=15, weight='bold')
    ax[1].xaxis.set_major_locator(ticker.FixedLocator(positions4))
    ax[1].xaxis.set_major_formatter(ticker.FixedFormatter(labels4))
    ax[1].yaxis.set_major_locator(ticker.FixedLocator(positions4))
    ax[1].yaxis.set_major_formatter(ticker.FixedFormatter(labels4))
    # add_colorbar(im01, label=label, format=FormatStrFormatter('%.1e'))

    # single pixel + disk + bicone
    im10 = ax[2].imshow(psf3_data[53:203, 52:202], origin='lower',
                        norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=psf3_data.max()), cmap='afmhot')
    ax[2].text(48, 125, 'F1500W', color='w', fontsize=15, weight='bold')
    ax[2].xaxis.set_major_locator(ticker.FixedLocator(positions4))
    ax[2].xaxis.set_major_formatter(ticker.FixedFormatter(labels4))
    ax[2].yaxis.set_major_locator(ticker.FixedLocator(positions4))
    ax[2].yaxis.set_major_formatter(ticker.FixedFormatter(labels4))
    # add_colorbar(im10, label=label, format=FormatStrFormatter('%.1e'))

    # complicated AGN
    im11 = ax[3].imshow(psf4_data[53:203, 52:202], origin='lower',
                        norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=psf4_data.max()), cmap='afmhot')
    ax[3].text(48, 125, 'F1800W', color='w', fontsize=15, weight='bold')
    ax[3].xaxis.set_major_locator(ticker.FixedLocator(positions4))
    ax[3].xaxis.set_major_formatter(ticker.FixedFormatter(labels4))
    ax[3].yaxis.set_major_locator(ticker.FixedLocator(positions4))
    ax[3].yaxis.set_major_formatter(ticker.FixedFormatter(labels4))
    # add_colorbar(im11, label=label, format=FormatStrFormatter('%.1e'))

    # complicated AGN
    im21 = ax[4].imshow(psf5_data[53:203, 52:202], origin='lower',
                        norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=psf5_data.max()), cmap='afmhot')
    ax[4].text(48, 125, 'F2100W', color='w', fontsize=15, weight='bold')
    ax[4].xaxis.set_major_locator(ticker.FixedLocator(positions4))
    ax[4].xaxis.set_major_formatter(ticker.FixedFormatter(labels4))
    ax[4].yaxis.set_major_locator(ticker.FixedLocator(positions4))
    ax[4].yaxis.set_major_formatter(ticker.FixedFormatter(labels4))
    # add_colorbar(im21, label=label, format=FormatStrFormatter('%.1e'))

    ax[0].tick_params(direction='in', color='w', length=6, axis='both', left=True, top=False, right=False, bottom=True,
                      labelleft=True, labeltop=False, labelright=False, labelbottom=True)
    ax[1].tick_params(direction='in', color='w', length=6, axis='both', left=True, top=False, right=False, bottom=True,
                      labelleft=True, labeltop=False, labelright=False, labelbottom=True)
    ax[2].tick_params(direction='in', color='w', length=6, axis='both', left=True, top=False, right=False, bottom=True,
                      labelleft=True, labeltop=False, labelright=False, labelbottom=True)
    ax[3].tick_params(direction='in', color='w', length=6, axis='both', left=True, top=False, right=False, bottom=True,
                      labelleft=True, labeltop=False, labelright=False, labelbottom=True)
    ax[4].tick_params(direction='in', color='w', length=6, axis='both', left=True, top=False, right=False, bottom=True,
                      labelleft=True, labeltop=False, labelright=False, labelbottom=True)

    fig.subplots_adjust(hspace=0, wspace=0)

    plt.suptitle('Reference PSFs', fontsize=20)
    fig.supxlabel('X (pixels)', fontsize=15)
    fig.supylabel('Y (pixels)', fontsize=15, x=0.007)
    plt.savefig(out_path + 'Model_Inputs_PSFs.png')
    plt.show()


