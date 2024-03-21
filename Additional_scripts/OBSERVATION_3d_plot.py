from Convenience_Functions import *

# 9 May 2023: Author: Amilcar Torres-Quijano
# Edited by: Mason Leist
# License: MIT License
#
# Copyright (c) 2021 Amílcar R. Torres-Quijano
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# display the 3D-signals and corresponding observed/deconvolved images for a single filter
single_3D = False

# display the 3D-signals and corresponding observed/deconvolved images for all five filters in a horizontal 4x5 grid
grid_3D_horizontal = True

# display the 3D-signals and corresponding observed/deconvolved images for all five filters in a vertical 5x4 grid
grid_3D_vertical = False

# display images peak to corresponding deconvolved filter-> for the grid_3D_horizontal and grid_3D_vertical plots ONLY
norm_indiv = True
# display images peak to the F2100W deconvolved image-> for the grid_3D_horizontal and grid_3D_vertical plots ONLY
norm_f2100w = False

if single_3D:
    # display the 3D-signals and corresponding observed/deconvolved images for a single filter
    # set the observed filter
    filter = 'F2100W'

    # set the box size based on the aperture fit size (radius = 17.8 pixels, diameter = ~36 pixels)
    box = 18

    # set the import path
    wcs_path = 'Images/JWST/5_ERS_Data/5_Deconvolution_results/NGC5728_obs2/WCS_test/'
    # set the output path
    out_path = 'Images/JWST/5_ERS_Data/5_Deconvolution_results/NGC5728_obs2/Kraken/3D_signals/'

    # Read in the image array FITS files
    jwst_miri = wcs_path + str(filter) + '_OBS_final.fits'
    deconvolved = wcs_path + str(filter) + '_DEC_final.fits'

    # obs centroid (y,x): 127,127
    # deconvolved centroid (y,x): 129,128
    # set the observed centroid
    ceny, cenx = 100, 100

    # read in the observed data
    input_file = get_pkg_data_filename(jwst_miri)
    hdu = fits.open(input_file)
    data = hdu[0].data

    # read in the deconvolved data
    dec_input_file = get_pkg_data_filename(deconvolved)
    dec_hdu = fits.open(dec_input_file)
    dec_data = dec_hdu[0].data

    # set the normalization scale parameters based on the deconvolved image peak
    if filter == 'F560W':
        scale = 15564.4
    elif filter == 'F1000W':
        scale = 5498.14
    elif filter == 'F1500W':
        scale = 46488.5
    elif filter == 'F1800W':
        scale = 49645.9
    else:
        scale = 61463.5

    # normalize the observed data array to the deconvolved data array: 0->1
    data = data/scale
    dec_data = dec_data/scale

    # trim the image array values for 3D-signal plotting
    data_3d = data[ceny - box:ceny + box, cenx - box:cenx + box]
    dec_data_3d = dec_data[ceny - box:ceny + box, cenx - box:cenx + box]

    # save the x,y coordinates for each 3D-signal plot
    obsY = range(data_3d.shape[0])
    obsX = range(data_3d.shape[1])
    decY = range(dec_data_3d.shape[0])
    decX = range(dec_data_3d.shape[1])

    # plot the data
    fig = plt.figure(figsize=(10, 8), constrained_layout = True)
    # set the RA and Dec labels for the images
    labels = ['-10', '-5', '0', '5', '10']
    positions = [10, 50, 100, 150, 190]
    # set the RA and Dec labels for the 3D-signals
    labels_3d = ['-2', '-1', '0', '1', '2']
    positions_3d = [1, 9, 18, 27, 35]

    # plot the observed image array with box focusing on central 36 pixels
    ax1 = fig.add_subplot(221)
    im1 = ax1.imshow(data, origin='lower',norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=1), cmap='RdYlBu_r')
    ax1.set_title(str(filter) + ' JWST/MIRI', fontsize = 15)
    ax1.hlines(82, 82, 118, color='k')
    ax1.hlines(118, 82, 118, color='k')
    ax1.vlines(82, 82, 118, color='k')
    ax1.vlines(118, 82, 118, color='k')
    add_colorbar(im1, label = 'Normalized Flux')
    # set the figure axis labels and titles
    ax1.set_xlabel('RA-offset (")')
    ax1.set_ylabel('DEC-offset (")')
    ax1.xaxis.set_major_locator(ticker.FixedLocator(positions))
    ax1.xaxis.set_major_formatter(ticker.FixedFormatter(labels))
    ax1.yaxis.set_major_locator(ticker.FixedLocator(positions))
    ax1.yaxis.set_major_formatter(ticker.FixedFormatter(labels))
    ax1.tick_params(direction='in', color='w', length=6, left=True, top=False, right=False, bottom=True, labelleft=True,labeltop=False, labelright=False, labelbottom=True)

    # plot the deconvolved image array with box focusing on central 36 pixels
    ax2 = fig.add_subplot(222)
    im2 = ax2.imshow(dec_data, origin='lower',norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=1), cmap='RdYlBu_r')
    ax2.set_title(str(filter) + ' Deconvolved', fontsize = 15)
    ax2.hlines(82, 82, 118, color='k')
    ax2.hlines(118, 82, 118, color='k')
    ax2.vlines(82, 82, 118, color='k')
    ax2.vlines(118, 82, 118, color='k')
    add_colorbar(im2, label = 'Normalized Flux')
    # set the figure axis labels and titles
    ax2.set_xlabel('RA-offset (")')
    ax2.set_ylabel('DEC-offset (")')
    ax2.xaxis.set_major_locator(ticker.FixedLocator(positions))
    ax2.xaxis.set_major_formatter(ticker.FixedFormatter(labels))
    ax2.yaxis.set_major_locator(ticker.FixedLocator(positions))
    ax2.yaxis.set_major_formatter(ticker.FixedFormatter(labels))
    ax2.tick_params(direction='in', color='w', length=6, left=True, top=False, right=False, bottom=True, labelleft=True, labeltop=False, labelright=False, labelbottom=True)
    # Set orientation compass: x,y
    # create the height of the arrow
    # size = [txt1-x, txt1-y, txt2-x, txt2-y,
    # arrow1-x, arrow1-y, arrow1-xlen, arrow1-ylen,
    # arrow2-x, arrow2-y, arrow2-xlen, arrow2-ylen]
    size = [177, 43, 141, 6,
            183, 10, 0, 25,
            183, 10, -25, 0
            ]
    # add text
    ax2.text(size[0], size[1], 'N', color='w', fontsize=12)
    ax2.text(size[2], size[3], 'E', color='w', fontsize=12)
    arrow1 = plt.arrow(size[4], size[5], size[6], size[7], width=1, color='w')
    # create the length of the arrow
    arrow2 = plt.arrow(size[8], size[9], size[10], size[11], width=1, color='w')

    # add arrows
    ax2.add_patch(arrow1)
    ax2.add_patch(arrow2)

    # plot the 3D-signal figures: observed
    ax3 = fig.add_subplot(223, projection='3d')
    obs_X, obs_Y = np.meshgrid(obsX,obsY)
    ax3.plot_surface(obs_X,obs_Y,data_3d,cmap='magma',norm=SymLogNorm(linthresh=100000, linscale=0.5, vmin=None, vmax=None, clip=False, base=10))
    #set the Z-scale
    ax3.set_zlim(0, 1)
    #Default Viewing Angle
    ax3.view_init(8,-130)
    # set the figure axis labels and titles
    ax3.set_xlabel('RA-offset (")')
    ax3.set_ylabel('DEC-offset (")')
    ax3.set_zlabel('Normalized Flux')
    ax3.xaxis.set_major_locator(ticker.FixedLocator(positions_3d))
    ax3.xaxis.set_major_formatter(ticker.FixedFormatter(labels_3d))
    ax3.yaxis.set_major_locator(ticker.FixedLocator(positions_3d))
    ax3.yaxis.set_major_formatter(ticker.FixedFormatter(labels_3d))

    # plot the 3D-signal figures: observed
    ax4 = fig.add_subplot(224, projection='3d')
    dec_X, dec_Y = np.meshgrid(decX,decY)
    ax4.plot_surface(dec_X,dec_Y,dec_data_3d,cmap='magma',norm=SymLogNorm(linthresh=100000, linscale=0.5, vmin=None, vmax=None, clip=False, base=10))
    #set the Z-scale
    ax4.set_zlim(0, 1)
    #Default Viewing Angle
    ax4.view_init(8,-130)
    # set the figure axis labels and titles
    ax4.set_xlabel('RA-offset (")')
    ax4.set_ylabel('DEC-offset (")')
    ax4.set_zlabel('Normalized Flux')
    ax4.xaxis.set_major_locator(ticker.FixedLocator(positions_3d))
    ax4.xaxis.set_major_formatter(ticker.FixedFormatter(labels_3d))
    ax4.yaxis.set_major_locator(ticker.FixedLocator(positions_3d))
    ax4.yaxis.set_major_formatter(ticker.FixedFormatter(labels_3d))

    plt.savefig(out_path + str(filter) + '_single_grid.png')
    plt.show()

if grid_3D_horizontal:
    # display the 3D-signals and corresponding observed/deconvolved images for a single filter
    # set the box size based on the aperture fit size (radius = 17.8 pixels, diameter = ~36 pixels)
    box = 18
    # set the import path
    wcs_path = 'Images/JWST/5_ERS_Data/5_Deconvolution_results/NGC5728_obs2/WCS_test/'
    # set the output path
    out_path = 'Images/JWST/5_ERS_Data/5_Deconvolution_results/NGC5728_obs2/Kraken/3D_signals/'
    # obs centroid (y,x): 127,127
    # deconvolved centroid (y,x): 129,128
    # set the observed centroid
    ceny, cenx = 100, 100
    # set the normalization scale parameters based on the deconvolved image peak
    # NOTE: only one of these can be set at a time

    # display the 3D-signals and corresponding observed/deconvolved images for all five filters
    filters = ['F560W', 'F1000W', 'F1500W', 'F1800W', 'F2100W']
    # set the deconvolved image peaks
    norm_factor = [15564.4, 6471.72, 46488.5, 49645.9, 61463.5]

    # read in the cropped + WCS rotated FITS files, append the images to an array
    obs_imgs = []
    dec_imgs = []
    obs_3d = []
    dec_3d = []
    obs_X = []
    obs_Y = []
    dec_X = []
    dec_Y = []
    for i in range(0, int(len(filters))):
        jwst_miri = wcs_path + str(filters[i]) + '_OBS_final.fits'
        deconvolved = wcs_path + str(filters[i]) + '_DEC_final.fits'

        # read in the observed data
        input_file = get_pkg_data_filename(jwst_miri)
        hdu = fits.open(input_file)
        data = hdu[0].data

        # read in the deconvolved data
        dec_input_file = get_pkg_data_filename(deconvolved)
        dec_hdu = fits.open(dec_input_file)
        dec_data = dec_hdu[0].data

        if norm_indiv:
            scale = norm_factor[i]
            # set the picture save name
            save_name = 'Grid_horizontal_peak_indiv.png'
            label = 'Normalized Flux'

        if norm_f2100w:
            scale = 61463.5
            # set the picture save name
            save_name = 'Grid_horizontal_peak_f2100w.png'
            label = 'Normalized Flux'

        # normalize the observed data array to the deconvolved data array: 0->1
        data = data / scale
        dec_data = dec_data / scale

        # trim the image array values for 3D-signal plotting
        data_3d = data[ceny - box:ceny + box, cenx - box:cenx + box]
        dec_data_3d = dec_data[ceny - box:ceny + box, cenx - box:cenx + box]

        # save the x,y coordinates for each 3D-signal plot
        obsY = range(data_3d.shape[0])
        obsX = range(data_3d.shape[1])
        decY = range(dec_data_3d.shape[0])
        decX = range(dec_data_3d.shape[1])

        # crop images
        data = data[81:119, 81:119]
        dec_data = dec_data[81:119, 81:119]

        # append the images to individual arrays
        obs_imgs.append(data)
        dec_imgs.append(dec_data)
        obs_3d.append(data_3d)
        dec_3d.append(dec_data_3d)
        obs_X.append(obsX)
        obs_Y.append(obsY)
        dec_X.append(decX)
        dec_Y.append(decY)

    # plot the data
    fig = plt.figure(figsize=(24,16), tight_layout = True)
    # set the RA and Dec labels for the images
    labels = ['-2', '-1', '0', '1', '2']
    positions = [1, 9, 18, 27, 35]
    # set the RA and Dec labels for the 3D-signals
    labels_3d = ['-2', '-1', '0', '1', '2']
    positions_3d = [1, 9, 18, 27, 35]

    for i in range(0, int(len(filters))):
        # plot the observed image array with box focusing on central 36 pixels
        # horizontal array
        ax1 = fig.add_subplot(4,5,i+1)
        im1 = ax1.imshow(obs_imgs[i], origin='lower', norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=1), cmap='RdYlBu_r')
        ax1.set_title(str(filters[i]), fontsize=20)

        # set the figure axis labels and titles
        ax1.set_xlabel('RA-offset (")')
        ax1.set_ylabel('DEC-offset (")')
        ax1.xaxis.set_major_locator(ticker.FixedLocator(positions))
        ax1.xaxis.set_major_formatter(ticker.FixedFormatter(labels))
        ax1.yaxis.set_major_locator(ticker.FixedLocator(positions))
        ax1.yaxis.set_major_formatter(ticker.FixedFormatter(labels))
        ax1.tick_params(direction='in', color='w', length=6, left=True, top=False, right=False, bottom=True,labelleft=True, labeltop=False, labelright=False, labelbottom=True)

        # plot the deconvolved image array with box focusing on central 36 pixels
        # horizontal array
        ax2 = fig.add_subplot(4,5,i+6)
        im2 = ax2.imshow(dec_imgs[i], origin='lower', norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=1),cmap='RdYlBu_r')

        if i == 4:
            add_colorbar(im1, label=label)
            add_colorbar(im2, label=label)

        # if i == 4:
        #     # Set orientation compass: x,y
        #     # create the height of the arrow
        #     # size = [txt1-x, txt1-y, txt2-x, txt2-y,
        #     # arrow1-x, arrow1-y, arrow1-xlen, arrow1-ylen,
        #     # arrow2-x, arrow2-y, arrow2-xlen, arrow2-ylen]
        #     size = [177, 43, 141, 6,
        #             183, 10, 0, 25,
        #             183, 10, -25, 0
        #             ]
        #     # add text
        #     ax2.text(size[0], size[1], 'N', color='w', fontsize=12)
        #     ax2.text(size[2], size[3], 'E', color='w', fontsize=12)
        #     arrow1 = plt.arrow(size[4], size[5], size[6], size[7], width=1, color='w')
        #     # create the length of the arrow
        #     arrow2 = plt.arrow(size[8], size[9], size[10], size[11], width=1, color='w')

        # set the figure axis labels and titles
        ax2.set_xlabel('RA-offset (")')
        ax2.set_ylabel('DEC-offset (")')
        ax2.xaxis.set_major_locator(ticker.FixedLocator(positions))
        ax2.xaxis.set_major_formatter(ticker.FixedFormatter(labels))
        ax2.yaxis.set_major_locator(ticker.FixedLocator(positions))
        ax2.yaxis.set_major_formatter(ticker.FixedFormatter(labels))
        ax2.tick_params(direction='in', color='w', length=6, left=True, top=False, right=False, bottom=True,
                        labelleft=True, labeltop=False, labelright=False, labelbottom=True)

        # plot the 3D-signal figures: observed
        # horizontal array
        ax3 = fig.add_subplot(4, 5, i + 11, projection='3d')
        obsx, obsy = np.meshgrid(obs_X[i], obs_Y[i])
        print(len(obsx))
        ax3.plot_surface(obsx, obsy, obs_3d[i], cmap='RdYlBu_r',norm=SymLogNorm(linthresh=100000, linscale=0.5, vmin=None, vmax=None, clip=False, base=10))
        # set the Z-scale
        ax3.set_zlim(0, 1)
        # Default Viewing Angle
        ax3.view_init(8, -130)
        # set the figure axis labels and titles
        ax3.set_xlabel('RA-offset (")')
        ax3.set_ylabel('DEC-offset (")')
        ax3.zaxis.set_rotate_label(False)
        ax3.set_zlabel(label, rotation=90)
        ax3.xaxis.set_major_locator(ticker.FixedLocator(positions_3d))
        ax3.xaxis.set_major_formatter(ticker.FixedFormatter(labels_3d))
        ax3.yaxis.set_major_locator(ticker.FixedLocator(positions_3d))
        ax3.yaxis.set_major_formatter(ticker.FixedFormatter(labels_3d))

        # plot the 3D-signal figures: observed
        # horizontal array
        ax4 = fig.add_subplot(4,5,i+16, projection='3d')
        decx, decy = np.meshgrid(dec_X[i], dec_Y[i])
        ax4.plot_surface(decx, decy, dec_3d[i], cmap='RdYlBu_r',norm=SymLogNorm(linthresh=100000, linscale=0.5, vmin=None, vmax=None, clip=False, base=10))
        # set the Z-scale
        ax4.set_zlim(0, 1)
        # Default Viewing Angle
        ax4.view_init(8, -130)
        # set the figure axis labels and titles
        ax4.set_xlabel('RA-offset (")')
        ax4.set_ylabel('DEC-offset (")')
        ax4.zaxis.set_rotate_label(False)
        ax4.set_zlabel(label, rotation=90)
        ax4.xaxis.set_major_locator(ticker.FixedLocator(positions_3d))
        ax4.xaxis.set_major_formatter(ticker.FixedFormatter(labels_3d))
        ax4.yaxis.set_major_locator(ticker.FixedLocator(positions_3d))
        ax4.yaxis.set_major_formatter(ticker.FixedFormatter(labels_3d))

        if i == 0:
            ax1.annotate('JWST/MIRIM', xy = (-1, 0.5), xytext = (-1,5), xycoords = ax1.yaxis.label,
                         textcoords = 'offset points', fontsize = 20, ha = 'right', va = 'center', rotation = 'vertical')
            ax2.annotate('Deconvolved', xy = (-1, 0.5), xytext = (-1,5), xycoords = ax2.yaxis.label,
                         textcoords = 'offset points', fontsize = 20, ha = 'right', va = 'center', rotation = 'vertical')
            ax3.annotate('JWST/MIRIM', xy = (-1, 0.5), xytext = (-1,5), xycoords = ax3.zaxis.label,
                         textcoords = 'offset points', fontsize = 20, ha = 'right', va = 'center', rotation = 'vertical')
            ax4.annotate('Deconvolved', xy = (-1, 0.5), xytext = (-1,5), xycoords = ax4.zaxis.label,
                         textcoords = 'offset points', fontsize = 20, ha = 'right', va = 'center', rotation = 'vertical')


    plt.savefig(out_path + save_name)
    plt.show()

if grid_3D_vertical:
    # display the 3D-signals and corresponding observed/deconvolved images for a single filter
    # set the box size based on the aperture fit size (radius = 17.8 pixels, diameter = ~36 pixels)
    box = 18
    # set the import path
    wcs_path = 'Images/JWST/5_ERS_Data/5_Deconvolution_results/NGC5728_obs2/WCS_test/'
    # set the output path
    out_path = 'Images/JWST/5_ERS_Data/5_Deconvolution_results/NGC5728_obs2/Kraken/3D_signals/'
    # obs centroid (y,x): 127,127
    # deconvolved centroid (y,x): 129,128
    # set the observed centroid
    ceny, cenx = 100, 100
    # set the normalization scale parameters based on the deconvolved image peak
    # NOTE: only one of these can be set at a time

    # display the 3D-signals and corresponding observed/deconvolved images for all five filters
    filters = ['F560W', 'F1000W', 'F1500W', 'F1800W', 'F2100W']
    # set the deconvolved image peaks
    norm_factor = [15564.4, 5498.14, 46488.5, 49645.9, 61463.5]

    # read in the cropped + WCS rotated FITS files, append the images to an array
    obs_imgs = []
    dec_imgs = []
    obs_3d = []
    dec_3d = []
    obs_X = []
    obs_Y = []
    dec_X = []
    dec_Y = []
    for i in range(0, int(len(filters))):
        jwst_miri = wcs_path + str(filters[i]) + '_OBS_final.fits'
        deconvolved = wcs_path + str(filters[i]) + '_DEC_final.fits'

        # read in the observed data
        input_file = get_pkg_data_filename(jwst_miri)
        hdu = fits.open(input_file)
        data = hdu[0].data

        # read in the deconvolved data
        dec_input_file = get_pkg_data_filename(deconvolved)
        dec_hdu = fits.open(dec_input_file)
        dec_data = dec_hdu[0].data
        #dec_data[np.isnan(dec_data)] = 0

        if norm_indiv:
            scale = norm_factor[i]
            # set the picture save name
            save_name = 'Grid_vertical_peak_indiv.png'
            label = 'Normalized Flux'

        if norm_f2100w:
            scale = 61463.5
            # set the picture save name
            save_name = 'Grid_vertical_peak_f2100w.png'
            label = 'Normalized Flux'

        # normalize the observed data array to the deconvolved data array: 0->1
        data = data / scale
        data[np.isnan(data)] = 0
        dec_data = dec_data / scale

        print(str(filters[i]) + ' JWST/MIRI scaled to F2100W peak:\t ' + str(data.max()))
        print(str(filters[i]) + ' Kraken scaled to F2100W peak:\t ' + str(dec_data.max()))
        print('--------------------------------------------------')

        # trim the image array values for 3D-signal plotting
        data_3d = data[ceny - box:ceny + box, cenx - box:cenx + box]
        dec_data_3d = dec_data[ceny - box:ceny + box, cenx - box:cenx + box]

        # save the x,y coordinates for each 3D-signal plot
        obsY = range(data_3d.shape[0])
        obsX = range(data_3d.shape[1])
        decY = range(dec_data_3d.shape[0])
        decX = range(dec_data_3d.shape[1])

        # append the images to individual arrays
        obs_imgs.append(data)
        dec_imgs.append(dec_data)
        obs_3d.append(data_3d)
        dec_3d.append(dec_data_3d)
        obs_X.append(obsX)
        obs_Y.append(obsY)
        dec_X.append(decX)
        dec_Y.append(decY)

    # plot the data
    fig = plt.figure(figsize=(20,21), constrained_layout = True)
    # set the RA and Dec labels for the images
    labels = ['-10', '-5', '0', '5', '10']
    positions = [10, 50, 100, 150, 190]
    # set the RA and Dec labels for the 3D-signals
    labels_3d = ['-2', '-1', '0', '1', '2']
    positions_3d = [1, 9, 18, 27, 35]

    count1 = 1
    count2 = 2
    count3 = 3
    count4 = 4
    for i in range(0, 5):

        # plot the observed image array with box focusing on central 36 pixels
        # horizontal array
        ax1 = fig.add_subplot(5,4,i+count1)
        im1 = ax1.imshow(obs_imgs[i], origin='lower', norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=1), cmap='RdYlBu_r')

        if i == 0:
            ax1.set_title('JWST/MIRI', fontsize=25)
            ax1.text(70, 180, str(filters[i]), color='w', fontsize=18, weight='bold')
        else:
            ax1.text(65, 180, str(filters[i]), color='w', fontsize=18, weight='bold')

        ax1.hlines(82, 82, 118, color='k')
        ax1.hlines(118, 82, 118, color='k')
        ax1.vlines(82, 82, 118, color='k')
        ax1.vlines(118, 82, 118, color='k')
        add_colorbar(im1, label=label)
        # set the figure axis labels and titles
        ax1.set_xlabel('RA-offset (")')
        ax1.set_ylabel('DEC-offset (")')
        ax1.xaxis.set_major_locator(ticker.FixedLocator(positions))
        ax1.xaxis.set_major_formatter(ticker.FixedFormatter(labels))
        ax1.yaxis.set_major_locator(ticker.FixedLocator(positions))
        ax1.yaxis.set_major_formatter(ticker.FixedFormatter(labels))
        ax1.tick_params(direction='in', color='w', length=6, left=True, top=False, right=False, bottom=True,labelleft=True, labeltop=False, labelright=False, labelbottom=True)


        # plot the 3D-signal figures: observed
        # horizontal array
        ax2 = fig.add_subplot(5, 4, i+count2, projection='3d')
        obsx, obsy = np.meshgrid(obs_X[i], obs_Y[i])
        ax2.plot_surface(obsx, obsy, obs_3d[i], cmap='magma',norm=SymLogNorm(linthresh=100000, linscale=0.5, vmin=None, vmax=None, clip=False, base=10))
        # set the Z-scale
        ax2.set_zlim(0, 1)
        # Default Viewing Angle
        ax2.view_init(8, -130)
        # set the figure axis labels and titles
        ax2.set_xlabel('RA-offset (")')
        ax2.set_ylabel('DEC-offset (")')
        ax2.set_zlabel(label)
        ax2.xaxis.set_major_locator(ticker.FixedLocator(positions_3d))
        ax2.xaxis.set_major_formatter(ticker.FixedFormatter(labels_3d))
        ax2.yaxis.set_major_locator(ticker.FixedLocator(positions_3d))
        ax2.yaxis.set_major_formatter(ticker.FixedFormatter(labels_3d))

        # plot the deconvolved image array with box focusing on central 36 pixels
        # horizontal array
        ax3 = fig.add_subplot(5,4,i+count3)
        im3 = ax3.imshow(dec_imgs[i], origin='lower', norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=1),cmap='RdYlBu_r')

        if i == 0:
            ax3.set_title('Deconvolved', fontsize=25)
            ax3.text(70, 180, str(filters[i]), color='w', fontsize=18, weight='bold')
        else:
            ax3.text(65, 180, str(filters[i]), color='w', fontsize=18, weight='bold')


        ax3.hlines(82, 82, 118, color='k')
        ax3.hlines(118, 82, 118, color='k')
        ax3.vlines(82, 82, 118, color='k')
        ax3.vlines(118, 82, 118, color='k')
        add_colorbar(im3, label=label)

        if i == 4:
            # Set orientation compass: x,y
            # create the height of the arrow
            # size = [txt1-x, txt1-y, txt2-x, txt2-y,
            # arrow1-x, arrow1-y, arrow1-xlen, arrow1-ylen,
            # arrow2-x, arrow2-y, arrow2-xlen, arrow2-ylen]
            size = [179, 43, 142, 7,
                    183, 10, 0, 25,
                    183, 10, -25, 0
                    ]
            # add text
            ax3.text(size[0], size[1], 'N', color='w', fontsize=12)
            ax3.text(size[2], size[3], 'E', color='w', fontsize=12)
            arrow1 = plt.arrow(size[4], size[5], size[6], size[7], width=1, color='w')
            # create the length of the arrow
            arrow2 = plt.arrow(size[8], size[9], size[10], size[11], width=1, color='w')

            # add arrows
            ax3.add_patch(arrow1)
            ax3.add_patch(arrow2)

        # set the figure axis labels and titles
        ax3.set_xlabel('RA-offset (")')
        ax3.set_ylabel('DEC-offset (")')
        ax3.xaxis.set_major_locator(ticker.FixedLocator(positions))
        ax3.xaxis.set_major_formatter(ticker.FixedFormatter(labels))
        ax3.yaxis.set_major_locator(ticker.FixedLocator(positions))
        ax3.yaxis.set_major_formatter(ticker.FixedFormatter(labels))
        ax3.tick_params(direction='in', color='w', length=6, left=True, top=False, right=False, bottom=True,
                        labelleft=True, labeltop=False, labelright=False, labelbottom=True)

        # plot the 3D-signal figures: observed
        # horizontal array
        ax4 = fig.add_subplot(5,4,i+count4, projection='3d')
        decx, decy = np.meshgrid(dec_X[i], dec_Y[i])
        ax4.plot_surface(decx, decy, dec_3d[i], cmap='magma',norm=SymLogNorm(linthresh=100000, linscale=0.5, vmin=None, vmax=None, clip=False, base=10))
        # set the Z-scale
        ax4.set_zlim(0, 1)
        # Default Viewing Angle
        ax4.view_init(8, -130)
        # set the figure axis labels and titles
        ax4.set_xlabel('RA-offset (")')
        ax4.set_ylabel('DEC-offset (")')
        ax4.set_zlabel(label)
        ax4.xaxis.set_major_locator(ticker.FixedLocator(positions_3d))
        ax4.xaxis.set_major_formatter(ticker.FixedFormatter(labels_3d))
        ax4.yaxis.set_major_locator(ticker.FixedLocator(positions_3d))
        ax4.yaxis.set_major_formatter(ticker.FixedFormatter(labels_3d))

        # update the count
        count1 += 3
        count2 += 3
        count3 += 3
        count4 += 3

    fig.subplots_adjust(hspace=0, wspace=0)
    plt.savefig(out_path + save_name)
    plt.show()

#################################################### PRACTICE CODE #####################################################

# def plot_3d_log(file_name, title, ext):
#     '''
#     :param file_name: user-defined path to input FITS file
#     :param title: user-defined title for image plot
#     :param ext: user-defined extension for data
#     :return:
#     3d-signal plot
#     '''
#
#     # read in the image data
#     input_file = get_pkg_data_filename(file_name)
#     hdu = fits.open(input_file)
#     data = hdu[ext].data
#
#     # trim the image array to focus on the central 30x30 x,y coordinates around the nucleus
#     #              Y        X
#     data = data[ceny-box:ceny+box, cenx-box:cenx+box]
#     ##THE STUFF FOR 3D PLOTS
#     x = range(data.shape[1])
#     y = range(data.shape[0])
#     fig = plt.figure()
#     ax = fig.add_subplot(111, projection='3d')
#     X, Y = np.meshgrid(x,y)
#     test = ax.plot_surface(X,Y,data,cmap='magma',norm=SymLogNorm(linthresh=100000, linscale=0.5, vmin=None, vmax=None, clip=False, base=10))
#     print(np.amax(data))
#     ax.set_title(title, y=0.9)
#     #set the Z-scale
#     # NOTE: each image Z-scale is set to the corresponding deconvolved image peak
#     ax.set_zlim(0, scale)
#     #Default Viewing Angle
#     ax.view_init(8,-130)
#
#     ax.set_xlabel('X Position (px)')
#     ax.set_ylabel('Y Position (px)')
#     ax.zaxis.set_major_formatter(FormatStrFormatter('%.0e'))
#     ax.set_zlabel('Counts')
#     plt.savefig(out_path + 'Deconvolved_' + str(filter) + '_log_3d.png')
#     plt.show()
#     return
#
# def plot_3d_log_wire(file_name, title, ext):
#     '''
#     :param file_name: user-defined path to input FITS file
#     :param title: user-defined title for image plot
#     :param ext: user-defined extension for data
#     :return:
#     3d-signal plot
#     '''
#     # read in the image array
#     input_file = get_pkg_data_filename(file_name)
#     hdu = fits.open(input_file)
#     data = hdu[ext].data
#
#     # trim the image array to focus on the central 30x30 x,y coordinates around the nucleus
#     #              Y        X
#     data = data[ceny-box:ceny+box, cenx-box:cenx+box]
#
#     ##THE STUFF FOR 3D PLOTS
#     x = range(data.shape[1])
#     y = range(data.shape[0])
#
#     fig = plt.figure()
#     ax = fig.add_subplot(111, projection='3d')
#     X, Y = np.meshgrid(x,y)
#     test = ax.plot_wireframe(X,Y,data,cmap='magma',norm=SymLogNorm(linthresh=100000, linscale=0.5, vmin=None, vmax=None, clip=False, base=10))
#     print(np.amax(data))
#     ax.set_title(title, y=0.9)
#
#     # NOTE: each image Z-scale is set to the corresponding deconvolved image peak
#     ax.set_zlim(0, scale)
#
#     #Default Viewing Angle
#     ax.view_init(8,-130)
#
#     ax.set_xlabel('X Position (px)')
#     ax.set_ylabel('Y Position (px)')
#     ax.set_zlabel('Counts')
#     ax.zaxis.set_major_formatter(FormatStrFormatter('%.0e'))
#     plt.savefig(out_path + 'Deconvolved_' + str(filter) + '_log_wire_3d.png')
#     plt.show()
#     return
#
# # JWST/MIRI
# #plot_3d_log(jwst_miri, str(filter)+' JWST/MIRI ', 1)
# #plot_3d_log_wire(jwst_miri, str(filter)+' JWST/MIRI ', 1)
#
# # deconvolved
# plot_3d_log(deconvolved, str(filter)+' Kraken deconvolved', 1)
# plot_3d_log_wire(deconvolved, str(filter)+' Kraken deconvolved', 1)
#
# # plot single observed/deconvolved images?
# image = True
# # plot 2x2 image grid?
# grid = True
#
# if image:
#     # read in the observed image data
#     img_file = get_pkg_data_filename(jwst_miri)
#     hdu = fits.open(img_file)
#     img= hdu[1].data
#
#     # read in the deconvolved image data
#     dec_file = get_pkg_data_filename(deconvolved)
#     hdu = fits.open(dec_file)
#     dec = hdu[1].data
#
#
#
#     # plot the observed image
#     fig, ax = plt.subplots(figsize=(6, 5), constrained_layout = True)
#     # obs centroid (y,x): 127,127
#     im1 = ax.imshow(img[27:227, 27:227], origin = 'lower',norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=scale), cmap='RdYlBu_r' )
#     ax.set_title(str(filter) + ' JWST/MIRI', fontsize = 15)
#     ax.hlines(85, 85, 115, color='k')
#     ax.hlines(115, 85, 115, color='k')
#     ax.vlines(85, 85, 115, color='k')
#     ax.vlines(115, 85, 115, color='k')
#     ax.set_xlabel('X Position (px)')
#     ax.set_ylabel('Y Position (px)')
#     # set the RA/DEC labels-x axis
#     ax.xaxis.set_major_locator(ticker.FixedLocator(positions))
#     ax.xaxis.set_major_formatter(ticker.FixedFormatter(labels))
#     ax.yaxis.set_major_locator(ticker.FixedLocator(positions))
#     ax.yaxis.set_major_formatter(ticker.FixedFormatter(labels))
#     ax.tick_params(direction='in', color='w', length=6, left=True, top=False, right=False, bottom=True, labelleft=True,labeltop=False, labelright=False, labelbottom=True)
#     add_colorbar(im1, format=FormatStrFormatter('%.e'))
#     plt.savefig(out_path + 'Observed_' + str(filter) + '.png')
#     plt.show()
#
#     # plot the deconvolved image
#     fig2, ax2 = plt.subplots(figsize=(6, 5), constrained_layout=True)
#     # deconvolved centroid (y,x): 129,128
#     im2 = ax2.imshow(dec[29:229, 28:228], origin='lower',norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=scale), cmap='RdYlBu_r')
#     ax2.set_title(str(filter) + ' Kraken deconvolved', fontsize = 15)
#     ax2.hlines(85, 85, 115, color='k')
#     ax2.hlines(115, 85, 115, color='k')
#     ax2.vlines(85, 85, 115, color='k')
#     ax2.vlines(115, 85, 115, color='k')
#     ax2.set_xlabel('X Position (px)')
#     ax2.set_ylabel('Y Position (px)')
#     # set the RA/DEC labels-x axis
#     ax2.xaxis.set_major_locator(ticker.FixedLocator(positions))
#     ax2.xaxis.set_major_formatter(ticker.FixedFormatter(labels))
#     ax2.yaxis.set_major_locator(ticker.FixedLocator(positions))
#     ax2.yaxis.set_major_formatter(ticker.FixedFormatter(labels))
#     ax2.tick_params(direction='in', color='w', length=6, left=True, top=False, right=False, bottom=True, labelleft=True,labeltop=False, labelright=False, labelbottom=True)
#     add_colorbar(im2, format=FormatStrFormatter('%.e'))
#     plt.savefig(out_path + 'Deconvolved_' + str(filter) + '.png')
#     plt.show()
#
# if grid:
#     # plot 2x2 image grid
#     # read in the .pngs
#     img = np.asarray(Image.open(out_path + 'Observed_'+str(filter)+'.png'))
#     dec = np.asarray(Image.open(out_path + 'Deconvolved_'+str(filter)+'.png'))
#     img_3d = np.asarray(Image.open(out_path + 'Observed_' + str(filter)+'_log_3d.png'))
#     dec_3d = np.asarray(Image.open(out_path + 'Deconvolved_' + str(filter)+'_log_3d.png'))
#
#     # crop the 3d-image arrays
#     img_3d = img_3d[70:550, 70:550]
#     dec_3d = dec_3d[70:550, 70:550]
#
#
#     fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(7,6), constrained_layout = True)
#     im1 = ax[0][0].imshow(img)
#     im2 = ax[0][1].imshow(dec)
#     im3 = ax[1][0].imshow(img_3d)
#     im4 = ax[1][1].imshow(dec_3d)
#
#     # hide the axis labels
#     ax[0][0].axis('off')
#     ax[0][1].axis('off')
#     #ax[1][0].axis('off')
#     #ax[1][1].axis('off')
#     ax[1][0].tick_params(direction='in', color='w', length=6, left=False, top=False, right=False, bottom=False, labelleft=False,labeltop=False, labelright=False, labelbottom=False)
#     ax[1][1].tick_params(direction='in', color='w', length=6, left=False, top=False, right=False, bottom=False, labelleft=False,labeltop=False, labelright=False, labelbottom=False)
#     plt.savefig(out_path + str(filter) + '_grid.png')
#     plt.show()
