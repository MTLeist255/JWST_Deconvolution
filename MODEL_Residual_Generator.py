# 2022 May 16: Code used to generate residual (background subtracted) MIRISim outputs
# NOTE: this is a shortcut for what the JWST pipeline should do

from Convenience_Functions import *

# Photometrically rescale the model image?
rescale = False

# set the simulation filter
filter = 'F560W'

# set the model name
idx = 3
model = ['Model_single_pixel',
        'Model_single_pixel_disk',
        'Model_single_pixel_disk_bicone',
        'Model_AGN_complicated'
        ]

# set the observed wavelength
if filter == 'F560W':
    obs_wave = 5.6
elif filter == 'F1000W':
    obs_wave = 10
elif filter == 'F1500W':
    obs_wave = 15
elif filter == 'F1800W':
    obs_wave = 18
elif filter == 'F2100W':
    obs_wave = 21

# pixel scale
pixel_scale = 0.1110

# define the import path for the image and background
im_path = 'Images/JWST/2_MIRISim_outputs/Deconvolution_testing/Raw_outputs/FITS/Ideal_Model/' + str(model[idx])
back_path = 'Images/JWST/2_MIRISim_outputs/Deconvolution_testing/Raw_outputs/background_FITS/'

# define the output path for the residual
out_dir = 'Images/JWST/2_MIRISim_outputs/Deconvolution_testing/Residual_outputs/FITS/Ideal_Model/'

# read in RAW MIRISim image
raw_image = get_pkg_data_filename(im_path + '/det_image_seq1_MIRIMAGE_'+str(filter)+'exp1.fits')
image_list = fits.open(raw_image)
image_data = image_list[1].data

# save the header info
hdr = image_list[0].header.copy()
# copy the SCI extension header information
hdr1 = image_list[1].header.copy()

# trim from 4x4 array to necessary 2x2 array with relevant data
# select integration number
image_data = image_data[0]
# select frame number
image = image_data[4]

# read in simulated background image, observed at same wavelength
sub = get_pkg_data_filename(back_path + '/det_image_seq1_MIRIMAGE_'+str(filter)+'exp1.fits')
sub_list = fits.open(sub)
# access the single frame data (0th integration - 5th Frame)
sub_data = sub_list[1].data
sub = sub_data[0]
sub1 = sub[4]

# Create residual image
# residual = image - background
residual = image - sub1

# mask negative pixel values and background pixel values
residual[residual<0] = 1

# Photometrically re-scale the image?
if rescale:
    # save the background subtracted image toa  dummy FITS file
    fits.writeto(im_path + '/Dummy_FITS.fits', residual, overwrite=True)
    dumb = get_pkg_data_filename(im_path + '/Dummy_FITS.fits')
    dumb_list = fits.open(dumb)
    # access the single frame data (0th integration - 5th Frame)
    dumb_data = dumb_list[0].data

    # measure aperture photometry -> used for scaling
    center_quad = centroid_quadratic(dumb_data, xpeak=128, ypeak=128, fit_boxsize=20)
    center = (center_quad[0], center_quad[1])

    # 8. measure photometry
    decon_flux = MAE_Aperture_Photometry(HDUlist=dumb,
                                         ext=0,
                                         radius=17.8,
                                         normalize=False,
                                         centroid=center,
                                         pad=0.5)

    print('Aperture photom. before scaling: ', decon_flux[0])

    # normalize  image
    residual = residual / residual.max()

    # manually tweak bad pixel values
    # residual[102:105,58:61]=0.0009

    # rescale the image
    scale = 8.68e3
    residual = residual * scale

    # save the rescaled image to dummy fits file, read in and measure photometry
    fits.writeto(im_path + '/Dummy_FITS2.fits', residual, overwrite=True)

# Update and save the FITS header info
hdr1['BUNIT'] = 'DN'
#hdr1['BUNIT'] = 'Jy'
hdr1[''] = 'and Background Subtracted info'
hdr1['PIXELSCL'] = (pixel_scale, 'pixel scale')
hdr1['WAVELEN'] = (obs_wave, 'Observed wavelength')
hdr1['FILTER'] = (filter, 'Observation filter')
hdr1['INSTRUM'] = ('MIRI', 'Observation instrument')
hdr1['BACKSUB'] = ('TRUE', 'single-frame MIRISim background subtracted')
hdr1['MDLHIST'] = ('PARAMS', 'Input models params were derived from the')
hdr1['MDLHIS2'] = ('_____', str(obs_wave)+ ' um observation of NGC 5728, Rosario et al.')

# set the individual model parameters in the FITS header
hdr1['PNTSRC'] = (7e6, 'Point source peak')
hdr1['DISK'] = ('disk', '2DGausian(3e4, 10, 10, 1, 1.5)')
hdr1['BICSIZ'] = (22, 'Ionization cone pixel size')
hdr1['BICSCL'] = ('250:1', 'Ionization cone scaling factor')
hdr1['GALAXY'] = (2e4, 'Normalized background galaxy factor')

if rescale:
    hdr1['SCLFCTR'] = (scale, 'Normalized image scaling factor')

# create new FITS header information for deconvolved image
new_hdul = fits.HDUList()
# append new header information -> single image
new_hdul.append(fits.ImageHDU(residual, header=hdr1))
new_hdul.append(fits.ImageHDU(header=hdr))

# save FITS file od optimum deconvolved image with correct FITS header information
new_hdul.writeto(out_dir + str(model[idx]) + '/SCALED_' + str(filter) + '_' + str(model[idx]) + '.fits', overwrite=True)

# close the image list
image_list.close()
sub_list.close()



