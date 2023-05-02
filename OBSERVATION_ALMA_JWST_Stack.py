# 7 Dec 2022: practice stacking ALMA NGC 5728 CO 2-1 image and JWST 10 um image
# try observed and deconvolved?

from Convenience_Functions import *

im_path = 'Images/JWST/5_ERS_Data/1_Raw_Data/ngc5728_MIRI_Imaging/'
RL_path = 'Images/JWST/4_Deconvolution/1_Richardson_Lucy/NGC5728/FITS/'
jy_path = 'Images/JWST/2_MIRISim_outputs/Deconvolution_testing/Residual_outputs/FITS/JWST_Observations/'
out_path = 'Images/JWST/5_ERS_Data/1_Raw_Data/ngc5728_MIRI_Imaging/'
K_path = 'Images/JWST/4_Deconvolution/5_Kraken/NGC5728/FITS/'

# read in JWST/MIRI data
blue = get_pkg_data_filename(im_path + 'NGC5728_F1000W_final.fits')
blue_list = fits.open(blue)
blue_im = blue_list[1].data

# read in kraken data
green = get_pkg_data_filename(K_path + 'NGC5728_F1000W_final_DET_SAMP.fits')
green_list = fits.open(green)
green_im = green_list[12].data

# overwrite the kraken FITS header info
prime_hdr = blue_list[0].header.copy()
sci_hdr = blue_list[1].header.copy()
hdu = fits.HDUList()
hdu.append(fits.PrimaryHDU(header=prime_hdr))
# 8b. Append the science header info and MIRISim image data to the extension 1
hdu.append(fits.ImageHDU(green_im, header=sci_hdr))
fits.writeto('wcs_test.fits', green_im, sci_hdr, overwrite=True)
green_list.flush()

green = get_pkg_data_filename('wcs_test.fits')
green_list = fits.open(green)
green_im = green_list[0].data

# read in reference ALMA image
red = get_pkg_data_filename('ngc5728_alma_full_pyspeckit_fit_single_gaussian_28-08-2018_physical_params.fits')
red_list = fits.open(red)
red_im = red_list[0].data

# reproject the JWST to the ALMA image orientation?
# OR Orientate the deconvolved image to the observed image orientation
blue_img = reproject_interp(green_list[0], red_list[0].header)
blue_img = blue_img[0]
fits.writeto('wcs_test.fits', blue_img, red_list[0].header, overwrite=True)

plt.imshow(blue_img, origin = 'lower', norm =ImageNormalize(stretch=LogStretch(), vmin=0, vmax = blue_im.max()))
plt.title('Reprojected coordinates?')
plt.show()

# determine centroids of each image and align accordingly
# trim ALMA image: center (x,y) = (254, 250)
red_im = red_im[185:315, 189:319]

# resize blue imge: center (x,y) = (256, 251)
# resize blue_img (kraken) center (x,y) = (385, 181)
blue_img = blue_img[116:246, 320:450]
fits.writeto('wcs_test.fits', blue_img, red_list[0].header, overwrite=True)

# try different method to plot RGB
img = np.zeros((red_im.shape[0], red_im.shape[1],3), dtype=float)
img[:,:,0] = img_scale.asinh(red_im, scale_min=0, scale_max = red_im.max())
img[:,:,1] = img_scale.asinh(blue_img, scale_min=0, scale_max = green_im.max())

fig, axes = plt.subplots(1, figsize=(7, 5), tight_layout=True)
axes.set_aspect('equal')
x1, y1 = blue_img.shape

# plot the combined image
plt.imshow(img, origin = 'lower', extent = (int(0), int(x1), int(0), int(y1)))
positions = [1, 32.5, 65, 97.5, 125]
labels = ['-7', '-3.5', '0', '3.5', '7']
axes.xaxis.set_major_locator(ticker.FixedLocator(positions))
axes.xaxis.set_major_formatter(ticker.FixedFormatter(labels))
axes.yaxis.set_major_locator(ticker.FixedLocator(positions))
axes.yaxis.set_major_formatter(ticker.FixedFormatter(labels))
mono = {'family' : 'monospace'}
# , fontdict = mono
plt.text(70,120,'ALMA/CO (2-1)', color = 'r', fontsize = 15)
plt.text(2,120,'JWST/MIRI 10 μm', color = 'lime', fontsize = 15)
plt.text(2,110,'Deconvolved', color = 'lime', fontsize = 15, fontdict = mono)
plt.xlabel('RA-Offset (")')
plt.ylabel('DEC-Offest (")')

# set compass params
size = [115.5, 18.5, 101.5, 4, 117, 5, 0, 10, 117, 5, -10, 0]
#Set orientation compass: x,y
# add text
plt.text(size[0], size[1], 'N', color='w')
plt.text(size[2], size[3], 'E', color='w')

# create the height of the arrow
arrow1 = plt.arrow(size[4], size[5], size[6], size[7], width=0.3, color='w')
# create the length of the arrow
arrow2 = plt.arrow(size[8], size[9], size[10], size[11], width=0.3, color='w')

# add arrows
axes.add_patch(arrow1)
axes.add_patch(arrow2)
plt.savefig('Stacked_Kraken_ALMA.png')
plt.show()