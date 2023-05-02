# 2022 May 13: Generate reference PSFs using WebbPSF
import webbpsf

from Convenience_Functions import *

# format the output PSFs for Kraken deconvolution?
Kraken = False

# Create a MIRI instance and calculate a PSF from the given filter
filter = 'F2100W'
miri = webbpsf.MIRI()
miri.filter = filter
out_dir = 'Images/JWST/5_ERS_Data/4_STAGE3_outputs/NGC5728_obs1/RL_Formatted/'
name = 'RL_reference_PSF_'+str(filter)+'_DET_DIST_FULL_original2.fits'
psf = miri.calc_psf(fov_pixels =130, oversample=4, display = False, outfile= out_dir + name)

# read in reference PSF
miri_psf = get_pkg_data_filename(out_dir + name)
psf_list = fits.open(miri_psf)
psf2 = psf_list[0].data
print(psf2.shape)
psf3 = psf_list[1].data
print(psf3.shape)
psf4 = psf_list[2].data
print(psf4.shape)
psf5 = psf_list[3].data
print(psf5.shape)
psf_list.close()

# Display the PSF
fig, ax = plt.subplots(nrows=1, ncols=4, figsize=(10,8), tight_layout = True)
ax[0].imshow(psf2, origin = 'lower', norm = ImageNormalize(stretch=LogStretch(), vmin=0, vmax=psf2.max()), cmap = 'turbo')
ax[0].set_title('index 0')
ax[1].imshow(psf3, origin = 'lower', norm = ImageNormalize(stretch=LogStretch(), vmin=0, vmax=psf3.max()), cmap = 'turbo')
ax[1].set_title('index 1')
ax[2].imshow(psf4, origin = 'lower', norm = ImageNormalize(stretch=LogStretch(), vmin=0, vmax=psf4.max()), cmap = 'turbo')
ax[2].set_title('index 2')
ax[3].imshow(psf5, origin = 'lower', norm = ImageNormalize(stretch=LogStretch(), vmin=0, vmax=psf5.max()), cmap = 'turbo')
ax[3].set_title('index 3')
plt.show()

# read in the generated PSFs and format to a 1024x1024 array for Kraken deconvolution
if Kraken:
    # read in the MIRI data
    image = get_pkg_data_filename(out_dir + 'Kraken_reference_PSF_'+str(filter)+'_DET_DIST_FULL.fits')
    image_list = fits.open(image)
    image_data = image_list[3].data

    # save the FITS header info
    prime_hdr = image_list[0].header.copy()
    sci_hdr = image_list[1].header.copy()

    # pad the image to a 1024x1024 image array
    print('Image shape before pad: ', image_data.shape)
    naxis = (1024, 1024)
    pad_data = np.zeros(naxis)
    pad_data[384:640, 384:640] += image_data
    print('Image shape after pad: ', pad_data.shape)

    # save the formatted image to a FITS file
    base_filename = 'Kraken_reference_PSF_'+str(filter)+'_DET_DIST_obs1.fits'
    outfile = os.path.join(out_dir, base_filename)
    hdu = fits.HDUList()

    sci_hdr['NAXIS'] = (1024, 'New image array x-axis length')
    sci_hdr['NAXIS'] = (1024, 'New image array y-axis length')

    # append the simulated image
    # append new header information -> single image
    hdu.append(fits.PrimaryHDU(pad_data, header=prime_hdr))
    hdu.append(fits.ImageHDU(header=sci_hdr))
    hdu.writeto(outfile, overwrite=True)
    image_list.flush()
