# 2022 May 13: Generate reference PSFs using WebbPSF

from Convenience_Functions import *

# format the output PSFs for Kraken deconvolution?
Kraken = False

# set the output directory and name
#out_dir = 'Images/JWST/5_ERS_Data/4_STAGE3_outputs/NGC4388/MIRI_Imaging/Deconvolution_Leist/Reference_PSFs/'
out_dir = 'Images/JWST/5_ERS_Data/5_Deconvolution_results/NGC4388/Deconvolution/PSFs/'

# Create a MIRI instance and calculate a PSF from the given filter
filter = 'F560W'
name = 'WebbPSF_'+str(filter)+'.fits'
miri = webbpsf.MIRI()
miri.filter = filter


# 15 Jan 2024: I should be able to generate PSFs close to the observation OTD but for some reason I can't? Come back and
# fix this
#miri.load_wss_opd_by_date('2023-07-07T13:33:30.016', plot=True)

# calculate PSF
psf = miri.calc_psf(fov_pixels =256, oversample=4, display = False, outfile= out_dir + name)

# read in reference PSF, save the detector-sampled PSF for decovolution
miri_psf = get_pkg_data_filename(out_dir + name)
psf_list = fits.open(miri_psf)
# save the FITS header info
prime_hdr = psf_list[0].header.copy()
sci_hdr = psf_list[1].header.copy()

# check the PSF imiage array sizes
oversamp_psf = psf_list[0].data
print(oversamp_psf.shape)
detsamp_psf = psf_list[1].data
print(detsamp_psf.shape)
overdist_psf = psf_list[2].data
print(overdist_psf.shape)
detdist_psf = psf_list[3].data
print(detdist_psf.shape)

# reshape the PSF to 294x295
#detdist_psf = detdist_psf[0:295, 1:295]
# save the formatted image to a FITS file
base_filename = str(filter) + '_PSF.fits'
outfile = os.path.join(out_dir, base_filename)
hdu = fits.HDUList()

#sci_hdr['NAXIS'] = (294, 'New image array x-axis length')
#sci_hdr['NAXIS'] = (295, 'New image array y-axis length')

# append the simulated image
# append new header information -> single image
hdu.append(fits.PrimaryHDU(detdist_psf, header=prime_hdr))
hdu.append(fits.ImageHDU(header=sci_hdr))
hdu.writeto(outfile, overwrite=True)
psf_list.flush()

# Display the PSF
fig, ax = plt.subplots(nrows=1, ncols=4, figsize=(10,3), constrained_layout = True)
ax[0].imshow(oversamp_psf, origin = 'lower', norm = ImageNormalize(stretch=LogStretch(), vmin=0, vmax=oversamp_psf.max()), cmap = 'RdYlBu_r')
ax[0].set_title('OVERSAMP (idx 0)')
ax[1].imshow(detsamp_psf, origin = 'lower', norm = ImageNormalize(stretch=LogStretch(), vmin=0, vmax=detsamp_psf.max()), cmap = 'RdYlBu_r')
ax[1].set_title('DET_SAMP (idx 1)')
ax[2].imshow(overdist_psf, origin = 'lower', norm = ImageNormalize(stretch=LogStretch(), vmin=0, vmax=overdist_psf.max()), cmap = 'RdYlBu_r')
ax[2].set_title('OVERDIST (idx 2)')
ax[3].imshow(detdist_psf, origin = 'lower', norm = ImageNormalize(stretch=LogStretch(), vmin=0, vmax=detdist_psf.max()), cmap = 'RdYlBu_r')
ax[3].set_title('DET_DIST (idx 3)')
fig.suptitle(filter, fontsize=25)
#plt.savefig(out_dir + 'WebbPSF_PSF_grid_'+filter+'.png')
plt.show()

# read in the generated PSFs and format to a 1024x1024 array for Kraken deconvolution
if Kraken:
    # read in the MIRI data
    image = get_pkg_data_filename(out_dir + name)
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
    base_filename = 'Kraken_' + name
    outfile = os.path.join(out_dir, base_filename)
    hdu = fits.HDUList()

    sci_hdr['NAXIS'] = (256, 'New image array x-axis length')
    sci_hdr['NAXIS'] = (256, 'New image array y-axis length')

    # append the simulated image
    # append new header information -> single image
    hdu.append(fits.PrimaryHDU(pad_data, header=prime_hdr))
    hdu.append(fits.ImageHDU(header=sci_hdr))
    hdu.writeto(outfile, overwrite=True)
    image_list.flush()


psf_list.close()