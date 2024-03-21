# 17 Jan 2024: code to read in observation data and:
# 1) saving CRPIX, NAXIS, TARG_RA/DEC from header
# 2) determine image centroid (update CRPIX as needed)
# 3) pad the observation data to a 1024x1024 array
# 4) generate DET-DIST webbPSF images (oversample = 4) to a 1024x1024 array
# 4b) generate PSF close to ODT
# 4c) generate PSF scaled + close to ODT

# NOTE: this code is necessary as both Richardson-Lucy and Kraken return the previous deconvolution iteration centered
# on the image array. If the centroid of the object is not at the center of the image array the resulting iteration
# will not be aligned with the PSF and the deconvolution fails

from Convenience_Functions import *

# generate a reference PSF?
generate_psf = False

# set the obs info
objname = 'NGC4388'
#define filter to generate
filterchoice = 4
filtername = ['F560W', 'F1000W', 'F1500W', 'F1800W', 'F2100W']

# set input/output paths
in_path = 'Images/JWST/5_ERS_Data/5_Deconvolution_results/'+objname+'/'
im_out_path = in_path + 'Deconvolution/Padded_images/'
imagefile = 'jw02064-o009_t005_miri_'+filtername[filterchoice]+'-sub256_astromtweaked_i2d.fits'

# read in observational data
hdulist = fits.open(in_path+imagefile)
science_image = hdulist[1].data
# save header info: TARG_RA, TARG_DEC, NAXIS1, NAXIS2
sci_hdr       = hdulist[0].header
wcs_header    = hdulist[1].header

# set the WCS coordinates
imwcs = WCS(wcs_header)
# save image array shape
#x,y = wcs_header['NAXIS1'], wcs_header['NAXIS2']
y,x = 180,180
# set the nuclear position from TARG_RA, TARG_DEC
nuclear_coords = (sci_hdr['TARG_RA'], sci_hdr['TARG_DEC'])
# update the wcs header with the RA-, DEC info
wcs_header['TARG_RA'] = (sci_hdr['TARG_RA'], 'Target RA')
wcs_header['TARG_DEC'] = (sci_hdr['TARG_DEC'], 'Target DEC')
# set the image pixel scale (assuming square pixels)
pix_scale = np.sqrt(imwcs.pixel_scale_matrix[0][0]**2 + imwcs.pixel_scale_matrix[0][1]**2)*3600.0
# determine the image centroid. NOTE: centroid is assumed to be image peak as AGN acts as an essential point source
peaky, peakx = np.unravel_index(np.argmax(science_image,axis=None), science_image.shape)

# scale down image size and remeasure peak
science_image = science_image[peaky-90:peaky+90, peakx-90:peakx+90]
peaky2, peakx2 = np.unravel_index(np.argmax(science_image,axis=None), science_image.shape)



# update with measured centroid- used to align PSF to whole pixel
# 17 Jan 2024 NOTE: center is off by 1-pixel
wcs_header['NEW_CP1'], wcs_header['NEW_CP2'] = peaky2, peakx2
print(peaky2, peakx2)
wcs_header['CRPNOTE'] = ('CRPIX-Update', 'pix position used to align PSF + pad img array')

# pad image array to a 1024x1024 grid
naxis = (1024,1024)
pad_im_arr = np.zeros(naxis)

# determine offsets
offy = 511-peaky2
offx = 511-peakx2
print(offy, offy+y)
# pad image array from NAXIS1, NAXIS2 -> 1024,1024
# NOTE: this offsets the obs y,x centroid FIRST from the array y,x centroid then adds the offset+(y,x) array length
pad_im_arr[offy:offy+y, offx:offx+x] = science_image
# mask negative pixels
pad_im_arr[pad_im_arr<0]=0

# save the image array size
newx, newy = pad_im_arr.shape
print('OBS image array shape (x,y): ', pad_im_arr.shape)

# update the header with the padded image array size
wcs_header['NAXIS1'] = (newx, 'new x-axis size')
wcs_header['NAXIS2'] = (newy, 'new y-axis size')

#base_filename = objname + '_' + filtername[filterchoice] + '_1024.fits'
base_filename = 'test.fits'
outfile = os.path.join(im_out_path, base_filename)
hdu = fits.HDUList()

# save the FITS file with header info
# pad the 1024,1024 array to idx 0 with the updated wcs header info
hdu.append(fits.PrimaryHDU(pad_im_arr, header=wcs_header))
hdu.writeto(outfile, overwrite=True)
hdulist.flush()

if generate_psf:
    psf_out_path = in_path + 'Deconvolution/PSFs/standard_PSF/'

    # generate the MIRI DET_DIST PSF, make sure the centroid is at 512,512 on a 1024,1024 array
    # Create a MIRI instance and calculate a PSF from the given filter
    miri = webbpsf.MIRI()
    miri.filter = filtername[filterchoice]

    # 15 Jan 2024: I should be able to generate PSFs close to the observation OTD but for some reason I can't? Come back and
    # fix this
    # miri.load_wss_opd_by_date('2023-07-07T13:33:30.016', plot=True)

    # calculate PSF
    out_name = 'WebbPSF_' + filtername[filterchoice] + '.fits'
    # NOTE: to keep DET_DIST pixel scales correct build the array sampled to the MIRI SUB256 (i.e. 256x256) then ADD
    # the image array to a larger (1024x1024) array of zeros
    psf = miri.calc_psf(fov_pixels=1024, oversample=4, display=False, outfile=psf_out_path + out_name)

    # read in the DET_DIST PSF and save to its own FITS file
    hdulist = fits.open(psf_out_path + out_name)
    # read in DET_DIST index PSF
    sci_image = hdulist[3].data
    # save PSF header info
    psf_sci_hdr = hdulist[0].header
    print('PSF image array shape (x,y): ', sci_image.shape)

    # save PSF
    #fits.writeto(psf_out_path + filtername[filterchoice] + '_MIRI_DET_DIST.fits', sci_image, psf_sci_hdr,overwrite=True)

