# 17 Jan 2024: source code used to generate reference PSFs for deconvolution
# NOTE: source code modified from David Rosario (https://github.com/vikalibrate/GATOS_Deconvolution/blob/main/deconvolution_prep.py)

from Convenience_Functions import *

# import reference image
objname = 'NGC4388'
#define filter to generate
filterchoice = 0
filtername = ['F560W', 'F1000W', 'F1500W', 'F1800W', 'F2100W']
# radii of 80%EE from MIRI JDox (July 2022)
ee80radii   = [0.422, 0.636, 0.932, 1.110, 1.276]

# set input/output paths
in_path = 'Images/JWST/5_ERS_Data/5_Deconvolution_results/NGC4388/'
out_path = in_path + 'Deconvolution/PSFs/'
imagefile = 'jw02064-o009_t005_miri_'+filtername[filterchoice]+'-sub256_astromtweaked_i2d.fits'

# read in image data
hdulist = fits.open(in_path+imagefile)
science_image = hdulist[1].data
# save header info: TARG_RA, TARG_DEC, NAXIS1, NAXIS2
sci_hdr       = hdulist[0].header
wcs_header    = hdulist[1].header

# set the WCS coordinates
imwcs = WCS(wcs_header)
if filterchoice >1:
    # pad CRPIX1 +10 to account for offset from actual centroid and centroid given in FITS header
    ceny, cenx = wcs_header['CRPIX1']+8, wcs_header['CRPIX2']
else:
    # save image centroid from WCS coordinates
    ceny, cenx = wcs_header['CRPIX1'], wcs_header['CRPIX2']

# save image array shape
x,y = wcs_header['NAXIS1'], wcs_header['NAXIS2']
# set the nuclear position from TARG_RA, TARG_DEC
nuclear_coords = (sci_hdr['TARG_RA'], sci_hdr['TARG_DEC'])
graft_coords = (sci_hdr['TARG_RA'], sci_hdr['TARG_DEC'])
# set the image pixel scale (assuming square pixels)
pix_scale = np.sqrt(imwcs.pixel_scale_matrix[0][0]**2 + imwcs.pixel_scale_matrix[0][1]**2)*3600.0

# initialise the PSF object
inst = webbpsf.MIRI()
inst.aperturename = 'MIRIM_SUB256'
inst.filter = filtername[filterchoice]
# WebbPSF scale is same as native MIRI imager scale to 0.3% accuracy, but force for full accuracy
inst.pixelscale = pix_scale

# The i2d image is larger than the original detector size, so the position of the nucleus will
# different here than on the detector. However, assuming the center corresponds to approximately
# the centre of the detector, the offset from the centre should still be good enough for PSF
# construction.
imrectsize=min(science_image.shape)
nuclear_pixel = imwcs.all_world2pix([nuclear_coords],0)

xcentoff = int(nuclear_pixel[0][0]+0.5-imrectsize/2.0)
ycentoff = int(nuclear_pixel[0][1]+0.5-imrectsize/2.0)
inst.options['source_offset_r'] = np.sqrt(xcentoff**2 + ycentoff**2)*pix_scale
inst.options['source_offset_theta'] = np.rad2deg(np.arctan2(ycentoff,xcentoff)) - 90.0

# Pad the image area to ensure that we can get a full PSF even if the nucleus is close to the edge
psfrectsize = int(np.ceil(10*ee80radii[filterchoice]/pix_scale)) # final size of calibration PSF image

# Ensure it is odd, so PSF is centered on a pixel
if (psfrectsize % 2) == 0:
    psfrectsize += 1

outname1 = 'WebbPSF_'+filtername[filterchoice] +'.fits'
modelpsfhdu = inst.calc_psf(oversample=4,fov_pixels=imrectsize+psfrectsize, outfile= out_path + outname1)

# read in and resize the DET_DIST PSF
hdulist = fits.open(out_path+outname1)
# read in DET_DIST index PSF
sci_image = hdulist[3].data
# save PSF header info
psf_sci_hdr       = hdulist[0].header

# Resize the image array to match the observed image array size and centroid positions
new_shape = (y,x)
resized_image = zoom(sci_image, (new_shape[0] / sci_image.shape[0], new_shape[1] / sci_image.shape[1]))

# Calculate the shift for centering: x
shift_x = round(cenx - new_shape[0] // 2)
# offset by 1/2 a pixel to account for shift rounding: y
shift_y = round((ceny - new_shape[1] // 2)-1)

# Shift the resized image array
shifted_image = np.roll(resized_image, (shift_y, shift_x), axis=(0, 1))

# new image array shape (y,x)
newy, newx = shifted_image.shape

# update NAXIS, CRPIX and save to FITS file
psf_sci_hdr['NAXIS1'] = (newx, 'x-axis length')
psf_sci_hdr['NAXIS2'] = (newy, 'y-axis length')
psf_sci_hdr['CRPIX1'] = (ceny, 'PSF y-centroid determined by '+objname+' '+filtername[filterchoice]+' obs')

if filterchoice>1:
    # annotate CRPIX1 manual correction in header
    psf_sci_hdr['NOTES-2'] = ('CRPIX-OFF', 'CRPIX1 value from obs header off by ~8-pixs')
    psf_sci_hdr['NOTES-3'] = ('CRPIX-OFF', 'CRPIX1 manually offset ~8-pixs to correct')

psf_sci_hdr['CRPIX2'] = (cenx, 'PSF x-centroid determined by '+objname+' '+filtername[filterchoice]+' obs')
psf_sci_hdr['NOTES-1'] = ('CRPIX', 'CRPIX values determined to nearest whole pixel')



fits.writeto(out_path + filtername[filterchoice] + '_MIRI_DET_DIST_1024.fits', shifted_image, psf_sci_hdr, overwrite=True)

# close the importer
hdulist.close()