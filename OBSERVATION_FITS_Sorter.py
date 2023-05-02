# 26 April 2023: Sort through the deconvolved FITS files and save the optimum deconvolved image to a single FITS file,
# based on our merit function convergance criteria

from Convenience_Functions import *

# set the optimum deconvolved iteration to save
# NOTE: the 0th-array contains the primary header, +1 to the index for the correct image
save = 66+1

# set the Directory
Dir = '5_Kraken'

# set the filter used
filter = 'F2100W'

# set the FITS import/output path
path = 'Images/JWST/4_Deconvolution/'+str(Dir)+'/NGC5728/obs2/FITS/'
out_path = 'Images/JWST/5_ERS_Data/4_Stage3_outputs/NGC5728_obs2/Single_DECONV_FITS/'

# set the observation name
name = 'RL_NGC5728_'+str(filter)+'_obs2_DET_SAMP_XY.fits'

# read in the deconvolution data array
decon = get_pkg_data_filename(path + name)
# get the necessary image data
decon_list = fits.open(decon)
# save the optimum deconvolved image

decon_data = decon_list[save].data

# copy the primary and science header
prime_hdr = decon_list[0].header.copy()
sci_hdr = decon_list[1].header.copy()

# save the FITS file and header info to a seperate file
base_filename = 'NGC5728_'+str(filter)+'_DECONV_RL_Leist.fits'
outfile = os.path.join(out_path, base_filename)
hdu = fits.HDUList()

# append the primary header
hdu.append(fits.PrimaryHDU(header=prime_hdr))
# append the original image and science header
# set the deconvolution method
method = 'Kraken'
# 14a. Append the primary header info to extension 1
sci_hdr['BUNIT'] = 'mJy/sr'
sci_hdr[''] = ' and Kraken MFBD Deconvolution Parameters'
sci_hdr['SAMP'] = ('DET_SAMP', 'PSF extension used for deconvolution')
sci_hdr['METHOD'] = ('Kraken MFBD', 'Deconvolution method')
sci_hdr['ITERATS'] = (str(save), 'Optimum deconvolution number')
hdu.append(fits.ImageHDU(decon_data, header=sci_hdr))
hdu.writeto(outfile, overwrite=True)
