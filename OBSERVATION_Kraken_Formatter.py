# 8 mar 2023: format the JWST/MIRI observations to the Kraken input format

from Convenience_Functions import *

resize_256 = False
resize_1024 = True

# set the input and output path
obj_name = 'NGC3227'
im_path = 'Images/JWST/5_ERS_Data/5_Deconvolution_results/'+str(obj_name)+'/'
out_path = im_path + 'Deconvolution/Padded_Images/'

# set the observation filter
filter = '1800'
name = 'jw02064-o007_t004_miri_f'+str(filter)+'w-sub256_astromtweaked_i2d'

# read in the MIRI data
image = get_pkg_data_filename(im_path + name + '.fits')
image_list = fits.open(image)
image_data = image_list[1].data

# save the FITS header info
prime_hdr = image_list[0].header.copy()
sci_hdr = image_list[1].header.copy()

# set the centroid (y,x) positions from the header info
ceny, cenx = sci_hdr['CRPIX1'], sci_hdr['CRPIX2']
# set the shape of the image array
x,y = sci_hdr['NAXIS1'], sci_hdr['NAXIS2']

if resize_1024:
    # Find image peak and recenter image to center of 1024x024 array based on this
    peaky, peakx = np.unravel_index(np.argmax(image_data, axis=None), image_data.shape)
    # update with measured centroid- used to align PSF to whole pixel
    # 17 Jan 2024 NOTE: center is off by 1-pixel
    sci_hdr['NEW_CP1'], sci_hdr['NEW_CP2'] = peaky, peakx
    print('Centroid from header: ', ceny, cenx)
    print('Centroid from image peak: ', peaky, peakx)
    sci_hdr['CRPNOTE'] = ('CRPIX-Update', 'peak pix position used to align PSF + pad img array')

    # pad image array to a 1024x1024 grid
    naxis = (1024, 1024)
    pad_im_arr = np.zeros(naxis)

    # determine offsets
    offy = 511 - peaky
    offx = 511 - peakx
    # pad image array from NAXIS1, NAXIS2 -> 1024,1024
    # NOTE: this offsets the obs y,x centroid FIRST from the array y,x centroid then adds the offset+(y,x) array length
    pad_im_arr[offy:offy + y, offx:offx + x] = image_data

    # mask negative pixels in the image array
    pad_im_arr[pad_im_arr<0] = 0

    # save the formatted image to a FITS file
    base_filename = str(obj_name) + '_F' + str(filter) + 'W_1024.fits'
    outfile = os.path.join(out_path, base_filename)
    hdu = fits.HDUList()

    sci_hdr['NAXIS'] = (1024, 'New image array x-axis length')
    sci_hdr['NAXIS'] = (1024, 'New image array y-axis length')

    # append the simulated image
    # append resized image + science header
    hdu.append(fits.PrimaryHDU(pad_im_arr, header=sci_hdr))
    # append the original image + primary header
    hdu.append(fits.ImageHDU(image_data, header=prime_hdr))
    hdu.writeto(outfile, overwrite=True)
    image_list.flush()


if resize_256:
    # crop the image array to 256x256
    print('Image shape before trim: ', image_data.shape)
    # reshape image array: (y,x)
    image_data = image_data[round(ceny) - 128:round(ceny) + 128, round(cenx) - 128:round(cenx) + 128]
    print('Image shape after trim: ', image_data.shape)

    # pad negative pixel values = 1
    image_data[image_data < 0] += 0
    pad_data = image_data

    # save the formatted image to a FITS file
    base_filename = 'PAD_' + str(obj_name) + '_F' + str(filter) + 'W.fits'
    outfile = os.path.join(out_path, base_filename)
    hdu = fits.HDUList()

    sci_hdr['NAXIS'] = (256, 'New image array x-axis length')
    sci_hdr['NAXIS'] = (256, 'New image array y-axis length')

    # append the simulated image
    # append new header information -> single image
    # append new header information -> single image
    hdu.append(fits.PrimaryHDU(pad_data, header=prime_hdr))
    hdu.append(fits.ImageHDU(header=sci_hdr))
    hdu.writeto(outfile, overwrite=True)
    image_list.flush()
