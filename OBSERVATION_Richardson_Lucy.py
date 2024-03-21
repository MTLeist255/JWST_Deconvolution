# 16 Jan 2024: simple code to test how well the function 'measure_merits' works

from Convenience_Functions import *

warnings.filterwarnings('ignore')

# deconvolve the image data
deconvolve = True
# measure the merit functions after deconvolution
# NOTES: (1) can be run independently from the above if statement, (2) must be set to TRUE for the remaining if statements to work
measure_merit_functions = False
# rotate the final observed and deconvolved images to N-E orientation and save
rotate = False
# save the merit functions to a .txt file and corresponding plot?
save_merits = False
# plot the final merit functions?
plot_merits = False

# set the total number of iterations
dec_iter = 100
skip=0

# set the object name
obj_name = ['NGC4388']
# set the reference PSF used for deconvolution
psf = ['observed_psf']
# set the observation filter
filter = ['F560W', 'F2100W']

# iterate through each object
for i in range(0, int(len(obj_name))):
    # keep track of progress
    print('Object: ', obj_name[i])
    # iterate through each reference PSF
    for j in range(0, int(len(psf))):
        # keep track of progress
        print('Reference PSF: ', psf[j])
        # iterate through each observed filter
        for k in range(0, int(len(filter))):
            # keep track of progress
            print('Filter: ', filter[k])

            # set the input/output paths
            global_path = 'Images/JWST/5_ERS_Data/5_Deconvolution_results/' + obj_name[i] + '/Deconvolution/'
            im_path = global_path + 'Padded_images/'
            psf_path = global_path + 'PSFs/' + psf[j] + '/'
            out_path = global_path + 'Results/Richardson_Lucy/' + psf[j] + '/'

            # read in the MIRI data
            miri_name = obj_name[i] + '_' + filter[k] + '_1024.fits'
            image = get_pkg_data_filename(im_path + miri_name)
            image_list = fits.open(image)
            # save the image data
            image_data = image_list[0].data
            # save the header info
            sci_hdr = image_list[0].header

            # read in the PSF data
            # set the PSF name
            if j == 0:
                # standard PSF
                #psf_name = filter[k] + '_MIRI_DET_DIST.fits'
                psf_name = obj_name[i] + '_' + filter[k] + '_norm_observedpsf_1024.fits'
            elif j == 1:
                # scaled PSF
                psf_name = obj_name[i] + '_' + filter[k] + '_scaledpsfmodel_1024.fits'
            else:
                # observed PSF
                psf_name = obj_name[i] + '_' + filter[k] + '_observedpsf_1024.fits'

            # get the PSF data
            psf_img = get_pkg_data_filename(psf_path + psf_name)
            psf_list = fits.open(psf_img)
            psf_data = psf_list[0].data

            # Richardson-Lucy deconvolve
            image_arr = []
            for l in range(0, dec_iter):
                if l == 0:
                    # set the 0th iteration as the original image
                    dec_data = image_data
                else:
                    # RL deconvolve and measure image stats
                    dec_data = richardson_lucy_np(image_data, psf_data, l)

                # append the image data
                image_arr.append(dec_data)

            # save all of the deconvolved images to a single FITS cube
            im_arr = np.array(image_arr)
            base_filename = obj_name[i] + '_' + filter[k] + '_DECONV.fits'
            outfile = os.path.join(out_path, base_filename)

            # update the FITS header
            sci_hdr['METHOD'] = ('Richardson-Lucy', 'Deconvolution method used')
            sci_hdr['REFPSF'] = (psf[j], 'Reference PSF used during deconvolution')
            sci_hdr['MAXITER'] = (dec_iter, 'Maximum number of deconvolution iterations')

            # save the FITS data
            hdu = fits.PrimaryHDU(im_arr, header=sci_hdr)
            hdu1 = fits.HDUList([hdu])
            hdu.writeto(outfile, overwrite=True)






