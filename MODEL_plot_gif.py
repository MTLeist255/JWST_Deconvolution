# 2022 July 9: code to arbitrarily save deconvolution images to create gifs

# import convenience functions

from Convenience_Functions import *

# set necessary file input/output paths
# set the model name
model = 'residual_MIRISim_single_pixel_5um'

# set the input path for the deconvolved image and merit function measurements
im_path = 'Images/JWST/4_Deconvolution/1_Richardson_Lucy/FITS/'

# set the output path for the deconvolved image and merit function measurements
mes_path = 'Images/JWST/4_Deconvolution/1_Richardson_Lucy/Measurements/'

# set the output path for all the images
out_path = 'Images/JWST/4_Deconvolution/1_Richardson_Lucy/Images/GIF/' + model + '/'

# set the aperture radial size
# set the optimum deconvolved number (determine din P_Merit_Functions.py)
optimum = 17.8

#*********************************************************************************************************************#
# Read in the observed image
image = get_pkg_data_filename(im_path + 'deconvolution_RL_' + model + '.fits')

# read in image data for plotting
image_list = fits.open(image)

# read in the merit function data
o_nint, o_flux, o_fwhm, o_strehl = np.loadtxt(mes_path + model + '/RL_Image_Plot_' + str(model) + '.txt', unpack = True)

# display/save each output in the range 0 -> optimum
total_im = []
for i in range(0,int(optimum+1)):
    image_data = image_list[i].data
    total_im.append(image_data)

    show_Image(image_data, title='Richardson-Lucy',
               normalize=True,
               text='Iteration: ' + str(i)
               + '\nFWHM: ' + str(o_fwhm[i]),
               cmap='viridis',
               print_params=False,
               display=False,
               save=out_path + '/GIF_n_' + str(i)
               )
