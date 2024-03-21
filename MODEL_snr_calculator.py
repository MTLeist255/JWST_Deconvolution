# 25 Aug 2023: simple code to calculate the cummulative SNR of a group of images and plot the results

from Convenience_Functions import *

# set the aperture radisu for photom
aper_rad =87.4

# plot the V1 and V2 image grids
dir1 = '2_MIRISim_outputs/Deconvolution_testing/Residual_outputs/'
dir2 = '4_Deconvolution/5_Kraken'
sub_dir = '26sep2023_fc_noscale_unmasked/'
dir = [dir1, dir2]

# set the output
out_path = 'Images/JWST/4_Deconvolution/Total_Comparisons/tests/'
out_name = 'Fig_V6.png'
mods = [0,1,2,3,4,5,6,7,8,9,10,11]

sim = []
k = []
k_iter = []
rl = []
rl_iter = []
for i in range(0, int(len(dir))):
    for j in range(0, int(len(mods))):
        if i > 0:
            im_path = 'Images/JWST/' + str(dir[i]) + '/FITS/Ideal_Model/Model_AGN_complicated/current_12Aug2023/V6/final_deconv/' + str(sub_dir)
            name2 = 'F2100W_V6_DECONV_'+str(mods[j])+'.fits'
        else:
            im_path = 'Images/JWST/' + str(dir[i]) + '/FITS/Ideal_Model/Model_AGN_complicated/current_12Aug2023/V6/' + str(sub_dir)
            name2 = 'F2100W_V6_model_'+str(mods[j])+'.fits'

        # import data
        miri = get_pkg_data_filename(im_path + name2)
        miri_list = fits.open(miri)
        miri_data = miri_list[0].data

        # measure the SNR
        snr = MAE_SNR_calculator(miri_data, extract=128, box_size=10, centroid=(128, 128), aper_rad=aper_rad)
        print('SNR: ', snr)

