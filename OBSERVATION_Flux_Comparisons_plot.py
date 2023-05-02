# 2022 AUG 23: simple code to compare flux measurements using different aperture sizes

from Convenience_Functions import *

model = 'residual_resized_NGC5728_F2100W_final'

# set the input path for the deconvolved image and merit function measurements
im_path = 'Images/JWST/4_Deconvolution/1_Richardson_Lucy/FITS/'

# set the output path for the deconvolved image and merit function measurements
mes_path = 'Images/JWST/4_Deconvolution/1_Richardson_Lucy/Measurements/'

# set the output path for all the images
out_path = 'Images/JWST/4_Deconvolution/1_Richardson_Lucy/Images/' + model

# set the out[ut directory based on reference PSF model
psf_model = '/WebbPSF'

# set the output directory to the aperture fit
fitA = 'AGN'
fitP = 'PSF'
fitG = 'Galaxy'

# read in the merit function data
nint1, flux1, fwhm1, strehl1 = np.loadtxt(mes_path + model + psf_model+ '/'+fitP+'/20220822/centroid_fixed_fwhm/RL_Image_Plot_' + str(model) + '.txt', unpack = True)
nint2, flux2, fwhm2, strehl2 = np.loadtxt(mes_path + model + psf_model+ '/'+fitP+'/20220822/centroid_1dg_fwhm/RL_Image_Plot_' + str(model) + '.txt', unpack = True)
nint3, flux3, fwhm3, strehl3 = np.loadtxt(mes_path + model + psf_model+ '/'+fitP+'/20220822/centroid_fixed_fwhm_scaled/RL_Image_Plot_' + str(model) + '.txt', unpack = True)
nint4, flux4, fwhm4, strehl4 = np.loadtxt(mes_path + model + psf_model+ '/'+fitP+'/20220822/centroid_1dg_fwhm_scaled/RL_Image_Plot_' + str(model) + '.txt', unpack = True)
nint5, flux5, fwhm5, strehl5 = np.loadtxt(mes_path + model + psf_model+ '/'+fitP+'/20220822/centroid_fixed/RL_Image_Plot_' + str(model) + '.txt', unpack = True)
nint6, flux6, fwhm6, strehl6 = np.loadtxt(mes_path + model + psf_model+ '/'+fitP+'/20220822/centroid_1dg/RL_Image_Plot_' + str(model) + '.txt', unpack = True)
nint7, flux7, fwhm7, strehl7 = np.loadtxt(mes_path + model + psf_model+ '/'+fitP+'/20220822/centroid_fixed+1/RL_Image_Plot_' + str(model) + '.txt', unpack = True)
nint8, flux8, fwhm8, strehl8 = np.loadtxt(mes_path + model + psf_model+ '/'+fitP+'/20220822/centroid_1dg+1/RL_Image_Plot_' + str(model) + '.txt', unpack = True)
nint9, flux9, fwhm9, strehl9 = np.loadtxt(mes_path + model + psf_model+ '/'+fitP+'/20220822/centroid_fixed+2/RL_Image_Plot_' + str(model) + '.txt', unpack = True)
nint10, flux10, fwhm10, strehl10 = np.loadtxt(mes_path + model + psf_model+ '/'+fitP+'/20220822/centroid_1dg+2/RL_Image_Plot_' + str(model) + '.txt', unpack = True)


n_int = []
flux_1 = []
flux_2 = []
flux_3 = []
flux_4 = []
flux_5 = []
flux_6 = []
flux_7 = []
flux_8 = []
flux_9 = []
flux_10 = []

for i in range(0, 50):
    # convert nint from float -> int
    nint_1 = nint1[i]
    nint_1 = int(nint_1)
    n_int.append(nint_1)

    # set rest as floats
    flux_1a = flux1[i]
    flux_2a = flux2[i]
    flux_3a = flux3[i]
    flux_4a = flux4[i]
    flux_5a = flux5[i]
    flux_6a = flux6[i]
    flux_7a = flux7[i]
    flux_8a = flux8[i]
    flux_9a = flux9[i]
    flux_10a = flux10[i]

    # convert Jy -> mJy
    flux_1a = flux_1a * 1000
    flux_2a = flux_2a * 1000
    flux_3a = flux_3a * 1000
    flux_4a = flux_4a * 1000
    flux_5a = flux_5a * 1000
    flux_6a = flux_6a * 1000
    flux_7a = flux_7a * 1000
    flux_8a = flux_8a * 1000
    flux_9a = flux_9a * 1000
    flux_10a = flux_10a * 1000

    fluxA = '{:.7f}'.format(flux_1a)
    fluxB = '{:.7f}'.format(flux_2a)
    fluxC = '{:.7f}'.format(flux_3a)
    fluxD = '{:.7f}'.format(flux_4a)
    fluxE = '{:.7f}'.format(flux_5a)
    fluxF = '{:.7f}'.format(flux_6a)
    fluxG = '{:.7f}'.format(flux_7a)
    fluxH = '{:.7f}'.format(flux_8a)
    fluxI = '{:.7f}'.format(flux_9a)
    fluxJ = '{:.7f}'.format(flux_10a)

    fluxA = float(fluxA)
    fluxB = float(fluxB)
    fluxC = float(fluxC)
    fluxD = float(fluxD)
    fluxE = float(fluxE)
    fluxF = float(fluxF)
    fluxG = float(fluxG)
    fluxH = float(fluxH)
    fluxI = float(fluxI)
    fluxJ = float(fluxJ)

    flux_1.append(fluxA)
    flux_2.append(fluxB)
    flux_3.append(fluxC)
    flux_4.append(fluxD)
    flux_5.append(fluxE)
    flux_6.append(fluxF)
    flux_7.append(fluxG)
    flux_8.append(fluxH)
    flux_9.append(fluxI)
    flux_10.append(fluxJ)


# Set the (x,y) data points
x = n_int
y1 = flux_1
y2 = flux_2
y3 = flux_3
y4 = flux_4
y5 = flux_5
y6 = flux_6
y7 = flux_7
y8 = flux_8
y9 = flux_9
y10 = flux_10


# plot the date
fig2 = plt.figure(figsize=(15, 6), tight_layout=True)
plt.subplot(121)
plt.plot(1,1,color='w',label = 'Aperture:')
plt.plot(x,y1, color = 'r', label = 'varies w/FWHM', linestyle = '-')
plt.plot(x,y3, color = 'g', label = 'scaled + varies w/ FWHM', linestyle = '--')
plt.plot(x,y5, color = 'b', label = 'fixed', linestyle = '-.')
plt.plot(x,y7, color = 'k', label = 'fixed+1', linestyle = ':')
plt.plot(x,y9, color = 'c', label = 'fixed+2', linestyle = 'dotted')
plt.title('centroid fixed')
plt.legend()

plt.subplot(122)
plt.plot(1,1,color='w',label = 'Aperture:')
plt.plot(x,y2, color = 'r', label = 'varies w/FWHM', linestyle = '-')
plt.plot(x,y4, color = 'g', label = 'scaled + varies w/ FWHM', linestyle = '--')
plt.plot(x,y6, color = 'b', label = 'fixed', linestyle = '-.')
plt.plot(x,y8, color = 'k', label = 'fixed+1', linestyle = ':')
plt.plot(x,y10, color = 'c', label = 'fixed+2', linestyle = 'dotted')
plt.legend()
plt.title('centroid_1dg')

plt.suptitle('21 μm NGC 5728: aperture and centroid comparisons', fontsize=15)
fig2.supxlabel('Iterations', fontsize=15)
fig2.supylabel('Flux Density (mJy)', fontsize=15)
plt.savefig(mes_path + model + psf_model+ '/'+fitP+'_flux_centroid_comparisons.png')
plt.show()