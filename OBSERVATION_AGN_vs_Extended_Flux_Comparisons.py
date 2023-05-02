# 14 Nov 2022: Simple code to measure the flux between the point source and extended emission in thr GATOS AGN sample

from Convenience_Functions import *

# read in data
im_path = 'Images/JWST/5_ERS_Data/1_Raw_data/ngc5728_MIRI_Imaging/'

# observed
image = get_pkg_data_filename(im_path + 'NGC5728_F560W_final.fits')
image_list = fits.open(image)
obs_img = image_list[0].data
# deconvolved
dec_image = get_pkg_data_filename(im_path + 'RL_Deconvolution/NGC5728_F1000W_DECONV_Leist.fits')
dec_image_list = fits.open(dec_image)
dec_img = dec_image_list[1].data

# set the FWHM value (pixel)for each image
dec_fwhm = 1.891

# measure centroid using a quadratic fit
guess = (151, 137)
box = 15
center_quad = centroid_quadratic(dec_img,
                                 xpeak=guess[0],
                                 ypeak=guess[1],
                                 fit_boxsize=box)
center = (center_quad[0], center_quad[1])

# measure the flux of the point source based on the aperture size
position = [(center[0]-0.2), (center[1]-0.2)]

aperture1 = photutils_CircularAperture(position, r=dec_fwhm)
phot_table = aperture_photometry(dec_img, aperture1, method='exact')
phot_table['aperture_sum'].info.format = '%.8g'
aperture_sum = phot_table['aperture_sum'][0]
aper_total = float(aperture_sum)

# measure the extended emission for each image
# Perform annulus photometry and remove the background from the total flux measured above
# define the inner and outer annulus radius
inner_rad = dec_fwhm
outer_rad = 12.153

# Background estimation
# Create the annulus aperture
position1 = [(center[0]), (center[1])]
# position1 = [x,y]
annulus_aperture = CircularAnnulus(position1, r_in=inner_rad, r_out=outer_rad)
phot_table3 = aperture_photometry(dec_img, annulus_aperture, method='exact')
phot_table3['aperture_sum'].info.format = '%.8g'
aperture_sum3 = phot_table3['aperture_sum'][0]
aper_total3 = float(aperture_sum3)

# Simple mean within a circular annulus
aperstats = ApertureStats(dec_img, annulus_aperture)

# mean local per-pixel background
bkg_mean = aperstats.mean

# The total background within the circular aperture is the mean local per-pixel background * the circular aperture
# area. Find the aperture area using the same area over which photometry was performed.
area = aperture1.area_overlap(dec_img)

# calculate the total backgroud within the circular aperture
total_bkg = bkg_mean * area

# calculate the background subtracted sum
final_phot = aper_total - total_bkg

# testing only
aperture2 = photutils_CircularAperture(position, r=outer_rad)
phot_table2 = aperture_photometry(dec_img, aperture2, method='exact')
phot_table2['aperture_sum'].info.format = '%.8g'
aperture_sum2 = phot_table2['aperture_sum'][0]
aper_total2 = float(aperture_sum2)

# subtract?
subtract = aper_total2 - aper_total

print('Point source flux (mJ/sr): ', aper_total)
print('Extended emission flux (mJy/sr): ', final_phot)
print('Extended emission flux [annulus - aperture] (mJy/sr): ', subtract)

# display the apertures to make sure the're fitting correctly
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 7), tight_layout=True)
aper_total = aper_total
final_phot = final_phot

# trim decimals
total = '{:.2e}'.format(aper_total)
sub = '{:.2e}'.format(subtract)

# set the image positions
x_positions = [106, 128.5, 151, 173.5, 196]
y_positions = [93, 115.5, 138, 160.5, 182]
labels = ['-10', '-5', '0', '5', '10']

ap_patches = aperture1.plot(color='r', lw=2)
ap_patches2 = aperture2.plot(color='w', lw=2, linestyle = '--')
im00 = plt.imshow(dec_img, origin='lower', cmap='viridis', norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=5300))
plt.title('5.6 μm Deconvolution', fontsize = 20)
ax.xaxis.set_major_locator(ticker.FixedLocator(x_positions))
ax.xaxis.set_major_formatter(ticker.FixedFormatter(labels))
ax.yaxis.set_major_locator(ticker.FixedLocator(y_positions))
ax.yaxis.set_major_formatter(ticker.FixedFormatter(labels))
ax.text(111, 103, 'Aperture flux (mJy/sr):  ' + str(total), color = 'r', fontsize = 18)
ax.text(111, 98, 'Annulus flux (mJy/sr):   ' + str(sub), color = 'w', fontsize = 18)
ax.set_xlim(106, 196)
ax.set_ylim(93, 183)
ax.set_ylabel('Y (")', fontsize = 15)
ax.set_xlabel('X (")', fontsize = 15)
add_colorbar(im00, label='Log-scaled Flux Density (mJy/sr)', format=FormatStrFormatter('%.1e'))
plt.savefig('5um_Deconvolution_Comparison.png')


plt.grid(False)
plt.show()

