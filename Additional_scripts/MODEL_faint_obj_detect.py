# 18 AUg 2023: simple code to test the detectability of a faint object within our model

from Convenience_Functions import *

# read in FITS image
# read in data
miri = True
dec = False

if miri:
    path = 'Images/JWST/2_MIRISim_outputs/Deconvolution_testing/Residual_outputs/FITS/Ideal_Model/Model_AGN_complicated/current_12Aug2023/V5/'
    name = 'F2100W_V5_model_offset_19'

if dec:
    path = 'Images/JWST/4_Deconvolution/5_Kraken/FITS/Ideal_Model/Model_AGN_complicated/current_12Aug2023/V5/final_deconv/'
    name = 'F2100W_V5_DECONV_offset_16'
raw_image = get_pkg_data_filename(path + name + '.fits')
image_list = fits.open(raw_image)
image = image_list[0].data
# crop image data
image = image[53:203, 53:203]

# find sources with peaks 3-sigma above the background
mean, median, std = sigma_clipped_stats(image, sigma=3.0)
threshold = median + (3.0 * std)
tbl = find_peaks(image, threshold, box_size = 25)
print('photutils.find_peaks')
print(tbl)
print()

# find sources using DAOStarFinder?
daofind = DAOStarFinder(fwhm=3.0, threshold=3.0*std)
sources = daofind(image-median)

print('DAOStarFinder')
print('ID | x_centroid | y_centroid | flux')
pos = []
for i in range(0, int(len(sources['id']))):
    src_id = sources['id'][i]
    src_x = sources['xcentroid'][i]
    src_y = sources['ycentroid'][i]
    src_flux = sources['flux'][i]
    print(src_id, src_x, src_y, src_flux)

    positions = np.transpose((sources[i]['xcentroid'], sources[i]['ycentroid']))
    pos.append(positions)

# plot detections
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(7,7), constrained_layout=True)
ax.imshow(image, origin='lower', norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=image.max()), cmap='RdYlBu_r')
aper = photutils_CircularAperture(pos, r=2.0)
aper.plot(color='k', lw=1.5, alpha=1)

# add plot pretty things
positions1 = [5, 25, 50, 75, 100, 125, 145]
labels = ['-75', '-50', '-25', '0', '25', '50', '75']
ax.xaxis.set_major_locator(ticker.FixedLocator(positions1))
ax.xaxis.set_major_formatter(ticker.FixedFormatter(labels))
ax.xaxis.set_major_locator(ticker.FixedLocator(positions1))
ax.xaxis.set_major_formatter(ticker.FixedFormatter(labels))
ax.yaxis.set_major_locator(ticker.FixedLocator(positions1))
ax.yaxis.set_major_formatter(ticker.FixedFormatter(labels))
ax.yaxis.set_major_locator(ticker.FixedLocator(positions1))
ax.yaxis.set_major_formatter(ticker.FixedFormatter(labels))
ax.set_ylabel('Y (pix)', fontsize=15)
ax.set_xlabel('X (pix)', fontsize=15)
ax.set_title('MIRISim', fontsize=15)
ax.hlines(5, 18, 27, color = 'w')
ax.text(5, 8, '9 pix ~ 1"', color = 'w', fontsize=20)
plt.savefig(path+ 'images/' + name + '.png')
plt.show()


