# 8 May 2022: This is a callable class containing various convenience functions for this program. For each code in this package
# start by calling: 'from Convenience_Functions import *'. This will load all imports and methods below.

# necessary imports
from math import log10, sqrt
from skimage.metrics import structural_similarity as ssim
from scipy.ndimage import zoom
from astropy.coordinates import SkyCoord
from matplotlib._color_data import CSS4_COLORS
from matplotlib.patches import Circle
from mpl_toolkits import axes_grid1
import os
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from astropy.io import fits
import numpy as np
from photutils.aperture import aperture_photometry
import matplotlib.pyplot as plt
from astropy.visualization import simple_norm, ImageNormalize, MinMaxInterval, LogStretch, SqrtStretch, AsinhStretch, SinhStretch, PowerStretch
from photutils.aperture import CircularAperture as photutils_CircularAperture
from photutils.aperture import SkyCircularAperture as photutils_SkyAperture
from photutils.aperture import CircularAnnulus
from astropy.stats import sigma_clipped_stats
from astropy.convolution import Gaussian2DKernel, convolve
from photutils.centroids import centroid_2dg, centroid_com, centroid_quadratic, centroid_1dg
from photutils import RectangularAperture
from astropy.io import fits
from matplotlib.ticker import FormatStrFormatter
import scipy.interpolate, scipy.ndimage
import matplotlib
from scipy.ndimage import shift
import matplotlib.ticker as ticker
from matplotlib.ticker import MaxNLocator
from astropy import units as u
from photutils.psf import subtract_psf
from matplotlib.colors import LogNorm, SymLogNorm, Normalize
from mpl_toolkits.axes_grid1 import ImageGrid, make_axes_locatable
from sklearn.metrics import mean_squared_error
from math import sqrt
from skimage import restoration
import warnings
from astropy.utils.exceptions import AstropyUserWarning
from poppy import measure_fwhm, display_profiles, radial_profile, measure_ee, display_ee, display_psf
from poppy import *
import webbpsf
from specutils import Spectrum1D
from specutils.fitting import fit_lines, fit_generic_continuum
from webbpsf import measure_strehl, JWInstrument
from astropy.utils.data import get_pkg_data_filename
from astropy.modeling.models import Gaussian2D, Moffat2D, Gaussian1D
from photutils.psf import create_matching_kernel, TopHatWindow
from astropy.modeling import fitting, models
from FITS_tools.cube_regrid import downsample_cube
from photutils.detection import DAOStarFinder
from numpy.fft import fftn, ifftn, fftshift
from astropy.convolution import convolve_fft
from scipy import ndimage
import scipy.fft as fft
from skimage.color import rgb2hsv, rgb2gray
from skimage import color, data, metrics, restoration
import skimage.measure as measure
from skimage.io import imread, imshow
import skimage
from scipy.signal import convolve2d
from photutils.aperture import ApertureStats
from astropy.table import Table
import numpy as np
from astropy.coordinates import Angle
from photutils.aperture import EllipticalAperture as photutils_EllipticalAperture
from astropy.wcs import WCS
from skimage import color, data, restoration
from reproject import reproject_interp
import matplotlib.colors as mcolors
from reproject import reproject_exact
from photutils.psf.matching import resize_psf
from reproject import reproject_interp
from reproject.mosaicking import reproject_and_coadd
from reproject.mosaicking import find_optimal_celestial_wcs
from astropy.visualization import simple_norm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from PIL import Image
from astropy import stats
from astropy.stats import SigmaClip
from photutils.background import Background2D, MedianBackground
from photutils.datasets import make_noise_image
import math
from astropy.stats import sigma_clipped_stats
from photutils.detection import find_peaks
import pandas as pd
import warnings
from matplotlib.ticker import AutoMinorLocator
from photutils.profiles import RadialProfile, CurveOfGrowth
import scipy.optimize as opt
import scipy.ndimage as ndi
import photutils
import scipy.stats as sci_stats
from photutils.segmentation import detect_threshold, detect_sources
import aplpy as ap
from astropy.nddata import Cutout2D

# code used for RGB creation
import img_scale

#****************************************** Deconvolution Algorithms  *************************************************#

# 4) 11 April 2022: Richardson-Lucy deconvolution method- developed by Brian Northan
def richardson_lucy_np(image, psf, num_iters):
    '''11 April 2022: Non-circulant Richardson-Lucy deconvolution algorithm- developed by Brian
    Northan'''
    otf = fftn(fftshift(psf))
    otf_ = np.conjugate(otf)
    estimate = image
    #estimate = np.ones(image.shape)/image.sum()

    for i in range(num_iters):
        # print(i)
        reblurred = ifftn(fftn(estimate) * otf)
        ratio = image / (reblurred + 1e-30)
        estimate = estimate * (ifftn(fftn(ratio) * otf_)).astype(float)

    return estimate

#****************************************** Merit functions ############################################################

def measure_merits(data, radius, boxsize, smooth_radius, center, convert = False, sanity_check = False, plot=False):
    '''16 Jan 2024: new compact form to measure merit functions (FWHM reduction and flux conservation. Method makes use
    of (1) photutils.centroids.centroid_quadratic, (2) photutils.CircularAperture, (3) photutils.aperture_photometry,
    and (4) photutils.ApertureStats

        Parameters
    -----------------------
    data:           np.array
                    2D-data array containing the image data
    radius:         int
                    phot_aper radius size (also used as the annulus inner radius size for the smoothness parameter)
    boxsize         float
                    Resize the image data (smoothness parameter)
                    NOTE: this is the size of the object that encloses 95% of the flux (measured BEFORE this function
                    is called)
    smooth_radius   float
                    The boxcar filter width (equal to the filter diffraction limit) (smoothness parameter)
    center:         list
                    x,y initial guess positions for the image centroid
    convert:        boolean
                    COnvert FWHM values from pixels -> arsecs, aperture flux values from mJy/sr -> Jy
    sanity_check:   boolean
                    print out each measured value is called
    fwhm_method:    str
                    17 Jan. 2024: FOR TESTING ONLY
                    Toggle between measuring the FWHM using a 1D or 2D-Gaussian fit
    plot:           boolean
                    plot the radial profile and S-index measurement
    -----------------------
        Returns
    FWHM, aperture flux, smoothness, and centroid position'''

    # set the centroid location
    # NOTE: measured BEFORE this method is called
    center=center

    # measure the FWHM using radial_profile_FWHM
    fwhmrad = radial_profile_FWHM(data, center=center, plot=plot)
    fwhm = fwhmrad[1]

    # measure the photometry using photutils.CircularAperture and photutils.aperture_photometry
    phot_aper = photutils_CircularAperture(center, r=radius)
    phot_table = aperture_photometry(data, phot_aper, method='exact')
    aper_flux = np.float64(phot_table['aperture_sum'][0])

    # measure the smoothness from Concentration, Asymmetry, and Smoothness (CAS) statistics
    S_param = smoothness(data, boxradius=radius, boxsize=boxsize, smooth_radius=smooth_radius, sanity_check=False,plot=plot)


    # 16 Jan 2024: add code for background subtraction?

    if convert:
        # return the FWHM in arcseconds
        fwhm = fwhmrad[0]

        # 16 Jan 2024: convert aperture flux from mJy/sr -> Jy
        # 19 Jan 2024: double-check conversion factor is correct
        factor = (0.1110 ** 2) * 2.34806E-05
        aper_flux = aper_flux * factor

    # scale the outputs to 5-decimal places and convert back to float
    fwhm = '%.8f' % fwhm
    fwhm = np.float64(fwhm)
    aper_flux = '%.8f' % aper_flux
    aper_flux = np.float64(aper_flux)
    S_param = '%.8f' % S_param
    S_param = np.float64(S_param)

    if sanity_check:
        print('User-defined box size: ', boxsize)
        print('User-defined initial centroid guess (x,y): ', center)
        print('Measured centroid (x,y): ', center[0], center[1])
        print('User-defined aperture radius sizes (fwhm, flux): ', radius[0], radius[1])

        if convert:
            print('Merit functions (FWHM ("), aperture flux (Jy)): ', fwhm, aper_flux)
        else:
            print('Merit functions (FWHM (pix), aperture flux (mJy/sr)): ', fwhm, aper_flux)

        print('Smoothness (S): ', S_param)

    # return the FWHM, aperture photometry, smoothness, and centroid
    return fwhm, aper_flux, S_param, center

def radial_profile_FWHM(data, center, plot = False):
    '''14 Feb 2024: measure the FWHM of an image array through a 1D-Gaussian fit of the radial profile


        Parameters
    -----------------------
    data:           np.array
                    2D-data array containing the image data
    radius:         list
                    fwhm_aper and phot_aper radius sizes
    center:         list
                    x,y initial guess positions for the image centroid
    boxsize:        int
                    Box size to search for centroid position
    plot            boolean
                    plot the radial profile and 1D-Gaussian fit
    -----------------------
        Returns
    FWHM ("), FWHM (pixels)'''

    x,y = data.shape
    # set the radius to calculate the number of bins for the radial profile
    radius = np.linspace(0.11, 256, 256)
    #radius = np.arange(256)

    center = center

    # measure the FWHM from the image array radial profile using photutils.profiles
    rp = RadialProfile(data, center, radius, error=None, mask=None, method='exact')

    # fit a 1D-Gaussian to the radial profile
    rp.gaussian_fit

    # return the fwhm values from the 1D-Gaussian fit (in pixels)
    pix_fwhm = rp.gaussian_fwhm

    # convert from pixels to arcseconds
    arc_fwhm = pix_fwhm * 0.11

    if plot:
        # plot the radial profile and 1D-Gaussian fit
        plt.plot(rp.radius, rp.profile, label='Radial Profile')
        plt.plot(rp.radius, rp.gaussian_profile, label='Gaussian Fit\n FWHM (") = ' + str('%.5f' % arc_fwhm) + '\nFWHM (pix) = ' + str('%.5f' % pix_fwhm))
        plt.legend()
        plt.show()

    return arc_fwhm, pix_fwhm

def COG(data, xycent, max_level, params = False):
    '''23 Feb 2024: code to generate a curb of growth plot and return the pixel radius where 99% of the object flux is contained.
    parameters
    __________
    data:       numpy.array
                2D-image array data
    xycent:     tuple
                tuple of x,y image array float centroid values
    params:     boolean
                prints the individual steps for testing
    outputs
    pixel radius containing 99% of the enclosed energy, curve of growth plot,
    '''
    y,x = data.shape

    # set the radial bins such that the boundary is the x-max lim
    cenx, ceny = xycent[0], xycent[1]
    radx = 2*round(int(x-cenx))

    # set the radius to calculate the number of bins for the encircled energy (EE)
    # 12 Mar 2024: UNCOMMENT FOR TESTING ONLY
    #radius = np.arange(1, radx/2)
    radius = np.arange(1, radx)

    # measure the Curve of Growth from the image array using photutils.profiles
    cog = CurveOfGrowth(data, xycent, radius, method='exact')

    # save the encircled energy profile
    EE = cog.profile

    # save the EE radius profile
    rad_profile = cog.radius

    # return the maximum energy value
    maxEE = max(EE)

    # find the aperture radius that corresponds to user-defined maximum energy value
    max99 = maxEE * max_level/100

    # save the radial levels were the EE is at 50%, 80% and 95%
    for level in [max99]:
        if (EE >= level).any():
            EElev = rad_profile[np.where(EE >= level)[0][0]]
            # Return the radial size at 99%
            if level < max99 + 0.1:
                max99_rad = EElev

    if params:
        print('(y,x) image shape: ', y,x)
        print('(y,x) image centroid: ', cenx, ceny)
        print('Maximum number of radial bins to compute the encircled energy: ', radx)
        print('Max encircled energy: ', maxEE)
        print('95% max encircled energy: ', max99)
        print('Aperture size to enclose 95% max encircled energy: ', max99_rad)

    return max99_rad, EE

def sky_smoothness(data, smooth_radius):
    """
    Compute the smoothness of the sky background. NOTE: Lotz et al. (2004) computed the background smoothness as a
    function of the aperture area used to measure the galaxy flux. Here we have modified this to measure the background
    as the sum of 4 10x10 pixel boxes along the corners of the image area following the conventions of
    Torres-Quijano et al. (2021).

    Parameters
    __________
    data            numpy.array
                    Array containing the 2D-FITS image data
    smooth_radius   float
                    The boxcar filter width (equal to the filter diffraction limit)
    Outputs
    ___________
    The sum of the sky background flux
    """

    y,x = data.shape
    extract = int(round(y/2))
    box_size = 10

    # select a region of sky background along the image edge
    # define the four box origins
    x1, y1 = 0, 0
    x2, y2 = extract * 2, 0
    x3, y3 = 0, extract * 2
    x4, y4 = extract * 2, extract * 2

    # define the regions to measure stats from
    box1 = data[y1:y1 + box_size, x1:x1 + box_size]
    box2 = data[y2:y2 + box_size, x2 - box_size: x2]
    box3 = data[y3 - box_size:y3, x3:x3 + box_size]
    box4 = data[y4 - box_size:y4, x4 - box_size:x4]

    # compute the boxcar filter size for each box
    #box1
    bkg_smooth1 = ndi.uniform_filter(box1, round(smooth_radius))
    bkg_diff1 = box1 - bkg_smooth1
    bkg_diff1[bkg_diff1 < 0] = 0.0
    # #box2
    bkg_smooth2 = ndi.uniform_filter(box2, round(smooth_radius))
    bkg_diff2 = box2 - bkg_smooth2
    bkg_diff2[bkg_diff2 < 0] = 0.0
    # #box3
    bkg_smooth3 = ndi.uniform_filter(box3, round(smooth_radius))
    bkg_diff3 = box3 - bkg_smooth3
    bkg_diff3[bkg_diff3 < 0] = 0.0
    # #box4
    bkg_smooth4 = ndi.uniform_filter(box4, round(smooth_radius))
    bkg_diff4 = box4 - bkg_smooth4
    bkg_diff4[bkg_diff4 < 0] = 0.0

    # compute the sums for each smoothed background difference
    sum1 = np.sum(bkg_diff1) / float(box1.size)
    sum2 = np.sum(bkg_diff2) / float(box2.size)
    sum3 = np.sum(bkg_diff3) / float(box3.size)
    sum4 = np.sum(bkg_diff4) / float(box4.size)

    # return the final value as the average of the four summations
    final_sum = (sum1+sum2+sum3+sum4)

    return final_sum

def smoothness(data, boxradius, boxsize, smooth_radius, sanity_check = False, plot=False):
    """
    Calculate the smoothness (a.k.a. clumpiness) as defined in eq. (11) from Lotz et al. (2004).

    NOTES:
        1) The original definition by Conselice (2003) includes an additional factor of 10.
        2) Lotz et al. (2004) computed the smoothness as a function of the Petrosian radius (rp). For deconvolution purposes
        the Petrosian radius is too small of an aperture to capture the orange-peel effect we want to measure, so instead
        we measured the flux within an annulus whose inner radius is determined by the size (pixel) to enclose the first
        Airy ring of the F2100W JWST MIRI/MIRIM PSF filter (see Leist et al. 2024 for more details) while the outer
        radius is set as the limit where 95% of the object flux is contained.
        3) Lotz et al. (2004) set the boxcar filter width (σ) equal to 0.25rp. We set the filter width equal to the
        diffraction limit for each filter.

    Parameters
    __________
    data            numpy.array
                    Array containing the 2D-FITS image data
    boxradius       float
                    Set the annulus inner radius size
    boxsize         float
                    Resize the image data
    smooth_radius   float
                    The boxcar filter width (equal to the filter diffraction limit)
    sanity_check    boolean
                    Print the individual step measurements for sanity checking
    plot            boolean
                    Plot the I(i,j), I_S(i,j), and residual images.
    Outputs
    ___________
    The smoothness (S) of an image
    """

    # Compute the image centroid
    # NOTE: assumed to be the brightest pixel in the image array
    imy1, imx1 = np.unravel_index(np.argmax(data, axis=None), data.shape)

    # reshape the data array to a size that contains 95% of the flux
    # Note: the x,y values are padded +/- 0.5 pixels to give a whole integer result
    resize_data = data[int(imy1-((boxsize/2)-0.5)):int(imy1+((boxsize/2)+0.5)), int(imx1-((boxsize/2)-0.5)):int(imx1+((boxsize/2)+0.5))]

    # Compute the image centroid of the reshaped image array
    imy2, imx2 = np.unravel_index(np.argmax(resize_data, axis=None), resize_data.shape)
    # determine the centroid from a quadratic fit using the image peak as the initial guess
    xycent = centroid_quadratic(resize_data, xpeak=imx2, ypeak=imy2, fit_boxsize=15)

    # Exclude central region during smoothness calculation:
    # annulus size set by the aperture radius
    r_in = boxradius
    r_out = (boxsize/2)-1
    ap = photutils.aperture.CircularAnnulus(xycent, r_in, r_out)

    # Note: ndi.uniform_filter will only work if the radius size is given as a whole integer. We therefore rounded the
    # radius size to the nearest whole integer.
    data_smooth = ndi.uniform_filter(resize_data, size=round(smooth_radius))

    # compute the residual image: data-data_smooth
    data_diff = resize_data - data_smooth
    # set negative pixels to zero
    data_diff[data_diff < 0] = 0.0

    # measure the photometry of the original and residual images
    ap_flux = ap.do_photometry(resize_data, method='exact')[0][0]
    ap_diff = ap.do_photometry(data_diff, method='exact')[0][0]

    # compute the sky smoothness value
    sky_bkg = sky_smoothness(data, smooth_radius)

    # compute the Smoothness parameter
    S = np.abs(ap_diff - sky_bkg) / ap_flux

    if sanity_check:
        # print the individual measurements for sanity checking
        print('Smoothness filter size: ', smooth_radius)
        print('Input image centroid (x,y): ', imx1, imy1)
        print('Image shape containing 95% of object flux (x,y): ', resize_data.shape)
        print('Resized image centroid (x,y): ', imx2, imy2)
        print('Annulus radius (in | out): ', r_in, ' | ', r_out)
        print('Annulus flux (resized image): ', ap_flux)
        print('Annulus flux (difference image): ', ap_diff)
        print('Difference sky background flux: ', sky_bkg)
        print('Smoothness (S): ', S, '\n')

    if plot:
        # plot the three images: (1) resized, (2) smoothed resized, (3) difference
        fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(14,6), constrained_layout=True)

        # take the resized image data shape
        y,x = resize_data.shape


        # resized image
        ax[0].imshow(resize_data, origin='lower', norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=resize_data.max()),cmap='turbo')
        ax[0].set_title('I(i,j)', fontsize=20)

        # smoothed image
        ax[1].imshow(data_smooth, origin='lower', norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=data_smooth.max()),cmap='turbo')
        ax[1].set_title('$I^{γ}_{S}$(i,j)', fontsize=20)
        ax[1].text(5, 5, 'γ = ' + str(smooth_radius) + ' (pixels)', color = 'w', fontsize=20)

        # differenced image
        ax[2].imshow(data_diff, origin='lower', norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=data_diff.max()),cmap='turbo')
        ax[2].set_title('I(i,j) - $I^{γ}_{S}$(i,j)', fontsize=20)
        # plot the annulus data?
        ap.plot(color='w', lw=1)

        # set the plot x,y labels
        positions = [0, y/2, int(round(y))-1]
        labels = ['-'+ str(y/2), '0', str(y/2)]

        for i in range(0,3):
            ax[i].xaxis.set_major_locator(ticker.FixedLocator(positions))
            ax[i].xaxis.set_major_formatter(ticker.FixedFormatter(labels))
            ax[i].yaxis.set_major_locator(ticker.FixedLocator(positions))
            ax[i].yaxis.set_major_formatter(ticker.FixedFormatter(labels))

            if i == 0:
                ax[i].tick_params(direction='in', color='w', length=6, axis='both', left=True, top=False,right=False, bottom=True, labelleft=True, labeltop=False, labelright=False,labelbottom=True)
            else:
                ax[i].tick_params(direction='in', color='w', length=6, axis='both', left=True, top=False, right=False,bottom=True, labelleft=False, labeltop=False, labelright=False, labelbottom=True)

        fig.supylabel('Y (pixels)', fontsize=20)
        fig.supxlabel('X (pixels)', fontsize=20)
        plt.savefig('S_index_plot.pdf')
        plt.show()

    return S

#****************************************** Centroiding algorithms  ***************************************************#

def centroid_DAO(HDUlist = None, ext = None, fwhm = None, threshold = None, sigma = None):
    '''30 Sep 2022: Determine the centroid of an image using DAOStarFinder
    Parameter
    __________
    data:       numpy.array
                2D-array containing the image data
    fwhm:       float
                fwhm guess
    threshold:  float
                value in which the threshold the image below
    sigma:      float
                set the sigma threshold
    Returns
    _______
    center:     numpy.array
                x,y centroid coordinates
    '''
    # Read in data
    if HDUlist == HDUlist:
        HDUlist = fits.open(HDUlist)
    else:
        raise ValueError('Input must be a FITS HDUlist')

    # read the user-defined extension
    if ext == ext:
        data = HDUlist[ext].data
    else:
        # assume extension is 0
        data = HDUlist[0].data

    if fwhm == None:
        raise ValueError('Please set the fwhm guess')
    else:
        fwhm = fwhm

    if threshold == None:
        raise ValueError('Please set the image threshold value')
    else:
        threshold = float(threshold)

    if sigma == None:
        raise ValueError('Please set the sigma value')
    else:
        sigma = float(sigma)


    # determine the stddev above which to threshold
    mean, median, std = sigma_clipped_stats(data, sigma = sigma)

    # calculate the centroid of the image by determing the local density maxima that have a peak amplitude
    # greather than the given threshold and have a size/shape similiar to the defined 2D-Gaussian kernel, w/
    # first guess at FWHM value
    center_find = DAOStarFinder(fwhm=fwhm, threshold=threshold * std, exclude_border=False)
    centerD = center_find(data - median)

    x = []
    y = []
    for col in centerD.colnames:
        centerD[col].info.format = '%.8g'

    for num in range(len(centerD)):
        if (centerD['peak'][num] == np.amax(centerD['peak'])):
            # return the x,y centroid position
            pos_x = float(centerD['xcentroid'][num])
            x.append(pos_x)
            pos_y = float(centerD['ycentroid'][num])
            y.append(pos_y)

    # append the x and y centroid position to a single array
    center = x + y

    return center

#****************************************** Radial profile tools ******************************************************#

def MAE_radial_profile(HDUlist = None, ext=None, binsize=None, max_radius = int,  normalize = False, stddev=False,
                       print_params = False, centroid=None):
    '''Compute a discrete radial profile evaluated on the provided binsize.

    Code taken/modified from pydatatut.pdf and poppy.utils.
    Copyright (C) 2010-2021 Association of Universities for Research in
    Astronomy (AURA)

    Redistribution and use in source and binary forms, with or without modification, are permitted provided that the
    following conditions are met:

    Redistributions of source code must retain the above copyright notice, this list of conditions and the following
    disclaimer.

    Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following
    disclaimer in the documentation and/or other materials provided with the distribution.

    The name of AURA and its representatives may not be used to endorse or promote products derived from this software
    without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY AURA “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL AURA BE
     LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
     LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
     HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
      OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

    Parameters
    __________
    HDUList:                    fits.HDUlist
                                FITS file containing image to display
    ext:                        int
                                FITS extension, default = 0
                                    binsize:                    float
                                Size of step for the profile. Default is pixel size.
    centroid:                   string
                                Option to chose between 'centroid_1dg' and 'DAO_StarFinder'
                                1) 'default'- use DAOStarFinder to find the objects centroid
                                2) '1dg'- fit a 1D-Gaussian to the object
    max_radius:                 int
                                Sets the size to measure the radial profile
    normalize:                  boolean
                                OPTIONAL: normalize the input data to a peak of 1
    stddev:                     bool
                                OPTIONAL: compute the standard deviation in each radial bin.
    print_params:               boolean
                                OPTIONAL: print the results of each step for testing purposes
    Returns
    ________________
    results: tuple
         radius, profile, encircle energy (EE). EE (50,80,95), EE_radius (50,80,95), x,y centroid and (optional) sttdev
    '''

    # Read in data
    if HDUlist == HDUlist:
        HDUlist = fits.open(HDUlist)
    else:
        raise ValueError('Input must be a FITS HDUlist')

    # read the user-defined extension
    if ext == ext:
        data = HDUlist[ext].data
    else:
        # assume extension is 0
        data = HDUlist[0].data

    # Set the max radius to measure out to
    if max_radius == None:
        print('ERROR 1: No max radius value defined.')
    else:
        max_radius = max_radius

    # OPTIONAL: normalize data
    if normalize:
        data = data / data.max()

    if centroid == '1dg':
        # find the center of a 1D-Gaussian fit to the data
        center = centroid_1dg(data)

    elif centroid == 'quad':
        center = centroid_quadratic(data)

    elif centroid == 'exact':
        # set the centroid as the center of the image array
        y, x = data.shape
        im_y = int(x / 2)
        im_x = int(x / 2)
        center = [(im_x), (im_y)]

    elif centroid == centroid:
        center = centroid

    else:
        raise ValueError('No centroiding method defined, options are: 1) fit, 2) exact')

    # save as float for centroid testing
    cenx = center[0]
    ceny = center[1]

    # convert to int for processing
    X1 = int(center[0])
    Y1 = int(center[1])

    # read the pixel scale and observed wavelength
    # pixel_scale = HDUlist[ext].header['PIXELSCL']
    # pixel_scale = float(pixel_scale)
    pixel_scale = 0.1110
    #
    # wavelength = HDUlist[ext].header['WAVELEN']
    # wavelength = float(wavelength)
    wavelength = 5.6

    # Set max radial size to measure
    data = data[Y1-max_radius:Y1+max_radius, X1-max_radius:X1+max_radius]

    # Find the centroid of the new scaled image
    # find the center of a 1D-Gaussian fit to the data
    center = centroid_quadratic(data, fit_boxsize=15)
    Y = int(center[1])
    X = int(center[0])

    # Set the image boundaries
    y, x = np.indices(data.shape)

    # Set the binsize for the size of the profile step
    if binsize == None:
        binsize = pixel_scale
    else:
        binsize = binsize

    # Determine the radii of all pixels
    rad_all = np.sqrt((x - X) ** 2 + (y - Y) ** 2) * pixel_scale / binsize

    # Get sorted indices
    sort_int = np.argsort(rad_all.flat)

    # Sorted radii
    rad_sort = rad_all.flat[sort_int]

    # Sort image values by radii
    im_sort = data.flat[sort_int]

    # Interger part of the radii
    rad_im = rad_sort.astype(int)

    # Find numbers in the same radius bin by looking for where the radius changes
    # value and determining the distance between these changes.
    # Assume all radii are present
    delta_rad = rad_im[1:] - rad_im[:-1]

    # Determine the location of the changed radius
    rad_ind = np.where(delta_rad)[0]

    # Number in the radius bin
    rad_num = rad_ind[1:] - rad_ind[:-1]

    # Cumulative sum to figure out sums for each radii bin
    cum_sum_im = np.cumsum(im_sort, dtype=float)

    # Sum for image values in radius bin
    tbin = cum_sum_im[rad_ind[1:]] - cum_sum_im[rad_ind[:-1]]

    # Create the radial profile
    radial_profile = tbin / rad_num

    # Pre-pend the initial element the above code misses
    rp_final = np.empty(len(radial_profile) + 1)
    if rad_ind[0] != 0:
        # if there are multiple elements in the center bin, average them
        rp_final[0] = cum_sum_im[rad_ind[0]] / (rad_ind[0] + 1)
    else:
        # Otherwise if there's just one value take it
        rp_final[0] = cum_sum_im[0]

    rp_final[1:] = radial_profile

    # These should be centered in the bins, so add a half [1/3(?)]
    radius = np.arange(len(rp_final)) * binsize + binsize * 0.3

    # Determine the Encircled Energy
    EE = cum_sum_im[rad_ind]

    # Return the encircled energy at the given radius
    # 50% EE
    fifty = EE.max() * 0.5
    # 80% EE
    eighty = EE.max() * 0.8
    # 95% EE
    ninety_five = EE.max() * 0.95
    # 100% EE
    EE_full = EE.max()
    # save these three EE values
    EE_rad = [fifty, eighty, ninety_five, EE_full]

    # save the radial levels were the EE is at 50%, 80% and 95%
    levels = []
    for level in [EE_rad[0], EE_rad[1], EE_rad[2], EE_rad[3]]:

        if (EE >= level).any():
            EElev = radius[np.where(EE >= level)[0][0]]

            # Return the radial size at 50%
            if level < EE_rad[0] + 1:
                levels.append(EElev)
            # return the radial size at 80%
            elif level < EE_rad[1] + 1:
                levels.append(EElev)
            elif level < EE_rad[2] + 1:
                # return the radial size at 95%
                    levels.append(EElev)
            else:
                # return the radial size at 100%
                levels.append(EElev)

    if print_params:
        print('Method: MAE_radial_profile')
        print('Method: MAE_measure_FWHM_gaussian')
        print('Observed image shape: ', data.shape)
        print('Pixel scale: ', pixel_scale)
        print('Observed wavelength: ', wavelength)
        print('Image centroid (x,y): ', cenx, ceny)
        print('Radii of all pixels: ', rad_all)
        print('Location of changed radius: ', rad_ind)
        print('Number in the radius bin: ', rad_num)
        print('Cummulative sum to figure out sums for each radii bin: ', cum_sum_im)
        print('Sum for image values in radius bin', tbin)
        print('radius measurements: ', radius)
        print('Radial profile plot: ', radial_profile)
        print('Enclosed energy plot: ', EE)
        print('Flux level values: ', EE_rad)
        print('Flux radii: ', levels)
        print('centroid (x,y): ', cenx, ceny)

    # Calculate the standard deviation
    if stddev == False:
        # return: radius, radial profile, enclosed energy, energy at 50%, 80%, 95% and 100%, radius at 50, 80, 95 and 100%
        # x-centroid position, y-centroid position, image plate scale, wavelength
        return (radius, rp_final, EE, EE_rad, levels,  cenx, ceny, pixel_scale, wavelength)

    elif stddev == True:
        stddevs = np.zeros_like(rp_final)
        r_pix = rad_all * binsize
        for i, rr in enumerate(radius):
            if i == 0:
                wg = np.where(rad_all < rr + binsize / 2)
            else:
                wg = np.where((r_pix >= (rr - binsize / 2)) & (r_pix < (rr + binsize / 2)))
            stddevs[i] = data[wg].std()

        # return: radius, radial profile, enclosed energy, energy at 50%, 80% and 95%, radius at 50, 80 and 95%
        # x-centroid position, y-centroid position, image plate scale, wavelength and standard deviation
        return (radius, rp_final, EE, EE_rad, levels, cenx, ceny, pixel_scale, wavelength, stddevs)

def MAE_display_profile(HDUlist=None, ext=None, max_radius = None, title = None, binsize = None, print_params = False,
                        plot = None, display = False, save = None, centroid = None):
    '''Produce two plots of PSF radial profile and encircled energy (EE). Option to
     plot just the EE plot or the radial and EE plot

    Code taken/modified from pydatatut.pdf and poppy.utils.
    Copyright (C) 2010-2021 Association of Universities for Research in
    Astronomy (AURA)

    Redistribution and use in source and binary forms, with or without modification, are permitted provided that the
    following conditions are met:

    Redistributions of source code must retain the above copyright notice, this list of conditions and the following
    disclaimer.

    Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following
    disclaimer in the documentation and/or other materials provided with the distribution.

    The name of AURA and its representatives may not be used to endorse or promote products derived from this software
    without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY AURA “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL AURA BE
     LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
     LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
     HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
      OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

     Parameters
     __________
     HDUList:                    fits.HDUlist
                                 FITS file containing image to display
     ext:                        int
                                 FITS extension, default = 0
     max_radius:                 int
                                 Sets the size to measure the radial profile
     title:                      string
                                 title of the output plot
     print_params:               boolean
                                OPTIONAL: print the results of each step for testing purposes
     plot:                       string
                                 Sets the plot to be displayed, options are:
                                 1) 'both'-   plot both the radial profile and encircled energy plots
                                 2) 'radial'- plot the radial profile only
                                 3) 'EE'-     plot EE only
    display:                    boolean
                                Toggle whether to display the image or not
    save:                       string
                                OPTIONAL: the path to save the output image in
     Returns
     __________
     Plot: 			2D-array
                     Image array containing either the radial profile, EE_plot or both
     '''

    # read in the data
    if HDUlist == HDUlist:
        HDUlist = HDUlist
    else:
        raise ValueError('Input must be a FITS HDUlist')

    # set the extension to the data
    if ext == ext:
        ext = ext
    else:
        raise ValueError('No extension to data given')

    # Set the max radius to measure out to
    if max_radius == max_radius:
        max_radius = max_radius
    else:
        raise ValueError('Max radius not set')

    # Set the binsize for the size of the profile step
    if binsize == binsize:
        binsize = binsize
    else:
        binsize = None

    # set the title
    if title == title:
        title = title

    if centroid == '1dg':
        centroid= centroid

    elif centroid == 'quad':
        centroid = centroid

    elif centroid == 'exact':
        centroid = centroid

    elif centroid == centroid:
        centroid = centroid

    # Measure the radial profile of the image
    radius, radial_profile, EE, EE_rad, levels, x, y, pixel_scale, wavelength = MAE_radial_profile(HDUlist = HDUlist,
                                                                                                          ext=ext,
                                                                                                          binsize=binsize,
                                                                                                          max_radius = max_radius,
                                                                                                   centroid = centroid)

    if print_params:
        print('Method: MAE_display_profile')
        print('Max radius: ', max_radius)
        print('User-defined title: ', title)
        print('NOTE: following params determine by MAE_radial_profile method')
        print('Radius measurements: ', radius)
        print('Radial profile plot: ', radial_profile)
        print('Enclosed energy plot: ', EE)
        print('Flux level values: ', EE_rad)
        print('Flux radii: ', levels)
        print('centroid (x,y): ', x, y)
        print('pixel scale: ', pixel_scale)
        print('Obs wavelength: ', wavelength)

    # plot user-selected outputs
    if plot == 'both':
        plt.figure(figsize=(10, 5))
        plt.clf()

        # Plot radial data
        plt.subplot(2, 1, 1)
        plt.title(title)
        plt.ylabel('Radial Profile')
        plt.semilogy(radius, radial_profile)

        # Measure FWHM
        fwhm = MAE_measure_FWHM(HDUlist = HDUlist, ext = ext, level=0.5, binsize=binsize, max_radius=max_radius, centroid = centroid)

        # Add FWHM values to radial plot
        plt.text(fwhm[0], radial_profile[0] * 0.5, 'FWHM=%.3f"' % fwhm[0])

        # Plot encircled energy
        plt.subplot(2, 1, 2)
        plt.semilogy(radius, EE, color='r')
        plt.xlabel('Radius [arcsec]')
        plt.ylabel('Enclosed Energy')

        # Plot different EE levels
        # set radial boundaries: 50%, 80% and 95%
        for level in [EE_rad[0], EE_rad[1], EE_rad[2]]:
            if (EE > level).any():
                EElev = radius[np.where(EE > level)[0][0]]
                levels = []
                # Return the radial size at 50%
                if level < EE_rad[0] + 1:
                    yoffset = 0
                    levels.append(0.5)
                # return the radial size at 80%
                elif level < EE_rad[1] + 1:
                    yoffset = 0
                    levels.append(0.8)
                # return the radial size at 95%
                else:
                    yoffset = -0.05
                    levels.append(0.95)

                levels = np.array(levels)

                # plt.text(radius+1, EE level)
                plt.text(EElev + 0.1, level + yoffset, 'EE=%2d%% at r=%.3f"' % (levels * 100, EElev))

    elif plot == 'radial':
        plt.figure(figsize=(10, 5))
        plt.clf()

        # Plot radial data
        plt.title(title)
        plt.ylabel('Radial Profile')
        plt.semilogy(radius, radial_profile)

        # Measure FWHM
        fwhm = MAE_measure_FWHM(HDUlist = HDUlist, ext = ext, level=0.5, binsize=binsize, max_radius=max_radius,centroid= centroid)

        # Add FWHM values to radial plot
        plt.text(fwhm[0], radial_profile[0] * 0.5, 'FWHM=%.3f"' % fwhm[0])

    elif plot == 'EE':
        plt.figure(figsize=(10, 5))
        plt.clf()
        plt.title(title)
        plt.semilogy(radius, EE, color='r')
        plt.xlabel('Radius [arcsec]')
        plt.ylabel('Enclosed Energy')

        # Plot different EE levels
        # set radial boundaries: 50%, 80% and 95%
        for level in [EE_rad[0], EE_rad[1], EE_rad[2]]:
            if (EE > level).any():
                EElev = radius[np.where(EE > level)[0][0]]
                levels = []
                # Return the radial size at 50%
                if level < EE_rad[0] + 1:
                    yoffset = 0
                    levels.append(0.5)
                # return the radial size at 80%
                elif level < EE_rad[1] + 1:
                    yoffset = 0
                    levels.append(0.8)
                # return the radial size at 95%
                else:
                    yoffset = -0.05
                    levels.append(0.95)

                levels = np.array(levels)
                plt.text(EElev + 0.1, level + yoffset, 'EE=%2d%% at r=%.3f"' % (levels * 100, EElev))

    else:
        raise ValueError('No plotting option selected')

    if save:
        plt.savefig(save + '.png')

    if display:
        plt.show()

#****************************************** FWHM measuring tools ******************************************************#

def MAE_measure_FWHM(HDUlist = None, ext = None, level=0.5, binsize=None, max_radius=int, stddev = False,
                     print_params = False, centroid=None):
    '''Measure the FWHM by interpolation of the radial profile

    Code taken/modified from pydatatut.pdf and poppy.utils.
    Copyright (C) 2010-2021 Association of Universities for Research in
    Astronomy (AURA)

    Redistribution and use in source and binary forms, with or without modification, are permitted provided that the
    following conditions are met:

    Redistributions of source code must retain the above copyright notice, this list of conditions and the following
    disclaimer.

    Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following
    disclaimer in the documentation and/or other materials provided with the distribution.

    The name of AURA and its representatives may not be used to endorse or promote products derived from this software
    without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY AURA “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL AURA BE
     LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
     LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
     HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
      OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

    Parameters
    __________
    HDUList:                    fits.HDUlist
                                FITS file containing image to display
    ext:                        int
                                FITS extension, default = 0
    binsize:                    float
                                Size of step for the profile. Default is pixel size.
    max_radius:                 int
                                Sets the size to measure the radial profile
    stddev:                     bool
                                Compute the standard deviation in each radial bin.
    print_params:               boolean
                                OPTIONAL: print the results of each step for testing purposes
    Returns
    __________
    FWHM: float
        The FWHM of the image, computed in both arcsec and pixel
    '''

    # read in the data
    if HDUlist == HDUlist:
        HDUlist = HDUlist
    else:
        raise ValueError('Input must be a FITS HDUlist')

    # set the extension to the data
    if ext == ext:
        ext = ext
    else:
        # assume extension is 0
        ext = ext

    # Set the binsize for the size of the profile step
    if binsize == binsize:
        binsize = binsize
    else:
        binsize = None

    # Set the max radius to measure out to
    if max_radius == max_radius:
        max_radius = max_radius
    else:
        raise ValueError('Max radius not set')

    # toggle whether to return stddev or not
    if stddev:
        stddev = stddev

    if centroid == '1dg':
        centroid= centroid

    elif centroid == 'quad':
        centroid = centroid

    elif centroid == 'exact':
        centroid = centroid

    elif centroid == centroid:
        centroid = centroid

    else:
        raise ValueError('No centroiding method defined, options are: 1) fit, 2) exact')


    # Measure the radial profile of the image
    if stddev:
        radius, radial_profile, EE, EE_rad, levels, x, y, pixel_scale, wavelength, stddev= MAE_radial_profile(HDUlist = HDUlist,
                                                                                                          ext=ext,
                                                                                                          binsize=binsize,
                                                                                                          max_radius = max_radius,
                                                                                                          stddev=stddev,
                                                                                                              centroid = centroid)
    else:
        radius, radial_profile, EE, EE_rad, levels, x, y, pixel_scale, wavelength = MAE_radial_profile(HDUlist = HDUlist,
                                                                                                          ext=ext,
                                                                                                          binsize=binsize,
                                                                                                          max_radius = max_radius,
                                                                                                       centroid=centroid)


    # Determine the peak of the profile
    rp_max = radial_profile.max()

    # Find the lower width limit
    wlower = np.where(radial_profile < rp_max * level)
    if len(wlower[0]) == 0:
        raise ValueError(
            "The supplied array's pixel values never go below {0:.2f} of its maximum, {1:.3g}. " +
            "Cannot measure FWHM.".format(level, rp_max))

    wmin = np.min(wlower[0])

    # Measure just past the half way mark
    winterp = np.arange(0, wmin + 2, dtype=int)[::-1]

    # Set the interpreter
    if len(winterp) < 6:
        kind = 'linear'
    else:
        kind = 'cubic'

    # 5) Set the interpreter to measure the FWHM
    interp_hw = scipy.interpolate.interp1d(radial_profile[winterp], radius[winterp], kind=kind)

    # 6) Save FWHM measure, and convert to pixels based on input pixel scale
    # a) Original (arcsec)
    fwhm_arcsec = 2 * interp_hw(rp_max * level)

    # b) Pixel scale
    fwhm_pixel = fwhm_arcsec / pixel_scale

    if print_params:
        print('Method: MAE_measure_FWHM')
        print('Max radius: ', max_radius)
        print('NOTE: following params determine by MAE_radial_profile method')
        print('Radius measurements: ', radius)
        print('Radial profile plot: ', radial_profile)
        print('Enclosed energy plot: ', EE)
        print('Flux level values: ', EE_rad)
        print('Flux radii: ', levels)
        print('centroid (x,y): ', x, y)
        print('pixel scale: ', pixel_scale)
        print('Obs wavelength: ', wavelength)
        print('Interpreter kind: ', kind)
        print('Radial profile peak: ', rp_max)
        print('FWHM (arcsec/pixel): ', fwhm_arcsec, fwhm_pixel)

    return fwhm_arcsec, fwhm_pixel

def MAE_measure_FWHM_gaussian(HDUlist=None, ext=None, threshold = None, print_params = False, plot = False, centroid=None, wave=None):
    '''2022 Mar 16: A modified version of the strehl measuring tool developed by Marshall Perrin at STScI.
    This code measures the FWHM by fitting a 1D Gaussian to the data.

    Code taken/modified from pydatatut.pdf and poppy.utils.
    Copyright (C) 2010-2021 Association of Universities for Research in
    Astronomy (AURA)

    Redistribution and use in source and binary forms, with or without modification, are permitted provided that the
    following conditions are met:

    Redistributions of source code must retain the above copyright notice, this list of conditions and the following
    disclaimer.

    Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following
    disclaimer in the documentation and/or other materials provided with the distribution.

    The name of AURA and its representatives may not be used to endorse or promote products derived from this software
    without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY AURA “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL AURA BE
     LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
     LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
     HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
      OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

     Parameters
     __________
     HDUList:                    fits.HDUlist
                                 FITS file containing image to display
     ext:                        int
                                 FITS extension, default = 0
    threshold:              float
                                  Value below/above which to measure the FWHM
    plot:                        boolean
                                   Plot the resultant Gaussian fit
    print_params:               boolean
                                OPTIONAL: print the results of each step for testing purposes
    Returns
    __________
    FWHM:                       float
                                 The FWHM of the image, computed in both arcsec and pixel
    centroid:                float
                                 The centroid position (x,y) of the object
    '''

    # Read in data
    if HDUlist == HDUlist:
        HDUlist = fits.open(HDUlist)
    else:
        raise ValueError('Input must be a FITS HDUlist')

    # read the user-defined extension
    if ext == ext:
        data = HDUlist[ext].data
    else:
        # assume extension is 0
        data = HDUlist[0].data

    # read the pixel scale and observed wavelength
    #pixel_scale = HDUlist[ext].header['PIXELSCL']
    #pixel_scale = float(pixel_scale)
    pixel_scale = 0.11

    #wavelength = HDUlist[ext].header['WAVELEN']
    #wavelength = float(wavelength)
    wavelength = wave

    # prepare an array with radius in arcsec
    y, x = np.indices(data.shape, dtype=float)

    if centroid == '1dg':
        # find the center of a 1D-Gaussian fit to the data
        center = centroid_1dg(data)

    elif centroid == 'quad':
        center = centroid_quadratic(data)

    elif centroid == 'exact':
        # set the centroid as the center of the image array
        y, x = data.shape
        im_y = int(x / 2)
        im_x = int(x / 2)
        center = [(im_x), (im_y)]
        print(center)

    elif centroid == centroid:
        center = centroid

    else:
        raise ValueError('No centroiding method defined, options are: 1) fit, 2) exact')

    # save as float for centroid testing
    cenx = center[0]
    ceny = center[1]

    # bound the fitting box: len(image_data(x)) - center_coordinate(x)
    # i.e. (x, y) = (180-90, 180-90) = (90,90)
    y -= float(ceny)
    x -= float(cenx)

    # set the radius in arcsec
    r = np.sqrt(x ** 2 + y ** 2) * pixel_scale

    # select pixels above the threshold
    # NOTE: image is normalized to a peak = 1 above
    wpeak = np.where(data > threshold)

    # determine the radius peak
    r_peak = r[wpeak]

    # determine the image data peak
    im_peak = data[wpeak]

    # Determine the best fit Gaussian parameters for normalized inputs
    # Determine the diffraction limit: diameter of JWST primary mirror
    D = 6.5
    diff_limit = ((1.22*(wavelength * 1e-6)) / D) * 206265

    # Determine the  standard deviation:
    stddev = diff_limit / 2 * np.sqrt(2*np.log(2))

    # Determine the amplitude:
    amplitude = 1 / (stddev*np.sqrt(2*np.pi))

    # Create 1D-Gaussian fit
    #g_init = models.Gaussian1D(amplitude=amplitude, mean=0, stddev=stddev)
    g_init = models.Gaussian1D()
    g_init.mean.fixed = True

    # Fit the Gaussian model to the radial and image peak
    fit_g = fitting.LevMarLSQFitter()
    g = fit_g(g_init, r_peak, im_peak)

    # convert the results from sigma parameter to FWHM
    # NOTE: astropy fitting does not constrain the stddev to be positive only for some reason, take the absolute value here
    fwhm = 2 * np.sqrt(2 * np.log(2)) * np.abs(g.stddev)

    # Save FWHM measure, and convert to pixels based on input pixel scale
    # Original (arcsec)
    fwhm_arcsec = fwhm

    # Pixel scale
    fwhm_pixel = fwhm_arcsec / pixel_scale

    if print_params:
        print('Method: MAE_measure_FWHM_gaussian')
        print('Observed image shape: ', data.shape)
        print('Pixel scale: ', pixel_scale)
        print('Observed wavelength: ', wavelength)
        print('Image centroid (x,y): ', cenx, ceny)
        print('Diffraction limit (used for Gaussian fit): ', diff_limit)
        print('Standard deviation (used for Gaussian fit): ', stddev)
        print('Amplitude (used for Gaussian fit): ', amplitude)
        print('fwhm (arcsec/pixel): ', fwhm_arcsec, fwhm_pixel)

    if plot:
        plt.loglog(r_peak, im_peak, linestyle='none', marker='o', alpha=0.5)
        rmin = r_peak[r_peak != 0].min()
        plotr = np.linspace(rmin, r_peak.max(), 30)

        plt.title('1D-Gaussian fit')
        plt.plot(plotr, g(plotr))
        plt.xlabel("Radius [arcsec]")
        plt.ylabel("Intensity relative to peak")

        plt.axhline(0.5, ls=":")
        plt.axvline(fwhm / 2, ls=':')
        plt.text(0.1, 0.2, 'FWHM={:.3f} arcsec'.format(fwhm_arcsec), transform=plt.gca().transAxes, )

        plt.gca().set_ylim(threshold * .5, 2)
        plt.show()

    return fwhm_arcsec, fwhm_pixel, cenx, ceny

def MAE_measure_FWHM_2Dgaussian(HDUlist = None, ext = None, plot=False, centroid=None, fit = None):
    '''2022 Mar 29: Fit a 2D-Gaussian to the data and measure the FWHM of this fit. Optionally, plot the input image
    data and the resultant fit. The Gaussian is fit by tacking creating a small box around the centroid of the image
    (used-defined) then fitting a Gaussian to this cut-out.

    Parameters
    __________
    data: numpy array
        Normalized 2D-array containing the image data
    pixel_scale: float
        pixel scale of the 2d-image array, can be found in FITS header
    plot: boolean
        option to return a side-by-side plot comparing the input image and reference 2D-Gaussian fit
    threshold: float
        Value above which the FWHM is measured
    filter: string
        enter the filter the observation was made in, used to calculate the theoretical FWHM for DAOStarFinder. Options include:
        F560W:  image observed at 5.6 μm
        F1000W: image observed at 10 μm
        F1500W: image observed at 15 μm
        F1800W: image observed at 18 μm
        F2100W: image observed at 21 μm
        2022 Mar 29: NOTE: find a more robust way of doing this
    centroid: string
        set the centroiding method to determine the centroid of the image. Options include
        1dg:        use centroid_1dg
        2dg:        use centroid_2dg
        com:        use centroid_com
        quadratic:  use centroid_quadratic
        DAO:        use DAOStarFinder
    fit: int
        set the bounding size of the box fit to the centroid

    Returns
    __________
    FWHM: float
        The FWHM of the image, computed in both arcsec and pixel
    centroid: float
        The centroid position (x,y) of the object
    fit: 2D numpy array
        save model fit for plotting outside of the method
    '''
    # Read in data
    # Read in data
    if HDUlist == HDUlist:
        HDUlist = fits.open(HDUlist)
    else:
        raise ValueError('Input must be a FITS HDUlist')

    # read the user-defined extension
    if ext == ext:
        data = HDUlist[ext].data
    else:
        # assume extension is 0
        data = HDUlist[0].data

    # read the pixel scale and observed wavelength
    pixel_scale = HDUlist[ext].header['PIXELSCL']
    pixel_scale = float(pixel_scale)

    wavelength = HDUlist[ext].header['WAVELEN']
    wavelength = float(wavelength)

    # prepare an array with radius in arcsec
    y, x = np.indices(data.shape, dtype=float)

    if centroid == '1dg':
        # find the center of a 1D-Gaussian fit to the data
        center = centroid_1dg(data)

    elif centroid == 'quad':
        center = centroid_quadratic(data)

    elif centroid == 'exact':
        # set the centroid as the center of the image array
        y, x = data.shape
        im_y = int(x / 2)
        im_x = int(x / 2)
        center = [(im_x), (im_y)]
        print(center)

    elif centroid == centroid:
        center = centroid

    else:
        raise ValueError('No centroiding method defined, options are: 1) fit, 2) exact')

    # determine size of the Gaussian fitting box
    if fit == None:
        raise Warning('WARNING: no fit given for 2D-Gaussian')
    else:
        bb = fit

    # save as float for centroid testing
    cenx = center[0]
    ceny = center[1]

    # bound the fitting box: len(image_data(x)) - center_coordinate(x)
    # i.e. (x, y) = (180-90, 180-90) = (90,90)
    y -= float(ceny)
    x -= float(cenx)

    # return the centroid locations
    Y = '{:.6f}'.format(center[1])
    X = '{:.6f}'.format(center[0])

    # bound the centroid for the 2D-fit later
    yc = int(center[1])
    xc = int(center[0])

    # Determine the best fit Gaussian parameters
    # set the size of the box which fits the 2D-Gaussian to the centroid
    # re-save the image shape
    box = data[yc - bb:yc + bb, xc - bb:xc + bb]
    #square[yc - bb:yc + bb, xc - bb:xc + bb] += box
    yp, xp = box.shape
    yg, xg = np.mgrid[:yp, :xp]

    # Determine the best fit Gaussian parameters for normalized inputs
    # Determine the diffraction limit: diameter of JWST primary mirror
    D = 6.5
    diff_limit = ((1.22*(wavelength * 1e-6)) / D) * 206265

    # Determine the  standard deviation:
    stddev = diff_limit / 2 * np.sqrt(2*np.log(2))

    # Determine the amplitude:
    amplitude = 1 / (stddev*np.sqrt(2*np.pi))

    # fit 2D-Gaussian
    model = models.Gaussian2D(amplitude = amplitude, x_mean = 0, y_mean = 0, x_stddev = stddev, y_stddev = stddev)
    #model = models.Gaussian2D()
    fit_f = fitting.LevMarLSQFitter()
    f = fit_f(model, xg, yg, box)

    # save 2D-fit to numpy array
    fit_A = np.array(f(xg, yg), dtype='>f8')

    # return FWHM measurements (x,y) in pixels
    x_fwhmP = f.x_fwhm
    y_fwhmP = f.y_fwhm

    # convert the FWHM measurements to arcsec/px
    x_fwhmA = x_fwhmP * pixel_scale
    y_fwhmA = y_fwhmP * pixel_scale

    # model fit
    # resizing functions to match fit size to image size
    yf = int(len(fit_A[1]) / 2)
    xf = int(len(fit_A[0]) / 2)
    #square[yc - yf:yc + yf, xc - xf:xc + xf] += fit_A
    #print(square.shape)

    if plot:
        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10, 8), squeeze=True)
        fig.tight_layout(pad=4.5)

        # normalize and log scale data and fit
        norm_im = simple_norm(data, 'log')
        norm_fit = simple_norm(fit_A, 'log')

        # input image
        axes[0].imshow(data, origin='lower', norm=norm_im)
        axes[0].set_xlabel('X Position (px)')
        axes[0].set_ylabel('Y Position (px)')
        axes[0].set_title('Input Image | filter: ' + filter)

        # set fit coordinates for text box
        y = int(len(fit_A[1])) - 50
        x = 5

        # format fwhm measurements
        arc = '{:.3f}'.format(x_fwhmA)
        pix = '{:.3f}'.format(x_fwhmP)

        axes[1].imshow(fit_A, origin='lower', norm=norm_fit)
        axes[1].set_xlabel('X Position (px)')
        axes[1].set_ylabel('Y Position (px)')
        axes[1].set_title('Model fit | fit size: ' + str(bb))
        axes[1].text(x,y,'\u0332'.join('FWHM ') + '\narcsec: ' + str(arc) + '\npixel: ' + str(pix), color = 'r',
                        fontsize = 10)

        plt.show()

    return x_fwhmP, y_fwhmP, x_fwhmA, y_fwhmA, X, Y, fit_A

# ****************************************** Flux measuring tools ******************************************************#

def MAE_Aperture_Photometry(HDUlist=None, ext=0, fwhm=None, radius = None, normalize = False, centroid = None,
                            print_params = False, plot=False, save= False, display = False, title = None, trim = None,
                            aperture = None, semi_major = None, semi_minor = None, theta = None, pad = None):
    '''8 May 2022: Aperture Photometry tool used to measure flux in
    2D_FITS image arrays.
    Parameters
    __________
    HDUList:                    fits.HDUlist
                                FITS file containing image to display
    ext:                        int
                                FITS extension, default = 0
    fwhm:                       string
                                Select the method to measure the FWHM of the image. The FWHM sets the radial size of the annulus, options are:
                                1) '1D-Gaussian'- measure the FWHM by fitting a 1D-Gaussian along the x-axis of the objects centroid
                                2) 'Radial'- measure the FWHM from the image radial profile
    mult_fwhm:                  string
                                Set the multiplier of the FWHM
    normalize:                  boolean
                                option to normalize the data to a peak pixel value of 1
    centroid:                   string
                                Set the method to determine objects centroid. Options are:
                                1) '1dg' - determine the centroid of the object using a 1D-Gaussian fit
                                2) 'quad'- determine the centroid of the object using a quadratic fit
                                3) 'exact'- set the objects centroid as the image array center
                                4) given - centroid given by user
    print_params:               boolean
                                OPTIONAL: print the outputs of each step for trouble-shooting
    plot:                       boolean
                                Option to plot the aperture and annulus fits
    save:                       boolean
                                OPTIONAL: save the aperture fit
    display:                    boolean
                                OPTIONAL: display the aperture fit
    trim                        int
                                OPTIONAL: trim the data array size by this value
    aperture:                   string
                                OPTIONAL: set the aperture used to measure photometry, default is CircularAperture
    semi_major:                 float
                                OPTIONAL: (only used when Elliptical aperture is defined) the semi-major axis of the
                                ellipse in pixels
    semi_minor:                 float
                                OPTIONAL: (only used when Elliptical aperture is defined) the semi-minor axis of the
                                ellipse in pixels
    theta:                      float
                                OPTIONAL: (only used when Elliptical aperture is defined) the rotation angle as an
                                angular quantity or value in radians (as a float) from the positive x axis. The rotation
                                angle increases counter clockwise
    pad:                        float
                                pad the centroid location for aperture fitting
    Returns
    __________
    Total flux:                 float
                                Total flux within the aperture
    Background subtracted flux: float
                                Total flux measured after subtraction of annulus flux from aperture
    Aperture radius:            float
                                radial size of the aperture (arcsec/pixel) to return 95% of the flux within the object
    centroid position:          float
                                Position of the centroid fit
    '''

    # Read in data
    if HDUlist == HDUlist:
        HDUlist1 = HDUlist
        HDUlist = fits.open(HDUlist)
    else:
        raise ValueError('Input must be a FITS HDUlist')


    # read the user-defined extension
    if ext == ext:
        data = HDUlist[ext].data
    else:
        # assume extension is 0
        data = HDUlist[0].data

    # option to normalize the data
    if normalize:
        data = data / data.max()
    else:
        data = data

    if centroid == '1dg':
        # save centroid option for FWHM measurement
        centroid_fwhm = centroid

        # find the center of a 1D-Gaussian fit to the data
        centroid = centroid_1dg(data)

    elif centroid == 'quad':
        # save centroid option for FWHM measurement
        centroid_fwhm = centroid

        centroid = centroid_quadratic(data)

    elif centroid == 'exact':
        # save centroid option for FWHM measurement
        centroid_fwhm = centroid

        # set the centroid as the center of the image array
        y, x = data.shape
        im_y = int(x / 2)
        im_x = int(x / 2)
        centroid = [(im_x),(im_y)]

    elif centroid == centroid:
        # append for aperture measurement
        centroid = centroid

    else:
        raise ValueError('No centroiding method defined, options are: 1) fit, 2) exact')

    if radius == None:
        raise ValueError('Please specify the multiplier of the FWHM')
    else:
        radius=radius

    # read the pixel scale and observed wavelength
    #pixel_scale = HDUlist[ext].header['PIXELSCL']
    #pixel_scale = float(pixel_scale)
    pixel_scale = 0.11
    #wavelength = HDUlist[ext].header['WAVELEN']
    #wavelength = float(wavelength)
    wavelength = 5
    # centroid set manually
    x = centroid[0]
    y = centroid[1]

    # find the total flux within the given aperture radius
    # Pad the centroid by a user-defined value? If no padding, pad = 0
    position = [(x-pad), (y-pad)]

    # convert to arcsec
    radius_arcsec = radius * pixel_scale

    # create a user-defined aperture, centered on the object centroid, and measure flux
    if aperture == 'Ellipse':
        aperture1 = EllipticalAperture(position, float(semi_major), float(semi_minor), float(theta))
    else:
        aperture1 = photutils_CircularAperture(position, r=radius)

    phot_table = aperture_photometry(data, aperture1, method = 'exact')
    phot_table['aperture_sum'].info.format = '%.8g'
    #aperture_sum = phot_table['aperture_sum'][0]/aperture1.area
    aperture_sum = phot_table['aperture_sum'][0]
    aper_total = float(aperture_sum)

    # Perform annulus photometry and remove the background from the total flux measured above
    # user-defined inner and outer annulus radius
    inner_rad = 5
    outer_rad = 15

    # Background estimation
    # Create the annulus aperture
    position1 = [x,y]
    annulus_aperture = CircularAnnulus(position1, r_in=inner_rad, r_out=outer_rad)

    # Simple mean within a circular annulus
    aperstats = ApertureStats(data, annulus_aperture)

    # mean local per-pixel background
    bkg_mean = aperstats.mean

    # The total background within the circular aperture is the mean local per-pixel background * the circular aperture
    # area. Find the aperture area using the same area over which photometry was performed.
    area = aperture1.area_overlap(data)

    # calculate the total backgroud within the circular aperture
    total_bkg = bkg_mean * area

    # calculate the background subtracted sum
    final_phot = aper_total - total_bkg

    # set the comparison value
    # total aperture flux
    compare_total = final_phot * 0.95

    # define new radius
    r = pixel_scale / 2

    # Finds the optimum radius to return 95% of the flux
    while True:
        aperture2 = photutils_CircularAperture(position, r)
        phot_table1 = aperture_photometry(data, aperture2, method = 'exact')
        phot_table1['aperture_sum'].info.format = '%.8g'
        aperture_sum1 = phot_table1['aperture_sum'][0]
        aper_95 = float(aperture_sum1)

        if aper_95 <= compare_total:
            # if 95% of the radius is not found, increase the aperture radius size by the image pixel scale
            r += pixel_scale
        else:
            break

    # Return the 95% aperture radial sizes in arcsec and pixels
    r_arcsec = r * pixel_scale
    r_pixel = r

    if title == None:
        title = 'Aperture photometry fit'
    else:
        title = title

    if print_params:
        print('Method: MAE_Aperture_Photometry')
        print('Observed image shape: ', data.shape)
        print('Pixel scale: ', pixel_scale)
        print('Observed wavelength: ', wavelength)
        print('Object FWHM (used in multiple to set the aperture radius): ', fwhm)
        print('centroid (x,y): ', x, y)
        print('Optimum aperture radius (arcsec/pixels): ', radius_arcsec, radius)
        print('Total flux: ', aper_total)
        print('Radius to return 95% flux (arcsec/pixel): ', r_arcsec, r_pixel)
        print('95% flux: ', aper_95)
        print('Annulus inner/outer radius (pixels): ', inner_rad, outer_rad)
        print('Background Estimation')
        print('Mean estimated background: ', bkg_mean)
        print('Aperture area: ', area)
        print('Estimated background: ', total_bkg)
        print('Background removed flux: ', final_phot)

    if plot:
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 7), tight_layout=True)
        aper_total = aper_total
        final_phot = final_phot

        # trim decimals
        total = '{:.3e}'.format(aper_total)
        BS_total = '{:.2e}'.format(final_phot)
        r_pixel = '{:.2e}'.format(r_pixel)

        # set the image positions
        positions = [0, 15, 25, 35, 49]
        labels = ['-25', '-10', '0', '10', '25']

        # set x,y limits on added text
        xlim = len(data[0]) - 10
        ylim = len(data[1]) - 10

        ap_patches = aperture1.plot(color='w', lw=1)
        #ap_patches = aperture1.plot(color='w', lw=1, label='Total photometry aperture\n (r=' + str(radius) + '*FWHM)')
        #ap_patches1 = aperture2.plot(color='r', lw=1, label='95% photometry aperture')
        #ap_patches3 = annulus_aperture.plot(color='y', lw=1, label='Background annulus')


        norm = ImageNormalize(stretch=LogStretch(), vmin=0, vmax = data.max())
        im = plt.imshow(data, origin='lower', cmap= 'RdYlBu_r', norm=norm)
        plt.title(title)
        ax.xaxis.set_major_locator(ticker.FixedLocator(positions))
        ax.xaxis.set_major_formatter(ticker.FixedFormatter(labels))
        ax.yaxis.set_major_locator(ticker.FixedLocator(positions))
        ax.yaxis.set_major_formatter(ticker.FixedFormatter(labels))
        plt.ylabel('Y (pixels)')
        plt.xlabel('X (pixels)')
        plt.grid(False)
        add_colorbar(im, label='Log-scaled Flux Density (mJy)', format=FormatStrFormatter('%.1e'))

    if save:
        plt.savefig(save + '.png')

    if display:
        plt.show()

    # return: total flux, background subtracted flux, 95% flux, the estimated background, FWHM, optimum aperture radius
    # (arcsec/pixel), 95% optimum aperture radius (arcsec/pixel), annulus radius (inner/outer)
    # centroid position (x,y)
    return aper_total, final_phot, aper_95, total_bkg, fwhm, radius_arcsec, radius, r_arcsec, r_pixel, inner_rad, outer_rad, x, y

#************************************ Strehl ratio measuring tools ****************************************************#

def MAE_calculate_PSF(pupil=None, effwave=None, pixel_scale = None, filter = None, over_sampling = None, outfile = None,
                  outdir = None, plot = False, print_params=False):
    '''2022 May 12: Method used to calculate the perfect diffraction-limited monochromatic point spread function (PSF).
    Parameters
    __________
    pupil:         string
                   User-defined telescope pupil image. Options are:
                   1) JWST
    effwave:       string
                   Observation wavelength
    pixel_scale:   string
                   Pixel scale of the observed image
    filter:        string
                   Filter used to make the observations
    over_sampling: int
                   Factor by which to over/under-sample the theoretical PSF
    outfile:       String
                   Path specifying where to save the trimmed calculated PSF FITS file
    outdir:        string
                   Set the ouput directory for the generated FITS file
    plot:          boolean
                   OPTIONAL: plot the generated PSF
    print_params:  boolean
                   OPTIONAL: print the results of each step for trouble-shooting and testing purposes
    Returns
    __________
    psf:           numpy.array
                   2D-numpy array containing the theoretical PSF image data
    '''

    # set the pupil image selection
    if pupil == None:
        raise ValueError('No pupil image specified. Option(s) include: 1) JWST, 2) ')
    elif pupil == pupil:
        pupil = pupil

    # Read in the pupil image
    # This code can be modified to read in other telescope pupil images as well
    if pupil == 'JWST':
        # Read-in the JWST pupil image
        pupil_im = get_pkg_data_filename('JWST_pupil.fits')
        HDUlist = fits.open(pupil_im)
        pupil = HDUlist[0].data
        print(pupil.shape)
        HDUlist.close()

        # set the pixel scale (in m)
        pupil_scale = HDUlist[0].header['PIXSCALE']
        pupil_scale = float(pupil_scale)

        # set the x-axis length: naxis1 -> x_axis
        naxis1 = HDUlist[0].header['NAXIS1']

        # set the y-axis length: n_axis2 -> y_axis
        naxis2 = HDUlist[0].header['NAXIS2']

    # Set the effective wavelength
    if effwave == None:
        raise ValueError('No effective wavelength given')
    else:
        effwave = float(effwave)

    # Set the effective wavelength
    if pixel_scale == None:
        raise ValueError('No effective wavelength given')
    else:
        or_pixel_scale = float(pixel_scale)

    # set the observation filter
    if filter == None:
        raise ValueError('No observation filter specified')
    else:
        filter = filter

    # Set the oversampling factor
    # NOTE: the following values are based on the MIRISim pixel scale (0.1110)
    # JWST PSF only: over_sampling = 4
    # JWST PSF observed by MIRI: over_sampling = 1
    if over_sampling == None:
        raise ValueError('No over-sampling factor given')
    else:
        over_sampling = int(over_sampling)

    # set the outfile directory path
    if outfile == None:
        raise ValueError('No output name given')
    else:
        outfile = outfile

    # set the outfile directory path
    if outdir == None:
        raise ValueError('No output directory given')
    else:
        outdir = outdir

    # Calculate the total extent of the pupil needed to Fast-Fourier Transform (FFT)
    # set conversion factor: radius -> arcsec ~206265
    rad2arcsec = 206264.8062409636

    # Find the extent (in m)
    extent_in_m = ((over_sampling * (effwave * 1e-6)) / or_pixel_scale) * rad2arcsec


    # Set the number of pixels for FFT
    npixfft = 2 * (0.5 * extent_in_m) / pupil_scale
    npixfft = int(npixfft)

    # Create a 2D-array of floats with the value npixfft
    exit_list = (npixfft, npixfft)
    exit_pupil = np.zeros(exit_list)

    # subtract (x_axis, y_axis) from exit_pupil to set pupil size
    exit_pupil[0:naxis1, 0:naxis2] = pupil

    # compute the psf
    compute_oversamp_psf = fftshift(np.abs(fftn(exit_pupil)**2))

    # Ensure the right number of pixels is given for re-binning
    n_impix = npixfft - (npixfft % over_sampling)

    # # re-bin the theoretical psf
    rebin_size = int(n_impix / over_sampling)
    psf = bin_ndarray(compute_oversamp_psf[0:n_impix, 0:n_impix], new_shape = (rebin_size, rebin_size), operation = 'sum')

    # return the sum of the elements in the array
    psf /= np.sum(psf)
    # Determine the size of the psf
    x, y = psf.shape
    # find the center of the image
    cenx, ceny = int(round(x/2)), int(round(y/2))
    # trim the PSF to a 256x256 array
    trim_psf = psf[ceny-512:ceny+512, cenx-512:cenx+512]
    # update the pixel_scale
    pixel_scale = or_pixel_scale / over_sampling

    # Save all necessary FITS header information and save the trimmed PSF to a FITS file
    save_fits = save_FITS(trim_psf, filename='calcPSF',
                             name=outfile,
                             output_dir=outdir,
                             pixel_scale=pixel_scale,
                             filter=filter,
                             wavelength=effwave,
                             oversample = over_sampling)

    if print_params:
        print('Method: MAE_calculate_PSF')
        print('Outputs in order\n')
        print('pupil scale: ', pupil_scale)
        print('x-axis length: ', naxis1)
        print('y-axis length: ', naxis2)
        print('effective wavelength: ', effwave)
        print('pixel scale: ', or_pixel_scale)
        print('over-sampling factor: ', over_sampling)
        print('extent_in_m [extent(m)]: ', extent_in_m)
        print('npixfft [number of pixels for fft]: ', npixfft)
        print('exit_pupil [shape before fft]: ', exit_pupil.shape)
        print('compute_oversamp_psf: ', compute_oversamp_psf.shape)
        print('n_impix: ', n_impix)
        print('Final re-binned psf shape: ', psf.shape)
        print('Trimmed PSF shape: ', trim_psf.shape)
        print('Over-sampled pixel scale: ', pixel_scale)

    if plot:
        # find the (x,y) limit of the trimmed image
        x, y = trim_psf.shape
        # plot the generated PSF image with relevant information
        norm = simple_norm(trim_psf, 'log')
        plt.figure(figsize = (6,5))
        im = plt.imshow(trim_psf, origin = 'lower', cmap = 'viridis', norm = norm)
        plt.title(str(effwave)+' μm theoretical PSF', fontsize = 20)
        plt.ylabel('Pixel column')
        plt.xlabel('Pixel row')
        plt.text(10, y - 15, 'Wavelength: ' + str(effwave) + ' μm', color='w', fontsize=10)
        plt.text(10, y - 30, 'Oversampling: ' + str(over_sampling), color='w', fontsize=10)
        plt.text(10, y - 45, 'Pixel scale: ' + str(pixel_scale), color='w', fontsize=10)
        plt.colorbar(im)
        plt.show()

    # return the theoretical psf
    return psf, trim_psf, save_fits

def calculate_Strehl(HDUlist=None, ext=0, calculate_PSF = None, pupil=None, over_sampling=None, filter = None,
                     normalize=False, centroid = None, psf_centroid = None, fwhm = None, mult_fwhm = None, outdir = None,
                     outfile = None, display_psf = False, plot=False,  print_params=False,  photom_rad = None,):
    '''2022 May 12: Method used to compute the Strehl ratio of an image, based on a calculated perfect diffraction-limited
    monochromatic point spread function (PSF). Code adapted from IDL to Python, base code provided by Marcos van Dam (Gemini).
    Parameters
    __________
    HDUlist:               fits.HDUlist
                           FITS file containing image to display
    ext:                   int
                           FITS extension, default = 0
    calculate_PSF:         string
                           Set the method used to calculate the theoretical reference PSF. Options are:
                           1) 'calcPSF'- Generates the theoretical PSF using the MAE_calculate_PSF method
                           2) 'WebbPSF'- Generates the theoretical PSF using WebbPSF
                           3) fits.HDUlist- user-entered FITS file containing the PSF
    pupil:                 string
                           Set the telescope pupil image to calculate the theoretical PSF. Options are (can be added):
                           1) JWST
    over_sampling:         int
                           Factor by which to over-sample the theoretical PSF
    filter:                string
                           Specify the filter used to make the observation, used for calculating theoretical PSF using
                           STScI method
    normalize:             boolean
                           OPTIONAL: normalize the data to a peak pixel value of 1
    centroid:              string
                           Set the method to determine objects centroid. Options are:
                           1) given- centroid coordinates given by user
                           2) '1dg'- determine the centroid of the object using a 1D-Gaussian fit
                           3) 'exact'- set the objects centroid as the image array center
    psf_centroid:          string
                           Set the method to determine PSFs centroid. Options are:
                           1) '1dg'- determine the centroid of the object using a 1D-Gaussian fit
                           2) 'exact'- set the objects centroid as the image array center
    fwhm:                  string
                            Used for the MAE_Aperture_Photometry tool
                            Select the method to measure the FWHM of the image. The FWHM sets the radial size of the annulus, options are:
                            1) '1D-Gaussian'- measure the FWHM by fitting a 1D-Gaussian along the x-axis of the objects centroid
                            2) 'Radial'- measure the FWHM from the image radial profile
    mult_fwhm:             string
                           Used for the MAE_Aperture_Photometry tool
                            Set the multiplier of the FWHM
    outdir:                string
                           Set the output directory for the theoretical PSF FITS file to be stored in
    outfile:               string
                           Set the name of the theoretical PSF
    display_psf:           boolean
                           OPTIONAL: display the PSF generated by WebbPSF
    plot:                  boolean
                           OPTIONAL: plot a side-by-side comparison of the reference image and the theoretical PSF
    print_params:          boolean
                           OPTIONAL: print the outputs of each step for trouble-shooting
    photom_rad:            string
                           OPTIONAL: Aperture photometry radial value, defined in arcsec
    Returns
    __________
    Strehl:               float
                          Strehl ratio of the observed image to the calculated PSF. Should be a float between 0.0 - 1.0
    '''

    # Read in image
    if HDUlist == HDUlist:
        # set HDUlist for aperture photometry only
        HDUlist1 = HDUlist

        # open the image data
        HDUlist = fits.open(HDUlist)
    else:
        raise ValueError('Input must be a FITS HDUlist')

    # read the user-defined extension
    if ext == ext:
        data = HDUlist[ext].data
    else:
        # assume extension is 0
        data = HDUlist[0].data

    # Set the pupil image for PSF calculations
    if pupil == None:
        raise ValueError('No telescope pupil defined')
    else:
        pupil2 = pupil

    # set the pixel scale
    pixel_scale = HDUlist[ext].header['PIXELSCL']
    pixel_scale = float(pixel_scale)

    # set the observed wavelength
    wavelength = HDUlist[ext].header['WAVELEN']
    wavelength = float(wavelength)

    # set the zise of the theoretical PSF, based on the axis size of the observed image
    axis = HDUlist[ext].header['NAXIS1']
    axis = int(axis)

    # set the over-sampling factor
    # NOTE: oversample = 1 -> sets the output image to MIRIs plate scale
    if over_sampling == None:
        raise ValueError('No over-sampling factor given')
    else:
        over_sampling = over_sampling

    # set the observation filter
    if filter == None:
        raise ValueError('No observation filter specified')
    else:
        filter = filter

    # Set the output name for the theoretical PSF
    if outfile == None:
        raise ValueError('No output name specified for the theoretical PSF')
    else:
        outfile = outfile

    # Set the output directory for the theoretical PSF
    if outdir == None:
        raise ValueError('No output path specified for the theoretical PSF')
    else:
        outdir = outdir

    # calculate theoretical PSF
    if calculate_PSF == None:
        raise ValueError('No method specified to calculate the theoretical PSF: Options are: 1) FFT, 2) WebbPSF')

    elif calculate_PSF == 'calcPSF':
        # Use the calculate_PSF method
        name = calculate_PSF + '_'
        if not display_psf:
            MAE_calculate_PSF(pupil=pupil2, effwave=wavelength, pixel_scale = pixel_scale, filter = filter,
                                    over_sampling=over_sampling, outfile = outfile, outdir=outdir, plot = False,
                                    print_params=False)

            # read in reference PSF
            miri_psf = get_pkg_data_filename(outdir + name + outfile + '.fits')
            psf_list = fits.open(miri_psf)
            psf = psf_list[0].data
            psf_list.close()

        else:
            MAE_calculate_PSF(pupil=pupil2, effwave=wavelength, pixel_scale = pixel_scale, filter = filter,
                                    over_sampling=over_sampling, outfile = outfile, outdir=outdir, plot = True,
                                    print_params=False)

            # read in reference PSF
            miri_psf = get_pkg_data_filename(outdir + name + outfile + '.fits')
            psf_list = fits.open(miri_psf)
            psf = psf_list[0].data
            psf_list.close()

    elif calculate_PSF == 'WebbPSF':
        # Use the STScI calculate PSF method
        name = calculate_PSF + '_'
        miri = webbpsf.MIRI()
        miri.filter = filter

        # calculate the theoretical PSF
        miri.calc_psf(fov_pixels=axis, oversample=over_sampling, display=False, outfile= outdir + name +
                                                                                                   outfile + '.fits')

        # read in reference PSF
        miri_psf = get_pkg_data_filename(outdir + name + outfile + '.fits')
        psf_list = fits.open(miri_psf)
        psf = psf_list[0].data
        psf_list.close()

        if display_psf:
            webbpsf.display_psf(psf)
            plt.show()

    elif calculate_PSF == calculate_PSF:
        # read in reference PSF
        miri_psf = get_pkg_data_filename(calculate_PSF)
        psf_list = fits.open(miri_psf)
        psf = psf_list[0].data
        psf_list.close()


    # OPTIONAL: normalize the image data before measuring photometry
    if normalize == True:
        data = data / data.max()
        # Normalize the psf data
        psf = psf / psf.max()
        # set the normalization param
        normalize = True

    else:
        normalize = False

    if centroid == '1dg':
        # set the centroid for apertue photometry
        center = centroid

        # find the center of a 1D-Gaussian fit to the data
        centroid = centroid_1dg(data)

    elif centroid == 'quad':
        # set the centroid for apertue photometry
        center = centroid

        # determine the centroid for peak position
        centroid = centroid_quadratic(data)

    elif centroid == centroid:
        # set the centroid for apertue photometry
        center = centroid

        centroid = centroid

    # determine the peak of the image and PSF
    # centroid set manually
    im_x = int(centroid[0])
    im_y = int(centroid[1])

    # find the centroid of the theoretical PSF based on the method used to calculate the theoretical psf
    if psf_centroid == '1dg':
        # find the centroid using the 'centroid_2dg' method
        psf_center = centroid_1dg(psf)

        # set centroid for aperture photom
        psf_centroid = psf_centroid

    elif psf_centroid == 'exact':
        # set the centroid as the center of the image array
        y,x = psf.shape
        psf_y = int(x/2)
        psf_x = int(y/2)
        psf_center = [(psf_x),(psf_y)]

        # set centroid for aperture photom
        psf_centroid = psf_centroid

    else:
        raise ValueError('No centroiding method defined, options are: 1) fit, 2) exact')

    # set the (x,y) coordinate positions
    psf_y = int(psf_center[0])
    psf_x = int(psf_center[1])

    # Find the image peak (I_peak) and the PSF peak (PSF_peak)
    # image peak
    I_peak = data[im_y, im_x]

    # PSF peak
    # set the peak by the centorid method
    PSF_peak = psf[psf_y, psf_x]


    # Measure the flux in the image and PSF
    # measure the FWHM of the image
    if fwhm == None:
        raise ValueError('Please specify the method used to measure the FWHM of the image. Options are: 2) Radial, 2) 1D-Gaussian')

    elif fwhm == '1D-Gaussian':
        fwhm = fwhm

    elif fwhm == 'Radial':
        fwhm = fwhm

    if mult_fwhm == None:
        raise ValueError('Please specify the multiplier of the FWHM')

    else:
        mult_fwhm = mult_fwhm

    # set the psf fwhm
    # 2022 08 11: return and find a more robust way of mearuring this
    if wavelength == 5.6:
        # FWHM values given from STScI: https://jwst-docs.stsci.edu/jwst-mid-infrared-instrument/miri-performance/miri-point-spread-functions
        psf_fwhm = 1.636
    elif wavelength == 10:
        psf_fwhm = 2.888
    elif wavelength == 15:
        psf_fwhm = 4.354
    elif wavelength == 18:
        psf_fwhm = 5.224
    elif wavelength == 21:
        psf_fwhm = 5.989

    # measure the flux in the image using the MAE_Aperture_Photometry tool
    flux_im = MAE_Aperture_Photometry(HDUlist=HDUlist1,
                                      ext=0,
                                      centroid=center,
                                      fwhm=fwhm,
                                      radius=mult_fwhm,
                                      normalize=normalize,
                                    )
    I_flux = flux_im[1]

    # measure the flux in the PSF using the MAE_Aperture_Photometry tool
    flux_psf = MAE_Aperture_Photometry(HDUlist=miri_psf,
                                       ext=0,
                                       centroid=psf_centroid,
                                       fwhm=psf_fwhm,
                                       radius=mult_fwhm,
                                       normalize=normalize,
                                    )

    PSF_flux = flux_psf[1]

    # Compute the intensity of the image [I(x=0) = image peak / image flux = I_peak / I_flux]
    I_intensity = (I_peak / I_flux)

    # Compute the intensity of the image [PSF(x=0) = PSF peak / PSF flux = PSF_peak / PSF_flux]
    PSF_intensity = (PSF_peak / PSF_flux)

    # Calculate the Strehl ratio
    # Strehl = Image intensity / PSF intensity = I(x=0) / PSF(x=0)
    Strehl = (I_intensity / PSF_intensity)

    if plot:
        # plot the side-by-side image comparisons with the measured Strehl ratio
        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10, 5))
        fig.suptitle('Strehl ratio = {:.5f}'.format(Strehl), fontsize=20)
        fig.supxlabel('Pixel row')
        fig.supylabel('Pixel column')

        # normalize the image
        norm = simple_norm(data, 'log')

        # normalize the psf
        normpsf = simple_norm(psf, 'log')

        # plot image
        # cmap = 'gist_heat'
        axes[0].imshow(data, origin='lower', cmap = 'viridis', norm=norm)
        axes[0].set_title('Observed Image')

        # plot PSF
        axes[1].imshow(psf, origin='lower', cmap = 'viridis', norm=normpsf)
        axes[1].set_title('Reference PSF')
        plt.show()

    if print_params:
        print('Method: calculate_Strehl')
        print('Observed image shape: ', data.shape)
        print('Theoretical PSF shape: ', psf.shape)
        print('Pixel scale: ', pixel_scale)
        print('Observed wavelength: ', wavelength)
        print('Pupil choice: ', pupil)
        print('Observation filter specified: ', filter)
        print('Specifed output directory: ', outdir + outfile + '.fits')
        print('Over-sampling factor: ', over_sampling)
        print('Image centroid (x,y): ', im_x, im_y)
        print('PSF centroid (x,y)', psf_x, psf_y)
        print('Image flux: ', '{:.5f}'.format(I_flux))
        print('PSF flux: ', '{:.5f}'.format(PSF_flux))
        print('Image peak: ', I_peak)
        print('PSF peak: ', PSF_peak)
        print('Image intensity (I_peak/I_flux): ', '{:.5f}'.format(I_intensity))
        print('PSF intensity (PSF_peak/PSF_flux): ', '{:.5f}'.format(PSF_intensity))
        print('Strehl ratio (I_intensity/PSF_intensity): ', '{:.5f}'.format(Strehl))

    # raise an error if the Strehl ratio is ever greater than 1 (testing only)
    if Strehl > 1.00000:
         raise ValueError('Strehl ratio is >1, Strehl ratio = ' + '{:.5f}'.format(Strehl))

    # return the measured Strehl ratio
    return Strehl, I_peak, I_flux, I_intensity, PSF_peak, PSF_flux, PSF_intensity

#*************************************** Save FITS files tools ********************************************************#

def save_FITS(data, filename=None, name=None, output_dir=None, pixel_scale=float, filter=None, wavelength=float,
              instrument = None, oversample = None):
    '''2022 May 9: Method to save 2D-data array into a FITS file for different input formats
    Parameters
    __________
    data:         2D-data array
                  Array containing image data
    filename:     string
                  User-defined filename, used for saving convention. Options are:
                  1) 'observation'-               saves the observed image into a FITS file
                  2) 'deconvolution'-             saves a single deconvolved image into a FITS file
                  3) 'deconvolution_iteratation'- saves each deconvolved image as a function of iteration into a single FITS file
                  4) None-                        saves a standard FITS file
    name:         string
                  Name of the FITS file
    output_dir:   string
                  Directory specified to save the FITS file into
    pixel_scale:  float
                  User-defined plate scale of the image
    filter:       string
                  Observation filter of the image
    wavelength:   float
                  Observed image wavelength
    instrument:   string
                  OPTIONAL: Specify the instrument the observation was made with
    oversample:   string
                  OPTIONAL: specify the over-sampling factor used to generate the theoretical PSF
    Returns
    __________
    FITS file:    HDUlist
                  Formatted FITS file containgin necessary headers and image data
'''

    suffix = '.fits'
    name = name
    output_dir = output_dir

    # Save the formatted observation FITS file
    if filename == 'observation':
        base_filename = filename
        outfile = os.path.join(output_dir, base_filename + '_' + name + suffix)
        hdu = fits.PrimaryHDU(data)

        # read in and set the image pixel scale
        if pixel_scale == None:
            raise ValueError('Pixel scale must be given')
        else:
            hdu.header['PIXELSCL'] = pixel_scale

        # read in and set the image observation filter
        if filter == None:
            raise ValueError('Filter must be given')
        else:
            hdu.header['FILTER'] = filter

        # read in and set the image observed wavelength
        if wavelength == None:
            raise ValueError('Observation wavelength must be given')
        else:
            hdu.header['WAVELEN'] = wavelength

        # read in and set the observation instrument
        if instrument == None:
            raise ValueError('Instrument used for observation must be given')
        else:
            hdu.header['INSTRUM'] = instrument

        # check if FITS file is created and, if it is, remove it
        if os.path.isfile(outfile):
            os.remove(outfile)

        hdu.writeto(outfile, overwrite=True)

    # Save the formatted residual (background-subtracted) observation FITS file
    elif filename == 'residual':
        base_filename = filename
        outfile = os.path.join(output_dir, base_filename + '_' + name + suffix)
        hdu = fits.PrimaryHDU(data)

        # read in and set the image pixel scale
        if pixel_scale == None:
            raise ValueError('Pixel scale must be given')
        else:
            hdu.header['PIXELSCL'] = pixel_scale

        # read in and set the image observation filter
        if filter == None:
            raise ValueError('Filter must be given')
        else:
            hdu.header['FILTER'] = filter

        # read in and set the image observed wavelength
        if wavelength == None:
            raise ValueError('Observation wavelength must be given')
        else:
            hdu.header['WAVELEN'] = wavelength

        # check if FITS file is created and, if it is, remove it
        if os.path.isfile(outfile):
            os.remove(outfile)

        # read in and set the observation instrument
        if instrument == None:
            raise ValueError('Instrument used for observation must be given')
        else:
            hdu.header['INSTRUM'] = instrument

        hdu.writeto(outfile, overwrite=True)

    # Save the pipeline reduced observation FITS file
    elif filename == 'pipeline':
        base_filename = filename
        outfile = os.path.join(output_dir, base_filename + '_' + name + suffix)
        hdu = fits.PrimaryHDU(data)

        # read in and set the image pixel scale
        if pixel_scale == None:
            raise ValueError('Pixel scale must be given')
        else:
            hdu.header['PIXELSCL'] = pixel_scale

        # read in and set the image observation filter
        if filter == None:
            raise ValueError('Filter must be given')
        else:
            hdu.header['FILTER'] = filter

        # read in and set the image observed wavelength
        if wavelength == None:
            raise ValueError('Observation wavelength must be given')
        else:
            hdu.header['WAVELEN'] = wavelength

        # check if FITS file is created and, if it is, remove it
        if os.path.isfile(outfile):
            os.remove(outfile)

        # read in and set the observation instrument
        if instrument == None:
            raise ValueError('Instrument used for observation must be given')
        else:
            hdu.header['INSTRUM'] = instrument

        hdu.writeto(outfile, overwrite=True)

    # Save the formatted deconvolved FITS file
    elif filename == 'deconvolution':
        base_filename = filename
        outfile = os.path.join(output_dir, base_filename + '_' + name + suffix)
        hdu = fits.HDUList()

        count = 0
        for i in data:

            if count == 0:

                # append the observed image
                hdu.append(fits.ImageHDU(data[count]))

                # read in and set the image pixel scale
                if pixel_scale == None:
                    raise ValueError('Pixel scale must be given')
                else:
                    hdu[count].header['PIXELSCL'] = pixel_scale

                # read in and set the image observation filter
                if filter == None:
                    raise ValueError('Filter must be given')
                else:
                    hdu[count].header['FILTER'] = str(filter)

                # read in and set the image observed wavelength
                if wavelength == None:
                    raise ValueError('Observation wavelength must be given')
                else:
                    hdu[count].header['WAVELEN'] = wavelength

                # read in and set the observation instrument
                if instrument == None:
                    raise ValueError('Instrument used for observation must be given')
                else:
                    hdu[count].header['INSTRUM'] = instrument

                    hdu.writeto(outfile, overwrite=True)
            else:
                # append each deconvolved image
                hdu.append(fits.ImageHDU(data[count]))

            count += 1

        hdu.writeto(outfile, overwrite=True)

    # Save the formatted calcPSF generated theoretical PSF FITS file
    elif filename == 'calcPSF':
        base_filename = filename
        outfile = os.path.join(output_dir, base_filename + '_' + name + suffix)
        hdu = fits.PrimaryHDU(data)

        # read in and set the image pixel scale
        if pixel_scale == None:
            raise ValueError('Pixel scale must be given')
        else:
            hdu.header['PIXELSCL'] = pixel_scale

        # read in and set the image observation filter
        if filter == None:
            raise ValueError('Filter must be given')
        else:
            hdu.header['FILTER'] = str(filter)

        # read in and set the image observed wavelength
        if wavelength == None:
            raise ValueError('Observation wavelength must be given')
        else:
            hdu.header['WAVELEN'] = wavelength

        # check if FITS file is created and, if it is, remove it
        if os.path.isfile(outfile):
            os.remove(outfile)

        # read in and set the observation instrument
        if oversample == None:
            raise ValueError('Over sampling factor used for theoretical PSF generation must be given')
        else:
            hdu.header['OVERSMP'] = oversample

        hdu.writeto(outfile, overwrite=True)

    else:
        outfile = os.path.join(output_dir, name + suffix)
        hdu = fits.PrimaryHDU(data)
        hdu.writeto(outfile, overwrite=True)

#******************************************* Various Plotting tools **************************************************#

def show_Image(data, title = None, text = None, xlim = None, ylim = None, normalize = False, cmap = None,
               labels = None, image_bounds = None, aperture=None, convert = None, compass = False,
               print_params = None, display = False, save = None):
    '''2022 May 24: Method used to display a single log-scaled image with the title, flux,
    FWHM and Strehl ratio measurements. Assumes image data is given in pixels, has option to convert to desired units.
    Parameters
    __________
    data:           numpy.array
                    2D-numpy array containing the image data
    title:          string
                    User-defined title of the image
    text:           string
                    User-defined image caption
    xlim:           int
                    set the x-limit for text plotting
    ylim:           int
                    set y-limit for text plotting
    normalize:      boolean
                    OPTIONAL: normalize the input data
    cmap:           string
                    User-defined colorband for the image
    labels          numpy.array
                    Array containing user-defined x,y,and colormap labels
    image_bounds:   numpy.array
                    OPTIONAL: array containing the vmin and vmax values for image plotting manually
    aperture:       numpy.array
                    OPTIONAL: create a circular aperture for plotting purposes only
                    aperture[0] = x-centroid coordinates
                    aperture[1] = y-centroid coordinates
                    aperture[2] = aperture radius
    display:        boolean
                    Toggle whether to display the image or not
    compass:        Boolean
                    OPTIONAL: add a N through E compass to the image
    convert:        float
                    OPTIONAl: Image is always displayed in pixels unless toggled
                    Converts pixels -> arcsec based on pixel scale given
    save:           string
                    OPTIONAL: the path to save the output image in
    print_params:   boolean
                    OPTIONAL: print the outputs of each step for trouble-shooting
    Returns
    __________
    image:  numpy.array
            2D-image displaying the flux, FWHM and Strehl measurements
            '''

    # find the (x,y) limit of the trimmed image
    x, y = data.shape

    # set the title
    if title == None:
        raise ValueError('Please specify a title, no title given')
    else:
        title = title

    # plot the generated PSF image with relevant information
    fig, ax = plt.subplots(1, figsize=(7, 5), tight_layout=True)
    ax.set_aspect('equal')

    # # # set the min/max values?
    if image_bounds == image_bounds:
        min = image_bounds[0]
        max = image_bounds[1]
    else:
        min = None
        max = None

    # Normalize the data for plotting
    if normalize == None:
        raise ValueError('Please specify normalization param')
    elif normalize == 'SimpleNorm':
        normO = simple_norm(data, 'log', clip=True)
    elif normalize == 'LogNorm':
        normO = LogNorm(vmin = min, vmax = max, clip=True)
    elif normalize == 'Normalize':
        normO = Normalize(vmin=min, vmax=max)
    elif normalize == 'ImNorm':
        normO = ImageNormalize(stretch=LogStretch(), vmin=min, vmax=max)
    elif normalize == 'AsinhNorm':
        normO = ImageNormalize(stretch=AsinhStretch(a=0.001), vmin=min, vmax=max)

    if cmap == None:
        raise ValueError('Please specify the cmap used. Main options are: 1) viridis, 2) gist-heat, 3) etc.')
    else:
        cmap = cmap

    if labels == None:
        labels = None
    else:
        x_label = labels[0]
        y_label = labels[1]
        c_label = labels[2]

    # toggle whatever units to convert image data too
    if convert == None:
        convert = 1
    else:
        convert = convert

    # 20220830: come back and make the x/y-limit setting more robust
    im = plt.imshow(data, origin='lower', cmap=cmap, norm=normO,
                        extent = (0*convert, x*convert, 0*convert, y*convert))
    plt.title(title, fontsize=20)
    plt.ylabel(y_label)
    plt.xlabel(x_label)
    plt.text(xlim*convert, ylim*convert, str(text), color='w', fontsize=10)
    add_colorbar(im, label = c_label, format = FormatStrFormatter('%.2e'))

    if x_label == 'RA Offset (")':
        # set for large image
        if x == 150:
            positions = [5*convert, 30*convert, 50*convert, 75*convert, 100*convert, 120*convert, 145*convert]
            labels = ['16', '8', '4', '0', '-4', '-8', '-16']
            ax.xaxis.set_major_locator(ticker.FixedLocator(positions))
            ax.xaxis.set_major_formatter(ticker.FixedFormatter(labels))
            ax.yaxis.set_major_locator(ticker.FixedLocator(positions))
            ax.yaxis.set_major_formatter(ticker.FixedFormatter(labels))
        elif x == 120:
            positions = [10*convert, 16.5*convert, 33*convert, 60*convert, 77*convert, 93.5*convert, 110*convert]
            labels = ['12', '6', '3', '0', '-3', '-6', '-12']
            ax.xaxis.set_major_locator(ticker.FixedLocator(positions))
            ax.xaxis.set_major_formatter(ticker.FixedFormatter(labels))
            ax.yaxis.set_major_locator(ticker.FixedLocator(positions))
            ax.yaxis.set_major_formatter(ticker.FixedFormatter(labels))
        elif x == 50:
            positions = [5*convert, 15*convert, 25*convert, 35*convert, 45*convert]
            labels = ['5', '2', '0', '-2', '-5']
            ax.xaxis.set_major_locator(ticker.FixedLocator(positions))
            ax.xaxis.set_major_formatter(ticker.FixedFormatter(labels))
            ax.yaxis.set_major_locator(ticker.FixedLocator(positions))
            ax.yaxis.set_major_formatter(ticker.FixedFormatter(labels))
    else:
        if x == 150:
            ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
            ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
        elif x == 120:
            ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
            ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
        elif x == 50:
            ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
            ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))



    # Set aperture fit params
    if aperture:
        # convert pixel values to arcsecond
        x = aperture[0]*convert
        y = aperture[1]*convert
        rad = aperture[2]*convert

        # Circle((x,y-centroids), aperture radius)
        circle = plt.Circle((x, y), rad, fill=False, color='white')
        ax.add_patch(circle)

    if compass:
        # set the x,y image data
        x1,y1 = data.shape

        # 20220901: hard code params for now -> return and make more robust later
        if labels[0] == 'X (")':
            if x1 == 150:
                # size = [txt1-x, txt1-y, txt2-x, txt2-y, arrow1-x, arrow1-y, arrow1-xlen, arrow1-ylen, arrow2-x, arrow2-y, arrow2-xlen, arrow2-ylen]
                size = [136, 8, 130.5, 26.5, 125, 15, 4, 8, 125, 15, 8, -4]
            if x1 == 120:
                # size = [txt1-x, txt1-y, txt2-x, txt2-y, arrow1-x, arrow1-y, arrow1-xlen, arrow1-ylen, arrow2-x, arrow2-y, arrow2-xlen, arrow2-ylen]
                size = [106, 8, 100.5, 26.5, 95, 15, 4, 8, 95, 15, 8, -4]
            else:
                # size = [txt1-x, txt1-y, txt2-x, txt2-y, arrow1-x, arrow1-y, arrow1-xlen, arrow1-ylen, arrow2-x, arrow2-y, arrow2-xlen, arrow2-ylen]
                size = [45, 1, 43, 10, 40, 5, 2, 3.5, 40, 5, 3.5, -2]

        else:
            # 50 -> x,y*1
            if x1 == 150:
                # size = [txt1-x, txt1-y, txt2-x, txt2-y, arrow1-x, arrow1-y, arrow1-xlen, arrow1-ylen, arrow2-x, arrow2-y, arrow2-xlen, arrow2-ylen]
                size = [138, 25, 123, 8, 140, 10, 0, 10, 140, 10, -10, 0]
            if x1 == 120:
                # size = [txt1-x, txt1-y, txt2-x, txt2-y, arrow1-x, arrow1-y, arrow1-xlen, arrow1-ylen, arrow2-x, arrow2-y, arrow2-xlen, arrow2-ylen]
                size = [108, 25, 93, 8, 110, 10, 0, 10, 110, 10, -10, 0]
            else:
                # size = [txt1-x, txt1-y, txt2-x, txt2-y, arrow1-x, arrow1-y, arrow1-xlen, arrow1-ylen, arrow2-x, arrow2-y, arrow2-xlen, arrow2-ylen]
                size = [44.2, 12, 37.5, 4.5, 45, 5, 0, 4, 45, 5, -4, 0]



        # add text
        plt.text(size[0]*convert, size[1]*convert, 'N', color = 'w')
        plt.text(size[2]*convert, size[3]*convert, 'E', color = 'w')

        # create the height of the arrow
        arrow1 = plt.arrow(size[4]*convert, size[5]*convert, size[6]*convert, size[7]*convert, width = 0.04, color = 'w')

        # create the length of the arrow
        arrow2 = plt.arrow(size[8]*convert, size[9]*convert, size[10]*convert, size[11]*convert, width = 0.04, color = 'w')

        # add arrows
        ax.add_patch(arrow1)
        ax.add_patch(arrow2)

    if print_params:
        print('Method: show_Image')
        print('Observed image shape: ', data.shape)
        print('User-specified title: ', title)
        print('User-defined image caption: ', str(text))
        print('User-specified cmap: ', cmap)

    if save:
        plt.savefig(save + '.png')

    if display:
        plt.show()

def show_Double_Image(data, data2, title = None, subtitle1 = None, subtitle2 = None, text1 = None, text2 = None,
                      xlim = None, ylim = None, normalize = False, cmap = None, labels = None, image_bounds = None,
                      aperture1 = None, aperture2 = None, convert = None, compass = False, print_params = None,
                      display = False, save = None):
    '''2022 June 12: Method used to display a side-by-side log-scaled image comparison with the title, flux,
    FWHM and Strehl ratio measurements. Assumes image data is given in pixels, has option to convert to desired units.
    Parameters
    __________
    data:           numpy.array
                    2D-numpy array containing the image data
    data2:          numpy.array
                    2D-numpy array containing the image data
    title:          string
                    User-defined title of the image
    subtitle1:      string
                    First user-defined sub-title
    subtitle2:      string
                    Second user-defined sub-title
    text1:          string
                    First user-defined image caption
    text2:          string
                    Second user-defined image caption
    xlim:           int
                    set the x-limit for text plotting
    ylim:           int
                    set y-limit for text plotting
    normalize:      boolean
                    OPTIONAL: normalize the input data
    cmap:           string
                    User-defined colorband for the image
    labels          numpy.array
                    Array containing user-defined x,y,and colormap labels
    image_bounds:   numpy.array
                    OPTIONAL: array containing the vmin and vmax values for image plotting manually
    aperture1:      numpy.array
                    OPTIONAL: create a circular aperture for plotting purposes only -> 1st-image
                    aperture2[0] = x-centroid coordinates
                    aperture2[1] = y-centroid coordinates
                    aperture2[2] = aperture radius
    aperture2:      numpy.array
                    OPTIONAL: create a circular aperture for plotting purposes only -> 1st-image
                    aperture1[0] = x-centroid coordinates
                    aperture1[1] = y-centroid coordinates
                    aperture1[2] = aperture radius
    convert:        float
                    OPTIONAl: Image is always displayed in pixels unless toggled
                    Converts pixels -> arcsec based on pixel scale given
    compass:        Boolean
                    OPTIONAL: add a N through E compass to the image
    display:        boolean
                    Toggle whether to display the image or not
    save:           string
                    OPTIONAL: the path to save the output image in
    print_params:   boolean
                    OPTIONAL: print the outputs of each step for trouble-shooting
    Returns
    __________
    image:  numpy.array
            2D-image displaying the flux, FWHM and Strehl measurements
            '''

    # find the (x,y) limit of the trimmed image
    x, y = data.shape

    # set the title
    if title == None:
        raise ValueError('Please specify a title, no title given')
    else:
        title = title

    # set the subtitles
    if subtitle1 == None:
        subtitle1 = 'Observed image'
    else:
        subtitle1 = subtitle1

    if subtitle2 == None:
        subtitle2 = 'Deconvolved image'
    else:
        subtitle2 = subtitle2

    # # # set the min/max values?
    if image_bounds == image_bounds:
        min1 = image_bounds[0]
        max1 = image_bounds[1]
        min2 = image_bounds[2]
        max2 = image_bounds[3]
    else:
        min1 = None
        max1 = None
        min2 = None
        max2 = None

    # Normalize the data for plotting
    if normalize == None:
        raise ValueError('Please specify normalization param')
    elif normalize == 'SimpleNorm':
        normO = simple_norm(data, 'log', clip=True)
        normD = simple_norm(data2, 'log', clip=True)
    elif normalize == 'LogNorm':
        normO = LogNorm(vmin=min1, vmax=max1)
        normD = LogNorm(vmin=min2, vmax=max2)
    elif normalize == 'Normalize':
        normO = Normalize(vmin=min1, vmax=max1)
        normD = Normalize(vmin=min2, vmax=max2)
    elif normalize == 'ImNorm':
        normO = ImageNormalize(stretch=LogStretch(), vmin=min1, vmax=max1)
        normD = ImageNormalize(stretch=LogStretch(), vmin=min2, vmax=max2)
    elif normalize == 'AsinhNorm':
        normO = ImageNormalize(stretch=AsinhStretch(a=0.001), vmin=min1, vmax=max1)
        normD = ImageNormalize(stretch=AsinhStretch(a=0.001), vmin=min2, vmax=max2)

    if labels == None:
        labels = None
    else:
        x_label = labels[0]
        y_label = labels[1]
        c_label = labels[2]

    if cmap == None:
        raise ValueError('Please specify the cmap used. Main options are: 1) viridis, 2) gist-heat, 3) etc.')
    else:
        cmap = cmap

    # set the conversion factor for the image data (if given)
    # toggle whatever units to convert image data too
    if convert == None:
        convert = 1
    else:
        convert = convert

    # 20220830: come back and make the x/y-limit setting more robust
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(15, 6), tight_layout=True)

    # Image 1
    im1 = ax1.imshow(data, origin='lower', cmap=cmap, norm=normO,
                     extent=(0 * convert, x * convert, 0 * convert, y * convert)
                     )
    ax1.set_title(subtitle1, fontsize=15)
    ax1.text(xlim * convert, ylim * convert, str(text1), color='w', fontsize=12)
    #ax1.text((xlim + 2) * convert, (ylim / 4.5) * convert, '1"', color='w', fontsize=15)
    #ax1.hlines(y=(ylim / 5) * convert, xmin=xlim * convert, xmax=(xlim + 9) * convert, linewidth=2, color='w')
    ax1.set_ylabel(y_label)
    ax1.set_xlabel(x_label)
    add_colorbar(im1, label=c_label, format=FormatStrFormatter('%.2e'))

    # Image 2
    im2 = ax2.imshow(data2, origin='lower', cmap=cmap, norm=normD,
                     extent=(0 * convert, x * convert, 0 * convert, y * convert)
                     )
    ax2.set_title(subtitle2, fontsize=15)
    ax2.text(xlim * convert, ylim * convert, str(text2), color='w', fontsize=12)
    ax2.set_ylabel(y_label)
    ax2.set_xlabel(x_label)
    add_colorbar(im2, label=c_label, format=FormatStrFormatter('%.2e'))

    if x_label == 'RA Offset (")':
        # set for large image
        if x == 150:
            positions = [5*convert, 30*convert, 50*convert, 75*convert, 100*convert, 120*convert, 145*convert]
            labels = ['16', '8', '4', '0', '-4', '-8', '-16']
            ax1.xaxis.set_major_locator(ticker.FixedLocator(positions))
            ax1.xaxis.set_major_formatter(ticker.FixedFormatter(labels))
            ax1.yaxis.set_major_locator(ticker.FixedLocator(positions))
            ax1.yaxis.set_major_formatter(ticker.FixedFormatter(labels))
            ax2.xaxis.set_major_locator(ticker.FixedLocator(positions))
            ax2.xaxis.set_major_formatter(ticker.FixedFormatter(labels))
            ax2.yaxis.set_major_locator(ticker.FixedLocator(positions))
            ax2.yaxis.set_major_formatter(ticker.FixedFormatter(labels))
        elif x == 120:
            positions = [10*convert, 16.5*convert, 33*convert, 60*convert, 77*convert, 93.5*convert, 110*convert]
            labels = ['12', '6', '3', '0', '-3', '-6', '-12']
            ax1.xaxis.set_major_locator(ticker.FixedLocator(positions))
            ax1.xaxis.set_major_formatter(ticker.FixedFormatter(labels))
            ax1.yaxis.set_major_locator(ticker.FixedLocator(positions))
            ax1.yaxis.set_major_formatter(ticker.FixedFormatter(labels))
            ax2.xaxis.set_major_locator(ticker.FixedLocator(positions))
            ax2.xaxis.set_major_formatter(ticker.FixedFormatter(labels))
            ax2.yaxis.set_major_locator(ticker.FixedLocator(positions))
            ax2.yaxis.set_major_formatter(ticker.FixedFormatter(labels))
        elif x == 50:
            positions = [5*convert, 15*convert, 25*convert, 35*convert, 45*convert]
            labels = ['5', '2', '0', '-2', '-5']
            ax1.xaxis.set_major_locator(ticker.FixedLocator(positions))
            ax1.xaxis.set_major_formatter(ticker.FixedFormatter(labels))
            ax1.yaxis.set_major_locator(ticker.FixedLocator(positions))
            ax1.yaxis.set_major_formatter(ticker.FixedFormatter(labels))
            ax2.xaxis.set_major_locator(ticker.FixedLocator(positions))
            ax2.xaxis.set_major_formatter(ticker.FixedFormatter(labels))
            ax2.yaxis.set_major_locator(ticker.FixedLocator(positions))
            ax2.yaxis.set_major_formatter(ticker.FixedFormatter(labels))
    else:
        if x == 150:
            ax1.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
            ax1.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
            ax2.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
            ax2.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
        elif x == 120:
            ax1.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
            ax1.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
            ax2.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
            ax2.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
        elif x == 50:
            ax1.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
            ax1.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
            ax2.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
            ax2.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))

    # Set aperture fit params-> image 1
    if aperture1:
        # convert pixel values to arcsecond
        x1 = aperture1[0] * convert
        y1 = aperture1[1] * convert
        rad1 = aperture1[2] * convert

        # Circle((x,y-centroids), aperture radius)
        ax1.add_patch(matplotlib.patches.Circle((x1, y1), rad1, fill=False, color='white'))

    # Set aperture fit params-> image 2
    if aperture2:
        # convert pixel values to arcsecond
        x2 = aperture2[0] * convert
        y2 = aperture2[1] * convert
        rad2 = aperture2[2] * convert

        # Circle((x,y-centroids), aperture radius)
        ax2.add_patch(matplotlib.patches.Circle((x2, y2), rad2, fill=False, color='white'))

    ax1.set_aspect('equal')
    ax2.set_aspect('equal')

    if compass:
        # set the x,y image data
        x1,y1 = data.shape

        # 20220901: hard code params for now -> return and make more robust later
        if labels[0] == 'X (")':
            if x1 == 150:
                # size = [txt1-x, txt1-y, txt2-x, txt2-y, arrow1-x, arrow1-y, arrow1-xlen, arrow1-ylen, arrow2-x, arrow2-y, arrow2-xlen, arrow2-ylen]
                size = [136, 8, 130.5, 26.5, 125, 15, 4, 8, 125, 15, 8, -4]
            if x1 == 120:
                # size = [txt1-x, txt1-y, txt2-x, txt2-y, arrow1-x, arrow1-y, arrow1-xlen, arrow1-ylen, arrow2-x, arrow2-y, arrow2-xlen, arrow2-ylen]
                size = [106, 8, 100.5, 26.5, 95, 15, 4, 8, 95, 15, 8, -4]
            else:
                # size = [txt1-x, txt1-y, txt2-x, txt2-y, arrow1-x, arrow1-y, arrow1-xlen, arrow1-ylen, arrow2-x, arrow2-y, arrow2-xlen, arrow2-ylen]
                size = [45, 1, 43, 10, 40, 5, 2, 3.5, 40, 5, 3.5, -2]
        else:
            # 50 -> x,y*1
            if x1 == 150:
                # size = [txt1-x, txt1-y, txt2-x, txt2-y, arrow1-x, arrow1-y, arrow1-xlen, arrow1-ylen, arrow2-x, arrow2-y, arrow2-xlen, arrow2-ylen]
                size = [138, 25, 123, 8, 140, 10, 0, 10, 140, 10, -10, 0]
            if x1 == 120:
                # size = [txt1-x, txt1-y, txt2-x, txt2-y, arrow1-x, arrow1-y, arrow1-xlen, arrow1-ylen, arrow2-x, arrow2-y, arrow2-xlen, arrow2-ylen]
                size = [108, 25, 93, 8, 110, 10, 0, 10, 110, 10, -10, 0]
            else:
                # size = [txt1-x, txt1-y, txt2-x, txt2-y, arrow1-x, arrow1-y, arrow1-xlen, arrow1-ylen, arrow2-x, arrow2-y, arrow2-xlen, arrow2-ylen]
                size = [44.2, 12, 37.5, 4.5, 45, 5, 0, 4, 45, 5, -4, 0]

        # add text
        plt.text(size[0]*convert, size[1]*convert, 'N', color = 'w')
        plt.text(size[2]*convert, size[3]*convert, 'E', color = 'w')

        # create the height of the arrow
        arrow1 = plt.arrow(size[4]*convert, size[5]*convert, size[6]*convert, size[7]*convert, width = 0.04, color = 'w')

        # create the length of the arrow
        arrow2 = plt.arrow(size[8]*convert, size[9]*convert, size[10]*convert, size[11]*convert, width = 0.04, color = 'w')

        # add arrows
        ax2.add_patch(arrow1)
        ax2.add_patch(arrow2)

    plt.suptitle(title, fontsize = 20)

    if print_params:
        print('Method: show_Image')
        print('Observed image shape: ', data.shape)
        print('User-specified title: ', title)
        print('1st-image caption: ', str(text1))
        print('2nd-image caption: ', str(text2))
        print('User-specified cmap: ', cmap)

    if save:
        plt.savefig(save + '.png')

    if display:
        plt.show()

def show_Grid10Image(data, wavelength=None, nint=None, flux=None, fwhm=None, strehl=None, cmap=None, print_params = False,
                     display=False, save=False):
    '''
    2022 May 24: Method used to plot a 2x5 grid of images.
    Parameters
    __________
    data:           numpy.array
                    2D-numpy array containing the image data
    wavelength:     float
                    wavelength the image was observed in
    nint:           int
                    Iteration number
    flux:           numpy.array
                    Flux measurement of the image
    fwhm:           numpy.array
                    FWHM measurement of the image
    strehl:         numpy.array
                    Strehl ratio measurement of the image
    cmap:           string
                    User-defined colorband for the image
    display:        boolean
                    Toggle whether to display the image or not
    save:           string
                    OPTIONAL: the path to save the output image in
    print_params:   boolean
                    OPTIONAL: print the outputs of each step for trouble-shooting
    Returns
    __________
    image:  numpy.array
            2D-image displaying the flux, FWHM and Strehl measurements
    '''
    # find the (x,y) limit of the trimmed image
    x, y = data[0].shape

    # set the observed wavelength
    if wavelength == None:
        raise ValueError('Please specify the observation wavelength, no wavelength given')
    else:
        if wavelength == 5.6:
            obs_wave = float(wavelength)
        else:
            obs_wave = int(wavelength)

    # set the iteration
    if nint == None:
        raise ValueError('Please specify the observation wavelength, no wavelength given')
    else:
        # pad the nint array by removing the last count
        rem_nint = len(nint)-1
        im_nint = nint[:rem_nint]

    # set the flux
    if flux == None:
        raise ValueError('Please enter the flux measured of this image')
    else:
        im_flux = flux

    # Set the FWHM
    if fwhm == None:
        raise ValueError('Please enter the FWHM (pixels) measured of this image')
    else:
        im_fwhm = fwhm

    # set the Strehl ratio
    if strehl == None:
        raise ValueError('Please enter the Strehl ratio measured of this image')
    else:
        im_strehl = strehl

    # set the cmap
    if cmap == None:
        raise ValueError('Please specify the cmap used. Main options are: 1) viridis, 2) gist-heat, 3) etc.')
    else:
        cmap = cmap

    # format outputs for plotting
    im_flux2 = []
    im_fwhm2 = []
    im_strehl2 = []
    for i in range(0, len(im_flux)):
        # flux
        im_flux1 = im_flux[i]
        im_flux1 = '{:.7}'.format(im_flux1)
        im_flux1 = str(im_flux1)
        im_flux2.append(im_flux1)

        # fwhm
        im_fwhm1 = im_fwhm[i]
        im_fwhm1 = '{:.5}'.format(im_fwhm1)
        im_fwhm1 = str(im_fwhm1)
        im_fwhm2.append(im_fwhm1)

        # strehl
        im_strehl1 = im_strehl[i]
        im_strehl1 = '{:.5}'.format(im_strehl1)
        im_strehl1 = str(im_strehl1)
        im_strehl2.append(im_strehl1)

    # normalize/log-scale images
    # the observation
    data2 = data[0]
    normD1 = simple_norm(data2, 'log')

    # deconvolved images
    for i in data:
        normD = simple_norm(i, 'log')

    # plot the observed image and the 1st 9 iterations
    fig = plt.figure(figsize=(12, 6))

    grid = ImageGrid(fig, 111,  # similar to subplot(111)
                     nrows_ncols=(2, 5),  # creates 9x9 grid of axes
                     axes_pad=0.5,  # pad between axes in inch.
                     share_all=True
                     )

    for ax, im, n_int, flux, fwhm, strehl in zip(grid, data, im_nint, im_flux2, im_fwhm2, im_strehl2):
        # append the observation first
        if n_int == 0:
            ax.imshow(im, origin='lower', norm=normD1, cmap=cmap)
            ax.set_title('NGC 5728' + str(obs_wave) + ' μm')
            plt.tight_layout(pad=1.5)
            ax.text(10, y - 20, 'Flux: ' + flux + ' (mJy)', color='w', fontsize=8)
            ax.text(10, y - 35, 'FWHM: ' + fwhm + '"', color='w', fontsize=8)
            ax.text(10, y - 50, 'Strehl: ' + strehl, color='w', fontsize=8)

        # now plot all the deconvolution outputs
        else:
            ax.imshow(im, origin='lower', norm=normD, cmap=cmap)
            ax.set_title('n-iteration: ' + str(n_int))
            plt.tight_layout(pad=1.5)
            ax.text(10, y - 20, 'Flux: ' + flux + ' (mJy)', color='w', fontsize=8)
            ax.text(10, y - 35, 'FWHM: ' + fwhm + '"', color='w', fontsize=8)
            ax.text(10, y - 50, 'Strehl: ' + strehl, color='w', fontsize=8)

        fig.suptitle('Richardson-Lucy Deconvolution', fontsize=20)
        fig.supxlabel('RA (")')
        fig.supylabel('DEC (")')

    if print_params:
        print('Method: show_Grid10Image')
        print('Observed image shape: ', data[0].shape)
        print('Observed image wavelength: ', obs_wave)
        print('Iteration number: ', nint)
        print('Observed image flux: ', im_flux)
        print('Observed image FWHM: ',im_fwhm)
        print('Observed image Strehl ratio: ',im_strehl)
        print('User-specified cmap: ', cmap)

    if save:
        plt.savefig(save + '.png')

    if display:
        plt.show()


def show_Grid25Image(data, wavelength=None, nint=None, flux=None, fwhm=None, strehl=None, cmap=None, print_params=False,
                     display=False, save=False):
    '''
    2022 May 24: Method used to plot a 5x5 grid of images.
    Parameters
    __________
    data:           numpy.array
                    2D-numpy array containing the image data
    wavelength:     float
                    wavelength the image was observed in
    nint:           int
                    Iteration number
    flux:           numpy.array
                    Flux measurement of the image
    fwhm:           numpy.array
                    FWHM measurement of the image
    strehl:         numpy.array
                    Strehl ratio measurement of the image
    cmap:           string
                    User-defined colorband for the image
    display:        boolean
                    OPTIONAL: Toggle whether to display the image or not
    save:           string
                    OPTIONAL: the path to save the output image in
    print_params:   boolean
                    OPTIONAL: print the outputs of each step for trouble-shooting
    Returns
    __________
    image:  numpy.array
            2D-image displaying the flux, FWHM and Strehl measurements
    '''
    # find the (x,y) limit of the image
    x, y = data[0].shape

    # set the observed wavelength
    if wavelength == None:
        raise ValueError('Please specify the observation wavelength, no wavelength given')
    else:
        if wavelength == 5.6:
            obs_wave = float(wavelength)
        else:
            obs_wave = int(wavelength)

    # set the iteration
    if nint == None:
        raise ValueError('No iteration measurements given, please specify the iteration numbers.')
    else:
        # pad the nint array by removing the last count
        rem_nint = len(nint)-1
        im_nint = nint[:rem_nint]

    # set the flux
    if flux == None:
        raise ValueError('Please enter the flux array measured of this image')
    else:
        im_flux = flux

    # Set the FWHM
    if fwhm == None:
        raise ValueError('Please enter the FWHM (pixels) array measured of this image')
    else:
        im_fwhm = fwhm

    # set the Strehl ratio
    if strehl == None:
        raise ValueError('Please enter the Strehl ratio array measured of this image')
    else:
        im_strehl = strehl

    # set the cmap
    if cmap == None:
        raise ValueError('Please specify the cmap used. Main options are: 1) viridis, 2) gist-heat, 3) etc.')
    else:
        cmap = cmap

    # format outputs for plotting
    im_flux2 = []
    im_fwhm2 = []
    im_strehl2 = []
    for i in range(0, len(im_flux)):
        # flux
        im_flux1 = im_flux[i]
        im_flux1 = '{:.7}'.format(im_flux1)
        im_flux1 = str(im_flux1)
        im_flux2.append(im_flux1)

        # fwhm
        im_fwhm1 = im_fwhm[i]
        im_fwhm1 = '{:.5}'.format(im_fwhm1)
        im_fwhm1 = str(im_fwhm1)
        im_fwhm2.append(im_fwhm1)

        # strehl
        im_strehl1 = im_strehl[i]
        im_strehl1 = '{:.5}'.format(im_strehl1)
        im_strehl1 = str(im_strehl1)
        im_strehl2.append(im_strehl1)

    # Split the data sets into 5x5 chunks
    image_grid = list()
    nint_grid = list()
    fwhm_grid = list()
    flux_grid = list()
    strehl_grid = list()

    chunk_size = 25
    for i in range(0, len(data), chunk_size):
        # append the image data in groups of 25
        image_grid.append(data[i:i + chunk_size])

        # append the iteration number in groups of 25
        nint_grid.append(im_nint[i:i + chunk_size])

        # append the flux measures in groups of 25
        flux_grid.append(im_flux2[i:i + chunk_size])

        # append the fwhm measures in groups of 25
        fwhm_grid.append(im_fwhm2[i:i + chunk_size])

        # append the strehl measures in groups of 25
        strehl_grid.append(im_strehl2[i:i + chunk_size])

    # normalize/log-scale images
    # the observation
    data2 = data[0]
    normD1 = simple_norm(data2, 'log')

    # the deconvolved images
    for i in data:
        normD = simple_norm(i, 'log')

    chunks = []
    for i in range(len(image_grid)):
        print('chunk #: ', i)
        fig = plt.figure(figsize=(14,14), tight_layout=True)
        grid = ImageGrid(fig, 111,  # similar to subplot(111)
                         nrows_ncols=(5, 5),  # creates 5x5 grid of axes
                         axes_pad=0.5,  # pad between axes in inch.
                         share_all=True
                         )

        for ax, im, n_int, flux, fwhm, strehl in zip(grid, image_grid[i], nint_grid[i], flux_grid[i], fwhm_grid[i], strehl_grid[i]):
            # plot the observed image and the 1st 25 iterations plus any subsequent grids
            # append the observation first
            if n_int == 0:
                ax.imshow(im, origin='lower', norm=normD1, cmap=cmap)
                ax.set_title('NGC 5728' + str(obs_wave) + ' μm')
                plt.tight_layout(pad=1.5)
                ax.text(10, y - 20, 'Counts: ' + flux + ' (mJy)', color='w', fontsize=8)
                ax.text(10, y - 35, 'FWHM: ' + fwhm + '"', color='w', fontsize=8)
                ax.text(10, y - 50, 'Strehl: ' + strehl, color='w', fontsize=8)

            # now plot all the deconvolution outputs
            else:
                ax.imshow(im, origin='lower', norm=normD, cmap=cmap)
                ax.set_title('n-iteration: ' + str(n_int))
                plt.tight_layout(pad=1.5)
                ax.text(10, y - 20, 'Counts: ' + flux + ' (mJy)', color='w', fontsize=8)
                ax.text(10, y - 35, 'FWHM: ' + fwhm + '"', color='w', fontsize=8)
                ax.text(10, y - 50, 'Strehl: ' + strehl, color='w', fontsize=8)

            fig.suptitle('Richardson-Lucy Deconvolution', fontsize=20)
            fig.supxlabel('RA (")')
            fig.supylabel('DEC (")')
        chunks.append(i)
        if save:
            plt.savefig(save + '_chunk_' + str(i) + '.png')

    if print_params:
        print('Method: show_Grid10Image')
        print('Observed image shape: ', data[0].shape)
        print('Observed image wavelength: ', obs_wave)
        print('Iteration number: ', nint)
        print('Observed image flux: ', im_flux)
        print('Observed image FWHM: ', im_fwhm)
        print('Observed image Strehl ratio: ', im_strehl)
        print('User-specified cmap: ', cmap)
        print('Number of chunks:', chunks)

    if display:
        plt.show()


def plot_dif_MeritFunctions(nint=None, flux=None, fwhm=None, print_params=False, display=False, plot=None,
                        save=False, optimum = None, title = None):
    ''' Plot the results of all three merit functions as a function od deconvolution iteration in a 3x1 subplot or
    as individual plots
    Parameters
    __________
    nint:            numpy.array
                     Array containing the number of iterations
    flux:            numpy.array
                     Array containing the flux measurements
    fwhm:            numpy.array
                     Array containing the fwhm measurements
    print_params:    boolean
                     OPTIONAL: print the outputs of each step for trouble-shooting
    display:         boolean
                     OPTIONAL: Toggle whether to display the plot or not
    plot:            string
                     OPTIONAL: choose to plot individual merit functions only. Options are:
                     1) None- Default, plots all the outputs
                     2) 'flux'- plot the results of the flux measurements only
                     3) 'FWHM'- plot the results of the FWHM measurements only
                     4) 'strehl'- plot the results of the Strehl measurements only
    save:            boolean
                     OPTIONAL: the path to save the output image in
    Optimum:         string
                     Set the optimum deconvolution number
    Returns
    __________
    plot:            numpy.array
                     2D-image displaying the flux, FWHM and Strehl measurements plot
    '''

    # set the iteration
    if nint == None:
        raise ValueError('No iteration measurements given, please specify the iteration numbers.')
    else:
        # pad the nint array by removing the last count
        im_nint = nint

    # set the flux
    if flux == None:
        raise ValueError('Please enter the flux array measured of this image')
    else:
        im_flux = flux

    # Set the FWHM
    if fwhm == None:
        raise ValueError('Please enter the FWHM (pixels) array measured of this image')
    else:
        im_fwhm = fwhm

    # Set the (x,y) data points
    x = im_nint
    left_y = im_flux
    cen_y = im_fwhm
    right_y = im_fwhm

    if plot == None:
        # plot the data
        fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(12, 6), tight_layout=True)

        # add trendline?
        z_fwhm = np.polyfit(x, cen_y, 2)
        p_fwhm = np.poly1d(z_fwhm)

        z_flux = np.polyfit(x, left_y, 2)
        p_flux = np.poly1d(z_flux)

        # Create the 1x3 grid
        # Plot the FWHM values
        ax[0].scatter(x, cen_y, color='b', marker = 'o')
        ax[0].plot(x, cen_y, color='r', ls = '--', label='FWHM')
        ax[0].set_title('FWHM vs Iteration', fontsize=10)
        ax[0].set_xlim(0, int(len(x)))
        ax[0].set_ylabel('FWHM (")')
        ax[0].axvline(x=optimum, ymin=0, ymax=len(left_y), color='g', ls='--', label='iteration ' + str(optimum))
        ax[0].yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
        ax[0].legend(loc='center right')

        # Plot the flux values
        ax[1].scatter(x, left_y, color='b', marker = 'o')
        ax[1].plot(x, left_y, color='r', ls='--', label = 'Flux')
        ax[1].set_title('Flux vs Iteration', fontsize=10)
        ax[1].set_xlim(0, int(len(x)))
        ax[1].set_ylabel('Flux Density (mJy)')
        ax[1].axvline(x=optimum, ymin=0, ymax=len(left_y), color='g', ls='--', label='iteration ' + str(optimum))
        ax[1].yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
        ax[1].legend(loc='center right')

        fig.supxlabel('Iterations')
        plt.suptitle(title, fontsize=15)

    elif plot == 'flux':
        # plot the data
        fig = plt.figure(figsize=(8, 6), tight_layout=True)

        plt.plot(x, left_y, label='Flux', color='b')
        plt.xlim(0, int(len(x)))
        plt.title('Flux vs Iteration', fontsize=12)
        plt.ylabel('Flux Density (Jy')
        plt.xlabel('Iterations')
        plt.axvline(x=optimum, ymin=0, ymax=len(left_y), color='r', ls='--', label='Flux limit\n iteration 43')
        plt.legend(loc='center right')

    elif plot == 'FWHM':
        # plot the data
        fig = plt.figure(figsize=(8, 6), tight_layout=True)

        plt.plot(x, cen_y, label='FWHM', color='b')
        plt.xlim(0, int(len(x)))
        plt.ylabel('FWHM (")')
        plt.xlabel('Iterations')
        plt.title('FWHM vs Iteration')
        plt.axvline(x=optimum, ymin=0, ymax=len(left_y), color='r', ls='--', label='FWHM limit\n iteration 30')
        plt.legend(loc='center right')

    elif plot == 'strehl':
        # plot the data
        fig = plt.figure(figsize=(8, 6), tight_layout=True)

        plt.plot(x, right_y, label='Strehl ratio', color='b')
        plt.title('Strehl ratio vs Iteration', fontsize=12)
        plt.xlim(0, int(len(x)))
        plt.ylabel('Strehl ratio')
        plt.legend(loc='upper right')

    if print_params:
        print('Method: plot_MeritFunctions')
        print('Iteration number: ', im_nint)
        print('Observed image flux: ', im_flux)
        print('Observed image FWHM: ', im_fwhm)

    if save:
        plt.savefig(save + '.png')

    if display:
        plt.show()

def plot_same_MeritFunctions(nint=None, flux=None, fwhm=None, strehl=None, print_params=False, display=False,
                        save=False,):
    ''' Plot the results of all three merit functions as a function of deconvolution iteration on the same plot
    Parameters
    __________
    nint:            numpy.array
                     Array containing the number of iterations
    flux:            numpy.array
                     Array containing the flux measurements
    fwhm:            numpy.array
                     Array containing the fwhm measurements
    strehl:          numpy.array
                     Array containing the Strehl ratio measurements
    print_params:    boolean
                     OPTIONAL: print the outputs of each step for trouble-shooting
    display:         boolean
                     OPTIONAL: Toggle whether to display the plot or not
    save:            boolean
                     OPTIONAL: the path to save the output image in
    Returns
    __________
    plot:            numpy.array
                     2D-image displaying the flux, FWHM and Strehl measurements plot
    '''

    # set the iteration
    if nint == None:
        raise ValueError('No iteration measurements given, please specify the iteration numbers.')
    else:
        # pad the nint array by removing the last count
        im_nint = nint

        # set the optimum deconvolution limit
        optimum = im_nint.pop(0)

    # set the flux
    if flux == None:
        raise ValueError('Please enter the flux array measured of this image')
    else:
        im_flux = flux
        y_lim = np.max(im_flux)

    # Set the FWHM
    if fwhm == None:
        raise ValueError('Please enter the FWHM (pixels) array measured of this image')
    else:
        im_fwhm = fwhm

    # set the Strehl ratio
    if strehl == None:
        raise ValueError('Please enter the Strehl ratio array measured of this image')
    else:
        im_strehl = strehl

    # plot the data
    fig = plt.figure(figsize=(8, 6), tight_layout=True)

    # Set the (x,y) data points
    x = im_nint
    left_y = im_flux
    cen_y = im_fwhm
    right_y = im_strehl

    # format outputs for plotting
    im_flux2 = []
    im_fwhm2 = []
    im_strehl2 = []
    for i in range(0, len(im_flux)):
        # flux
        im_flux1 = im_flux[i]
        im_flux1 = '{:.5}'.format(im_flux1)
        im_flux1 = str(im_flux1)
        im_flux2.append(im_flux1)

        # fwhm
        im_fwhm1 = im_fwhm[i]
        im_fwhm1 = '{:.4}'.format(im_fwhm1)
        im_fwhm1 = str(im_fwhm1)
        im_fwhm2.append(im_fwhm1)

        # strehl
        im_strehl1 = im_strehl[i]
        im_strehl1 = '{:.4}'.format(im_strehl1)
        im_strehl1 = str(im_strehl1)
        im_strehl2.append(im_strehl1)

    # Create the plot
    plt.plot(x, left_y, label='Counts', color='b')
    plt.plot(x, cen_y, label='FWHM', color='r')
    plt.plot(x, right_y, label='Strehl ratio', color='g')
    plt.text(optimum+5, y_lim-10, '{:s}'.format('\u0332'.join('Optimum Iteration: ' + str(optimum))) +
             '\nFlux: ' + str(im_flux2[optimum]) + ' (Jy)'
              + '\nFWHM: ' + str(im_fwhm2[optimum]) + '"'
             + '\nStrehl: ' + str(im_strehl2[optimum]), color = 'k', fontsize=10)
    plt.axvline(x=optimum, ymin=0, ymax = len(left_y), color='k', ls = '--')
    plt.title('Merit functions vs Iteration', fontsize=12)
    plt.xlim(0, int(len(x)))
    plt.xlabel('Iterations')
    plt.ylabel('Merit function')
    plt.legend(loc='center right')

    if print_params:
        print('Method: plot_MeritFunctions')
        print('Iteration number: ', im_nint)
        print('Observed image flux: ', im_flux)
        print('Observed image FWHM: ', im_fwhm)
        print('Observed image Strehl ratio: ', im_strehl)

    if save:
        plt.savefig(save + '.png')

    if display:
        plt.show()

#******************************************** Additional Tools *********************************************************#

def NormalizeData(data):
    '''Funcion to normalize your data between a range of 0 - 1'''
    return (data - np.min(data)) / (np.max(data) - np.min(data))

def bin_ndarray(ndarray, new_shape, operation='sum'):
    """
    2022 May 23: Code from https://gist.github.com/derricw/95eab740e1b08b78c03f
    Bins an ndarray in all axes based on the target shape, by summing or
        averaging.
    Number of output dimensions must match number of input dimensions.
    Example
    -------
    >>> m = np.arange(0,100,1).reshape((10,10))
    >>> n = bin_ndarray(m, new_shape=(5,5), operation='sum')
    >>> print(n)
    [[ 22  30  38  46  54]
     [102 110 118 126 134]
     [182 190 198 206 214]
     [262 270 278 286 294]
     [342 350 358 366 374]]
    """
    if not operation.lower() in ['sum', 'mean', 'average', 'avg']:
        raise ValueError("Operation {} not supported.".format(operation))
    if ndarray.ndim != len(new_shape):
        raise ValueError("Shape mismatch: {} -> {}".format(ndarray.shape,
                                                           new_shape))
    compression_pairs = [(d, c//d) for d, c in zip(new_shape,
                                                   ndarray.shape)]
    flattened = [l for p in compression_pairs for l in p]
    ndarray = ndarray.reshape(flattened)
    for i in range(len(new_shape)):
        if operation.lower() == "sum":
            ndarray = ndarray.sum(-1*(i+1))
        elif operation.lower() in ["mean", "average", "avg"]:
            ndarray = ndarray.mean(-1*(i+1))
    return ndarray

def add_colorbar(im, aspect=20, pad_fraction=0, **kwargs):
    """Add a vertical color bar to an image plot. Function
     to fit colorbar to image size"""
    divider = axes_grid1.make_axes_locatable(im.axes)
    width = axes_grid1.axes_size.AxesY(im.axes, aspect=1./aspect)
    pad = axes_grid1.axes_size.Fraction(pad_fraction, width)
    current_ax = plt.gca()
    cax = divider.append_axes("right", size=width, pad=pad)
    plt.sca(current_ax)
    return im.axes.figure.colorbar(im, cax=cax, **kwargs)

def calculate_difference(a, b):
    '''2022 Oct 27: Simple code to calculate percent difference'''

    diff = (np.abs(b-a)/a) * 100
    diff = float(diff)

    # scale the outputs to 5-decimal places and convert back to float
    diff = '%.8f' % diff
    diff = np.float64(diff)

    return diff

def MAE_SNR_calculator(data, extract = None, box_size = None, centroid = None, aper_rad = None):
    '''Code derived shamelessly from Amilcar Torres'Quijano on 22 July 2023. Measure the SNR in 4 background regions of
    the image, then compare the combined noise (median/mean?) to the source.
    https://github.com/Locuan/canaricam-pipelines/blob/main/snr_calculator.py
        Parameters
    __________
    data:            np.array
                     Imaging data
    extract:         int
                     define box size for data extraction
    box_size:        int
                     define box size for the background measurements
    centroid:        list
                     Initial x,y centroid guess
    aper_rad:        float
                     aperture radius size to measure photometry
    Returns
    __________
    snr:             int
                     signal-to-noise ratio of the data'''

    # set the data
    data = data

    # define box size for data extraction
    if extract == None:
        print('extract param missing')
    else:
        extract = extract

    # define box size for the background measurements
    if box_size == None:
        print('box_size param missing')
    else:
        box_size = box_size

    # set the centroid
    if centroid == None:
        print('centroid input missing')
    else:
        centroid = centroid
        position = [centroid[0], centroid[1]]
        box = 15
        center_quad = centroid_quadratic(data,
                                         xpeak=position[0],
                                         ypeak=position[1],
                                         fit_boxsize=box)
        center = (center_quad[0], center_quad[1])

    # set the aperture radius
    if aper_rad == None:
        print('aper_rad input missing')
    else:
        aper_rad = aper_rad

    # define the four box origins
    x1, y1 = 0, 0
    x2, y2 = extract * 2, 0
    x3, y3 = 0, extract * 2
    x4, y4 = extract * 2, extract * 2

    # define the regions to measure stats from
    box1 = data[y1:y1 + box_size, x1:x1 + box_size]
    box2 = data[y2:y2 + box_size, x2 - box_size: x2]
    box3 = data[y3 - box_size:y3, x3:x3 + box_size]
    box4 = data[y4 - box_size:y4, x4 - box_size:x4]

    # measure stats for each box
    box1_mean, box1_median, box1_std = round(np.mean(box1)), round(np.median(box1)), round(np.std(box1))
    box2_mean, box2_median, box2_std = round(np.mean(box2)), round(np.median(box2)), round(np.std(box2))
    box3_mean, box3_median, box3_std = round(np.mean(box3)), round(np.median(box3)), round(np.std(box3))
    box4_mean, box4_median, box4_std = round(np.mean(box4)), round(np.median(box4)), round(np.std(box4))
    # add all boxes to an array
    full_box = np.array([])
    # box1
    for j in range(len(box1)):
        for i in range(len(box1)):
            full_box = np.append(full_box, box1[j][i])
    # box2
    for j in range(len(box2)):
        for i in range(len(box2)):
            full_box = np.append(full_box, box2[j][i])
    # box3
    for j in range(len(box3)):
        for i in range(len(box3)):
            full_box = np.append(full_box, box3[j][i])
    # box4
    for j in range(len(box4)):
        for i in range(len(box4)):
            full_box = np.append(full_box, box4[j][i])

    # measure photometry of image data source
    # Create the Aperture and Annulus for the Standard Star
    aperture = photutils_CircularAperture(center, r=aper_rad)
    # Define the values for the Aperture Mask/Data
    # #The data portion being the array that corresponds to aperture pixels
    aperture_masks = aperture.to_mask(method='center')
    aperture_data = aperture_masks.multiply(data)
    # Store the values for the counts of the aperture and annulus
    aperture_counts = 0
    # Store the number of pixels for the aperture and annulus
    aperture_pixels = 0
    # Store the values within the aperture and annulus
    aperture_values = np.array([])
    # Statistics for the full box
    full_box_mean, full_box_median, full_box_noise = round(np.mean(full_box)), round(np.median(full_box)), np.std(full_box)

    # Add the median to the whole image
    for j in range(len(data)):
        for i in range(len(data[j])):
            # data[j][i] = data[j][i] - annulus_mean
            data[j][i] = data[j][i] - full_box_median
        # Aperture of whole image with mean
        aperture_data = aperture_masks.multiply(data)

    # Get the aperture counts
    for j in range(len(aperture_data)):
        for i in range(len(aperture_data)):
            if aperture_data[j][i] != 0:
                aperture_pixels += 1
                aperture_counts += aperture_data[j][i]
                aperture_values = np.append(aperture_values, aperture_data[j][i])

    #aperture_counts = round(aperture_counts)
    # Statistics for the Aperture
    #aperture_mean = round(np.mean(aperture_values))
    #aperture_median = round(np.median(aperture_values))
    #aperture_std = round(np.std(aperture_values))
    #signal_to_noise = round(aperture_counts/full_box_noise)
    signal_to_noise = aperture_counts / (full_box_noise * math.sqrt(aperture_pixels))

    # DELETE ME LATER-> convert counts -> mJy/sr?
    #cf = 0.42818
    #aperture_counts = aperture_counts*cf
    #full_box_noise = full_box_noise * cf
    #print()
    #print('apeture counts = ', aperture_counts)
    #print(aperture_pixels)
    #print('aperture noise = ', full_box_noise)
    #print('SNR = ', aperture_counts / (full_box_noise * math.sqrt(aperture_pixels)))
    #print()


    return signal_to_noise

def convert_units(input_file, output_file):
    '''19 Jan 2024: simple code toread in a .txt file, convert units from angstrom to um and e/photon to Jy'''
    with open(input_file, 'r') as f:
        lines = f.readlines()

    converted_data = []
    for line in lines:
        data = line.strip().split()
        if len(data) == 2:
            wavelength_angstrom = float(data[0])
            flux_e_photons = float(data[1])

            # Convert units
            wavelength_um = wavelength_angstrom / 10000  # Conversion from angstrom to micrometers
            flux_jansky = flux_e_photons * 1.98644586e-9  # Conversion from e/photons to jansky

            converted_data.append((wavelength_um, flux_jansky))

    with open(output_file, 'w') as f:
        for data_point in converted_data:
            f.write(f"{data_point[0]} {data_point[1]}\n")

# define CAS measureables
def calculate_cas_parameters(image):
    x,y = image.shape
    # Concentration (C)
    total_intensity = np.sum(image)
    sorted_intensity = np.sort(image.flatten())[::-1]
    concentration = sorted_intensity[:int(0.2 * len(sorted_intensity))].sum() / total_intensity

    # Asymmetry (A)
    mirrored_image = np.fliplr(image)
    asymmetry = np.abs(image - mirrored_image).sum() / total_intensity

    # Smoothness (S)
    smoothness = ssim(image, mirrored_image, data_range=x)

    return concentration, asymmetry, smoothness

def aperture_mean_nomask(ap, data, **kwargs):
    """
    NOTE: all code below referenced from the statmorph source code and modified for my testing purposes.
    If you want the code that works you should visit: https://github.com/vrodgom/statmorph/tree/master

    Calculate the mean flux of an image for a given photutils
    aperture object. Note that we do not use ``_aperture_area``
    here. Instead, we divide by the full area of the
    aperture, regardless of masked and out-of-range pixels.
    This avoids problems when the aperture is larger than the
    region of interest.
    """
    return ap.do_photometry(data, **kwargs)[0][0] / ap.area

def petro_resize(data, petro_extent_cas):
    ''''''
    # compute the image centroid position from a quadratic fit
    # guess the centroid position based on the image peak
    im_y, im_x = np.unravel_index(np.argmax(data, axis=None), data.shape)
    # determine the centroid from a quadratic fit using the image peak as the initial guess
    xycent = centroid_quadratic(data, xpeak=im_x, ypeak=im_y, fit_boxsize=15)

    # compute the Petrosian radius of the image, focused on the centroid position
    rad = rpetro_circ_generic(data, xycent)

    # set the bounding limits from the image based on the image centroid
    # Note: array slicing only deals in whole integers
    ceny, cenx = int(round(xycent[0])), int(round(xycent[1]))
    limit = int(round(petro_extent_cas*rad))

    # resize the image array to 1.5x the Petrosian radius
    # Note: petro_extent_cas is a global value = 1.5
    data = data[ceny-limit:ceny+limit, cenx-limit:cenx+limit]

    return data

def petrosian_function_circ(data, r, center, annulus_width, eta):
    """
    Helper function to calculate the circular Petrosian radius.

    For a given radius ``r``, return the ratio of the mean flux
    over a circular annulus divided by the mean flux within the
    circle, minus "eta" (eq. 4 from Lotz et al. 2004). The root
    of this function is the Petrosian radius.
    """

    # set the inner/outer annulus width and measure the flux within the annulus
    r_in = r - 0.5 * annulus_width
    r_out = r + 0.5 * annulus_width
    print('inner/outer radius: ', r_in, r_out)
    print('center: ', center)
    print('radius: ', r)
    circ_annulus = photutils.aperture.CircularAnnulus(center, r_in, r_out)

    # measure aperture flux
    circ_aperture = photutils.aperture.CircularAperture(center, r)

    # Force mean fluxes to be positive:
    circ_annulus_mean_flux = np.abs(aperture_mean_nomask(circ_annulus, data, method='exact'))
    circ_aperture_mean_flux = np.abs(aperture_mean_nomask(circ_aperture, data, method='exact'))

    print('annulus mean flux | aperture mean flux: ', circ_annulus_mean_flux, circ_aperture_mean_flux)

    if circ_aperture_mean_flux == 0:
        warnings.warn('[rpetro_circ] Mean flux is zero.', AstropyUserWarning)
        # If flux within annulus is also zero (e.g. beyond the image
        # boundaries), return zero. Otherwise return 1.0:
        ratio = float(circ_annulus_mean_flux != 0)
    else:
        ratio = circ_annulus_mean_flux / circ_aperture_mean_flux

    # return the ratio
    # NOTE: eta is a global parameter = 0.2
    # 1 March 2024 NOTE: return value as a whole interger because scipy.optimize.brentq does not work with float values?
    return int(round(np.abs(ratio-eta)))

def rpetro_circ_generic(data, center, annulus_width):
    """
    Compute the Petrosian radius for concentric circular apertures.

    Notes
    -----
    The so-called "curve of growth" is not always monotonic,
    e.g., when there is a bright, unlabeled and unmasked
    secondary source in the image, so we cannot just apply a
    root-finding algorithm over the full interval.
    Instead, we proceed in two stages: first we do a coarse,
    brute-force search for an appropriate interval (that
    contains a root), and then we apply the root-finder.

    """

    # Find appropriate range for root finder
    npoints = 100
    r_inner = annulus_width

    # compute the outer radius from the diagonal distance of the box cutoff to be measured
    nx, ny = data.shape
    #r_outer = np.sqrt(nx**2 + ny**2)
    r_outer = nx-1

    # compute the radius to measure the Petrosian radius
    assert r_inner < r_outer
    r_min, r_max = None, None
    for r in np.linspace(r_inner, r_outer, npoints):
        # 1 Mar 2024: play around with data or crop_data -> not sure which is correct here
        curval = petrosian_function_circ(data, r, center)
        if curval == 0:
            warnings.warn('[rpetro_circ] Found rpetro by chance?',AstropyUserWarning)
            return r
        elif curval > 0:  # we have not reached rpetro yet
            #print(curval, r)
            r_min = r
        else:  # we are beyond rpetro
            if r_min is None:
                warnings.warn('rpetro_circ < annulus_width! ' +'Setting rpetro_circ = annulus_width.', AstropyUserWarning)
                return r_inner
            r_max = r
            break

    assert r_min is not None
    if r_max is None:
        warnings.warn('[rpetro_circ] rpetro too large.', AstropyUserWarning)
        r_max = r_outer

    #print(r_min, r_max)
    # compute the petrosian radius
    rpetro_circ = opt.brentq(petrosian_function_circ(data, r, center), r_min, r_max, args=(center,), xtol=1e-6)

    return rpetro_circ
