# 2022 Oct 26: Base code to set the merit function parameters of each model/observation as a function of simulated wavelength

from Convenience_Functions import *

# toggle between: (1) centroid testing, (2) FWHM testing, (3) flux testing
centroid = False
fwhm = False

# toggle plots
plot = True
flux_decon = plot
flux_PS = plot
flux_bicone = plot
flux_galaxy = plot
aper21 = plot

# pad the image (by pixel)?
pad = 1

# set universal params (OPTIONAL)
fluxDEC = 17.8
fluxPS = 9.5
fluxMAJOR = 17
fluxMINOR = 15

# Set the filter
filter = 'F2100W'

# set the model type
model_type = 'Observation_model'

# set the model title
model = 'Model_AGN_complicated'

# set the import path for the observed image
im_path = 'Images/JWST/2_MIRISim_outputs/Deconvolution_testing/Residual_outputs/FITS/' +str(model_type) + '/' + str(model)

# set the output path for the deconvolved image FITS file
out_path = 'Images/JWST/4_Deconvolution/1_Richardson_Lucy/Measurements/' +str(model_type) + '/' + str(model)

# read in the image data
image = get_pkg_data_filename(im_path + '/' + str(filter) + '_' + str(model) + '.fits')

# read in image data for deconvolution and saving
image_list = fits.open(image)
image_data = image_list[1].data


# test centroiding
if centroid:
    box = 5.0
    while True:
        if box <= 70.1:
            # determine the image centroid
            guess = (130,132)

            # find the centroid
            center_quad = centroid_quadratic(image_data, xpeak=guess[0], ypeak=guess[1], fit_boxsize=box)
            center = (center_quad[0], center_quad[1])
            print(box, center)
            box += 0.1
            box = '{:.3}'.format(box)
            box = float(box)
        else:
            break

# set the optimum centroid params for the flux/fwhm test
guess = (130,132)

# box size (pixels)
box = 21.5
# find the centroid
center_quad = centroid_quadratic(image_data, xpeak=guess[0], ypeak=guess[1], fit_boxsize=box)
center = (center_quad[0], center_quad[1])
#print(box, center)

if fwhm:
    iter = 0.01
    while True:
        if iter < 3.5:
            # measure the FWHM and vary threshold param of Gaussian fit
            # measure FWHM
            # 2) Fitting a 1D_Gaussian method
            fwhm = MAE_measure_FWHM_gaussian(HDUlist=image,
                                             ext=1,
                                             threshold=iter,
                                             plot=False,
                                             print_params=False,
                                             centroid=center
                                             )

            print(iter, fwhm[1])
            iter += 0.01
            iter = '{:.3}'.format(iter)
            iter = float(iter)
        else:
            break

# measure the flux of the first airy ring -> deconvolution testing
if flux_decon:
    count = fluxDEC
    while True:
        if count <= fluxDEC:
            # measure photometry of the point source -> display aperture fit
            # measure the flux
            flux = MAE_Aperture_Photometry(HDUlist=image,
                                           ext=1,
                                           radius=count,
                                           normalize=False,
                                           centroid=center,
                                           plot=True,
                                           print_params=False,
                                           save=out_path + '/' + str(filter) + '/aper_' + str(count),
                                           display=True,
                                           title= str(filter) + ' ' + str(model),
                                           trim=25,
                                           pad = pad
                                           )

            print(count, flux[0])

            count += 0.1
            count = '{:.3}'.format(count)
            count = float(count)
        else:
            break

# measure the flux of the point source -> flux change testing
if flux_PS:
    count = fluxPS
    while True:
        if count <= fluxPS:
            # measure photometry of the point source -> display aperture fit
            # measure the flux
            flux = MAE_Aperture_Photometry(HDUlist=image,
                                           ext=1,
                                           radius=count,
                                           normalize=False,
                                           centroid=center,
                                           plot=True,
                                           print_params=False,
                                           save=out_path + '/' + str(filter) + '/Photometry_Tests/Flux_Conservation/PS_aper_' + str(count),
                                           display=True,
                                           title= str(filter) + ' ' + str(model),
                                           trim=25,
                                           pad = pad
                                           )

            print(count, flux[0])

            count += 0.1
            count = '{:.3}'.format(count)
            count = float(count)
        else:
            break

# measure the flux of the bicone -> flux change testing
if flux_bicone:
    # semi-major = width
    semi_major = fluxMAJOR
    # semi-minor = height
    semi_minor = fluxMINOR
    theta = 55.0
    while True:
        if theta <= 55.0:
            # measure photometry of the bicone -> display aperture fit
            # measure the flux
            flux = MAE_Aperture_Photometry(HDUlist=image,
                                           ext=1,
                                           radius=1,
                                           normalize=False,
                                           centroid=center,
                                           plot=True,
                                           print_params=False,
                                           save=out_path + '/' + str(filter) + '/Photometry_Tests/Flux_Conservation/bicone_aper_' + str(semi_major) + '_' + str(semi_minor),
                                           display=True,
                                           title=str(filter) + ' ' + str(model),
                                           trim=50,
                                           aperture='Ellipse',
                                           semi_major=semi_major,
                                           semi_minor=semi_minor,
                                           theta=55.0,
                                           pad = pad
                                           )

            print(semi_major, semi_minor, flux[0])

            #major += 1.0
            #major = '{:.3}'.format(major)
            #major = float(major)
            theta += 0.1
        else:
            break

# measure the flux of the galaxy -> flux change testing
if flux_galaxy:
    count = 75.0
    while True:
        if count <= 75.0:
            # measure photometry of the galaxy -> display aperture fit
            # measure the flux
            flux = MAE_Aperture_Photometry(HDUlist=image,
                                           ext=1,
                                           radius=count,
                                           normalize=False,
                                           centroid=center,
                                           plot=True,
                                           print_params=False,
                                           save=out_path + '/' + str(filter) + '/Photometry_Tests/Flux_Conservation/galaxy_aper_' + str(count),
                                           display=True,
                                           title=str(filter) + ' ' + str(model),
                                           trim=100,
                                           pad = pad
                                           )

            print(count, flux[0])

            count += 1.0
            count = '{:.3}'.format(count)
            count = float(count)
        else:
            break

# measure the flux of the galaxy -> flux change testing
if aper21:
    count = 17.8
    while True:
        if count <= 17.8:
            # measure photometry of the galaxy -> display aperture fit
            # measure the flux
            flux = MAE_Aperture_Photometry(HDUlist=image,
                                           ext=1,
                                           radius=count,
                                           normalize=False,
                                           centroid=center,
                                           plot=True,
                                           print_params=False,
                                           save=out_path + '/' + str(filter) + '/Photometry_Tests/21um_Aper/21um_aper_' + str(count),
                                           display=True,
                                           title=str(filter) + ' ' + str(model),
                                           trim=100,
                                           pad = pad
                                           )

            print(count, flux[0])

            count += 1.0
            count = '{:.3}'.format(count)
            count = float(count)
        else:
            break



