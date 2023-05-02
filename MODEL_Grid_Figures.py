# 2022 Oct 27: code to generate various gridded figures
# 1) Merit function outputs for each deconvolution algorithm
# 2) 5x6 image grid comparing the MIRISim, RL, IWFT, Sparse, AIDA, and Kraken images
# 3) Display a single merit function grid for individual deconvolution results

from Convenience_Functions import *

# set the deconvolution directory
Dir = '7_sparse_condat_vu'

# set the model index:
idx = 3

# toggle which grid to generate
# display single merit function grid
merit_single = False
# display grid of merit functions: 2x5 (1) Richardson-Lucy, (2) Wiener-Hunt, (3) ME, (4) AIDA, (5) Kraken
merit_grid = True

# grid the flux values as difference between instead of absolutes?
delta = False
# normalize the counts measurements?
normalize = True

# include diffraction limits in merit function plots?
# NOTE: these need to be input by hand for both the merit_single and merit_grid methods
diff_lim = True

# display image grid
img = True
# trim image sizes?
large = True
small = False

# set the iteration limit (NOTE: set +10 above longest iteration)
# NOTE: if iterations < 50 -> set limit = 50
single_xlim = 100
# set the buffer for the diffraction limit tags
single_buffer = 15

# set the iteration limit as a function of technique (NOTE: set +10 above longest iteration for each)
# (1) Richardson-Lucy, (2) IWFT, (3) SParse, (4) AIDA, (5) Kraken
# NOTE: if iterations < 50 -> set limit = 50
multi_xlim = [110,110,110,110,110]
# set the buffer for the diffraction limit tags as a function of deconvolution technique
multi_buffer =[30,30,30,30,30]

# set the input params
# set the simulated wavelength
wave = [5.6, 10.0, 15.0, 18.0, 21.0]
# model type
model_type = 'Ideal_Model'
# model title: (1) single pixel, (2) single pixel + disk, (3) single pixel + disk + ionization cone, (4) complicated AGN
model = ['Model_single_pixel', 'Model_single_pixel_disk', 'Model_single_pixel_disk_bicone', 'Model_AGN_complicated']
# merit function time
title = ['Toy Unresolved Point Source Model', 'Toy Polar Dust and Unresolved Point Source Model',
         'Toy Extended Dust, Polar Dust, and Unresolved Point Source Model', 'Toy AGN Model']
# set the observation filter
filter = ['F560W', 'F1000W', 'F1500W', 'F1800W', 'F2100W']
# set image cmap
#cmap = ['turbo', 'inferno', 'coolwarm', 'viridis', 'afmhot']
cmap = ['afmhot']
# set the colorbar label
label = 'Log-scaled Image Peak (counts/pixel)'
# set the image sampling
sample = 'DET_SAMP'


# set Richardson-Lucy deconvolution iterations: (1) single pixel, (2) single pixel + disk, (3) single pixel + disk + ionization cone, (4) complicated AGN
RL_opt5 = [35, 38, 42, 48]
RL_opt10 = [80, 68, 50, 88]
RL_opt15 = [145, 109, 70, 83]
RL_opt18 = [151, 126, 92, 97]
RL_opt21 = [154, 111, 87, 77]

# set IWFT deconvolution iterations: (1) single pixel, (2) single pixel + disk, (3) single pixel + disk + ionization cone, (4) complicated AGN
# NOTE: for IWFT results add +2 to the selected number
# NOTE: correct index = index +2
# idx0 = header, idx1 = mirisim, etc..
WH_opt5 = [80,54,40,93]
WH_opt10 = [120,56,47,78]
WH_opt15 = [78,85,64,91]
WH_opt18 = [99,102,48,98]
WH_opt21 = [100,94,61,89]

# set Sparse deconvolution iterations: (1) single pixel, (2) single pixel + disk, (3) single pixel + disk + ionization cone, (4) complicated AGN
ME_opt5 = [26, 25, 78, 82]
ME_opt10 = [30, 34, 73, 83]
ME_opt15 = [21, 60, 111, 89]
ME_opt18 = [26, 28, 100, 95]
ME_opt21 = [17, 30, 98, 85]

# set AIDA deconvolution iterations: (1) single pixel, (2) single pixel + disk, (3) single pixel + disk + ionization cone, (4) complicated AGN
AIDA_opt5 = [18, 13, 42, 27]
AIDA_opt10 = [54, 20, 24, 40]
AIDA_opt15 = [34, 33, 42, 41]
AIDA_opt18 = [51, 13, 45, 33]
AIDA_opt21 = [41, 14, 37, 36]

# set Kraken deconvolution iterations: (1) single pixel, (2) single pixel + disk, (3) single pixel + disk + ionization cone, (4) complicated AGN
kraken_opt5 = [17, 13, 22, 15]
kraken_opt10 = [21, 13, 33, 27]
kraken_opt15 = [25, 47, 18, 19]
kraken_opt18 = [29, 20, 23, 31]
kraken_opt21 = [37, 20, 23, 13]

# if called: plot grid of merit functions
if merit_grid:
    # set the output path
    out_path = 'Images/JWST/4_Deconvolution/Total_Comparisons'

    # set the directories and iterate through them: RL, IWFT, Sparse, AIDA, Kraken
    directories = ['1_Richardson_Lucy', '6_IWFT', '7_sparse_condat_vu', '4_AIDA', '5_Kraken']

    # set the output file name convention
    save_merit = out_path + '/' + str(model_type) + '/Merit/SCALED_' + str(model[idx])

    nintA = []
    flux5 = []
    flux10 = []
    flux15 = []
    flux18 = []
    flux21 = []
    fwhm5 = []
    fwhm10 = []
    fwhm15 = []
    fwhm18 = []
    fwhm21 = []

    # dummy index
    # append flux values to dummy index
    dummyA = []
    dummyB = []
    dummyC = []
    dummyD = []
    dummyE = []

    for i in range(0, int(len(directories))):
        # set the output path for the merit function measurements
        mes_path = 'Images/JWST/4_Deconvolution/' + str(directories[i]) + '/Measurements/' + str(model_type) + '/' + str(model[idx])

        # read in the merit function data for each simulated waveband
        nint5, flux_5, fwhm_5 = np.loadtxt(mes_path + '/F560W/RL_Image_Plot_' + str(model[idx]) + '_' + str(sample)+'_SCALED.txt',unpack=True)
        nint10, flux_10, fwhm_10 = np.loadtxt(mes_path + '/F1000W/RL_Image_Plot_' + str(model[idx]) + '_' + str(sample)+'_SCALED.txt', unpack=True)
        nint15, flux_15, fwhm_15 = np.loadtxt(mes_path + '/F1500W/RL_Image_Plot_' + str(model[idx]) + '_' + str(sample)+'_SCALED.txt', unpack=True)
        nint18, flux_18, fwhm_18 = np.loadtxt(mes_path + '/F1800W/RL_Image_Plot_' + str(model[idx]) + '_' + str(sample)+'_SCALED.txt', unpack=True)
        nint21, flux_21, fwhm_21 = np.loadtxt(mes_path + '/F2100W/RL_Image_Plot_' + str(model[idx]) + '_' + str(sample)+'_SCALED.txt', unpack=True)

        # normalize the flux counts?
        if normalize:
            A = flux_5
            B = flux_10
            C = flux_15
            D = flux_18
            E = flux_21

            flux_5 = []
            flux_10 = []
            flux_15 = []
            flux_18 = []
            flux_21 = []
            for x in range(0, multi_xlim[0]):
                # normalized to the MIRISim value?
                normA = A[x] / A[0]
                normB = B[x] / B[0]
                normC = C[x] / C[0]
                normD = D[x] / D[0]
                normE = E[x] / E[0]

                flux_5.append(normA)
                flux_10.append(normB)
                flux_15.append(normC)
                flux_18.append(normD)
                flux_21.append(normE)

        nintA.append(nint5)
        flux5.append(flux_5)
        flux10.append(flux_10)
        flux15.append(flux_15)
        flux18.append(flux_18)
        flux21.append(flux_21)
        fwhm5.append(fwhm_5)
        fwhm10.append(fwhm_10)
        fwhm15.append(fwhm_15)
        fwhm18.append(fwhm_18)
        fwhm21.append(fwhm_21)

        nintAa = []
        flux_5 = []
        flux_10 = []
        flux_15 = []
        flux_18 = []
        flux_21 = []
        fwhm_5 = []
        fwhm_10 = []
        fwhm_15 = []
        fwhm_18 = []
        fwhm_21 = []
        for i in range(0, int(len(flux5))):
            ninta = []
            flux_5a = []
            flux_10a = []
            flux_15a = []
            flux_18a = []
            flux_21a = []
            fwhm_5a = []
            fwhm_10a = []
            fwhm_15a = []
            fwhm_18a = []
            fwhm_21a = []
            for j in range(0, multi_xlim[i]):

                # convert nint from float -> int
                nint1 = nintA[i][j]
                nint1 = int(nint1)

                # read in FWHM data
                fwhmA = fwhm5[i][j]
                fwhmB = fwhm10[i][j]
                fwhmC = fwhm15[i][j]
                fwhmD = fwhm18[i][j]
                fwhmE = fwhm21[i][j]

                # trim float values
                Afwhm = '{:.3f}'.format(fwhmA)
                Bfwhm = '{:.3f}'.format(fwhmB)
                Cfwhm = '{:.3f}'.format(fwhmC)
                Dfwhm = '{:.3f}'.format(fwhmD)
                Efwhm = '{:.3f}'.format(fwhmE)

                # append everything
                fwhm_5a.append(fwhmA)
                fwhm_10a.append(fwhmB)
                fwhm_15a.append(fwhmC)
                fwhm_18a.append(fwhmD)
                fwhm_21a.append(fwhmE)

                # read in flux data
                fluxA = flux5[i][j]
                fluxB = flux10[i][j]
                fluxC = flux15[i][j]
                fluxD = flux18[i][j]
                fluxE = flux21[i][j]

                # toggle between absolute flux values or delta flux values
                if delta:
                    if j == 0:
                        # append current flux values to dummy index
                        dummyA.append(fluxA)
                        dummyB.append(fluxB)
                        dummyC.append(fluxC)
                        dummyD.append(fluxD)
                        dummyE.append(fluxE)
                        fluxA = 0
                        fluxB = 0
                        fluxC = 0
                        fluxD = 0
                        fluxE = 0
                    else:
                        # append current flux values to dummy index
                        dummyA.append(fluxA)
                        dummyB.append(fluxB)
                        dummyC.append(fluxC)
                        dummyD.append(fluxD)
                        dummyE.append(fluxE)

                        # calculate delta values
                        ratA = np.abs(fluxA/dummyA[j-1])
                        ratB = np.abs(fluxB/dummyB[j-1])
                        ratC = np.abs(fluxC/dummyC[j-1])
                        ratD = np.abs(fluxD/dummyD[j-1])
                        ratE = np.abs(fluxE/dummyE[j-1])

                        fluxA = ratA
                        fluxB = ratB
                        fluxC = ratC
                        fluxD = ratD
                        fluxE = ratE

                # trim float value for plotting merit functions
                Aflux = '{:.3f}'.format(fluxA)
                Bflux = '{:.3f}'.format(fluxB)
                Cflux = '{:.3f}'.format(fluxC)
                Dflux = '{:.3f}'.format(fluxD)
                Eflux = '{:.3f}'.format(fluxE)

                # append everything
                ninta.append(nint1)
                flux_5a.append(fluxA)
                flux_10a.append(fluxB)
                flux_15a.append(fluxC)
                flux_18a.append(fluxD)
                flux_21a.append(fluxE)

            # append each full list here
            nintAa.append(ninta)
            flux_5.append(flux_5a)
            flux_10.append(flux_10a)
            flux_15.append(flux_15a)
            flux_18.append(flux_18a)
            flux_21.append(flux_21a)
            fwhm_5.append(fwhm_5a)
            fwhm_10.append(fwhm_10a)
            fwhm_15.append(fwhm_15a)
            fwhm_18.append(fwhm_18a)
            fwhm_21.append(fwhm_21a)

    # append the iterations number: 5.6 -> 21 um
    RL_iter = [RL_opt5[idx], RL_opt10[idx], RL_opt15[idx], RL_opt18[idx], RL_opt21[idx]]
    WH_iter = [WH_opt5[idx]-2, WH_opt10[idx]-2, WH_opt15[idx]-2, WH_opt18[idx]-2, WH_opt21[idx]-2]
    ME_iter = [ME_opt5[idx], ME_opt10[idx], ME_opt15[idx], ME_opt18[idx], ME_opt21[idx]]
    AIDA_iter = [AIDA_opt5[idx], AIDA_opt10[idx], AIDA_opt15[idx], AIDA_opt18[idx], AIDA_opt21[idx]]
    kraken_iter = [kraken_opt5[idx], kraken_opt10[idx], kraken_opt15[idx], kraken_opt18[idx], kraken_opt21[idx]]

    # plot the data
    # 2022 Oct 28-> features to add later (1) plot diffraction limit as a function of wavelength, (2) plot ratios (obs/dec?)
    # Create the 1x2 grid
    fig, ax = plt.subplots(nrows=5, ncols=2, figsize=(14, 20), tight_layout = True)
    nint0 = range(0, multi_xlim[0])
    nint1 = range(0, multi_xlim[1])
    nint2 = range(0, multi_xlim[2])
    nint3 = range(0, multi_xlim[3])
    nint4 = range(0, multi_xlim[4])


    if normalize == True:
        idx2 = [0, 1, 2, 3, 4]
        for i in range(0, int(len(idx2))):
            ax[i][1].set_ylabel('Normalized Aperture Counts', fontsize=15)
            ax[i][1].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

        # # set the ylims
        if idx == 0:
            idx2 = [0, 1, 2, 3, 4]
            for i in range(0, int(len(idx2))):
                ax[i][1].set_ylim(0.55, 1.4)
        elif idx == 1:
            idx2 = [0, 1, 2, 3, 4]
            for i in range(0, int(len(idx2))):
                ax[i][1].set_ylim(0.7, 1.4)
        elif idx == 2:
            idx2 = [0, 1, 2, 3, 4]
            for i in range(0, int(len(idx2))):
                ax[i][1].set_ylim(0.85, 1.35)
        else:
            idx2 = [0, 1, 2, 3, 4]
            for i in range(0, int(len(idx2))):
                ax[i][1].set_ylim(0.9, 1.45)

        save_merit = save_merit + '_norm'
    elif delta == True:
        idx2 = [0, 1, 2, 3, 4]
        for i in range(0, int(len(idx2))):
            ax[i][1].set_ylabel('ΔAperture Counts', fontsize=15)
            ax[i][1].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        save_merit = save_merit + '_delta'
    else:
        idx2 = [0, 1, 2, 3, 4]
        for i in range(0, int(len(idx2))):
            ax[i][1].set_ylabel('Aperture Counts', fontsize=15)
            ax[i][1].yaxis.set_major_formatter(FormatStrFormatter('%.2e'))
            ax[i][1].set_ylim(1e4,4e7)
        save_merit = save_merit


    # Plot the FWHM values
    # plt.subplot(121)
    # richardson-Lucy
    # 5 um
    ax[0][0].plot(nint1, fwhm_5[0], color='b',  label='F560W   (' + str(RL_iter[0]) + ' iterations)', linewidth=1)
    ax[0][0].plot(RL_iter[0], fwhm_5[0][RL_iter[0]], color='b', marker='X', markersize=15)
    # 10 um
    ax[0][0].plot(nint1, fwhm_10[0], color='g',  label='F1000W (' + str(RL_iter[1]) + ' iterations)', linewidth=1)
    ax[0][0].plot(RL_iter[1], fwhm_10[0][RL_iter[1]], color='g', marker='X', markersize=15)
    # 15 um
    ax[0][0].plot(nint1, fwhm_15[0], color='y',  label='F1500W (' + str(RL_iter[2]) + ' iterations)', linewidth=1)
    ax[0][0].plot(RL_iter[2], fwhm_15[0][RL_iter[2]], color='y', marker='X', markersize=15)
    # 18 um
    ax[0][0].plot(nint1, fwhm_18[0], color='darkorange',  label='F1800W (' + str(RL_iter[3]) + ' iterations)', linewidth=1)
    ax[0][0].plot(RL_iter[3], fwhm_18[0][RL_iter[3]], color='darkorange', marker='X', markersize=15)
    # 21 um
    ax[0][0].plot(nint1, fwhm_21[0], color='r',  label='F2100W (' + str(RL_iter[4]) + ' iterations)', linewidth=1)
    ax[0][0].plot(RL_iter[4], fwhm_21[0][RL_iter[4]], color='r', marker='X', markersize=15)
    ax[0][0].set_xlim(0, multi_xlim[1])
    ax[0][0].set_ylim(0, int(max(fwhm_21[0]))+2)
    ax[0][0].set_ylabel('FWHM (pixel)', fontsize = 15)
    ax[0][0].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax[0][0].text(5, int(max(fwhm_21[0])), 'Richardson-Lucy', fontsize = 20)
    ax[0][0].legend(loc='upper right')
    ax[0][0].set_title('FWHM Limit', fontsize = 15)

    # wiener-hunt
    # 5 um
    ax[1][0].plot(nint1, fwhm_5[1], color='b',  label='F560W   (' + str(WH_iter[0]) + ' iterations)', linewidth=1)
    #ax[1][0].plot(WH_iter[0], fwhm_5[1][WH_iter[0]], color='b', marker='X', markersize=15)
    # 10 um
    ax[1][0].plot(nint1, fwhm_10[1], color='g',  label='F1000W (' + str(WH_iter[1]) + ' iterations)', linewidth=1)
    ax[1][0].plot(WH_iter[1], fwhm_10[1][WH_iter[1]], color='g', marker='X', markersize=15)
    # 15 um
    ax[1][0].plot(nint1, fwhm_15[1], color='y',  label='F1500W (' + str(WH_iter[2]) + ' iterations)', linewidth=1)
    ax[1][0].plot(WH_iter[2], fwhm_15[1][WH_iter[2]], color='y', marker='X', markersize=15)
    # 18 um
    ax[1][0].plot(nint1, fwhm_18[1], color='darkorange',  label='F1800W (' + str(WH_iter[3]) + ' iterations)', linewidth=1)
    ax[1][0].plot(WH_iter[3], fwhm_18[1][WH_iter[3]], color='darkorange', marker='X', markersize=15)
    # 21 um
    ax[1][0].plot(nint1, fwhm_21[1], color='r',  label='F2100W (' + str(WH_iter[4]) + ' iterations)', linewidth=1)
    ax[1][0].plot(WH_iter[4], fwhm_21[1][WH_iter[4]], color='r', marker='X', markersize=15)
    ax[1][0].set_xlim(0, multi_xlim[1])
    ax[1][0].set_ylim(0, int(max(fwhm_21[0]))+2)
    ax[1][0].set_ylabel('FWHM (pixel)', fontsize = 15)
    ax[1][0].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax[1][0].text(5, int(max(fwhm_21[0])), 'IWFT', fontsize = 20)
    ax[1][0].legend(loc='upper right')

    # multiscale maximum entropy
    # 5 um
    ax[2][0].plot(nint2, fwhm_5[2], color='b',  label='F560W   (' + str(ME_iter[0]) + ' iterations)', linewidth=1)
    ax[2][0].plot(ME_iter[0], fwhm_5[2][ME_iter[0]], color='b', marker='X', markersize=15)
    # 10 um
    ax[2][0].plot(nint2, fwhm_10[2], color='g',  label='F1000W (' + str(ME_iter[1]) + ' iterations)', linewidth=1)
    ax[2][0].plot(ME_iter[1], fwhm_10[2][ME_iter[1]], color='g', marker='X', markersize=15)
    # 15 um
    ax[2][0].plot(nint2, fwhm_15[2], color='y',  label='F1500W (' + str(ME_iter[2]) + ' iterations)', linewidth=1)
    ax[2][0].plot(ME_iter[2], fwhm_15[2][ME_iter[2]], color='y', marker='X', markersize=15)
    # 18 um
    ax[2][0].plot(nint2, fwhm_18[2], color='darkorange',  label='F1800W (' + str(ME_iter[3]) + ' iterations)', linewidth=1)
    ax[2][0].plot(ME_iter[3], fwhm_18[2][ME_iter[3]], color='darkorange', marker='X', markersize=15)
    # 21 um
    ax[2][0].plot(nint2, fwhm_21[2], color='r',  label='F2100W (' + str(ME_iter[4]) + ' iterations)', linewidth=1)
    ax[2][0].plot(ME_iter[4], fwhm_21[2][ME_iter[4]], color='r', marker='X', markersize=15)
    ax[2][0].set_xlim(0, multi_xlim[2])
    ax[2][0].set_ylim(0, int(max(fwhm_21[0]))+2)
    ax[2][0].set_ylabel('FWHM (pixel)', fontsize = 15)
    ax[2][0].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax[2][0].text(5, int(max(fwhm_21[0])), 'Sparse', fontsize = 20)
    ax[2][0].legend(loc='upper right')

    # AIDA
    # 5 um
    ax[3][0].plot(nint3, fwhm_5[3], color='b',  label='F560W   (' + str(AIDA_iter[0]) + ' iterations)', linewidth=1)
    ax[3][0].plot(AIDA_iter[0], fwhm_5[3][AIDA_iter[0]], color='b', marker='X', markersize=15)
    # 10 um
    ax[3][0].plot(nint3, fwhm_10[3], color='g',  label='F1000W (' + str(AIDA_iter[1]) + ' iterations)', linewidth=1)
    ax[3][0].plot(AIDA_iter[1], fwhm_10[3][AIDA_iter[1]], color='g', marker='X', markersize=15)
    # 15 um
    ax[3][0].plot(nint3, fwhm_15[3], color='y',  label='F1500W (' + str(AIDA_iter[2]) + ' iterations)', linewidth=1)
    ax[3][0].plot(AIDA_iter[2], fwhm_15[3][AIDA_iter[2]], color='y', marker='X', markersize=15)
    # 18 um
    ax[3][0].plot(nint3, fwhm_18[3], color='darkorange',  label='F1800W (' + str(AIDA_iter[3]) + ' iterations)', linewidth=1)
    ax[3][0].plot(AIDA_iter[3], fwhm_18[3][AIDA_iter[3]], color='darkorange', marker='X', markersize=15)
    # 21 um
    ax[3][0].plot(nint3, fwhm_21[3], color='r',  label='F2100W (' + str(AIDA_iter[4]) + ' iterations)', linewidth=1)
    ax[3][0].plot(AIDA_iter[4], fwhm_21[3][AIDA_iter[4]], color='r', marker='X', markersize=15)
    ax[3][0].set_xlim(0, multi_xlim[3])
    ax[3][0].set_ylim(0, int(max(fwhm_21[0]))+2)
    ax[3][0].set_ylabel('FWHM (pixel)', fontsize = 15)
    ax[3][0].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax[3][0].text(5, int(max(fwhm_21[0])), 'AIDA', fontsize = 20)
    ax[3][0].legend(loc='upper right')

    # kraken
    # 5 um
    ax[4][0].plot(nint4, fwhm_5[4], color='b',  label='F560W   (' + str(kraken_iter[0]) + ' iterations)', linewidth=1)
    ax[4][0].plot(kraken_iter[0], fwhm_5[4][kraken_iter[0]], color='b', marker='X', markersize=15)
    # 10 um
    ax[4][0].plot(nint4, fwhm_10[4], color='g',  label='F1000W (' + str(kraken_iter[1]) + ' iterations)', linewidth=1)
    ax[4][0].plot(kraken_iter[1], fwhm_10[4][kraken_iter[1]], color='g', marker='X', markersize=15)
    # 15 um
    ax[4][0].plot(nint4, fwhm_15[4], color='y',  label='F1500W (' + str(kraken_iter[2]) + ' iterations)', linewidth=1)
    ax[4][0].plot(kraken_iter[2], fwhm_15[4][kraken_iter[2]], color='y', marker='X', markersize=15)
    # 18 um
    ax[4][0].plot(nint4, fwhm_18[4], color='darkorange',  label='F1800W (' + str(kraken_iter[3]) + ' iterations)', linewidth=1)
    ax[4][0].plot(kraken_iter[3], fwhm_18[4][kraken_iter[3]], color='darkorange', marker='X', markersize=15)
    # 21 um
    ax[4][0].plot(nint4, fwhm_21[4], color='r',  label='F2100W (' + str(kraken_iter[4]) + ' iterations)', linewidth=1)
    ax[4][0].plot(kraken_iter[4], fwhm_21[4][kraken_iter[4]], color='r', marker='X', markersize=15)
    ax[4][0].set_xlim(0, multi_xlim[4])
    ax[4][0].set_ylim(0, int(max(fwhm_21[0]))+2)
    ax[4][0].set_ylabel('FWHM (pixel)', fontsize = 15)
    ax[4][0].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax[4][0].text(5, int(max(fwhm_21[0])), 'Kraken', fontsize = 20)
    ax[4][0].legend(loc='upper right')

    if diff_lim:
        save_merit = save_merit + '_diff_lim'

        # richardson-lucy
        ax[0][0].hlines(1.766, 2, 7, color='k', linestyle='--')
        ax[0][0].text(10, 1.766, 'F560W PSF FWHM', color='k')
        ax[0][0].hlines(2.891, 3,8, color='k', linestyle = '--')
        ax[0][0].text(11, 2.891, 'F1000W PSF FWHM', color='k')
        ax[0][0].hlines(4.39, 4, 9, color = 'k', linestyle = '--')
        ax[0][0].text(12, 4.39, 'F1500W PSF FWHM', color='k')
        ax[0][0].hlines(5.353, 3, 8, color='k', linestyle = '--')
        ax[0][0].text(11, 5.353, 'F1800W PSF FWHM', color='k')
        ax[0][0].hlines(5.996, 5, 10, color='k', linestyle = '--')
        ax[0][0].text(13, 5.996, 'F2100W PSF FWHM', color='k')

        # Wiener-Hunt
        ax[1][0].hlines(1.878, 40, 45, color='k', linestyle='--')
        ax[1][0].text(48, 1.878, 'F560W PSF FWHM', color='k')
        ax[1][0].hlines(2.975, 30,35, color='k', linestyle = '--')
        ax[1][0].text(38, 2.975, 'F1000W PSF FWHM', color='k')
        ax[1][0].hlines(4.422, 39, 44, color = 'k', linestyle = '--')
        ax[1][0].text(47, 4.422, 'F1500W PSF FWHM', color='k')
        ax[1][0].hlines(5.35, 36, 41, color='k', linestyle = '--')
        ax[1][0].text(44, 5.35, 'F1800W PSF FWHM', color='k')
        ax[1][0].hlines(6.101, 42,47, color='k', linestyle = '--')
        ax[1][0].text(50, 6.101, 'F2100W PSF FWHM', color='k')

        # Max. Entropy
        ax[2][0].hlines(1.868, 28, 33, color='k', linestyle='--')
        ax[2][0].text(36, 1.868, 'F560W PSF FWHM', color='k')
        ax[2][0].hlines(2.966, 24, 29, color='k', linestyle = '--')
        ax[2][0].text(32, 2.966, 'F1000W PSF FWHM', color='k')
        ax[2][0].hlines(4.433, 29, 32, color = 'k', linestyle = '--')
        ax[2][0].text(35, 4.433, 'F1500W PSF FWHM', color='k')
        ax[2][0].hlines(5.365, 26, 31, color='k', linestyle = '--')
        ax[2][0].text(34, 5.365, 'F1800W PSF FWHM', color='k')
        ax[2][0].hlines(6.121, 37,42, color='k', linestyle = '--')
        ax[2][0].text(45, 6.121, 'F2100W PSF FWHM', color='k')

        # AIDA
        ax[3][0].hlines(1.875, 30, 35, color='k', linestyle='--')
        ax[3][0].text(38, 1.875, 'F560W PSF FWHM', color='k')
        ax[3][0].hlines(2.779, 50, 55, color='k', linestyle = '--')
        ax[3][0].text(58, 2.779, 'F1000W PSF FWHM', color='k')
        ax[3][0].hlines(4.427, 12, 17, color = 'k', linestyle = '--')
        ax[3][0].text(20, 4.427, 'F1500W PSF FWHM', color='k')
        ax[3][0].hlines(5.366, 6, 11, color='k', linestyle = '--')
        ax[3][0].text(14, 5.366, 'F1800W PSF FWHM', color='k')
        ax[3][0].hlines(6.11, 9, 12, color='k', linestyle = '--')
        ax[3][0].text(15, 6.11, 'F2100W PSF FWHM', color='k')

        # Kraken
        diff_lim_new = [1.882, 2.982, 4.436, 5.373, 6.127]
        ax[4][0].hlines(1.748, 1, 6, color='k', linestyle='--')
        ax[4][0].text(9, 1.748, 'F560W PSF FWHM', color='k')
        ax[4][0].hlines(2.67, 2,7, color='k', linestyle = '--')
        ax[4][0].text(10, 2.67, 'F1000W PSF FWHM', color='k')
        ax[4][0].hlines(4.162, 1, 6, color = 'k', linestyle = '--')
        ax[4][0].text(9, 4.162, 'F1500W PSF FWHM', color='k')
        ax[4][0].hlines(4.823, 0, 5, color='k', linestyle = '--')
        ax[4][0].text(8, 4.823, 'F1800W PSF FWHM', color='k')
        ax[4][0].hlines(5.563, 0,5, color='k', linestyle = '--')
        ax[4][0].text(8, 5.563, 'F2100W PSF FWHM', color='k')

    # richardson-Lucy
    # 5 um
    ax[0][1].plot(nint0, flux_5[0], color='b', label='F560W   (' + str(RL_iter[0]) + ' iterations)', linewidth=1)
    # ax[0][1].plot(RL_iter[0], flux_5[0][RL_iter[0]], color='b', marker='X', markersize=15)
    # 10 um
    ax[0][1].plot(nint0, flux_10[0], color='g', label='F1000W (' + str(RL_iter[1]) + ' iterations)', linewidth=1)
    # ax[0][1].plot(RL_iter[1], flux_10[0][RL_iter[1]], color='g', marker='X', markersize=15)
    # 15 um
    ax[0][1].plot(nint0, flux_15[0], color='y', label='F1500W (' + str(RL_iter[2]) + ' iterations)', linewidth=1)
    # ax[0][1].plot(RL_iter[2], flux_15[0][RL_iter[2]], color='y', marker='X', markersize=15)
    # 18 um
    ax[0][1].plot(nint0, flux_18[0], color='darkorange', label='F1800W (' + str(RL_iter[3]) + ' iterations)',
                      linewidth=1)
    # ax[0][1].plot(RL_iter[3], flux_18[0][RL_iter[3]], color='darkorange', marker='X', markersize=15)
    # 21 um
    ax[0][1].plot(nint0, flux_21[0], color='r', label='F2100W (' + str(RL_iter[4]) + ' iterations)', linewidth=1)
    # ax[0][1].plot(RL_iter[4], flux_21[0][RL_iter[4]], color='r', marker='X', markersize=15)
    ax[0][1].set_title('Counts Conservation', fontsize=15)
    ax[0][1].set_xlim(0, multi_xlim[0])

    # wiener-hunt
    # 5 um
    ax[1][1].plot(nint1, flux_5[1], color='b', label='F560W   (' + str(WH_iter[0]) + ' iterations)', linewidth=1)
    ax[1][1].plot(WH_iter[0], flux_5[1][WH_iter[0]], color='b', marker='X', markersize=15)
    # 10 um
    ax[1][1].plot(nint1, flux_10[1], color='g', label='F1000W (' + str(WH_iter[1]) + ' iterations)', linewidth=1)
    #ax[1][1].plot(WH_iter[1], flux_10[1][WH_iter[1]], color='g', marker='X', markersize=15)
    # 15 um
    ax[1][1].plot(nint1, flux_15[1], color='y', label='F1500W (' + str(WH_iter[2]) + ' iterations)', linewidth=1)
    #ax[1][1].plot(WH_iter[2], flux_15[1][WH_iter[2]], color='y', marker='X', markersize=15)
    # 18 um
    ax[1][1].plot(nint1, flux_18[1], color='darkorange', label='F1800W (' + str(WH_iter[3]) + ' iterations)', linewidth=1)
    #ax[1][1].plot(WH_iter[3], flux_18[1][WH_iter[3]], color='darkorange', marker='X', markersize=15)
    # 21 um
    ax[1][1].plot(nint1, flux_21[1], color='r', label='F2100W (' + str(WH_iter[4]) + ' iterations)', linewidth=1)
    #ax[1][1].plot(WH_iter[4], flux_21[1][WH_iter[4]], color='r', marker='X', markersize=15)
    ax[1][1].set_xlim(0, multi_xlim[1])

    # Sparse Condat-Vu
    # 5 um
    ax[2][1].plot(nint2, flux_5[2], color='b', label='F560W   (' + str(ME_iter[0]) + ' iterations)', linewidth=1)
    #ax[2][1].plot(ME_iter[0], flux_5[2][ME_iter[0]], color='b', marker='X', markersize=15)
    # 10 um
    ax[2][1].plot(nint2, flux_10[2], color='g', label='F1000W (' + str(ME_iter[1]) + ' iterations)', linewidth=1)
    # ax[2][1].plot(ME_iter[1], flux_10[2][ME_iter[1]], color='g', marker='X', markersize=15)
    # 15 um
    ax[2][1].plot(nint2, flux_15[2], color='y', label='F1500W (' + str(ME_iter[2]) + ' iterations)', linewidth=1)
    # ax[2][1].plot(ME_iter[2], flux_15[2][ME_iter[2]], color='y', marker='X', markersize=15)
    # 18 um
    ax[2][1].plot(nint2, flux_18[2], color='darkorange', label='F1800W (' + str(ME_iter[3]) + ' iterations)',
                      linewidth=1)
    # ax[2][1].plot(ME_iter[3], flux_18[2][ME_iter[3]], color='darkorange', marker='X', markersize=15)
    # 21 um
    ax[2][1].plot(nint2, flux_21[2], color='r', label='F2100W (' + str(ME_iter[4]) + ' iterations)', linewidth=1)
    # ax[2][1].plot(ME_iter[4], flux_21[2][ME_iter[4]], color='r', marker='X', markersize=15)
    ax[2][1].set_xlim(0, multi_xlim[2])

    # AIDA
    # 5 um
    ax[3][1].plot(nint3, flux_5[3], color='b', label='F560W   (' + str(AIDA_iter[0]) + ' iterations)', linewidth=1)
    # ax[3][1].plot(AIDA_iter[0], flux_5[3][AIDA_iter[0]], color='b', marker='X', markersize=15)
    # 10 um
    ax[3][1].plot(nint3, flux_10[3], color='g', label='F1000W (' + str(AIDA_iter[1]) + ' iterations)', linewidth=1)
    # ax[3][1].plot(AIDA_iter[1], flux_10[3][AIDA_iter[1]], color='g', marker='X', markersize=15)
    # 15 um
    ax[3][1].plot(nint3, flux_15[3], color='y', label='F1500W (' + str(AIDA_iter[2]) + ' iterations)', linewidth=1)
    # ax[3][1].plot(AIDA_iter[2], flux_15[3][AIDA_iter[2]], color='y', marker='X', markersize=15)
    # 18 um
    ax[3][1].plot(nint3, flux_18[3], color='darkorange', label='F1800W (' + str(AIDA_iter[3]) + ' iterations)',
                      linewidth=1)
    # ax[3][1].plot(AIDA_iter[3], flux_18[3][AIDA_iter[3]], color='darkorange', marker='X', markersize=15)
    # 21 um
    ax[3][1].plot(nint3, flux_21[3], color='r', label='F2100W (' + str(AIDA_iter[4]) + ' iterations)', linewidth=1)
    # ax[3][1].plot(AIDA_iter[4], flux_21[3][AIDA_iter[4]], color='r', marker='X', markersize=15)
    ax[3][1].set_xlim(0, multi_xlim[3])

    # kraken
    # 5 um
    ax[4][1].plot(nint4, flux_5[4], color='b', label='F560W   (' + str(kraken_iter[0]) + ' iterations)',linewidth=1)
    # ax[4][1].plot(kraken_iter[0], flux_5[4][kraken_iter[0]], color='b', marker='X', markersize=15)
    # 10 um
    ax[4][1].plot(nint4, flux_10[4], color='g', label='F1000W (' + str(kraken_iter[1]) + ' iterations)',linewidth=1)
    # ax[4][1].plot(kraken_iter[1], flux_10[4][kraken_iter[1]], color='g', marker='X', markersize=15)
    # 15 um
    ax[4][1].plot(nint4, flux_15[4], color='y', label='F1500W (' + str(kraken_iter[2]) + ' iterations)',
                      linewidth=1)
    # ax[4][1].plot(kraken_iter[2], flux_15[4][kraken_iter[2]], color='r', marker='X', markersize=15)
    # 18 um
    ax[4][1].plot(nint4, flux_18[4], color='darkorange', label='F1800W (' + str(kraken_iter[3]) + ' iterations)',
                      linewidth=1)
    # ax[4][1].plot(kraken_iter[3], flux_18[4][kraken_iter[3]], color='c', marker='X', markersize=15)
    # 21 um
    ax[4][1].plot(nint4, flux_21[4], color='r', label='F2100W (' + str(kraken_iter[4]) + ' iterations)',
                      linewidth=1)
    # ax[4][1].plot(kraken_iter[4], flux_21[4][kraken_iter[4]], color='m', marker='X', markersize=15)
    ax[4][1].set_xlim(0, multi_xlim[4])

    fig.supxlabel('Iterations', fontsize=15)
    plt.suptitle(title[idx], fontsize=25, y=0.99)
    plt.savefig(save_merit + '.png')
    plt.show()

# if called: plot the image grid
if img:
    # set the path for the reference psf
    psf_path = 'Images/JWST/1_Input_Models/PSFs/'
    # set the output path
    out_path = 'Images/JWST/4_Deconvolution/Total_Comparisons/'+str(model_type)+'/Images/'
    # set the output file name convention
    save_img = out_path + str(model[idx]) +'/SCALED_Image_Grid_' + str(model[idx])

    # set the import path for the Richardson-Lucy deconvolved image FITS file
    RL_fits_path = 'Images/JWST/4_Deconvolution/1_Richardson_Lucy/FITS/' + str(model_type) + '/' + str(model[idx])
    # set the import path for the Wiener-Hunt deconvolved image FITS file
    WH_fits_path = 'Images/JWST/4_Deconvolution/6_IWFT/FITS/' + str(model_type) + '/' + str(model[idx])
    # set the import path for the SParse Condat-Vu deconvolved image FITS file
    ME_fits_path = 'Images/JWST/4_Deconvolution/7_sparse_condat_vu/FITS/' + str(model_type) + '/' + str(model[idx])
    # set the import path for the AIDA deconvolved image FITS file
    AIDA_fits_path = 'Images/JWST/4_Deconvolution/4_AIDA/FITS/' + str(model_type) + '/' + str(model[idx])
    # set the import path for the kraken MFBD deconvolved image FITS file
    kraken_fits_path = 'Images/JWST/4_Deconvolution/5_Kraken/FITS/' + str(model_type) + '/' + str(model[idx])

    # MIRISim observations
    image5 = get_pkg_data_filename(RL_fits_path + '/SCALED_F560W_' + str(model[idx]) + '_DET_SAMP.fits')
    image10 = get_pkg_data_filename(RL_fits_path + '/SCALED_F1000W_' + str(model[idx]) + '_DET_SAMP.fits')
    image15 = get_pkg_data_filename(RL_fits_path + '/SCALED_F1500W_' + str(model[idx]) + '_DET_SAMP.fits')
    image18 = get_pkg_data_filename(RL_fits_path + '/SCALED_F1800W_' + str(model[idx]) + '_DET_SAMP.fits')
    image21 = get_pkg_data_filename(RL_fits_path + '/SCALED_F2100W_' + str(model[idx]) + '_DET_SAMP.fits')

    # Read in the Wiener-Hunt deconvolved image data
    WHimg5 = get_pkg_data_filename(WH_fits_path + '/SCALED_F560W_' + str(model[idx]) + '_DET_SAMP.fits')
    WHimg10 = get_pkg_data_filename(WH_fits_path + '/SCALED_F1000W_' + str(model[idx]) + '_DET_SAMP.fits')
    WHimg15 = get_pkg_data_filename(WH_fits_path + '/SCALED_F1500W_' + str(model[idx]) + '_DET_SAMP.fits')
    WHimg18 = get_pkg_data_filename(WH_fits_path + '/SCALED_F1800W_' + str(model[idx]) + '_DET_SAMP.fits')
    WHimg21 = get_pkg_data_filename(WH_fits_path + '/SCALED_F2100W_' + str(model[idx]) + '_DET_SAMP.fits')

    # Read in the MEM deconvolved image data
    MEimg5 = get_pkg_data_filename(ME_fits_path + '/SCALED_F560W_' + str(model[idx]) + '_DET_SAMP.fits')
    MEimg10 = get_pkg_data_filename(ME_fits_path + '/SCALED_F1000W_' + str(model[idx]) + '_DET_SAMP.fits')
    MEimg15 = get_pkg_data_filename(ME_fits_path + '/SCALED_F1500W_' + str(model[idx]) + '_DET_SAMP.fits')
    MEimg18 = get_pkg_data_filename(ME_fits_path + '/SCALED_F1800W_' + str(model[idx]) + '_DET_SAMP.fits')
    MEimg21 = get_pkg_data_filename(ME_fits_path + '/SCALED_F2100W_' + str(model[idx]) + '_DET_SAMP.fits')

    # # Read in the AIDA deconvolved image data
    AIDAimg5 = get_pkg_data_filename(AIDA_fits_path + '/SCALED_F560W_' + str(model[idx]) + '_DET_SAMP.fits')
    AIDAimg10 = get_pkg_data_filename(AIDA_fits_path + '/SCALED_F1000W_' + str(model[idx]) + '_DET_SAMP.fits')
    AIDAimg15 = get_pkg_data_filename(AIDA_fits_path + '/SCALED_F1500W_' + str(model[idx]) + '_DET_SAMP.fits')
    AIDAimg18 = get_pkg_data_filename(AIDA_fits_path + '/SCALED_F1800W_' + str(model[idx]) + '_DET_SAMP.fits')
    AIDAimg21 = get_pkg_data_filename(AIDA_fits_path + '/SCALED_F2100W_' + str(model[idx]) + '_DET_SAMP.fits')

    # Read in the Kraken MFBD deconvolved image data
    kraken_img5 = get_pkg_data_filename(kraken_fits_path + '/SCALED_F560W_' + str(model[idx]) + '_DET_SAMP.fits')
    kraken_img10 = get_pkg_data_filename(kraken_fits_path + '/SCALED_F1000W_' + str(model[idx]) + '_DET_SAMP.fits')
    kraken_img15 = get_pkg_data_filename(kraken_fits_path + '/SCALED_F1500W_' + str(model[idx]) + '_DET_SAMP.fits')
    kraken_img18 = get_pkg_data_filename(kraken_fits_path + '/SCALED_F1800W_' + str(model[idx]) + '_DET_SAMP.fits')
    kraken_img21 = get_pkg_data_filename(kraken_fits_path + '/SCALED_F2100W_' + str(model[idx]) + '_DET_SAMP.fits')

    # get the necessary RL image data
    obs_list5 = fits.open(image5)
    obs_list10 = fits.open(image10)
    obs_list15 = fits.open(image15)
    obs_list18 = fits.open(image18)
    obs_list21 = fits.open(image21)

    # get the necessary WH image data
    WH_list5 = fits.open(WHimg5)
    WH_list10 = fits.open(WHimg10)
    WH_list15 = fits.open(WHimg15)
    WH_list18 = fits.open(WHimg18)
    WH_list21 = fits.open(WHimg21)

    # get the necessary Sparse Condat-Vu image data
    ME_list5 = fits.open(MEimg5)
    ME_list10 = fits.open(MEimg10)
    ME_list15 = fits.open(MEimg15)
    ME_list18 = fits.open(MEimg18)
    ME_list21 = fits.open(MEimg21)

    # get the necessary AIDA image data
    AIDA_list5 = fits.open(AIDAimg5)
    AIDA_list10 = fits.open(AIDAimg10)
    AIDA_list15 = fits.open(AIDAimg15)
    AIDA_list18 = fits.open(AIDAimg18)
    AIDA_list21 = fits.open(AIDAimg21)

    # get the necessary kraken image data
    kraken_list5 = fits.open(kraken_img5)
    kraken_list10 = fits.open(kraken_img10)
    kraken_list15 = fits.open(kraken_img15)
    kraken_list18 = fits.open(kraken_img18)
    kraken_list21 = fits.open(kraken_img21)

    # set the MIRISim image
    obs5 = obs_list5[1].data
    obs10 = obs_list10[1].data
    obs15 = obs_list15[1].data
    obs18 = obs_list18[1].data
    obs21 = obs_list21[1].data

    # Set the optimum Richardson-Lucy deconvolved image
    RL5 = obs_list5[RL_opt5[idx]].data
    RL10 = obs_list10[RL_opt10[idx]].data
    RL15 = obs_list15[RL_opt15[idx]].data
    RL18 = obs_list18[RL_opt18[idx]].data
    RL21 = obs_list21[RL_opt21[idx]].data

    # Set the optimum Wiener-Hunt deconvolved image
    WH5 = WH_list5[WH_opt5[idx]].data
    WH10 = WH_list10[WH_opt10[idx]].data
    WH15 = WH_list15[WH_opt15[idx]].data
    WH18 = WH_list18[WH_opt18[idx]].data
    WH21 = WH_list21[WH_opt21[idx]].data

    # Set the optimum MEM deconvolved image <- DEFINETLY CHANGE THIS
    ME5 = ME_list5[ME_opt5[idx]].data
    ME10 = ME_list10[ME_opt10[idx]].data
    ME15 = ME_list15[ME_opt15[idx]].data
    ME18 = ME_list18[ME_opt18[idx]].data
    ME21 = ME_list21[ME_opt21[idx]].data

    # Set the optimum AIDA deconvolved image
    AIDA5 = AIDA_list5[AIDA_opt5[idx]].data
    AIDA10 = AIDA_list10[AIDA_opt10[idx]].data
    AIDA15 = AIDA_list15[AIDA_opt15[idx]].data
    AIDA18 = AIDA_list18[AIDA_opt18[idx]].data
    AIDA21 = AIDA_list21[AIDA_opt21[idx]].data

    # Set the optimum Kraken MFBD deconvolved image
    kraken5 = kraken_list5[kraken_opt5[idx]].data
    kraken10 = kraken_list10[kraken_opt10[idx]].data
    kraken15 = kraken_list15[kraken_opt15[idx]].data
    kraken18 = kraken_list18[kraken_opt18[idx]].data
    kraken21 = kraken_list21[kraken_opt21[idx]].data

    # append the iterations number: 5.6 -> 21 um
    RL_iter = [RL_opt5[idx], RL_opt10[idx], RL_opt15[idx], RL_opt18[idx], RL_opt21[idx]]
    # NOTE: the correct IWFT index is -2 from the value listed above. The merit functions and FITS cube indices are different
    WH_iter = [WH_opt5[idx]-2, WH_opt10[idx]-2, WH_opt15[idx]-2, WH_opt18[idx]-2, WH_opt21[idx]-2]
    ME_iter = [ME_opt5[idx], ME_opt10[idx], ME_opt15[idx], ME_opt18[idx], ME_opt21[idx]]
    AIDA_iter = [AIDA_opt5[idx], AIDA_opt10[idx], AIDA_opt15[idx], AIDA_opt18[idx], AIDA_opt21[idx]]
    kraken_iter = [kraken_opt5[idx], kraken_opt10[idx], kraken_opt15[idx], kraken_opt18[idx], kraken_opt21[idx]]

    # set text locations according to image sizes
    xlim = 45
    ylim = 220
    size = 'Full'

    # set the x/y position labels
    positions = [0, 64, 128, 192, 255]
    labels = ['-128', '-64', '0', '64', '128']

    if large:
        # arbitraly trim the image data to a 150x150 image grid
        # MIRISim images
        obs5 = obs5[27:227, 29:229]
        obs10 = obs10[27:227, 29:229]
        obs15 = obs15[27:227, 29:229]
        obs18 = obs18[27:227, 29:229]
        obs21 = obs21[27:227, 29:229]

        # Richardson-Lucy images
        RL5 = RL5[27:227, 29:229]
        RL10 = RL10[27:227, 29:229]
        RL15 = RL15[27:227, 29:229]
        RL18 = RL18[27:227, 29:229]
        RL21 = RL21[27:227, 29:229]

        # # Wiener-Hunt images
        WH5 = WH5[27:227, 29:229]
        WH10 = WH10[27:227, 29:229]
        WH15 = WH15[27:227, 29:229]
        WH18 = WH18[27:227, 29:229]
        WH21 = WH21[27:227, 29:229]

        # ME images
        ME5 = ME5[27:227, 29:229]
        ME10 = ME10[27:227, 29:229]
        ME15 = ME15[27:227, 29:229]
        ME18 = ME18[27:227, 29:229]
        ME21 = ME21[27:227, 29:229]

        # AIDA images
        AIDA5 = AIDA5[27:227, 29:229]
        AIDA10 = AIDA10[27:227, 29:229]
        AIDA15 = AIDA15[27:227, 29:229]
        AIDA18 = AIDA18[27:227, 29:229]
        AIDA21 = AIDA21[27:227, 29:229]

        # kraken images
        kraken5 = kraken5[27:227, 29:229]
        kraken10 = kraken10[27:227, 29:229]
        kraken15 = kraken15[27:227, 29:229]
        kraken18 = kraken18[27:227, 29:229]
        kraken21 = kraken21[27:227, 29:229]

        # set the text areas
        xlim = 43
        ylim = 175
        size = 'Large'

        # set the x/y position labels
        labels = ['-100', '-50', '0', '50', '100']
        positions = [10, 50, 100, 150, 190]

    if small:
        # arbitraly trim the image data to a 50x50 image grid
        # MIRISim images
        obs5 = obs5[103:153, 103:153]
        obs10 = obs10[103:153, 103:153]
        obs15 = obs15[103:153, 103:153]
        obs18 = obs18[103:153, 103:153]
        obs21 = obs21[103:153, 103:153]

        # Richardson-Lucy images
        RL5 = RL5[103:153, 103:153]
        RL10 = RL10[103:153, 103:153]
        RL15 = RL15[103:153, 103:153]
        RL18 = RL18[103:153, 103:153]
        RL21 = RL21[103:153, 103:153]

        # Wiener-Hunt images
        WH5 = WH5[103:153, 103:153]
        WH10 = WH10[103:153, 103:153]
        WH15 = WH15[103:153, 103:153]
        WH18 = WH18[103:153, 103:153]
        WH21 = WH21[103:153, 103:153]

        # ME images
        ME5 = ME5[103:153, 103:153]
        ME10 = ME10[103:153, 103:153]
        ME15 = ME15[103:153, 103:153]
        ME18 = ME18[103:153, 103:153]
        ME21 = ME21[103:153, 103:153]

        # AIDA images
        AIDA5 = AIDA5[103:153, 103:153]
        AIDA10 = AIDA10[103:153, 103:153]
        AIDA15 = AIDA15[103:153, 103:153]
        AIDA18 = AIDA18[103:153, 103:153]
        AIDA21 = AIDA21[103:153, 103:153]

        # Kraken MFBD images
        kraken5 = kraken5[103:153, 103:153]
        kraken10 = kraken10[103:153, 103:153]
        kraken15 = kraken15[103:153, 103:153]
        kraken18 = kraken18[103:153, 103:153]
        kraken21 = kraken21[103:153, 103:153]

        # set the text areas
        xlim = 7
        ylim = 43
        size = 'Small'

        # set the x/y position labels
        positions = [0, 15, 25, 35, 49]
        labels = ['-25', '-10', '0', '10', '25']

    # iterate through each cmap -> plotting each individually
    for j in range(0, int(len(cmap))):
        # plot 7x5 grid image
        fig, ax = plt.subplots(nrows=6, ncols=5, figsize=(25, 25), tight_layout=True)

        # set the images max scale based on the brightest image of the set per wavelength
        max5 = WH5.max()
        max10 = WH10.max()
        max15 = WH21.max()
        max18 = WH21.max()
        max21 = WH21.max()

        # row 1: MIRISim observations
        # 5.6 um
        im10 = ax[0][0].imshow(obs5, origin='lower', norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=max5),cmap=str(cmap[j]))
        ax[0][0].set_ylabel('MIRISim', fontsize=25)
        ax[0][0].set_title('F560W', fontsize=25)
        ax[0][0].xaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[0][0].xaxis.set_major_formatter(ticker.FixedFormatter(labels))
        ax[0][0].yaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[0][0].yaxis.set_major_formatter(ticker.FixedFormatter(labels))
        add_colorbar(im10, label=label, format=FormatStrFormatter('%.1e'))

        # 10 um
        im11 = ax[0][1].imshow(obs10, origin='lower',norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=max10),cmap=str(cmap[j]))
        ax[0][1].set_title('F1000W', fontsize=25)
        ax[0][1].xaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[0][1].xaxis.set_major_formatter(ticker.FixedFormatter(labels))
        ax[0][1].yaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[0][1].yaxis.set_major_formatter(ticker.FixedFormatter(labels))
        add_colorbar(im11, label=label, format=FormatStrFormatter('%.1e'))

        # 15 um
        im12 = ax[0][2].imshow(obs15, origin='lower',norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=max15),cmap=str(cmap[j]))
        ax[0][2].set_title('F1500W', fontsize=25)
        ax[0][2].xaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[0][2].xaxis.set_major_formatter(ticker.FixedFormatter(labels))
        ax[0][2].yaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[0][2].yaxis.set_major_formatter(ticker.FixedFormatter(labels))
        add_colorbar(im12, label=label, format=FormatStrFormatter('%.1e'))

        # 18 um
        im13 = ax[0][3].imshow(obs18, origin='lower',norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=max18),cmap=str(cmap[j]))
        ax[0][3].set_title('F1800W', fontsize=25)
        ax[0][3].xaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[0][3].xaxis.set_major_formatter(ticker.FixedFormatter(labels))
        ax[0][3].yaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[0][3].yaxis.set_major_formatter(ticker.FixedFormatter(labels))
        add_colorbar(im13, label=label, format=FormatStrFormatter('%.1e'))

        # 21 um
        im14 = ax[0][4].imshow(obs21, origin='lower',norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=max21),cmap=str(cmap[j]))
        ax[0][4].set_title('F2100W', fontsize=25)
        ax[0][4].xaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[0][4].xaxis.set_major_formatter(ticker.FixedFormatter(labels))
        ax[0][4].yaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[0][4].yaxis.set_major_formatter(ticker.FixedFormatter(labels))
        add_colorbar(im14, label=label, format=FormatStrFormatter('%.1e'))

        # row 2: Richardson-Lucy deconvolved images
        # 5.6 um
        im20 = ax[1][0].imshow(RL5, origin='lower', norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=max5),cmap=str(cmap[j]))
        ax[1][0].set_ylabel('Richardson-Lucy', fontsize=25)
        ax[1][0].text(xlim, ylim, 'Iterations: ' + str(RL_iter[0]), color='w', fontsize=20)
        ax[1][0].xaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[1][0].xaxis.set_major_formatter(ticker.FixedFormatter(labels))
        ax[1][0].yaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[1][0].yaxis.set_major_formatter(ticker.FixedFormatter(labels))
        add_colorbar(im20, label=label, format=FormatStrFormatter('%.1e'))

        # 10 um
        im21 = ax[1][1].imshow(RL10, origin='lower', norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=max10),cmap=str(cmap[j]))
        ax[1][1].text(xlim, ylim, 'Iterations: ' + str(RL_iter[1]), color='w', fontsize=20)
        ax[1][1].xaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[1][1].xaxis.set_major_formatter(ticker.FixedFormatter(labels))
        ax[1][1].yaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[1][1].yaxis.set_major_formatter(ticker.FixedFormatter(labels))
        add_colorbar(im21, label=label, format=FormatStrFormatter('%.1e'))

        # 15 um
        im22 = ax[1][2].imshow(RL15, origin='lower', norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=max15),cmap=str(cmap[j]))
        ax[1][2].text(xlim, ylim, 'Iterations: ' + str(RL_iter[2]), color='w', fontsize=20)
        ax[1][2].xaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[1][2].xaxis.set_major_formatter(ticker.FixedFormatter(labels))
        ax[1][2].yaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[1][2].yaxis.set_major_formatter(ticker.FixedFormatter(labels))
        add_colorbar(im22, label=label, format=FormatStrFormatter('%.1e'))

        # 18 um
        im23 = ax[1][3].imshow(RL18, origin='lower', norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=max18),cmap=str(cmap[j]))
        ax[1][3].text(xlim, ylim, 'Iterations: ' + str(RL_iter[3]), color='w', fontsize=20)
        ax[1][3].xaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[1][3].xaxis.set_major_formatter(ticker.FixedFormatter(labels))
        ax[1][3].yaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[1][3].yaxis.set_major_formatter(ticker.FixedFormatter(labels))
        add_colorbar(im23, label=label, format=FormatStrFormatter('%.1e'))

        # 21 um
        im24 = ax[1][4].imshow(RL21, origin='lower', norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=max21),cmap=str(cmap[j]))
        ax[1][4].text(xlim, ylim, 'Iterations: ' + str(RL_iter[4]), color='w', fontsize=20)
        ax[1][4].xaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[1][4].xaxis.set_major_formatter(ticker.FixedFormatter(labels))
        ax[1][4].yaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[1][4].yaxis.set_major_formatter(ticker.FixedFormatter(labels))
        add_colorbar(im24, label=label, format=FormatStrFormatter('%.1e'))

        # # row 3: Wiener-Hunt deconvolved images
        # # 5.6 um
        im30 = ax[2][0].imshow(WH5, origin='lower', norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=max5),cmap=str(cmap[j]))
        ax[2][0].set_ylabel('IWFT', fontsize=25)
        ax[2][0].text(xlim, ylim, 'Iterations: ' + str(WH_iter[0]), color='w', fontsize=20)
        ax[2][0].xaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[2][0].xaxis.set_major_formatter(ticker.FixedFormatter(labels))
        ax[2][0].yaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[2][0].yaxis.set_major_formatter(ticker.FixedFormatter(labels))
        add_colorbar(im30, label=label, format=FormatStrFormatter('%.1e'))

        # # # # 10 um
        im31 = ax[2][1].imshow(WH10, origin='lower', norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=max10),cmap=str(cmap[j]))
        ax[2][1].text(xlim, ylim, 'Iterations: ' + str(WH_iter[1]), color='w', fontsize=20)
        ax[2][1].xaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[2][1].xaxis.set_major_formatter(ticker.FixedFormatter(labels))
        ax[2][1].yaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[2][1].yaxis.set_major_formatter(ticker.FixedFormatter(labels))
        add_colorbar(im31, label=label, format=FormatStrFormatter('%.1e'))
        #
        # # # # # 15 um
        im32 = ax[2][2].imshow(WH15, origin='lower', norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=max15), cmap=str(cmap[j]))
        ax[2][2].text(xlim, ylim, 'Iterations: ' + str(WH_iter[2]), color='w', fontsize=20)
        ax[2][2].xaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[2][2].xaxis.set_major_formatter(ticker.FixedFormatter(labels))
        ax[2][2].yaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[2][2].yaxis.set_major_formatter(ticker.FixedFormatter(labels))
        add_colorbar(im32, label=label, format=FormatStrFormatter('%.1e'))
        #
        # # # # # 18 um
        im33 = ax[2][3].imshow(WH18, origin='lower', norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=max18),cmap=str(cmap[j]))
        ax[2][3].text(xlim, ylim, 'Iterations: ' + str(WH_iter[3]), color='w', fontsize=20)
        ax[2][3].xaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[2][3].xaxis.set_major_formatter(ticker.FixedFormatter(labels))
        ax[2][3].yaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[2][3].yaxis.set_major_formatter(ticker.FixedFormatter(labels))
        add_colorbar(im33, label=label, format=FormatStrFormatter('%.1e'))
        #
        # # # # 21 um
        im34 = ax[2][4].imshow(WH21, origin='lower', norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=max21),cmap=str(cmap[j]))
        ax[2][4].text(xlim, ylim, 'Iterations: ' + str(WH_iter[4]), color='w', fontsize=20)
        ax[2][4].xaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[2][4].xaxis.set_major_formatter(ticker.FixedFormatter(labels))
        ax[2][4].yaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[2][4].yaxis.set_major_formatter(ticker.FixedFormatter(labels))
        add_colorbar(im34, label=label, format=FormatStrFormatter('%.1e'))

        # # row 4: MEM deconvolved images
        im40 = ax[3][0].imshow(ME5, origin='lower', norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=max5),cmap=str(cmap[j]))
        ax[3][0].set_ylabel('Sparse', fontsize=25)
        ax[3][0].text(xlim, ylim, 'Iterations: ' + str(ME_iter[0]), color='w', fontsize=20)
        ax[3][0].xaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[3][0].xaxis.set_major_formatter(ticker.FixedFormatter(labels))
        ax[3][0].yaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[3][0].yaxis.set_major_formatter(ticker.FixedFormatter(labels))
        add_colorbar(im40, label=label, format=FormatStrFormatter('%.1e'))

        # 10 um
        im41 = ax[3][1].imshow(ME10, origin='lower', norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=max10),cmap=str(cmap[j]))
        ax[3][1].text(xlim, ylim, 'Iterations: ' + str(ME_iter[1]), color='w', fontsize=20)
        ax[3][1].xaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[3][1].xaxis.set_major_formatter(ticker.FixedFormatter(labels))
        ax[3][1].yaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[3][1].yaxis.set_major_formatter(ticker.FixedFormatter(labels))
        add_colorbar(im41, label=label, format=FormatStrFormatter('%.1e'))

        # 15 um
        im42 = ax[3][2].imshow(ME15, origin='lower', norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=max15),cmap=str(cmap[j]))
        ax[3][2].text(xlim, ylim, 'Iterations: ' + str(ME_iter[2]), color='w', fontsize=20)
        ax[3][2].xaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[3][2].xaxis.set_major_formatter(ticker.FixedFormatter(labels))
        ax[3][2].yaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[3][2].yaxis.set_major_formatter(ticker.FixedFormatter(labels))
        add_colorbar(im42, label=label, format=FormatStrFormatter('%.1e'))

        # 18 um
        im43 = ax[3][3].imshow(ME18, origin='lower', norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=max18),cmap=str(cmap[j]))
        ax[3][3].text(xlim, ylim, 'Iterations: ' + str(ME_iter[3]), color='w', fontsize=20)
        ax[3][3].xaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[3][3].xaxis.set_major_formatter(ticker.FixedFormatter(labels))
        ax[3][3].yaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[3][3].yaxis.set_major_formatter(ticker.FixedFormatter(labels))
        add_colorbar(im43, label=label, format=FormatStrFormatter('%.1e'))

        # 21 um
        im44 = ax[3][4].imshow(ME21, origin='lower', norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=max21),cmap=str(cmap[j]))
        ax[3][4].text(xlim, ylim, 'Iterations: ' + str(ME_iter[4]), color='w', fontsize=20)
        ax[3][4].xaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[3][4].xaxis.set_major_formatter(ticker.FixedFormatter(labels))
        ax[3][4].yaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[3][4].yaxis.set_major_formatter(ticker.FixedFormatter(labels))
        add_colorbar(im44, label=label, format=FormatStrFormatter('%.1e'))

        # # row 5: AIDA deconvolved images
        im50 = ax[4][0].imshow(AIDA5, origin='lower', norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=max5),cmap=str(cmap[j]))
        ax[4][0].set_ylabel('AIDA', fontsize=25)
        ax[4][0].text(xlim, ylim, 'Iterations: ' + str(AIDA_iter[0]), color='w', fontsize=20)
        ax[4][0].xaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[4][0].xaxis.set_major_formatter(ticker.FixedFormatter(labels))
        ax[4][0].yaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[4][0].yaxis.set_major_formatter(ticker.FixedFormatter(labels))
        add_colorbar(im50, label=label, format=FormatStrFormatter('%.1e'))

        # 10 um
        im51 = ax[4][1].imshow(AIDA10, origin='lower', norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=max10),cmap=str(cmap[j]))
        ax[4][1].text(xlim, ylim, 'Iterations: ' + str(AIDA_iter[1]), color='w', fontsize=20)
        ax[4][1].xaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[4][1].xaxis.set_major_formatter(ticker.FixedFormatter(labels))
        ax[4][1].yaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[4][1].yaxis.set_major_formatter(ticker.FixedFormatter(labels))
        add_colorbar(im51, label=label, format=FormatStrFormatter('%.1e'))

        # 15 um
        im52 = ax[4][2].imshow(AIDA15, origin='lower', norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=max15),cmap=str(cmap[j]))
        ax[4][2].text(xlim, ylim, 'Iterations: ' + str(AIDA_iter[2]), color='w', fontsize=20)
        ax[4][2].xaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[4][2].xaxis.set_major_formatter(ticker.FixedFormatter(labels))
        ax[4][2].yaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[4][2].yaxis.set_major_formatter(ticker.FixedFormatter(labels))
        add_colorbar(im52, label=label, format=FormatStrFormatter('%.1e'))

        # 18 um
        im53 = ax[4][3].imshow(AIDA18, origin='lower', norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=max18),cmap=str(cmap[j]))
        ax[4][3].text(xlim, ylim, 'Iterations: ' + str(AIDA_iter[3]), color='w', fontsize=20)
        ax[4][3].xaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[4][3].xaxis.set_major_formatter(ticker.FixedFormatter(labels))
        ax[4][3].yaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[4][3].yaxis.set_major_formatter(ticker.FixedFormatter(labels))
        add_colorbar(im53, label=label, format=FormatStrFormatter('%.1e'))

        # 21 um
        im54 = ax[4][4].imshow(AIDA21, origin='lower', norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=max21),cmap=str(cmap[j]))
        ax[4][4].text(xlim, ylim, 'Iterations: ' + str(AIDA_iter[4]), color='w', fontsize=20)
        ax[4][4].xaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[4][4].xaxis.set_major_formatter(ticker.FixedFormatter(labels))
        ax[4][4].yaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[4][4].yaxis.set_major_formatter(ticker.FixedFormatter(labels))
        add_colorbar(im54, label=label, format=FormatStrFormatter('%.1e'))

        # # # # # row 6: Kraken-MFBD deconvolved images
        # # # # # 5.6 um
        im60 = ax[5][0].imshow(kraken5, origin='lower', norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=max5), cmap=str(cmap[j]))
        ax[5][0].set_ylabel('Kraken', fontsize=25)
        ax[5][0].text(xlim, ylim, 'Iterations: ' + str(kraken_iter[0]), color='w', fontsize=20)
        ax[5][0].xaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[5][0].xaxis.set_major_formatter(ticker.FixedFormatter(labels))
        ax[5][0].yaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[5][0].yaxis.set_major_formatter(ticker.FixedFormatter(labels))
        add_colorbar(im60, label=label, format=FormatStrFormatter('%.1e'))
        # 10 um
        im61 = ax[5][1].imshow(kraken10, origin='lower', norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=max10),cmap=str(cmap[j]))
        ax[5][1].text(xlim, ylim, 'Iterations: ' + str(kraken_iter[1]), color='w', fontsize=20)
        ax[5][1].xaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[5][1].xaxis.set_major_formatter(ticker.FixedFormatter(labels))
        ax[5][1].yaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[5][1].yaxis.set_major_formatter(ticker.FixedFormatter(labels))
        add_colorbar(im61, label=label, format=FormatStrFormatter('%.1e'))
        # 15 um
        im62 = ax[5][2].imshow(kraken15, origin='lower', norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=max15),cmap=str(cmap[j]))
        ax[5][2].text(xlim, ylim, 'Iterations: ' + str(kraken_iter[2]), color='w', fontsize=20)
        ax[5][2].xaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[5][2].xaxis.set_major_formatter(ticker.FixedFormatter(labels))
        ax[5][2].yaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[5][2].yaxis.set_major_formatter(ticker.FixedFormatter(labels))
        add_colorbar(im62, label=label, format=FormatStrFormatter('%.1e'))
        # 18 um
        im63 = ax[5][3].imshow(kraken18, origin='lower', norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=max18),cmap=str(cmap[j]))
        ax[5][3].text(xlim, ylim, 'Iterations: ' + str(kraken_iter[3]), color='w', fontsize=20)
        ax[5][3].xaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[5][3].xaxis.set_major_formatter(ticker.FixedFormatter(labels))
        ax[5][3].yaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[5][3].yaxis.set_major_formatter(ticker.FixedFormatter(labels))
        add_colorbar(im63, label=label, format=FormatStrFormatter('%.1e'))
        # 21 um
        im64 = ax[5][4].imshow(kraken21, origin='lower', norm=ImageNormalize(stretch=LogStretch(), vmin=0, vmax=max21),cmap=str(cmap[j]))
        ax[5][4].text(xlim, ylim, 'Iterations: ' + str(kraken_iter[4]), color='w', fontsize=20)
        ax[5][4].xaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[5][4].xaxis.set_major_formatter(ticker.FixedFormatter(labels))
        ax[5][4].yaxis.set_major_locator(ticker.FixedLocator(positions))
        ax[5][4].yaxis.set_major_formatter(ticker.FixedFormatter(labels))
        add_colorbar(im64, label=label, format=FormatStrFormatter('%.1e'))

        fig.suptitle(title[idx], fontsize=40, y=0.99)
        fig.supxlabel('X (pixels)', fontsize=35, y=0.005, x = 0.51)
        fig.supylabel('Y (pixels)', fontsize=35, x=0.01)
        plt.savefig(save_img + '_SCALED_' + size + '_' + str(cmap[j]) + '.png')
        plt.show()

    # close all image importers
    obs_list5.close()
    obs_list10.close()
    obs_list15.close()
    obs_list18.close()
    obs_list21.close()

# if called: plot the single merit function measurements
if merit_single:

    # set the output path for the merit function measurements
    mes_path = 'Images/JWST/4_Deconvolution/'+str(Dir)+'/Measurements/' + str(model_type) + '/' + str(model[idx])

    # set the output path
    out_path = 'Images/JWST/4_Deconvolution/'+str(Dir)+'/Images/' + str(model_type) + '/' + str(model[idx]) + '/Total_Comparisons'
    # set the save file name
    save_merit = out_path + '/SCALED_Merit_Function_' + str(model[idx]) + '_' + str(sample)

    # read in the merit function data for each waveband
    nint5, flux5, fwhm5 = np.loadtxt(mes_path + '/F560W/RL_Image_Plot_' + str(model[idx]) + '_'+str(sample)+'_SCALED.txt', unpack=True)
    nint10, flux10, fwhm10 = np.loadtxt(mes_path + '/F1000W/RL_Image_Plot_' + str(model[idx]) + '_'+str(sample)+'_SCALED.txt', unpack=True)
    nint15, flux15, fwhm15 = np.loadtxt(mes_path + '/F1500W/RL_Image_Plot_' + str(model[idx]) + '_'+str(sample)+'_SCALED.txt', unpack=True)
    nint18, flux18, fwhm18 = np.loadtxt(mes_path + '/F1800W/RL_Image_Plot_' + str(model[idx]) + '_'+str(sample)+'_SCALED.txt', unpack=True)
    nint21, flux21, fwhm21 = np.loadtxt(mes_path + '/F2100W/RL_Image_Plot_' + str(model[idx]) + '_'+str(sample)+'_SCALED.txt', unpack=True)

    # sent the merit function values
    if Dir == '1_Richardson_Lucy':
        merit_iter = [RL_opt5[idx], RL_opt10[idx], RL_opt15[idx], RL_opt18[idx], RL_opt21[idx]]
    elif Dir == '2_Wiener_Hunt':
        merit_iter = [WH_opt5[idx], WH_opt10[idx], WH_opt15[idx], WH_opt18[idx], WH_opt21[idx]]
    elif Dir == '3_Max_Entropy':
        merit_iter = [ME_opt5[idx], ME_opt10[idx], ME_opt15[idx], ME_opt18[idx], ME_opt21[idx]]
    elif Dir == '4_AIDA':
        merit_iter = [AIDA_opt5[idx], AIDA_opt10[idx], AIDA_opt15[idx], AIDA_opt18[idx], AIDA_opt21[idx]]
    elif Dir == '5_Kraken':
        merit_iter = [kraken_opt5[idx], kraken_opt10[idx], kraken_opt15[idx], kraken_opt18[idx], kraken_opt21[idx]]
    elif Dir == '7_sparse_condat_vu':
        merit_iter = [ME_opt5[idx], ME_opt10[idx], ME_opt15[idx], ME_opt18[idx], ME_opt21[idx]]

    nint = []
    flux_5 = []
    flux_10 = []
    flux_15 = []
    flux_18 = []
    flux_21 = []
    fwhm_5 = []
    fwhm_10 = []
    fwhm_15 = []
    fwhm_18 = []
    fwhm_21 = []
    plot_flux_5 = []
    plot_flux_10 = []
    plot_flux_15 = []
    plot_flux_18 = []
    plot_flux_21 = []
    plot_fwhm_5 = []
    plot_fwhm_10 = []
    plot_fwhm_15 = []
    plot_fwhm_18 = []
    plot_fwhm_21 = []
    for i in range(0, single_xlim):
        # convert nint from float -> int
        nint1 = nint5[i]
        nint1 = int(nint1)

        # read in flux data
        fluxA = flux5[i]
        fluxB = flux10[i]
        fluxC = flux15[i]
        fluxD = flux18[i]
        fluxE = flux21[i]

        # read in FWHM data
        fwhmA = fwhm5[i]
        fwhmB = fwhm10[i]
        fwhmC = fwhm15[i]
        fwhmD = fwhm18[i]
        fwhmE = fwhm21[i]

        # trim float value for plotting merit functions
        Aflux = '{:.3f}'.format(fluxA)
        Bflux = '{:.3f}'.format(fluxB)
        Cflux = '{:.3f}'.format(fluxC)
        Dflux = '{:.3f}'.format(fluxD)
        Eflux = '{:.3f}'.format(fluxE)
        Afwhm = '{:.3f}'.format(fwhmA)
        Bfwhm = '{:.3f}'.format(fwhmB)
        Cfwhm = '{:.3f}'.format(fwhmC)
        Dfwhm = '{:.3f}'.format(fwhmD)
        Efwhm = '{:.3f}'.format(fwhmE)

        # trim float value for plotting images
        Aflux_plot = '{:.2f}'.format(fluxA)
        Bflux_plot = '{:.2f}'.format(fluxB)
        Cflux_plot = '{:.2f}'.format(fluxC)
        Dflux_plot = '{:.2f}'.format(fluxD)
        Eflux_plot = '{:.2f}'.format(fluxE)
        Afwhm_plot = '{:.2f}'.format(fwhmA)
        Bfwhm_plot = '{:.2f}'.format(fwhmB)
        Cfwhm_plot = '{:.2f}'.format(fwhmC)
        Dfwhm_plot = '{:.2f}'.format(fwhmD)
        Efwhm_plot = '{:.2f}'.format(fwhmE)

        # append everything
        nint.append(nint1)
        flux_5.append(fluxA)
        flux_10.append(fluxB)
        flux_15.append(fluxC)
        flux_18.append(fluxD)
        flux_21.append(fluxE)
        fwhm_5.append(fwhmA)
        fwhm_10.append(fwhmB)
        fwhm_15.append(fwhmC)
        fwhm_18.append(fwhmD)
        fwhm_21.append(fwhmE)
        plot_flux_5.append(Aflux_plot)
        plot_flux_10.append(Bflux_plot)
        plot_flux_15.append(Cflux_plot)
        plot_flux_18.append(Dflux_plot)
        plot_flux_21.append(Eflux_plot)
        plot_fwhm_5.append(Afwhm_plot)
        plot_fwhm_10.append(Bfwhm_plot)
        plot_fwhm_15.append(Cfwhm_plot)
        plot_fwhm_18.append(Dfwhm_plot)
        plot_fwhm_21.append(Efwhm_plot)

    # plot the data
    # Create the 1x2 grid
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(12, 6), tight_layout=True)

    # Plot the FWHM values
    # plt.subplot(121)
    # 5 um
    ax[0].plot(nint, fwhm_5, color='b', linewidth=1)
    ax[0].plot(merit_iter[0], fwhm_5[merit_iter[0]], color='b', marker='X', markersize=15)
    # 10 um
    ax[0].plot(nint, fwhm_10, color='g', linewidth=1)
    ax[0].plot(merit_iter[1], fwhm_10[merit_iter[1]], color='g', marker='X', markersize=15)
    # 15 um
    ax[0].plot(nint, fwhm_15, color='y', linewidth=1)
    ax[0].plot(merit_iter[2], fwhm_15[merit_iter[2]], color='y', marker='X', markersize=15)
    # 18 um
    ax[0].plot(nint, fwhm_18, color='darkorange', linewidth=1)
    ax[0].plot(merit_iter[3], fwhm_18[merit_iter[3]], color='darkorange', marker='X', markersize=15)
    # 21 um
    ax[0].plot(nint, fwhm_21, color='r', linewidth=1)
    ax[0].plot(merit_iter[4], fwhm_21[merit_iter[4]], color='r', marker='X', markersize=15)

    ax[0].set_title('FWHM', fontsize=13)
    ax[0].set_xlim(0, single_xlim)
    ax[0].set_ylim(0, int(max(fwhm_21))+2)
    ax[0].set_ylabel('FWHM (pixel)')
    ax[0].text(5,int(max(fwhm_21))+1, 'Max. Entropy', fontsize = 20)
    ax[0].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

    if diff_lim:
        # set the calculated diffraction limits
        diff_lim = [1.953, 3.487, 5.232, 6.278, 7.324]
        ax[0].hlines(diff_lim[0], 0, int(len(fwhm_5)), color = 'b', linestyle = '--')
        ax[0].text(int(single_xlim-single_buffer), diff_lim[0] + 0.1, 'F560W Diff. Limit', color='b')
        ax[0].hlines(diff_lim[1], 0, int(len(fwhm_5)), color='g', linestyle = '--')
        ax[0].text(int(single_xlim-single_buffer), diff_lim[1] + 0.1, 'F1000W Diff. Limit', color='g')
        ax[0].hlines(diff_lim[2], 0, int(len(fwhm_5)), color = 'y', linestyle = '--')
        ax[0].text(int(single_xlim-single_buffer), diff_lim[2] + 0.1, 'F1500W Diff. Limit', color='y')
        ax[0].hlines(diff_lim[3], 0, int(len(fwhm_5)), color='darkorange', linestyle = '--')
        ax[0].text(int(single_xlim-single_buffer), diff_lim[3] + 0.1, 'F1800W Diff. Limit', color='darkorange')
        ax[0].hlines(diff_lim[4], 0, int(len(fwhm_5)), color = 'r', linestyle = '--')
        ax[0].text(int(single_xlim-single_buffer), diff_lim[4] + 0.1, 'F2100W Diff. Limit', color='r')

        save_merit = save_merit + '_diff_lim'

    # Plot the flux values
    # plt.subplot(122)
    # 5 um
    ax[1].plot(nint, flux_5, color='b',  linewidth=1,  label='F560W (' + str(merit_iter[0]) + ' iterations)')
    #ax[1].plot(merit_iter[0], flux_5[merit_iter[0]], color='b', marker='X', markersize=15)
    # 10 um
    ax[1].plot(nint, flux_10, color='g',  linewidth=1,  label='F1000W (' + str(merit_iter[1]) + ' iterations)')
    #ax[1].plot(merit_iter[1], flux_10[merit_iter[1]], color='g', marker='X', markersize=15)
    # 15 um
    ax[1].plot(nint, flux_15, color='y',  linewidth=1,  label='F1500W (' + str(merit_iter[2]) + ' iterations)')
    #ax[1].plot(merit_iter[2], flux_15[merit_iter[2]], color='y', marker='X', markersize=15)
    # 18 um
    ax[1].plot(nint, flux_18, color='darkorange',  linewidth=1,  label='F1800W (' + str(merit_iter[3]) + ' iterations)')
    #ax[1].plot(merit_iter[3], flux_18[merit_iter[3]], color='darkorange', marker='X', markersize=15)
    # 21 um
    ax[1].plot(nint, flux_21, color='r',  linewidth=1,  label='F2100W (' + str(merit_iter[4]) + ' iterations)')
    #ax[1].plot(merit_iter[4], flux_21[merit_iter[4]], color='r', marker='X', markersize=15)

    ax[1].set_title('Flux', fontsize=13)
    ax[1].set_xlim(0, single_xlim)
    ax[1].set_ylabel('Flux Density (counts/pixel)')
    ax[1].yaxis.set_major_formatter(FormatStrFormatter('%.2e'))
    ax[1].legend(loc='lower right')

    fig.supxlabel('Iterations')
    plt.suptitle(title[idx], fontsize=18)
    plt.savefig(save_merit + '.png')
    plt.show()