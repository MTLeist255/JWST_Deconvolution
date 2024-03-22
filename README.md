# JWST_Deconvolution Overview
Applications of different deconvolution algorithms to a simulated toy AGN model and JWST/MIRI Imaging. Source code (provided as is) for Leist et al. 2024. The Jupyter notebook gives an example of the deconvolution methodology followed. In it, the user can:

	1) Read in the original and PSF FITS image files
	2) Non-circulant Richardson-Lucy deconvolve the image data
	3) Measure merit functions and determine convergence criteria

This process is repeatable for the other deconvolution algorithms used in Leist et al. 2024. We encourage any interested party to read the algorithms seminal papers, download the relevant source code, and try for themselves.

**Kraken:** Contact author for permission to use code (email: [douglas.hope@protonmail.com](douglas.hope@protonmail.com)) <br>
        [Hope, D. A., Jefferies, S. M., Causi, G. L., et al. 2022, ApJ, 926, 88](https://iopscience.iop.org/article/10.3847/1538-4357/ac2df3)
        
**non-circulant Richardson-Lucy:** Code developed by Brian Northan ([source code](https://github.com/clij/clij2-fft)) <br>
[Richardson, W. H. 1972, JOSA, 62, 55](https://opg.optica.org/josa/fulltext.cfm?uri=josa-62-1-55&id=54565) <br>
[Lucy, L. B. 1974, AJ, 79, 745](https://ui.adsabs.harvard.edu/abs/1974AJ.....79..745L/abstract) <br>
[Ingaramo, M. C., York, A. G., Hoogendoorn, E., et al. 2014, ChemPhysChem, 15, 794](https://chemistry-europe.onlinelibrary.wiley.com/doi/10.1002/cphc.201300831) <br>
First introduced in the [2014 Grand Challenge on Deconvolution](https://bigwww.epfl.ch/deconvolution/challenge2013/index.html?p=doc_math_rl)

**Adaptive Imaging Deconvolution Algorithm (AIDA)**: Code developed by Erik Hom ([source code](https://github.com/erikhom/aida)) <br>
[Hom, E. F. Y., Marchis, F., Lee, T. K., et al. 2007, JOSAA, 24, 1580](https://opg.optica.org/josaa/fulltext.cfm?uri=josaa-24-6-1580&id=134611)

**Sparse regularization with the Condat-Vu (SCV) algorithm:** Code developed by "a community of people, among which the CEA Neurospin UNATI and CEA CosmoStat laboratories, in particular A. Grigis, J.-L. Starck, P. Ciuciu, and S. Farrens" ([source code](https://github.com/CEA-COSMIC/pysap-astro)) <br>
[Farrens, S., Ngolè Mboula, F. M., & Starck, J. L. 2017, A&A, 601, A66](https://www.aanda.org/articles/aa/full_html/2017/05/aa29709-16/aa29709-16.html)

**Iterative Weiner Filtering and Thresholding (IWFT)**: Contact author for permission to use code (email: [sroubekf@utia.cas.cz](sroubekf@utia.cas.cz)) <br>
[Šroubek, F., Kerepecký, T., & Kamenický, J. 2019, in 2019 27th European Signal Processing Conf. (EUSIPCO) (New York: IEEE), 1](https://ieeexplore.ieee.org/document/8903114)

# Installation
To run the Jupyter notebook we recommend installing Conda then creating a Conda environment with Python 3.9:

Conda create -n JWST_deconvolution python=3.9 <br>
Conda activate JWST_deconvolution

Then pip installing the following dependancies:

Pip install Jupyter <br>
	    webbpsf <br>
	    astropy <br>
	    'photutils[all]' <br>
	    numpy <br>
	    matplotlib <br>
	    reproject <br>
	    poppy <br>
	    scipy <br>
     
# Citing
If you use this code for a scientific publication, please cite the following article:

[Leist, M. T., Packham, C., Rosario, D. J. V., et al. 2024, AJ, 167, 96](https://iopscience.iop.org/article/10.3847/1538-3881/ad1886)
 
# Licensing
Licensed under MIT license.
