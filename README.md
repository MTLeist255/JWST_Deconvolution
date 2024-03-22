# JWST_Deconvolution
Applications of different deconvolution algorithms to a simulated toy AGN model and JWST/MIRI Imaging. Source code (provided as is) for Leist et al. 2024. The Jupyter notebook gives an example of the deconvolution methodology followed. In it, the user can:

	1) Read in the original and PSF FITS image files
	2) Non-circulant Richardson-Lucy (Richardson 1972; Lucy 1974; XXX et al. 2014) deconvolve the image data
	3) Measure merit functions and determine convergence criteria

This process is repeatable for the other deconvolution algorithms used in Leist et al. 2024. We encourage any interested party to read the algorithms seminal papers, download the relevant source code, and try for themselves.

**Kraken:**: Contact author for permission to use code (give Doug's email) <br>
        [Hope, D. A., Jefferies, S. M., Causi, G. L., et al. 2022, ApJ, 926, 88](https://iopscience.iop.org/article/10.3847/1538-4357/ac2df3)
        
**non-circulant Richardson-Lucy:** Code developed by Brian Northan ([source code](https://github.com/clij/clij2-fft)) <br>
Richardson 1972 <br>
Lucy 1974 <br>
XXX et al. 2014? <br>
First introduced in the 2014 Grand Challenge on Deconvolution: [code link](https://bigwww.epfl.ch/deconvolution/challenge2013/index.html?p=doc_math_rl)

[Adaptive Imaging Deconvolution Algorithm (AIDA)](https://github.com/erikhom/aida) <br>
Hom et al. 2007

[Sparse regularization with the Condat-Vu (SCV) algorithm](https://github.com/CEA-COSMIC/pysap-astro) <br>
Farrens et al. 2017

Iterative Weiner Filtering and Thresholding (IWFT): Contact author for permission to use code (give XXX email) <br>
XXX et al. 2019

# Installation
To run the Jupyter notebook we recommend installing Conda then creating a Condo environment with Python 3.9:

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

	citation, link (Leist et al. 2024)
 
# Licensing
Licensed under MIT license.
