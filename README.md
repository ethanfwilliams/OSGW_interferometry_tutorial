# Supporting code for OSGW interferometry

Williams, E.F., et al. (submitted) "Surface gravity wave interferometry and ocean current monitoring with ocean-bottom DAS," <i>EarthArXiv</i>

This repository contains example scripts to compute OSGW cross-correlations, measure dispersion, and invert for current speed. The original code for the paper (to produce and process over 45000x215 cross-correlations) was written in a (messy but fast) combination of Python, Fortran, and CUDA Fortran. Here, I have re-written the key steps in simple Python in order to process a small demonstration set of cross-correlation pairs. 

