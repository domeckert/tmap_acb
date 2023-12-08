# tmap_acb
Thermodynamic map generation from X-ray multiband images

This is a set of Python scripts that allows the user to generate thermodynamic maps for galaxy clusters and groups from a set of X-ray images in various energy bands. The tool uses the adaptive circular binning (ACB) algorithm to accumulate a target number of photons around each pixel, determine the pseudo-spectrum of the source in the corresponding region and fit it with a spectral template generated using the APEC code implemented in XSPEC.

- _templates.py_ : Generate template count rates from the APEC model in a set of energy bands

- _skybkg__spec.py_ : Compute template count rates from a sky background model

- _tmap__acb.py_ : Thermodynamic map generation

Running each of the tools without arguments provides a short description of their usage
