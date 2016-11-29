# README

## Description
The program determines regular, i.e. constant within each exposure, velocity off-sets between quasar/stellar absorption line spectra.
The latter should be normalised and be output files of [UVES_popler](https://github.com/MTMurphy77/UVES_popler) software by Michael Murphy.

## Installation
The software was written to work with `python 2.7`. 
To get it to work, include the path to the `voffset` executable file into your system PATH.
That is it, you can run it from any directory on your machine.

**Python modules used:**

- Barak
- astropy
- numpy
- scipy
- matplotlib (if you use `--plot`)
- sys


## Usage

`voffset <list of spectra> --an anchor [--option]`

where:
```
<list of spectra> - a list of paths to the fits files with spectra to cross-correlate with the reference spectrum,

anchor - a path to the file with the reference spectrum (should follow `--an`).
```
Options:
```
--fwhm - folowed with a path to the file with Gaussian FWHM kernels (in km/s) for each of the spectra including the reference one 
(if not specified, 3-pixel width is used instead). Example of the file with Gaussian kernels can be found in the `example` folder.

--plot - create a pdf file with an off-set distribution bar chart.
```
The program creates a file with two columns: file names of the specified spectra, velocity off-sets in km/s.
The first line corresponds to the reference spectrum.

## How it works

I colloquially call the process of finding the off-sets as a cross-correlation. However, it is not precisely what voffset does. 
A simple cross-correlation procedure is limited to a pixel size 
(~1-3 km/s for typical spectra from high resolution spectrographs of HIRES/Keck and UVES/VLT), 
while typical off-sets have sub-pixel values. 
I thus implemented a method similar to one described by [Evans & Murphy (2013)](https://arxiv.org/pdf/1310.5703v1.pdf).

First, the spectra are convolved with a Gaussian of constant FWHM in velocity space using [Barak](https://github.com/nhmc/Barak) package by Neil Crighton. 
The latter can be specified using `--fwhm` option. 
The suggestion is to use resolution (FWHM) of the telescope. 
If no file is specified, the program will convolve the spectra with 3-pixel width instead. 
Then, flux and error arrays of the reference exposure (specified with `--an`) are interpolated with cubic spline to get continuous functions 
of wavelength and the off-set. 
The velocity off-set between each examined exposure and the reference one is given by the minimum 
of a chi^2-function with respect to the off-set. The program uses only valid (unaffected by CRs, etc.) pixels covered 
by both exposures. It determines the chi^2 minimums for all the specified spectra and writes them into a "voffset_result.dat" file. 
The program knows 2 ways of minimising: numerical and "direct" (you can choose either, see 179th line of the `voffset` file). 
The default one is the latter where chi^2 is directly calculated in a range of +/-5 km/s with a step of 0.001 km/s 
(you can tweak it in the code) and then the minimum value is taken out.
