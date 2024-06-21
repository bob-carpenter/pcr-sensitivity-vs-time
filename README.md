# SARS-Cov-2 PCR test sensitivity versus time

This repository contains the reproducible code for the analysis of Covid PCR test sensitivity versus time and the source for the poster presented at ISBA 2024.  The data included is from patients hospitalized with Covid in the UK.

This is joint work with Thomas Ward (UK Health Security Agency).

## Generating the poster

Change directories to `pcr-sensitivity-vs-time/latex` and run `pdflatex` on `isba2024.tex`.  There is no cross-referencing or bibliography so running once is sufficient.

## Generating the analyses

Change directories to `pcr-sensitivity-vs-time/R` and run the script `fit.R`.  This will generate pdf files in the same directory and will not overwrite the image data in the LaTeX directory.

## Data Source

The data is from Mallet et al., 2020, [At what times during infection is SARS-CoV-2 detectable and no longer detectable using RT-PCR-based tests? A systematic review of individual participant data](https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01810-8). _BMC Medicine_.



