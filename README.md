# ERAP2 shapes the HLA-A29 peptidome
Scripts and data used in analyses for “ERAP2 facilitates a cryptic binding motif of  HLA-A29”. Preprint url;

This repository contains:

# Data

# Data sets used for the analyses in the study

*Biological replicate 1*

1. MS of peptide elution from immunopurification by the human monoclonal antibody DK1G8 (anti-HLA-A29)
https://doi.org/10.5281/zenodo.3833697

2. MS of peptide elution from immunopurification by the pan-class I antibody W6/32 (anti-HLA-ABC)
https://doi.org/10.5281/zenodo.3833778

*Biological replicate 2*

3. MS of peptide elution from immunopurification by the human monoclonal antibody DK1G8 (anti-HLA-A29)
https://doi.org/10.5281/zenodo.3833784

4. MS of peptide elution from immunopurification by the pan-class I antibody W6/32 (anti-HLA-ABC)
https://doi.org/10.5281/zenodo.3833791

# R Scripts

**1.Peptide abundance.R**

An R script based on the workflow for differential expression analysis using limma R package following Kammers and associates. Kammers, K., Cole, R. N., Tiengwe, C. & Ruczinski, I. Detecting Significant Changes in Protein Abundance. EuPA open proteomics 7, 11–19 (2015). 
adapted from: http://www.biostat.jhsph.edu/~kkammers/software/eupa/R_guide.html


**2. NMDS.R**
An R script written by Sarkizova et al Nature Biotechnology volume 38, pages199–209(2020) and adapted for non-metric multidimensional scaling of peptidome data in this study. This script requires the 948 9-mers from HLA-A29 which is loaded in the script ("HLA_A29_948_9mers.RData)


