# Miller-et-al_ProteinPolarity
Accompanies Miller SG, Hoh M, Ebmeier CC, Tay, J, Ahn NG, Cooperative polarization of MCAM/CD146 and ERM family proteins in melanoma. Contains R and MATLAB code for the 3P (protein polarity by percentile) image analysis pipeline, kymograph analysis, and BioID proteomics analysis. Repository for Supplemental Data containing processed BioID proteomics results and imaging data comprised of manual annotation, machine learning classifications, and single-cell feature values from 3P pipeline analyses.

# License and Copyright Information
Copyright (C) 2023  University of Colorado Boulder.

Author contributions: Suzannah Miller wrote the code contained in this repository. The concept and code for applying a percentile threshold to identify bright regions within a cell were adapted from code developed by Jian Wei Tay.

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. A copy of the license is provided in the LICENSE file.

The corresponding author is Natalie Ahn, University of Colorado Boulder.

# 3P Pipeline Materials
The 3P pipeline is a method for single-cell analysis of fluorescence microscopy data which identifies cells with polarized proteins and characterizes the localization of polarized proteins. Image segmentation and feature extraction is performed in MATLAB. Additional features are extracted using CellProfiler (Stirling et al., BMC Bioinformatics 22, 433). Features are combined and pre-processed for ML in R. Images are sampled and blinded for manual labeling, and then training of a machine learning model for classifying polarized cells is implemented in R using the caret package (Kuhn, 2019). Additionally, the code is included for analysis of the colocalization features and generation of the figures in Miller _et al_. 

The implementation of the pipeline is described extensively in the Supplemental Materials for Miller _et al_.

The "3P-Pipeline" folder is divided into three sections. 

The "Feature-Selection" folder contains code, input, and output data used during optimization of the pipeline to select features for machine learning. This was used to generate Supplemental Figure S7A-C in Miller _et al_.

The "Final-pipeline" folder contains all code for implementation of the pipeline, including MATLAB code for image processing and feature extraction, the CellProfiler pipeline, and R Markdown files for machine learning using iterative builing of a training set. The input and output data and the classification models for the dataset shown in Miller _et al_. Figure 7 and 9 (classification of polarized MCAM, MSN, EZR, and p-ERM) are included in the folder to allow the user to run the R Markdown file "MachineLearning_MCAM-ERM.Rmd". This file also contains the statistical analyses of the results and analysis of the co-polarization and colocalization of pairs of polarized proteins. 

The "siRNA-ERM-Depletion" folder contains the code for analyzing the dataset shown in Figure 5B,D, Figure 6, and Supplemental Figures S7D-S9 in Miller _et al_. 

# Additional Content
The "BioID" folder contains the R Markdown file for performing the analysis of the MCAM BioID data.

The "Kymograph" folder contains the MATLAB scripts for generating kymographs from the live-cell imaging data of WM239a cell lines expressing MCAM-GFP, LifeAct-mTagBFP, and either MSN-mCherry, EZR-mCherry, or RDX-mCherry.

The "Supplemental-Data-Files" is a repository for supplemental data which accompanies Miller et al. It contains Supplemental Table S1 (BioID peptide quantification) and the consolidated feature data, ML classifications, and ML labeling for 3P pipeline analyses. Note that the "3P-Pipeline" folder includes additional input and output data for running the 3P analyses.

# Software Requirements
Analyses were run using MATLAB 2020a including the Image Processing Toolbox (version 11.1) and the Statistics and Machine Learning Toolbox (version 11.7), CellProfiler 4.2.1 or earlier (https://cellprofiler.org/), and R version R4.1.0 or later. Nikon .nd2 image files were imported into MATLAB using The Open Microscopy Environment Bio-Formats toolbox (https://docs.openmicroscopy.org/bio-formats/)https://docs.openmicroscopy.org/bio-formats/).


