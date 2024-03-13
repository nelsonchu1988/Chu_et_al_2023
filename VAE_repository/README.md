# Code Submission for Variational Autoencoder Analysis of Bone Marrow

This repository contains the code for the variational autoencoder (VAE) analysis of bone marrow, including datasets and reproduction of VAE-related figures in the manuscript (i.e., figures 3c-3j, 4k-4o, 4q, and 5k-5o, 5q).

## System Requirements

- MATLAB R2022a or later version for Microsoft Windows, macOS, and Linux operating systems.

The following toolboxes need to be installed in MATLAB:
- Statistics and Machine Learning Toolbox
- Deep Learning Toolbox 

## Installation

Please follow these steps:
1. Download this repository.
2. Unzip the downloaded file into your desired directory.
3. Launch MATLAB.
4. Navigate to this repository folder located in the directory where you unzipped the file.

## Demo

Follow these steps:
1. Launch MATLAB and navigate to this repository folder located in the directory where you unzipped the file.
2. Open the script file 'main.m' and run the script.

This will initiate the analysis to:
1. Load the dataset (fig3-CXCL12, fig4-CD146, or fig5-CD271; switch between them in line 12 of main.m).
2. Load its corresponding pretrained VAE model.
3. Generate the figures in its corresponding output folder (replicate_fig3_CXCL12, replicate_fig4_CD146, or replicate_fig4_CD271). 

The runtime is typically within 3 minutes.
