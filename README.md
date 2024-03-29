# Code for analysis of the Dipole Estimator used in Secrest et. al. (2021) and Secrest et. al. (2022)

## Introduction
This repository contains the Julia and Python code used in analysis of the fit_dipole routine of the Healpy Python library. This is an ongoing project under the guidance of Prof. Sarkar (Oxford), Dr. Rameez (TIFR) and Prof. Singh (IISER Mohali). Each branch contains different files with the highest numbered branch containing the latest (and most optimised) iteration of the process. Relevant information about each branch is as follows

### Master Branch
This is the primary branch of the project. The Simulation.jl file contains different Julia functions utilised for processing in V1 and V2 branches. The Simulation_v2.jl file contains different Julia functions utilised for processing in V3. For scrutinizing the code, this branch is the one that should be utilised.

### V1
This branch is the earliest version of the code. Here, the datasets were generated using Simulation notebook. Then, the fit_dipole routine was utilised for fitting the dipole to the simulated data in the FitDipole Notebook. The results were then used to calculate the bias of the estimator for different catalogue sizes in Simulation2, creating a scatter plot of dipole magnitude vs offset from CMB dipole in Simulation 3 and finally, the dipoles were converted to galactic coordinates in Equatorial2Galactic and then visualised in Dipole Visualisation notebook. For using this version, it is suggested that 3 different folders 'pixel';'galacticpixel' and 'supergalacticpixel' be created in the directory from where the code is being processed.

### V2
This branch in the second version of the code. Here, the fit_dipole function of Healpy has been adapted in Julia. Then, the dataset generation, dipole calculation, bias estimation and scatter plot generation were all implemented using Process.jl. Finally, dipole visulization in galactic coordinates was done in Visualization.py.

### V3
This branch in the third version of the code. Here, the pixelised boosting function has been adapted in Julia. Then, the dataset generation, dipole calculation, bias estimation and scatter plot generation were all implemented using Process_v2.jl. Finally, dipole visulization and histogram of dipole amplitude vs offset were done in Visualization_v2.py.