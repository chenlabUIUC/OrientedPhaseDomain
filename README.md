# OrientedPhaseDomain

This code package is developed to map the strain, oriented phase domain, and pair distribution function for the spinel MnO2 cathode particles with inserted ions.
Current version: 1.0
Date: 7/2022
Test: Matlab 2016b on Microsoft Windows 10 using an Intel(R) Core(TM) i5-8250U CPU @ 1.60GHz processor
For more information about the project, algorithms, and related publications please refer to the Chen Group website.

# Reference

# Getting started

(A) For mapping the strain and domain structures in cathode particles, visit folder "MapStrainAndDomain"

The matrixes of reciprocal lattice vectors are obtained by analyzing the 4D-STEM datasets using the “imToolBox” software package. For the instructions on how to analyze the raw 4D-STEM datasets and how to obtain the reciprocal lattice vectors, please refer to: https://github.com/flysteven/imToolBox

An example of matrixes of reciprocal lattice vectors are provided ((g220).xlsx and (g002).xlsx), as well as the corresponding mask of the cathode particle. The following procedures illustrate how to map strain and domain distribution based on the matrixes of reciprocal lattice vectors. The following steps should be executed only once:

1. Open "calculate_strain_g220_g002.m" and run it.
2. Open "CreatPremask.m" and run it. Import ""mask-Pre.tif" in ImageJ and manually create a mask by delienating the particle shape.
3. Open "MapStrainforXandY.m" and run it. Save strain maps.
4. Open "CalculateTetragonalityHistogram.m" and run it. Histogram will be generated.
5. Open "FitTetragonalityHist.m" and run it. Use the "cftool" for the fitting of the histogram. Set the constraints on the fitting parameters as needed to obtain a good fitting.
6. Open "MapDomain.m". Set the threshold from the fitting results to distinguish the [100]t and [111]t. An example of thresholds is provided in the code. Run it. Domain map will be genrated.
7. Open ""CalculateDomainStatistics.m". Run it. A summary of domain number, domain area, and particle area will be calculated.

Note: "cmapStrain.mat" and "cmapDomain.mat" files are for the color scales of the strain and domain maps, respectively.

(B) For the calculation of radial distribution function and pair correlation function

For the radial distribution function: Examples of the domain center of mass coordinates are provided as excel files. The coordinates of the domain are in the unit of pixel (pixel size 2 nm). Each excel file contains the domain center of mass coordinates for the tetragonal domains within one cathode particle. Open "DomainRadialDistribution_averageFiles.m" and run it. The code calculates the radial distribution function for the same type of domain and different types of domains considering all the cathode particles.

For the pair correlation function: One example of the domain center of mass coordinates for a cathode particleis provided as an excel file. The coordinates of the domain are in the unit of pixel (pixel size 2 nm). Open "DomainRadialDistribution_corrLen.m" and run it. It will calculate the pair correlation function for the cathode particle.

(C) For the Mn white line ratio analysis and EELS mapping

Open "Hyperspy EELS Mn White Line Ratio.ipynb" using JupyterNotebook. Instructions on how to use the code are included in the beginning of the code and in the comments. We used Pearson method for the calculations of Mn white line ratio (L3/L2). 
Run the codes as instructed in the notebook. An example of the low-loss and core-loss EELS data for a pristine cathode NP is provided for running the code.
The code for generating the map of Mn white line ratio (L3/L2) is provided, together with an example data and the corresponding mask of the cathode particle.
