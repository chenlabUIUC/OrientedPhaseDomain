# OrientedPhaseDomain

This code package is developed to map the strain, oriented phase domain, and pair distribution function for the spinel MnO2 cathode particles with inserted ions.
Current version: 1.0
Date: 7/2022
Test: Matlab 2016b on Microsoft Windows 10 using an Intel(R) Core(TM) i5-8250U CPU @ 1.60GHz processor
For more information about the project, algorithms, and related publications please refer to the Chen Group website.

# Reference

# Getting started

(A) For mapping the strain and domain structures in cathode particles, visit folder "MapStrainAndDomain"

The matrixes of reciprocal lattice vectors are obtained by analyzing the 4D-STEM datasets using the “imToolBox” software package. For the instructions on how to analyze the raw 4D-STEM datasets, please refer to: https://github.com/flysteven/imToolBox

An example of matrixes of reciprocal lattice vectors are provided ((g220).xlsx and (g002).xlsx), as well as the corresponding mask of the cathode particle. The following procedures illustrate how to map strain and domain distribution based on the matrixes of reciprocal lattice vectors. The following steps should be executed only once:

1. Open "calculate_strain_g220_g002.m" and run it.
2. Open "MapStrainforXandY.m" and run it. Strain maps will be generated. Save strain maps as needed.
3. Open "CalculateTetragonalityHistogram.m" and run it. Histogram of tetragonality will be generated and saved.
4. Open "FitTetragonalityHist.m" and run it. An example of fitting results will be shown. Alternatively, one can choose to do fitting of the histogram using the "cftool". In "cftool", set the constraints on the fitting parameters as needed to obtain a good fitting. 
5. Open "MapDomain.m". Set the threshold from the fitting results to distinguish the [100]t and [111]t. An example of thresholds is provided in the code. Run it. Domain map will be genrated. Save strain maps as needed.
6. Open ""CalculateDomainStatistics.m". Run it. A summary of domain number, domain area, and particle area will be calculated and saved.

Note: "cmapStrain.mat" and "cmapDomain.mat" files are for the color scales of the strain and domain maps, respectively. "createFit.m" contains the example of fitting.

(B) For the calculation of radial distribution function and pair correlation function, visit folders "RadialDistributionFunction" and "PairCorrelationFunction"

For the radial distribution function: Open "DomainRadialDistribution_averageFiles.m" and run it. The code calculates the radial distribution function for the same type of domain and different types of domains considering all the cathode particles. Examples of the domain center of mass coordinates of each cathode particle are provided as excel files. Each excel file contains the domain center of mass coordinates for one cathode particle. The coordinates of the domain are in the unit of pixel (pixel size 2 nm). 

For the pair correlation function: Open "DomainRadialDistribution_corrLen.m" and run it. It will calculate the pair correlation function for the cathode particle. One example of the domain center of mass coordinates for a cathode particleis provided as an excel file. The coordinates of the domain are in the unit of pixel (pixel size 2 nm).

(C) For the Mn white line ratio analysis and EELS mapping, visit folder "EELSmapping"

For the installation of the JupyterNotebook and the setup of the environment, please refer to the "Instructions.docx" file in the "installation" folder.

In the "AnalysisEELSWhiteLineRatio" folder, open "Hyperspy EELS Mn White Line Ratio.ipynb" using JupyterNotebook. Instructions on how to use the code are included in the beginning of the code and in the comments. Run the codes as instructed in the notebook. We used Pearson method for the calculations of Mn white line ratio (L3/L2). To apply Pearson method, it is necessary to first run the codes for "Curve Fitting Method" to determine curve center. An example of the low-loss and core-loss EELS data for a pristine cathode NP is provided for running the code.

The code for generating the map of Mn white line ratio (L3/L2), together with an example data and the corresponding mask of the cathode particle, is also provided in the "GenerateMap" folder.
