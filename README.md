# OrientedPhaseDomain

This code package is developed to map the strain, oriented phase domain, and pair distribution function for the spinel MnO2 cathode particles during Mg ion insertion.
Current version: 1.0
Date: 7/2022
Test: Matlab 2016b on Microsoft Windows 10 using an Intel(R) Core(TM) i5-8250U CPU @ 1.60GHz processor
For more information about the project, algorithms, and related publications please refer to the Chen Group website.

# Reference

# Getting started
(A) For mapping the strain and domain structures in cathode particles
The matrixes of g-vectors are obtained from the 4D-STEM data using the “imToolBox” software. An example of matrix of g-vectors are provided ((g220).xlsx and (g002).xlsx). An example of mask is provided. The following steps should be executed only once (installation):
(Matlab version: R2016a)

1. Open "calculate_strain_g220_g002.m" and run it.
2. Open "CreatPremask.m" and run it. Import ""mask-Pre.tif" in ImageJ and manually create a mask by delienating the particle shape.
3. Open "MapStrainforXandY.m" and run it. Save strain maps.
4. Open "CalculateTetragonalityHistogram.m" and run it. Histogram will be generated.
5. Open "FitTetragonalityHist.m" and run it. Use the "cftool" for the fitting of the histogram. Set the constraints on the fitting parameters as needed to obtain a good fitting.
6. Open "MapDomain.m". Set the threshold from the fitting results to distinguish the [100]t and [111]t. An example of thresholds is provided in the code. Run it. Save domain map.

Note: "cmapStrain.mat" and "cmapDomain.mat" files are for the color scales of the strain and domain maps, respectively.

(B) For the calculation of pair distribution function and pair correlation function
For the pair distribution function: Examples of the domain center of mass coordinates are provided (in the unit of pixel, pixel size 2 nm). Execute the code "DomainRadialDistribution_averageFiles.m" will calculate and plot pair distribution function.
For the pair correlation function: Examples of the domain center of mass coordinates are provided (in the unit of pixel, pixel size 2 nm). Execute the code "DomainRadialDistribution_corrLen.m" will calculate the pair correlation function.
