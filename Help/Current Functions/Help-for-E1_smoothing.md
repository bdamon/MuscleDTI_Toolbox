# Help for [<i>E1_smoothing</i>](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Preprocessing-Functions/E1_smoothing.m), v. 1.0.0.

## Introduction

This help file contains information about
1) [Purpose of the Program](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Current%20Functions/Help-for-E1_smoothing.md#1-Purpose)
2) [Usage of the Program](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Current%20Functions/Help-for-E1_smoothing.md#2-Usage)
3) [Syntax](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Current%20Functions/Help-for-E1_smoothing.md#3-Syntax)
4) [Example Code](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Current%20Functions/Help-for-E1_smoothing.md#4-Example-Code)
5) [Acknowledgements](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Current%20Functions/Help-for-E1_smoothing.md#5-Acknowledgements)

## 1. Purpose 
This function provides a fast, automated, and robust discretized spline smoothing for data of arbitrary dimension. For skeletal muscle DTI applications, the function is used to smooth the muscle's principal diffusion eigenvector field. 

## 2. Usage
The user inputs the first eigenvector map of the muscle as a 4D array with [x,y,z,3] entries. The array corresponds to the 3D distribution of each spatial coordinate of the principal diffusion eigenvector field. 

## 3. Syntax
[E1map_smooth] = E1_smoothing(E1_map, input_smoothing_parameter) 

The input arguments are:
* <i>data:</i> [e1map] 4D volume containing the principal diffusion eigenvector field of the muscle (x,y,z,3)
* <i>data:</i> [input_smoothing_parameter] optional, level of smoothing can also be prespecified by the user (not defined by the algorithm), with a higher input smoothing parameter resulting in more smoothing. 
  
The output arguments are:
* <i>E1_map_smooth:</i> (x,y,z,3) smoothed

## 4. Example code

## 5. Acknowledgements
This function depends on the smoothn function created by Damien Garcia and is available in the MATLAB File exchange (https://www.mathworks.com/matlabcentral/fileexchange/25634-smoothn) and in the BiomeCardio website (https://www.biomecardio.com/matlab/smoothn_doc.html).

REFERENCES
Pineda Guzman, Roberto A.; Lockard, Carly A.; Zhou, Xingyu; Damon, Bruce M. "Improved DTI-Based Skeletal Muscle Architecture Estimation via Diffusion Weighted Image Denoising and First Eigenvector Map Smoothing." NMR in Biomedicine, vol. 38, no. 9, pp. e70100, July 2025.  
Garcia, Damien. “Robust Smoothing of Gridded Data in One and Higher Dimensions with Missing Values.” Computational Statistics & Data Analysis, vol. 54, no. 4, Elsevier BV, Apr. 2010, pp. 1167–78, doi:10.1016/j.csda.2009.09.020.
Garcia, Damien. “A Fast All-in-One Method for Automated Post-Processing of PIV Data.” Experiments in Fluids, vol. 50, no. 5, Springer Science and Business Media LLC, Oct. 2010, pp. 1247–59, doi:10.1007/s00348-010-0985-y.
