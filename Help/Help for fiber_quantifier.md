# Help for the function <i>fiber_quantifier</i>, v. 0.1.x

## Introduction

This help file contains information about
1) [Usage of the program](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help%20for%20fiber_quantifier.md#1-usage)
2) [Syntax](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help%20for%20fiber_quantifier.md#2-Syntax)
3) [Input Arguments](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help%20for%20fiber_quantifier.md#3-Input-Arguments)
4) [Output Arguments](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help%20for%20fiber_quantifier.md#4-Output-Arguments)
5) [Acknowledgements](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help%20for%20fiber_quantifier.md#5-Acknowledgements)
6) [Example Code](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help%20for%20fiber_quantifier.md#6-Example-Code)

## 1. Usage

The function fiber_quantifier is used to calculate the muscle architectural parameters pennation angle, fiber tract length, and curvature in the MuscleDTI_Toolbox. Calculations are only made for fiber tracts having 10 or more points. The method for pennation measurements is described in Lansdown et al, J Appl Physiol 2007, with slight modifications to improve computational efficiency. The method for curvature measurements is described in Damon et al, Magn Reson Imaging 2012. Best results are obtained with polynomial-fitted fiber tracts, calculated using fiber_smoother.

## 2. Syntax

[angle_list, distance_list, curvature_list, fiber_all_mm, n_points, apo_area] = ...
   fiber_quantifier(fiber_all, roi_mesh, fq_options)

## 3. Input Arguments
 
* fiber_all: A 4D matrix containing the fiber tract points, with units of pixels (in X and Y) or slice number (in Z).  

* roi_mesh: The mesh reconstruction of the aponeurosis that was used as the seed surface for fiber tracking, output from define_roi.

* fq_options: A user-defined structure containing the following fields:

   .dwi_res: a three element vector with the FOV, (assumed to be the same for the x and y directions), in-plane matrix size, and the slice thickness. The FOV and slice thickness must be specified in mm.

   .filt_kernel: Pennation angle is calculated as complement to the angle formed by position vectors along the tract and the normal vector to the seed point of the tract. To do so, two lines tangent to the point are calculated. When determining the lines, their slopes are median-filtered over the NxN window specified in filt_kernel; enter this as a single odd integer, e.g. fq_options.filt_kernel=N.  Note that median filtering results in a loss of edge information and that the number of lost rows and columns increases with the size of the filter kernel.

## 4. Output Arguments
* angle_list: The pennation angles, for fiber tracking point numbers 2-end. Pennation angles are reported in degrees.

* distance_list: Distances along the fiber tracts, for fiber tracking point numbers 1-end. Distances are reported in mm.

* curvature_list: The curvature values at eadh fiber tracking point, for fiber tracking point numbers starting at 2 and ending 3 points before the tract's end. The latter is to avoid abrupt changes in tract position.  Curvature values are reported in m-1.

* fiber_all_mm: The fiber tract points converted to units of mm.

* n_points: A 3D matrix (rows x columns x 3) containing the number of points used to quantify length, pennation, and curvature in each tract

* apo_area: A 2D matrix (rows x columns) containing the amount of apoeurosis area represented by each fiber tract

## 5. Acknowledgements
 People: Zhaohua Ding, Adam Anderson, Anneriet Heemskerk
 
 Grant support: NIH/NIAMS R01 AR050101, NIH/NIAMS R01 AR073831

## 6. Example Code
