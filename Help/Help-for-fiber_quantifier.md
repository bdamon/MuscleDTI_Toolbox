# Help for the function <i>fiber_quantifier</i>, v. 0.1.x

## Introduction

This help file contains information about
1) [Usage of the program](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_quantifier.md#1-usage)
2) [Syntax](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_quantifier.md#2-Syntax)
3) [Input Arguments](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_quantifier.md#3-Input-Arguments)
4) [Output Arguments](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_quantifier.md#4-Output-Arguments)
5) [Acknowledgements](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_quantifier.md#5-Acknowledgements)
6) [Example Code](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_quantifier.md#6-Example-Code)

## 1. Purpose
 
The function <i>fiber_quantifier</i> is used to calculate the muscle architectural parameters pennation angle, fiber tract length, and curvature in the MuscleDTI_Toolbox. 

[Back to the top](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_quantifier.md)

## 2. Usage
Calculations are only made for fiber tracts having six or more points. Information about each measurement follows:

* <i>Fiber tract length</i>: this is measured by summing the inter-point distances along the tract.

* <i>Pennation</i>: The method for pennation measurements is essentially as described in [Lansdown et al, J Appl Physiol 2007](https://pubmed.ncbi.nlm.nih.gov/17446411/). The approach to measuring pennation angle traditionally used in ultrasound imaging is to manually specify the line tangent to a muscle fascicle at the point of its insertion into the aponeurosis and a second line tangent to the aponeurosis. The angle formed by these lines is measured. In Lansdown et al., this concept was extended into 3D space by defining the plane tangent to the seed point and its normal vector. Also, position vectors between the seed point and points along thetract were defined.  Pennation angle was defined as the complement to the angel formed by the normal vector and the position vectors.  In the toolbox, this is modified slightly to calculated the normal vector as the cross product between two tangent lines to the seed point (one lying in the row direction of the aponeurosis mesh and one lying in the column direction). This change improves computational efficiency. 

* <i>Curvature</i>: The method for curvature measurements is described in [Damon et al, Magn Reson Imaging 2012](https://pubmed.ncbi.nlm.nih.gov/22503094/). Briefly, these use a discrete implementation of the Frenet-Serret equations. Specifically, the curvature K is defined in 

     dT/ds = K N
     
  where T is the tangent line to points along the curve, s is the step length between points, and N is the normal vector. In <i>fiber_quantifier</i>, K is calculated by multiplying each side of this equation by the Moore-Penrose pseudoinverse matrix of N.

For curvature, the best results are obtained with polynomial-fitted fiber tracts, calculated using <i>fiber_fitter</i>. Fiber tract length and pennation angle are unaffected by polynomial fitting.

[Back to the top](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_quantifier.md)

## 3. Syntax

[angle_list, distance_list, curvature_list, fiber_all_mm, n_points, apo_area] = ...
   fiber_quantifier(fiber_all, roi_mesh, fq_options);

The input arguments are:
 
* <i>fiber_all</i>: A 4D matrix containing the fiber tract points, with units of pixels (in X and Y) or slice number (in Z). This matrix could be substituted with fitted_fiber_all (the output of [<i>fiber_fitter</i>](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_fitter.md)

* <i>roi_mesh</i>: The mesh reconstruction of the aponeurosis that was used as the seed surface for fiber tracking, output from <i>define_roi</i>.

* <i>fq_options</i>: A user-defined structure containing the following fields:

   <i>.dwi_res</i>: a three element vector with the FOV, (assumed to be the same for the x and y directions), in-plane matrix size, and the slice thickness. The FOV and slice thickness must be specified in mm.

   <i>.filt_kernel</i>: As noted above, the pennation angle calculation relies on the measurement of two tangent lines along the roi_mesh. When determining these lines, their slopes are median-filtered over an NxN window as a way to mitigate the effects of digitization errors. The size of the median filter is specified in filt_kernel; enter this as a single odd integer, e.g. fq_options.filt_kernel = N.  Note that median filtering results in a loss of edge information; the number of rows and columns that will be lost increases with the size of the filter kernel.

The output arguments are:
* <i>angle_list</i>: The pennation angles, for fiber tracking point numbers 2-end. Pennation angles are reported in degrees.

* <i>distance_list</i>: Distances along the fiber tracts, for fiber tracking point numbers 1-end. Distances are reported in mm.

* <i>curvature_list</i>: The curvature values at each fiber tracking point, for fiber tracking point numbers starting at 2 and ending 3 points before the tract's end. The latter is to avoid abrupt changes in tract position.  Curvature values are reported in 1/m.

* <i>fiber_all_mm</i>: The fiber tract points converted to units of mm.

* <i>n_points</i>: A 3D matrix (rows x columns x 3) containing the number of points used to quantify length, pennation, and curvature in each tract

* <i>apo_area</i>: A 2D matrix (rows x columns) containing the amount of apoeurosis area represented by each fiber tract

[Back to the top](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_quantifier.md)

## 4. Example Code
The code below will measure fiber tract length, pennation angle, and curvature in fitted fiber tracts derived from DTI data have a matrix size of 192x192, a field of view of 192x192 mm, and a slice thickness of 7 mm; normal vectors to the seed points onteh aponeurosis mesh will be calcualted after smoothing the sloces over a 5x5 kernel:

%% Set options:

fq_options.dwi_res = [192 192 7];            %images have 192x192 FOV, 192x192 matrix, and 7 mm slice thickness

fq_options.filt_kernel = 5;                  %kernel size

%% call the function:

angle_list, distance_list, curvature_list, fiber_all_mm, n_points, apo_area] = ...

   fiber_quantifier(fiber_all, roi_mesh, fq_options);

[Back to the top](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_quantifier.md)
## 5. Acknowledgements
 People: Zhaohua Ding, Adam Anderson, Anneriet Heemskerk
 
 Grant support: NIH/NIAMS R01 AR050101, NIH/NIAMS R01 AR073831

[Back to the top](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_quantifier.md)
