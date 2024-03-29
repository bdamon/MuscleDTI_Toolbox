# Help for the function [<i>fiber_smoother</i>](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Tractography-Functions/fiber_smoother.m), v. 1.0.0

## Introduction

This help file contains information about
1) [Purpose of the Program](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_smoother.md#1-Purpose)
2) [Usage of the Program](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_smoother.md#2-Usage)
3) [Syntax](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_smoother.md#3-Syntax)
4) [Example Code](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_smoother.md#4-Example-Code)
5) [Acknowledgements](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_smoother.md#5-Acknowledgements)

## 1. Purpose
The function <i>fiber_smoother</i> is used to smooth and increase the spatial resolution along fiber tracts by performing a polynomial fit. 

[Back to the Top](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_smoother.md)

## 2. Usage

The user inputs the fiber tracts generated by [<i>fiber_track</i>](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_track.md) and a structure with options for implementing the polynomial fitting routines.  The [row column slice] positions are separately fitted to N<sup>th</sup>-order polynomials as functions of voxel-distance along the tract. The user selects the polynomial order separately for the row, column, and slice positions. The function returns the smoothed fiber tracts and matrices containing the polynomial coefficients for each point. 

To improve fitting and ensure that the fitted tract continues to originate from the seed point, the seed point is subtracted from the fiber tract prior to fitting. Then the <i>polyfit</i> function is used to fit the remaining row, column, and slice positions to polynomial functions. Finally, the <i>polyval</i> function is used to solve the polynomials at interpolation distances of interpolation_step.

This procedure is modified from [previous work](https://pubmed.ncbi.nlm.nih.gov/22503094/) to: 1) Fit the tract positions as functions of distance rather than point number, as required for tracking algorithms that use variable step sizes, such as FACT; and 2) Allow user selection of the polynomial order, including different polynomial orders for the row/column/slice positions.

[Back to the Top](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_smoother.md)

## 3. Syntax
[smoothed_fiber_all, pcoeff_r, pcoeff_c, pcoeff_s, n_points_fitted] = fiber_smoother(fiber_all, fs_options);

### Input Arguments

* <i>fiber_all</i>: The original fiber tracts, output from fiber_track

* <i>fs_options</i>: A structure containing the following fields:

   <i>.interpolation_step</i>: An interpolation interval for the fitted fiber tract, in units of pixels.  For example, setting interpolation_step to 0.25 would interpolate the fiber tract at intervals of 0.25 pixels.

   <i>.p_order</i></i>: A three-element vector containing the polynomial orders N<sub>R</sub>, N<sub>C</sub>, and N<sub>S</sub> to use when fitting the tracts; they should be specified in the sequence row-column-slice.

### Output Arguments

* <i>fitted_fiber_all</i>: A 4D matrix containing the smoothed fiber tracts

* <i>pcoeff_r, pcoeff_c, pcoeff_s</i>: Matrices of the polynomial coefficients for the tracts' row, column, and slice positions as functions of voxel distance, with dimensions of N<sub>R,A</sub> x N<sub>C,A</sub> x (n+1) 

* <i>n_points_smoothed</i>: The number of points in the fitted tracts.

[Back to the Top](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_smoother.md)

## 4. Example Code
Given:

1.	An anatomical image with variable name anat_image and having matrix size 192x192x44, field of view 192x192 mm, and slice thickness 7 mm;

2.	The aponeurosis mesh stored in a variable called roi_mesh; 

3.	DTI images having matrix size 192x192x44, field of view 192x192 mm, and slice thickness 7 mm; and

4.	The fiber tracts generated by fiber_track 

the following code will allow the user to 

1.	Smooth the fiber tracts using 3rd, 3rd, and 2nd order polynomials for the row, column, and slice directions, respectively; and

2.	Visualize the outcome, with the mesh displayed in gray and the fiber tracts displayed in red.

% Set smoothing options

fs_options.interpolation_step = 0.5;

fs_options.p_order = [3 3 2];
 
% Call the function:

[smoothed_fiber_all, pcoeff_r, pcoeff_c, pcoeff_s, n_points_fitter] = fiber_smoother(fiber_all, fs_options);

%Visualize the outcome:

plot_options.plot_mesh = 1; %do plot the aponeurosis mesh

plot_options.plot_mask = 0; %don't plot the mask

plot_options.plot_fibers = 1; %do plot the fibers

plot_options.anat_dims = [192 7]; %FOV and slice thickness of the images to be displayed, in mm

plot_options.anat_slices = 14:10:44; %display slices 14, 24, 34, and 44 

plot_options.mesh_size = [192 192]; %in-plane matrix size of the images used to create the mesh

plot_options.mesh_dims = [192 7]; %FOV and slice thickness of the images used to create the mesh, in mm

plot_options.mesh_color = [.75 .75 .75]; %make the mesh light gray

plot_options.fiber_color = [.8 .2 .2];

plot_options.dti_size = [192 192];

plot_options.dti_dims= [192 7];

smoothed_fiber_figure = fiber_visualizer(anat_image, plot_options, roi_mesh, mask, smoothed_fiber_all);

[Back to the Top](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_smoother.md)

## 5. Acknowledgements
People: Zhaohua Ding, Anneriet Heemskerk

Grant support: NIH/NIAMS R01 AR050101, NIH/NIAMS R01 AR073831

[Back to the Top](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_smoother.md)

