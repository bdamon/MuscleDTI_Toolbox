# Help for the function <i>fiber_fitter</i>, v. 0.1.x

## Introduction

This help file contains information about
1) [Usage of the program](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_fitter.md#1-usage)
2) [Syntax](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_fitter.md#2-Syntax)
3) [Input Arguments](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_fitter.md#3-Input-Arguments)
4) [Output Arguments](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_fitter.md#4-Output-Arguments)
5) [Acknowledgements](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_fitter.md#5-Acknowledgements)
6) [Example Code](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_fitter.md#6-Example-Code)

## 1. Usage
The function <i>fiber_fitter</i> is used to smooth fiber tracts and increase the spatial resolution of fiber tracts generated using the MuscleDTI_Toolbox. The x, y, and z positions are separately fitted to Nth order polynomials as functions of distance along the tract and are uniformly solved at interpolation distances of interpolation_step. The user selects the polynomial order.  The coefficients are estimated using the MATLAB function <i>polyfit</i> and smoothed tracts are generated using the function <i>polyval</i>. This procedure is modified from Damon et al, Magn Reson Imaging, 2012 in two ways. 
* The tract positions are fitted as functions of distance along the fiber tract rather than point number. This is required for tracking algorithms that use variable step sizes, such as FACT.  
* <i>fiber_fitter</i> allows the user to select the polynomial order separately for the X, Y,and Z coordinates. 

In the current version, it is the user's responsibility to select the polynomial order appropriately.  One way to do this would be to plot the points for an individual fiber tract and then plot the fitted results over that. Sample code for doing this is included below. However, future development will included automated selection of the most appropriate polynomial order.

## 2. Syntax
[angle_list, distance_list, curvature_list, fiber_all_mm, n_points, apo_area] = ...

fiber_quantifier(fiber_all, roi_mesh, fq_options);

## 3. Input Arugments

* <i>fiber_all</i>: The original fiber tracts, output from fiber_track

* <i>ff_options</i>: A structure containing the following fields:

   <i>.interpolation_step</i>: An interpolation interval for the fitted fiber tract, in units of pixels.  For example, setting interpolation_step to 0.25 would interpolate the fiber tract at intervals of 0.25 pixels.

   <i>.p_order</i></i>: A 3 element vector containing the polynomial orders, [Nx Ny Nz], to use when fitting the tracts. If a single number is entered (i.e., ff_options.p_order=N), then this order is used to fit the points in the X, Y and Z directions.

## 4. Output Arguments

* <i>fitted_fiber_all</i>: The fiber tracts following Nth order polynomial fitting

* <i>pcoeff_x</i>: A matrix of the polynomial coefficients for the tracts' x positions. The user is referred to help for the <i>polyfit</i> function to understand how to interpret these coefficients.

* <i>pcoeff_y</i>: A matrix of the polynomial coefficients for the tracts' y positions 

* <i>pcoeff_z</i>: A matrix of the polynomial coefficients for the tracts' z positions 

## 5. Acknowledgements
People: Zhaohua Ding, Anneriet Heemskerk

Grant support: NIH/NIAMS R01 AR050101, NIH/NIAMS R01 AR073831

## 6. Example Code
