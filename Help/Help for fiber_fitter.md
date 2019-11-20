# Help for the function <i>fiber_fitter</i>, v. 0.1.x

## Introduction

This help file contains information about
1) [Usage of the program](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help%20for%20fiber_fitter.md#1-usage)
2) [Syntax](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help%20for%20fiber_fitter.md#2-Syntax)
3) [Input Arguments](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help%20for%20fiber_fitter.md#3-Input-Arguments)
4) [Output Arguments](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help%20for%20fiber_fitter.md#4-Output-Arguments)
5) [Acknowledgements](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help%20for%20fiber_fitter.md#5-Acknowledgements)
6) [Example Code](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help%20for%20fiber_fitter.md#6-Example-Code)

## 1. Usage
The function <i>fiber_fitter</i> is used to smooth fiber tracts and increase the spatial resolution of fiber tracts generated using the MuscleDTI_Toolbox. The x, y, and z positions are separately fitted to Nth order polynomials as functions of distance along the tract and are uniformly solved at interpolation distances of interpolation_step. The user selects the polynomial order.  This procedure is modified from Damon et al, Magn Reson Imaging, 2012 to: 
   1) fit the tract positions as functions of distance rather than point number. This is required for tracking algorithms that use variable step sizes, such as FACT.
   2) allow selection of the polynomial order, including different polynomial orders for the X, Y, and X points.  
Future development will included automated selection of the most appropriate polynomial order.

## 2. Syntax
[angle_list, distance_list, curvature_list, fiber_all_mm, n_points, apo_area] = ...
fiber_quantifier(fiber_all, roi_mesh, fq_options);

## 3. Input Arugments

* fiber_all: The original fiber tracts, output from fiber_track

* ff_options: A structure containing the following fields:

   .interpolation_step: an interpolation interval for the fitted fiber tract, in units of pixels.  For example, setting interpolation_step to 0.25 would interpolate the fiber tract at intervals of 0.25 pixels.

   <i>.p_order</i>: A 3 element vector containing the polynomial orders, [Nx Ny Nz], to use when fitting the tracts. If a single number is entered (i.e., ff_options.p_order=N), then this order is used to fit the points in the X, Y and Z directions.

## 4. Output Arguments

* fitted_fiber_all: the fiber tracts following Nth order polynomial fitting

* pcoeff_x: a matrix of the polynomial coefficients for the tracts' x positions 

* pcoeff_y: a matrix of the polynomial coefficients for the tracts' y positions 

* pcoeff_z: a matrix of the polynomial coefficients for the tracts' z positions 

## 5. Acknowledgements
People: Zhaohua Ding, Anneriet Heemskerk

Grant support: NIH/NIAMS R01 AR050101, NIH/NIAMS R01 AR073831

## 6. Example Code
