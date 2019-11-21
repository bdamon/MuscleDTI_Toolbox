# Help for the function <i>fiber_selector</i>, v. 0.1.x

## Introduction

This help file contains information about
1) [Usage of the program](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_selector.md#1-usage)
2) [Syntax](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_selector.md#2-Syntax)
3) [Input Arguments](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_selector.md#3-Input-Arguments)
4) [Output Arguments](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_selector.md#4-Output-Arguments)
5) [Acknowledgements](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_selector.md#5-Acknowledgements)
6) [Example Code](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_selector.md#6-Example-Code)

## 1. Usage
The function fiber_selector is used to sample from a set of quantified fiber tracts geerated using the MuscleDTI_Toolbox in a way that uniformly and optimally characterizes the muscle's architectural properties.Two stages of selection are possible. The first, required, step is to implement an updated version of the quality algorithm described in Heemskerk et al, Magn Reson Med, 2008. Specifically, the fiber tracts are selected for having:

   * Monotonically increasing values in the Z direction. This prevents errors due to overfitting in the Z direction; 
   
   * A minimum length (in mm);
   
   * A range of acceptable pennation angles (in degrees);
   
   * A maximum curvature value (in m<sup>-1</sup>); and
   
   * Values within the 95% confidence intervals for length, pennation angle, and curvature of the surrounding 24 tracts.
  
The use of length, pennation, and curvature criteria require the user to use their knowledge of the expected patterns of muscle geometry to supply values that are reasonable but will not inappropriately bias the results. These selection criteria, as well as the number of tracts that were rejected because of these criteria, should be included in the Methods sections of publications.
  
The second step of the selection process is optional; it is engaged by including a field .sampling_density in the <i>fs_options</i> structure described below. The desirability of using the second-stage process is derived from the fact that the roi_mesh is required to have a fixed number of rows and columns throughout; but the aponeurosis i>in vivo</i> varies in width.  Thus, the density of initially propagated fiber tracts, in 1/mm<sup>2</sup>, varies throughout the mesh. Also, the first-stage selection process rejects some tracking results.  To account for these issues, the fiber tracts that pass the first stage of the selection process can be selected so that they occur at a uniform spatial frequency in the final dataset. The user sets this frequency in fs_options.sampling_density.
  
## 2. Syntax
 [final_fibers, final_curvature, final_angle, final_distance, qual_mask, num_tracked, mean_fiber_props, mean_apo_props] = ...
 
   fiber_selector(fitted_fiber_all, angle_list, distance_list, curvature_list, n_points, apo_area, roi_flag, fs_options)

## 3. Input Arguments
 * <i>fitted_fiber_all</i>: the smoothed fiber tracts

 * <i>angle_list</i>: Pennation angles for smoothed fiber tracts, derived from <i>fiber_quantifier</i>;

 * <i>distance_list</i>: Distance measures for smoothed fiber tracts, derived from <i>fiber_quantifier</i>;

 * <i>curvature_list</i>: Curvature measures for smoothed fiber tracts, derived from <i>fiber_quantifier</i>;

 * <i>n_points</i>: The number of points quantified per fiber tract, derived from <i>fiber_quantifier</i>;

 * <i>apo_area</i>: A matrix indicating the amount of aponeurosis area associated with each fiber tract, derived from <i>fiber_quantifier</i>;

 * <i>roi_flag</i>: A mask indicating fiber tracts that propagated at least one point in the function fiber_track, derived from <i>fiber_track</i>;

 * <i>fs_options</i>: A structure containing user-specified criteria for selecting the tracts:
 
     <i>.min_distance</i>: minimum distance for selected tracts, in mm
     
     <i>.min_pennation</i>: minimum pennation angle, in degrees 
     
     <i>.max_pennation</i>: maximum pennation angle, in degrees 
     
     <i>.max_curvature</i>: maximum curvature, in m<sup>-1</sup>
     
     <i>.sampling_density</i>: 

## 4. Output Arguments
 * <i>final_fibers</i>: the fiber tracts that passed all selection criteria

 * <i>final_curvature</i>: pointwise measurements of curvature for the final tracts.

 * <i>final_angle</i>: pointwise measurements of pennation angle for the final tracts.

 * <i>final_distance</i>: pointwise measurements of cumulative distance for the final tracts

 * <i>qual_mask</i>: a 3D matrix of the same row x column size as the roi_mesh, with 6 layers corresponding to each stage of the selection process. In each layer, ones correspond to retained fibers and zeros correspond to rejected fibers

 * <i>num_tracked</i>: the number of fibers for each of the following steps:
   1) the number of potential fiber tracts;
   2) the number of these tracts generated by fiber_track;
   3-7) the number of fiber tracts that met criteria 1-5 above, respectively.

 * <i>mean_fiber_props</i>: A 3D matrix (rows x columns x 5) containing the mean curvature, pennation, and length values along each of the tracts; the amount of aponeurosis area represented by each tract; and the number of points in each tract

 * <i>mean_apo_props</i>: A 1 x 3 vector containing the whole-muscle mean values for curvature, pennation, and fiber tract length, weighted by the amount of aponeurosis area represented by each tract.

## 5. Acknowledgements
 People: Zhaohua Ding, Anneriet Heemskerk
 
 Grant support: NIH/NIAMS R01 AR050101, NIH/NIAMS R01 AR073831

## 6. Example Code
