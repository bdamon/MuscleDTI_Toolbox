# Help for the function [<i>fiber_goodness</i>](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Tractography-Functions/fiber_goodness.m), v. 1.0

## Introduction

This help file contains information about
1) [Purpose of the program](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_goodness.md#1-purpose)
2) [Usage of the program](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_goodness.md#2-usage)
3) [Syntax](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_goodness.md#3-Syntax)
4) [Example Code](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_goodness.md#4-Example-Code)
5) [Acknowledgements](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_goodness.md#5-Acknowledgements)

## 1. Purpose
The function fiber_goodness is used to assess the goodness of fiber tract data, and reject outlying fiber tract data, in the MuscleDTI_Toolbox. An optional, second-stage selection process allow the user to uniformaly sample the aponeurosis mesh.

[Back to the Top](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_goodness.md)

## 2. Usage
The quality algorithm described in [Heemskerk et al, 2008](https://pubmed.ncbi.nlm.nih.gov/19161166/) is implemented, but updated to account for the inclusion of curvature in the architectural computations. Specifically, the fibers are selected for having:
 1) Monotonically increasing values in the Z direction. In fiber_track, the Z component of the first eigenvector is forced to be positive; so the Z components of the unfitted tracts must increase monotonically. However, negative steps in the Z direction could result from overfitting fiber tract smoothing process. Fiber tracts with monotonically increasing Z values are indicated by entering a value of 1 at the corresponding [row column] indices in the 1st level of the 3rd dimension of qual_mask.  The number of tracts that meet this criterion are calculated and added to the 1<sup>st</sup> index of the vector <i>num_tracked</i>. Tracts that meet this criterion are advanced to the next level of selection
 2) A minimum length, in mm. This is specified by the user as a field in the structure <i>fs_options</i>. Fiber tracts that meet this criterion are indicated by entering a value of 1 at the corresponding [row column] indices in the 2<sup>nd</sup> level of the 3<sup>rd</sup> dimension of qual_mask. The total number of tracts that meet this criterion are calculated and added to the 2nd index into the vector num_tracked. Tracts that meet this criterion are advanced to the next level of selection.
 3) An acceptable pennation angle, in degrees.  This range is defined by the user in the structure <i>fs_options</i>. Fiber tracts that meet this criterion are indicated by entering a value of 1 at the corresponding [row column] indices in the 3rd level of the 3<sup>rd</sup> dimension of <i>qual_mask</i>. The total number of tracts that meet this criterion are calculated and added to the 3rd index into the vector <i>num_tracked</i>. Tracts that meet this criterion are advanced to the next level of selection. 
 4) An acceptable curvature value, in m<sup>-1</sup>. This range is defined by the user in the structure <i>fs_options</i>. Fiber tracts that meet this criterion are indicated by entering a value of 1 at the corresponding [row column] indices in the 4<sup>th</sup> level of the 3<sup>rd</sup> dimension of <i>qual_mask</i>. The total number of tracts that meet this criterion are calculated and added to the 4<sup>th</sup> index into the vector <i>num_tracked</i>. Tracts that meet this criterion are advanced to the next level of selection.
 5) A length, pennation angle, and curvature that lies within the 95% confidence interval for length, pennation angle, and curvature set by the surrounding 24 tracts. Fiber tracts that meet this criterion are indicated by entering a value of 1 at the corresponding [row column] indices in the 5<sup>th</sup> level of the 3rd dimension of <i>qual_mask</i>. The total number of tracts that meet this criterion are calculated and added to the 5<sup>th</sup> index into the vector <i>num_tracked</i>. Tracts that meet this criterion are advanced to the next step of analysis.

The use of length, pennation, and curvature criteria require the user to use their knowledge of the expected patterns of muscle geometry to supply values that are reasonable but will not inappropriately bias the results. These selection criteria, as well as the number of tracts that were rejected because of these criteria, should be included in the Methods sections of publications.
  
An optional final of the selection process is to sample fiber tracts across the aponeurosis at a uniform spatial frequency.  Because the aponeurosis changes in width over the superior-inferior direction, the distance between seed points in the <i>roi_mesh</i> matrix varies. This would bias a simple average of the whole-muscle architectural properties toward narrower regions of the aponeurosis, where the seed points occur at higher sampling frequency. 

To avoid this problem, the user may include a field called sampling_frequency in the <i>fs_options</i> structure; if present, this produces an approximately uniform spatial sampling of the fiber tracts. To do so, <i>fiber_selector</i> first identifies all tracts that fall within the boundaries of each sampling interval. Then, the median values for length, pennation angle, and curvature are calculated. For each tract, a similarity index S is calculated that compares the length, mean pennation angle, and mean curvature of the Tth tract and the local median. The tract with the minimum value of S is taken as the most typical tract in the sampling interval and is used in further steps. If sampling_frequency is not included in <i>fs_options</i>, then this sampling does not occur and the tracts defined in step 5 are used in further steps. 

In either case, the preserved fiber tracts are stored in the matrix <i>final_fibers</i>; their structural properties are stored in the matrices <i>final_curvature</i>, <i>final_angle</i>, and <i>final_distance</i>; and the whole-muscle mean values for length, pennation angle, and curvature are stored in the matrix <i>mean_apo_props</i>.

[Back to the Top](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_goodness.md)

## 3. Syntax
[final_fibers, final_curvature, final_angle, final_distance, qual_mask, num_tracked, mean_fiber_props, mean_apo_props] = 

fiber_goodness(fiber_all, angle_list, distance_list, curvature_list, n_points, roi_flag, apo_area, roi_mesh, fg_options)

The input arguments are:
 * <i>fiber_all</i>: The fiber tracts from which selection will be made.  The matrix could be the output of <i>fiber_track</i> or the smoothed fiber tracts output from <i>fiber_smoother</i>

 * <i>angle_list, distance_list, curvature_list, n_points, apo_area</i>: The outputs of <i>fiber_quantifier</i>.

 * <i>roi_flaga</i>: The outputs of <i>fiber_quantifier</i>.

 * <i>roi_mesh</i>: A mask indicating fiber tracts that propagated at least one point, output from <i>fiber_track</i>;

 * <i>fg_options</i>: A structure containing user-specified criteria for selecting the tracts:
 
     <i>.dwi_res</i>: A three-element vector containing the field of view (mm), matrix size, and slice thickness (mm) of the diffusion-weighted images
 
     <i>.min_distance</i>: minimum distance for selected tracts, in mm
     
     <i>.min_pennation</i>: minimum pennation angle, in degrees 
     
     <i>.max_pennation</i>: maximum pennation angle, in degrees 
     
     <i>.max_curvature</i>: maximum curvature, in m<sup>-1</sup>
     
     <i>.sampling_density</i> (optional): The spatial frequency for uniform sampling of the aponeurosis mesh, in mm<sup>-1</sup>

The output arguments are:
 * <i>final_fibers</i>: the fiber tracts that passed all selection criteria

 * <i>final_curvature</i>: pointwise measurements of curvature for the final tracts.

 * <i>final_angle</i>: pointwise measurements of pennation angle for the final tracts.

 * <i>final_distance</i>: pointwise measurements of cumulative distance for the final tracts

 * <i>qual_mask</i>: a 3D matrix of the same row x column size as the <i>roi_mesh</i>, with 6 layers corresponding to each stage of the selection process. In each layer, ones correspond to retained fibers and zeros correspond to rejected fibers

 * <i>num_tracked</i>: the number of fibers for each of the following steps:
   1) the number of these tracts generated by <i>fiber_track</i>;
   2) the number of these tracts that were quantified by <i>fiber_quantifer</i>; and
   3-7) the number of fiber tracts that met criteria 1-5 above, respectively.

 * <i>mean_fiber_props</i>: A 3D matrix (rows x columns x 5) containing the mean curvature, pennation, and length values along each of the tracts; the amount of aponeurosis area represented by each tract; and the number of points in each tract.

 * <i>mean_apo_props</i>: A 1 x 3 vector containing the whole-muscle mean values for curvature, pennation, and fiber tract length, weighted by the amount of aponeurosis area represented by each tract.

[Back to the Top](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_goodness.md)

## 4. Example Code
Given:

1.	DTI images having matrix size 19219244, field of view 192192 mm, and slice thickness 7 mm; and

2.	The aponeurosis mesh stored in a variable called roi_mesh; 

3.	The fiber tracts generated by fiber_track 

4.	The pennation, length, curvature, and other data generated by fiber_quantifier 

the following code will allow the user to 

1.	Assess the goodness of the results and reject results that do not match the specified criteria (minimum length = 35 mm, range of pennation angles =0-40⁰, maximum curvature =40 m-1).

2.	Uniformly sample the aponeurosis at a spatial frequency of 1 tract/10 mm2.

% Set fiber goodness options

fg_options.dwi_res=[192 192 7]; % FOV (mm), matrix size, and dZ (mm) of DTI images 

fg_options.min_distance = 35; % fiber tracts must be 35 mm long 

fg_options.min_pennation = 0; % pennation angles must be positive 

fg_options.max_pennation = 40; % pennation angles must be smaller than 40 degrees

fg_options.max_curvature = 40; % curvatures must be smaller than 40 m^-1

fg_options.sampling_frequency= 10; % uniformly sample the aponeurosis with 1 tract/10 mm^2


% Call the function:

[final_fibers, final_curvature, final_angle, final_distance, qual_mask, num_tracked, mean_fiber_props, mean_apo_props] = 

fiber_goodness(fiber_all, angle_list, distance_list, curvature_list, n_points, roi_flag, apo_area, roi_mesh, fg_options) 


[Back to the Top](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_goodness.md)

## 5. Acknowledgements
 People: Zhaohua Ding, Anneriet Heemskerk
 
 Grant support: NIH/NIAMS R01 AR050101, NIH/NIAMS R01 AR073831

[Back to the Top](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_goodness.md)
