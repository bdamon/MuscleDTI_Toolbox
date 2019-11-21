# Help for the function <i>fiber_track</i>, v. 0.1.x

## Introduction

This help file contains information about
1) [Usage of the program](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_track.md#1-usage)
2) [Syntax](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_track.md#2-Syntax)
3) [Input Arguments](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_track.md#3-Input-Arguments)
4) [Output Arguments](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_track.md#4-Output-Arguments)
5) [Acknowledgements](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_track.md#5-Acknowledgements)
6) [Example Code](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_track.md#6-Example-Code)


## 1. Usage

The function <i>fiber_track</i> is used to fiber-track a muscle DTI dataset in the MuscleDTI_Toolbox. 
   
The required inputs include a 5D matrix containing the diffusion tensor; the mask delimiting the muscle of interest, the mesh reconstruction of the aponeurosis of muscle fiber insertion, and a structure called ft_options.  This structure allows the user to set fiber tracking options and define important properties about the images. These settings are described in detail below (Output Arguments).

Fibers are tracked from the mesh according to the selected propogation algorithm until they reach the edge of the mask or meet another stop criterion.  Stop criteria are set in the ft_options structure. See the description of the input arguments for additional information on these variables. 
   
The outputs include the fiber tracts, several variables describing the outcomes of the tracking, and selected data about the tracts.

## 2. Syntax

[fiber_all, roi_flag, stop_list, fiber_len, fa_all, md_all] = ...

   fiber_track(tensor_m, mask, roi_mesh, ft_options, plot_options, anat_image);

## 3. Input Arguments

* <i>tensor_m</i>: A 5D matrix containing rows, columns, slices, and the 3x3 diffusion tensor, calculated from pre-processing steps.

* <i>mask</i>: The mask delimiting the muscle to be fiber-tracked. It could be the output of define_muscle or any other image analysis program that creates a binary mask of the same size as the imaging data.   

* <i>roi_mesh</i>: The roi mesh, output from define_roi.  

* <i>ft_options</i>: A structure containing the following fields:

  <i>.ref_frame</i>: The frame of reference in which the diffusion directions are specified. For example, set ft_options.ref_frame='LPS'; if the left, posterior, and superior anatomical positions are (+).

  <i>.image_orient</i>: The orientation of the images. Specify the anatomical positions at the north and east edges of the image as A (anterior) or P (posterior) and right (R) or left (L).  Input the result as a 2-element string variable (ft_options.image_orient='RA', 'AL', etc.).

  <i>.mesh_dist</i>: The number of pixels to shift the mesh into the muscle, prior to fiber tracking. This can be a (+) or (-) number, depending on the desired direction of the shift.  This is most helpful in unipennate muscles in cases where it's necessary to shift the mesh by a small amount to ensure that the seed points lies within the muscle of interest.

  <i>.depth_ratio</i>: The ratio of slice thickness/in-plane resolution. Note that the function assumes equal in-plane voxel dimensions.

  <i>.prop_alg</i>: A string variable that specifies the method for determining the direction of fiber tract propagation. The available options include 'euler', 'tnsrln', 'rk4', and 'fact'.  These algorithms are described below:

  * <i>Euler</i>: The observed diffusion tensor D is diagonalized at the current fiber tracking point to find the direction indicated by the first eigenvector. Euler integration occurs as this direction is multiplied by the step size and added to the current point. The user must specify the step size in the field ft_options.step_size.
  
  * <i>Tensorlines</i>: The tensorlines algorithm (Lazar et al, Human Brain Mapping, 2003). The tensorlines algorithm combines streamline and tensor deflection components when calculating the propagation direction. The tensor deflection component can be weighted between deflected and non-deflected terms using the parameter w_punct. w_punct can vary from 0 to 1, where 1 provides full weighting for the deflected component. To set w_punct, create a field called ft_options.w_punct. The user must also specify the value of the largest eigenvalue throughout the muscle. To do so, create a field called ft_options.eigvalue1_max and set it to the largest eigenvalue observed in the muscle of interest.  The tensorline direction is multiplied by the step size and added to the current point. The user must specify the step size in the field ft_options.step_size.
  
  * <i>4th order Runge-Kutta integration of the first eigenvector</i>: The observed diffusion tensor D is diagonalized at the current fiber tracking point to find the direction indicated by the first eigenvector (i.e., the slope); this is repeated at three future points by following this direction. These directions are averaged per the Runge-Kutta method to set the direction of fiber tract propagation. This direction is multiplied by the step size and added to the current point. The user must specify the step size in the field ft_options.step_size. Note that if the 2nd, 3rd, or 4th order points fall outside of the mask, the Runge-Kutta method would normally fail; therefore the propagation algorithm automatically changes to Euler integration to avoid premature tract termination.
  
  * <i>FACT</i>: The FACT algorithm, as described by Mori et al (Ann Neurol, 1999). FACT is similar to Euler integration, except that the direction is changed as soon as the tract enters a new voxel. The step size does not need to be included in the ft_options structure, as FACT necessarily involves variable step sizes.

  <i>.step_size</i>: The fiber-tracking step size, in pixels. A step size of 1 reflects the voxel width. This is not used for FACT, but is required for other algorithms.

  <i>.term_mthd</i>: A string variable that specifies the method for determining whether or not to terminate a fiber tract. Any fiber tracking point that falls outside of the image mask will terminate the tract. Other available options include 'bin1' and 'bin2'. FACT uses its own criterion.

  * <i>bin1</i>: The angle and FA data from the current fiber tracking point are used to decide whether or not to terminate the tract. The allowable range of FA values is set in the field fa_thresh (see below). The angle criterion use the angle formed by two fiber tracking steps. One of these is always the current step.  The other step can lie N points back (see below). 
  
  * <i>bin2</i>: The angle and FA criteria from the current fiber tracking point are combined with those of a preceding point. The user sets the FA criteria as for bin1. If the two consecutive points have a disallowed FA value, then the tract terminates. For the angle criterion, the step angles are calculated 1) between the current step and its immediate predecessor and 2) between the current step and one looking back N points (N must be > 1). If both steps exceed the angle threshold, then the tract terminates. This option provides greater tolerance for errors in individual voxels than bin1.
  
  <i>.angle_thrsh</i>: A two-element vector containing the angle threshold in degrees and the number of look-back steps (used as described under term_mthd)

  <i>.fa_thrsh</i>: a two-element vector containing the lower and upper bounds of allowable FA values (used as described under term_mthd)
  
  The FACT algorithm uses its own method for tract termination. Thus, when the propogation algorithm is set to FACT, the user does not need to define term_mthd; however the user must create fields in ft_options called ft_options.r_crit and ft_options.num_fact_voxels. These are used to terminate tracts based on local variability in the first eigenvector. The reader is referred to the Mori paper to learn more about these parameters.

  <i>.r_crit</i>: A scalar quantity ranging from 0-1 that defines the allowable level of local variability in the direction of the first eigenvector.
  
  <i>.num_fact_voxels</i>: The number of local neighboring voxels to include in the calculation of r.
 
* <i>plot_options</i>: If specified, this calls the fiber_visualizer function to plot the fiber, mask, and roi mesh.
 
* <i>anat_image</i>: The structural images, of the same size as the DTI images.  These are required only if the user wishes to plot the fiber tracts.

## 4. Output Arguments
* <i>fiber_all</i>: The fiber tracts, with units of pixels. The rows and columns correspond to locations on the roi_mesh. Dimension 3 gives point numbers on the tract, and the fourth dimension has row, column, and slice coordinates.

* <i>roi_flag</i>: A matrix indicating the presence of fibers that propagated at least 1 point

* <i>stop_list</i>: A matrix containing the reason for fiber tract termination (4=mask, 3=curvature, 2=FA, 1=R (for FACT only))

* <i>fiber_len</i>: The length, in points, of each fiber tract. 

* <i>fa_all</i>: The pointwise FA values on each fiber tract.

* <i>md_all</i>: The pointwise mean diffusivities along each tract

## 5. Acknowledgements

 People: Zhaohua Ding, Adam Anderson, Amanda Buck, Anneriet Heemskerk, and Justin Montenegro
 
 Grant support: NIH/NIAMS R01 AR050101, NIH/NIAMS R01 AR073831

## 6. Example Code
