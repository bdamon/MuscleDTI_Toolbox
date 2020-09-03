# Help for the function <i>fiber_track</i>, v. 0.1.x

## Introduction

This help file contains information about
1) [Purpose of the Program](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_track.md#1-purpose)
2) [Usage of the Program](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_track.md#2-usage)
3) [Tracking Algorithms](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_track.md#3-tracking-algorithms)
4) [Syntax](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_track.md#4-Syntax)
5) [Example Code](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_track.md#5-Example-Code)
6) [Acknowledgements](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_track.md#6-Acknowledgements)


## 1. Purpose

The function <i>fiber_track</i> is used to fiber-track a muscle DTI dataset. 
   
## 2. Usage
The required inputs include a 5D matrix hold the diffusion tensor at each voxel and [row column slice] dimensions matching those of the DTMRI data; the muscle mask, output from define_muscle or other program; the aponeurosis mesh, output from define_roi; and a structure defining the fiber-tracking options.  This structure allows the user to set options such as the tracking algorithm, step size, laboratory frame of reference, image orientation, and tract termination method. 

Fibers are tracked from the aponeurosis mesh according to the selected propagation algorithm. Each fiber tract is propagated until it reaches the edge of the mask or meets another stop criterion, such as an excessive inter-segment angle or an out-of-bounds value for fractional anisotropy (FA).  See the description of the input arguments for additional information on these criteria. 

The outputs include the fiber tracts, variables describing the outcomes of the tracking, and selected data about the tracts.

## 3. Tracking Algorithms

### Consideration for Laboratory Frame of Reference and Image Orientation
When determining the direction of fiber tract propagation, the frame of reference for the diffusion-encoding directions, MATLAB’s use of row/column indexing, and the image orientation must be considered. For example, consider an image dataset that uses an LPS frame of reference (the left, posterior, and superior directions are the +X, +Y, and +Z directions of the laboratory frame of reference) and oriented such that the top edge of the image is the anatomical anterior direction and the right edge of the image is the anatomical left direction. The first eigenvector of the diffusion tensor is:

ε<sub>1</sub> =  |ε<sub>X</sub> ε<sub>Y</sub> ε<sub>Z</sub>|<sup>T</sup>,
    
with the superscript T indicating vector transposition and the subscripts X, Y, and Z respectively indicating the X, Y, and Z components of ε<sub>1</sub>, specified within the LPS frame of reference. With this frame of reference and image orientation, the +X and +Y directions correspond to increasing column and row indices, respectively, of the image matrix.  Because MATLAB’s convention is that row and column values are specified first and second, respectively, when indexing a matrix, ε<sub>1</sub> must be converted to [row column slice] indexing as follows

ε<sub>1</sub>' =  |ε<sub>R</sub> ε<sub>C</sub> ε<sub>S</sub>|<sup>T</sup> =  |ε<sub>Y</sub> ε<sub>X</sub> ε<sub>Z</sub>|<sup>T</sup>
 
with the subscripts R, C, and S here reflecting the row, column, and slice directions. 

The commonly used medical image formats (DICOM, NIFTII, etc.) may use different frames of reference from each other, and some conversion processes rotate images from what was prescribed on the scanner. Therefore, several different frames of reference and image orientation can be specified.  All currently allowed combinations are listed in the table below.

<table>
  <tr>
    <th>Image Orientation</th>
    <th>Laboratory Frame of Reference</th>
     <th>Conversion to ε<sub>1</sub>'</th>
  </tr>
  <tr>
    <td>January</td>
    <td>$100</td>
     <td> <</td>
  </tr>
</table>

## 4. Syntax

[fiber_all, roi_flag, stop_list, fiber_len, fa_all, md_all] = fiber_track(tensor_m, mask, roi_mesh, ft_options, plot_options, anat_image);

The input arguments are:

* <i>tensor_m</i>: A 5D matrix with the first-third dimensions matching the [rows columns slice] size of the DTI images and the fourth and fifth dimensions holding the 33 diffusion tensor at each voxel

* <i>mask</i>: The mask delimiting the muscle to be fiber-tracked.    

* <i>roi_mesh</i>: The mesh reconstruction of the aponeurosis of muscle fiber insertion, output from define_roi.  

* <i>ft_options</i>: A structure containing the following fields:

  <i>ref_frame</i>: The frame of reference in which the diffusion directions are specified. For example, set ft_options.ref_frame='LPS'; if the left, posterior, and superior anatomical positions are (+).

  <i>image_orient</i>: The orientation of the images. Specify the anatomical directions at the top and right edges of the image as A (anterior) or P (posterior) and right (R) or left (L).  Input the result as a 2-element string variable (ft_options.image_orient='RA', 'AL', etc.).

  <i>mesh_dist</i>: The number of pixels to shift the mesh into the muscle, prior to fiber tracking. This can be a (+) or (-) number, depending on the desired direction of the shift.  

  <i>depth_ratio</i>: The ratio of slice thickness/in-plane resolution. Note that the function assumes equal in-plane voxel dimensions.

  <i>prop_alg</i>: A string variable that specifies the method for determining the direction of fiber tract propagation. The available options include 'euler', 'rk4', and 'fact'.  The available options include 1) <i>euler</i>: Diagonalization of the observed diffusion tensor at the current fiber tracking point, followed by magnitude-sorting and Euler integration of the first eigenvector. The user must specify the step size in the field ft_options.step_size. 2) <i>rk4</i>: Diagonalization of the observed diffusion tensor, diagonalization to find the firsteigenvector, and integration according to a 4th-order Runge-Kutta method. The user must specify the step size in the field ft_options.step_size. 3) <i>fact</i>: an implementation of the FACT algorithm.

  <i>step_size</i>: The Euler and 4th-order Runge-Kutta methods require the user to set the fiber-tracking step size, in pixels. A step size of 1 reflects the voxel width.

  <i>term_mthd</i>: A string variable that specifies the method for determining whether to terminate a fiber tract or not. Any fiber tracking point that falls outside of the image mask will terminate the tract. Other criteria using the inter-segment angle and the FA, processed according to either of two algorithms. 1) bin1: Inter- segment angle and FA are used as binary criteria to decide whether to terminate the tract. The angle used is the angle formed by two fiber tracking steps. The user can decide whether to calculate this angle between the current step and its immediate predecessor (1-back) or between the current step and a step that looks back by M points. Using the lookback option allows a tract to correct its direction following an initially aberrant result. When bin1 is used the tract terminates if either the angle or FA value is disallowed for a single point. 2) bin2: At each point, the inter-segment angle and FA data are treated as for bin1, but two consecutive points must fail in order to terminate the tract.
  
The FA and inter-segment angle criteria are set in fields called <i>ft_options.angle_thrsh</i> and <ift_options.fa_thrsh</i>. <i>angle_thrsh</i> is a two-element vector containing the angle threshold in degrees and the number of look-back steps. <i>fa_thrsh</i> is a two-element vector containing the lower and upper bounds of allowable FA values

The FACT algorithm uses its own method for tract termination. Thus, when the propagation algorithm is set to FACT, the user does not need to define term_mthd; however the user must create fields in ft_options called ft_options.r_crit and ft_options.num_fact_voxels. <i>r_crit</i> is a scalar quantity ranging from 0-1 that defines the allowable level of local variability in the direction of the first eigenvector. <i>num_fact_voxels</i> is the number of local neighboring voxels to include in the calculation of <i>r</i>. These are used to terminate tracts based on local variability in the first eigenvector. 
 
The following input arguments are optional and are required only if the user wishes to plot the fiber tracts:
* <i>plot_options</i>: If specified, this calls the <i>fiber_visualizer</i> function to plot the fiber, mask, and roi mesh.
 
* <i>anat_image</i>: The structural images.

The output arguments are:
* <i>fiber_all</i>: The fiber tracts, with units of pixels. The rows and columns correspond to locations on the roi_mesh. Dimension 3 gives point numbers on the tract, and the fourth dimension has row, column, and slice coordinates.

* <i>roi_flag</i>: A matrix indicating the presence of fibers that propagated at least 1 point

* <i>stop_list</i>: A matrix containing the reason for fiber tract termination (4=mask, 3=curvature, 2=FA, or R (the last for FACT only))

* <i>fiber_len</i>: The length, in points, of each fiber tract. 

* <i>fa_all</i>: The pointwise FA values on each fiber tract.

* <i>md_all</i>: The pointwise mean diffusivities along each tract

## 5. Example Code
### Example 1
Given:
1.	An anatomical image with variable name anat_image and having matrix size 192x192x44, field of view 192x192 mm, and slice thickness 7 mm and 

2.	DTI images having matrix size 19219244, field of view 192x192 mm, and slice thickness 7 mm and with the diffusion tensor data stored in a matrix called tensor_m;

3.	The muscle mask stored in a variable called mask and the aponeurosis mesh stored in a variable called roi_mesh;

4.	All imaging data acquired with a laboratory frame of reference having left/anterior/superior as the positive X, Y, and Z directions and an image orientation with the north and east sides being the anatomical right and anterior directions; 

5.	No intention to shift the aponeurosis mesh from its original location; and

6.	Tracking options of: 4th-order Runge-Kutta integration, 0.5 voxel step size, the BIN2 tract termination algorithm with stop criteria of angle, 25⁰ and a two-point look back and FA bounds of 0.1-0.4, 

the following code will allow the user to:

1.	Generate the fiber tracts and

2.	Visualize the outcome.

% Tracking options:
ft_options.ref_frame = ‘LAS'; %left-anterior-superior directions are +X, +Y, +Z

ft_options.image_orient = ‘RA'; %image north is right side, image east is anterior 

ft_options.mesh_dist = 0; %don’t shift the mesh

ft_options.prop_algo = ‘rk4'; %Runge-Kutta

ft_options.step_size = 0.5; %0.5 pixel width step

ft_options.term_mthd = ‘bin2'; %BIN2 stop algorithm

ft_options.angle_thrsh = [25 2]; %>=25 degree inter-segment angle disallowed; look back two points

ft_options.fa_thrsh = [.1 .4]; %0.1<FA<0.4 range of FA values allowed

ft_options.depth_ratio=7; %ratio of ST/in-plane resolution of reconstructed images

% Set visualization options
plot_options.anat_dims = [192 7]; %FOV and slice thickness of the images to be displayed, in mm

plot_options.anat_slices = 14:10:44; %display slices 14, 24, 34, and 44 

plot_options.plot_mesh = 1; %do plot an aponeurosis mesh

plot_options.plot_mask = 0; %don’t plot the mask

plot_options.plot_fibers = 1; %do plot any fiber tracts

plot_options.mesh_size = [192 192]; %rows x columns of the images used to generate the mesh

plot_options.mesh_dims = [192 7]; %FOV and ST of the images used to create the mesh

plot_options.mesh_color = [0.75 0.75 0.75]; %make the mesh light gray

plot_options.mask_size = [192 192]; %rows x columns of the images used to generate the mask

plot_options.mask_dims = [192 7]; %FOV and ST of the images used to create the mask

plot_options.mask_color = [1 0 0]; %make the mask a red, semi-transparent overlay

plot_options.fiber_color = [.8 .2 .2]; %make the fibers red

plot_options.dti_size = [192 192]; %rows x columns of the DTI data

plot_options.dti_dims = [192 7]; %FOV and ST of the DTI data
 
% Call the function:
[fiber_all, roi_flag, stop_list, fiber_len, fa_all, md_all] = fiber_track...
      (tensor_m, mask, roi_mesh, ft_options, plot_options, anat_image);

 
4.5.2.	Example 2
1.	As for Example 1, except using FACT with RCrit equal to 0.8 and examined over 20 voxels.
the following code will allow the user to 
1.	Generate the fiber tracts and
2.	Visualize the fiber tracts and the aponeurosis mesh.

% Tracking options:
ft_options.ref_frame = ‘LAS'; %left-anterior-superior directions are +X, +Y, +Z

ft_options.image_orient = ‘RA'; %image north is right side, image east is anterior 

ft_options.mesh_dist = 0; %don’t shift the mesh

ft_options.prop_algo = 'fact'; %use FACT

ft_options.num_fact_voxels = 20; %calculate R over 20 neighboring voxels

ft_options.r_crit = 0.8; %Rcrit = 0.8

ft_options.depth_ratio = 7; %ST/in-plane resolution

% Set visualization options
plot_options.anat_dims = [192 7]; %FOV and slice thickness of the images to be displayed, in mm

plot_options.anat_slices = 14:10:44; %display slices 14, 24, 34, and 44 

plot_options.plot_mesh = 1; %do plot an aponeurosis mesh

plot_options.plot_mask = 0; %don’t plot the mask

plot_options.plot_fibers = 1; %do plot fiber tracts

plot_options.mesh_size = [192 192]; %rows x columns of the images used to generate the mesh

plot_options.mesh_dims = [192 7]; %FOV and ST of the images used to create the mesh

plot_options.mesh_color = [0.75 0.75 0.75]; %make the mesh light gray

plot_options.mask_size = [192 192];  %rows x columns of the images used to generate the mask

plot_options.mask_dims = [192 7]; %FOV and ST of the images used to create the mask

plot_options.mask_color = [1 0 0]; %make the mask a red, semi-transparent overlay

plot_options.fiber_color = [.8 .2 .2]; %make the fibers red

plot_options.dti_size = [192 192]; %rows x columns of the DTI data

plot_options.dti_dims = [192 7]; %FOV and ST of the DTI data

% Call the function:
[fiber_all, roi_flag, stop_list, fiber_len, fa_all, md_all] = fiber_track ...
      (tensor_m, mask, roi_mesh, ft_options, plot_options, anat_image);

## 6. Acknowledgements

 People: Zhaohua Ding, Adam Anderson, Amanda Buck, Anneriet Heemskerk, and Justin Montenegro
 
 Grant support: NIH/NIAMS R01 AR050101, NIH/NIAMS R01 AR073831
