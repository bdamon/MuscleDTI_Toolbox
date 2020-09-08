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
    <td>AL</td>
    <td>LPS</td>
     <td>|ε<sub>R</sub> ε<sub>C</sub> ε<sub>S</sub>|<sup>T</sup> =  |ε<sub>Y</sub> ε<sub>X</sub> ε<sub>Z</sub>|<sup>T</sup></td>
  </tr>
  <tr>
     <td>AL</td>
     <td>LAS</td>
     <td>|ε<sub>R</sub> ε<sub>C</sub> ε<sub>S</sub>|<sup>T</sup> =  |-ε<sub>Y</sub> ε<sub>X</sub> ε<sub>Z</sub>|</sup></td>
   </tr>
   <tr>
     <td>AL</td>
     <td>RAS</td>
     <td>|ε<sub>R</sub> ε<sub>C</sub> ε<sub>S</sub>|<sup>T</sup> =  |-ε<sub>Y</sub> -ε<sub>X</sub> ε<sub>Z</sub>|</sup></td>
   </tr>
   <tr>
    <td>RA</td>
    <td>LPS</td>
     <td>|ε<sub>R</sub> ε<sub>C</sub> ε<sub>S</sub>|<sup>T</sup> =  |ε<sub>X</sub> -ε<sub>Y</sub> ε<sub>Z</sub>|<sup>T</sup></td>
  </tr>
  <tr>
     <td>RA</td>
     <td>LAS</td>
     <td>|ε<sub>R</sub> ε<sub>C</sub> ε<sub>S</sub>|<sup>T</sup> =  |ε<sub>X</sub> ε<sub>Y</sub> ε<sub>Z</sub>|</sup></td>
   </tr>
   <tr>
     <td>RA</td>
     <td>RAS</td>
     <td>|ε<sub>R</sub> ε<sub>C</sub> ε<sub>S</sub>|<sup>T</sup> =  |-ε<sub>X</sub> ε<sub>Y</sub> ε<sub>Z</sub>|</sup></td>
   </tr>
</table>

### Fiber Tract Propagation and Termination
Tracts are initiated for every [row column] coordinate on the aponeurosis mesh using two nested for loops. Inside the inner loop, the seed point is determined; its location within the muscle mask is verified; the initial fiber-tracking step is found; a while loop is used to initiate and propagate a tract according to the selected algorithm until stop criteria are met; and the points are recorded in a matrix called fiber_all. fiber_all has dimensions of N<sub>R</sub>xN<sub>C,A</sub>xN<sub>P,Max</sub>x3, where N<sub>P,Max</sub> is the maximum number of points in any fiber tract; the fourth dimension holds the [row column slice] coordinates of each fiber-tracking point

The steps are:

Step 1 – Locate the seed point on the aponeurosis mesh: for each location on the aponeurosis mesh, the [row column slice] coordinates are used to form the seed point, P<sub>1</sub>.

Step 2 – Record the seed point in the fiber tract matrix or continue to the next location on the mesh: The location of P<sub>1</sub> inside the muscle mask is verified. If P<sub>1</sub> falls outside of the muscle mask, the tract does not propagate, a value of 4 is recorded at the [row column] location in the variable <i>stop_list</i>, and the loops continue to the next seed point. If P<sub>1</sub> falls within the mask, the seed point is added to fiber_all. P<sub>1</sub> is stored at index N<sub>P</sub> = 1, and the function proceeds to Step 3.

Step 3: Determine the initial fiber-tracking step, ∆S: Three options for determining ∆S are available: Euler integration, 4th-order Runge-Kutta integration, and FACT.  Although the present discussion is about the initial fiber-tracking step, to maintain generality we refer to the current fiber-tracking point as P<sub>n</sub>, the current first eigenvector as ε<sub>1,n)</sub>', and the current fiber-tracking step as ∆S<sub>n</sub>. The three tract propagation methods are implemented as follows:
* Euler: The initial direction of tract propagation is determined by rounding the [row column slice] coordinates of P_n, using them as indices into tensor_m, and retrieving D. D is diagonalized using the eig function; the eigenvalues are magnitude-sorted and ε_1 is identified. If ε<sub>1,Z</sub><0, ε<sub>1</sub> is multiplied by -1 so that tracts always propagate in ascending slice order. Then ε<sub>1</sub> is converted to ε<sub>1</sub>' as described above. Further, the slice-wise component of ε<sub>1</sub>', ε<sub>1,S</sub>', is adjusted to account for the aspect ratio of the voxel (the ratio of its thickness to its width) to determine a step direction:

   ε<sub>1,S'</sub>' = ε<sub>1,S</sub>' / (ΔZ/ΔX)
   
   where ΔZ/ΔX is the ratio of the slice thickness to the in-plane resolution. Then, the fiber-tracking step ΔS<sub>n,</sub> is calculated as 
   
   ΔS<sub>n</sub> = h * ε<sub>1</sub>'
   
   where h is the step size, expressed as a fraction of the voxel width, ∆X; it typically has values of 0.5-1.0 voxel width. As described below, h is set by the user.

* 4th-order Runge-Kutta Integration: 4th-order Runge-Kutta integration: Initially, ε<sub>1,n</sub>' is calculated as described above for Euler integration. Consistent with the 4th-order Runge-Kutta method, this step is modified as follows:
   
   ∆S<sub>n,0</sub> = h/2 * ε<sub>1,n</sub>

   with ∆S analogous to k in the Runge-Kutta method. A temporary second point, P<sub>n+1</sub>', is calculated as:
   
   P<sub>n+1</sub>'= P<sub>n</sub> + ∆S<sub>n,0</sub>

   The [row column slice] coordinates of P<sub>n+1</sub>' are rounded and used to look up D from the tensor_m matrix at P<sub>n+1</sub>'. D is diagonalized, the eigenvalues are magnitude-sorted, and ε<sub>1</sub>' and ∆S<sub>n</sub> are found at P<sub>n+1</sub>' ((ε<sub>1,n+1</sub>' and ∆S<sub>n+1</sub>, respectively). The step from P<sub>n+1</sub>' is calculated as:

   ∆S<sub>n+1</sub> = h/2 * ε<sub>1,n+1</sub>'

   A temporary third point is calculated as:

   P<sub>n+2</sub>'= P<sub>n+1</sub> + ∆S<sub>n+1</sub>

   The procedure is repeated for the third and fourth points, with:

   ∆S<sub>n+2</sub> = h * ε<sub>1,n+2</sub>'  and P<sub>n+3</sub>'= P<sub>n+2</sub> + ∆S<sub>n+2</sub>

   ε<sub>1,n+3</sub>' is then found. 
   
   The eigenvectors describing the muscle fiber orientations cannot be averaged directly.  Instead, for each point, the corresponding eigenvector is used to form a dyadic tensor as

   E<sub>D,n</sub> =  ε<sub>1,n</sub>' * ε<sub>1,n</sub><sup>T</sup>, 

   etc., where the superscript T indicates vector transposition.  Finally, the dyadic tensors are averaged using:

   1/6 ( E<sub>D,n0</sub> + 2 * E<sub>D,n+1</sub> + 2 * E<sub>D,n+ 2</sub> + E<sub>D,n+3</sub>)

   The mean dyadic tensor is diagonalized and its principal eigenvector, ε<sub>1</sub>', is found. This is used to calculate ∆S as described above.

* FACT: FACT follows the direction of the first eigenvector while the tract remains within a given voxel, but allows the tract direction to change as soon as the tract enters a new voxel.  To implement FACT, ε<sub>1</sub>' is initially determined as described for Euler integration. The tract directions in X, Y, and Z are calculated for values of h= {0.01,0.02,0.03,… 1} times the voxel width. The smallest value resulting in a new row, column, or slice coordinate is taken as the step size h, and ∆S is calculated as described above.

Regardless of the propagation method, the next point, P<sub>n+1</sub>, is calculated as P<sub>n+1</sub> =P<sub>n</sub> +∆S<sub>n</sub>. If it falls within the muscle mask, it is added to the <i>fiber_all</i> matrix

Step 4 – Propagate the tract: From point P<sub>n+1</sub>, the step ∆S<sub>n+1</sub> and the next fiber tract point P<sub>n+2</sub> are calculated according to the selected propagation algorithm, as described above. For Runge-Kutta integration, if the 2nd, 3rd, or 4th order points fall outside of the muscle mask, the propagation algorithm automatically changes to Euler integration until the tract is terminated. This is to avoid obtaining erroneous estimates of the fiber orientation from voxels outside of the muscle of interest.

Before being added to the <i>fiber_all</i> matrix, several termination criteria are applied. First, the location of P<sub>n+1</sub> within the muscle mask is verified; if not, tract propagation stops and a value of 4 is written into the <i>stop_list</i> matrix.  In addition, either of several algorithms may be applied and used to terminate tract propagation.
* BIN1: Two binary criteria are applied. First, the FA must fall within the bounds set by the user in <i>ft_options</i>. Also, the angle formed by the current and a previous fiber tracking step is calculated as:

ψ=cos<sup>-1</sup> (ε<sub>1,n</sup> ∙ ε<sub>1,n-p</sub>)

where p is the number of steps over which to look back and ∙ indicates the vector dot product. ψ must be smaller than the value specified by the user. When BIN1 is used, the tract terminates if either ψ or the FA value is disallowed for a single point. If the tract stops because the FA criterion failed, a value of 2 is written into the <i>stop_list</i> matrix; if the tract stops because the angle criterion failed, a value of 3 is written into the <i>stop_list</i> matrix.

* BIN2: The FA and ψ criteria are applied as described above, except that the ψ value must have been disallowed for two successive points. 

* FACT: The FACT algorithm for terminating fiber tract propagation was described by Mori et al. Briefly, the inner products are calculated between ε<sub>1</sub> in the current voxel and each of its 26 neighbors.  They are magnitude-sorted in descending order and then averaged over a user-specified number of voxels, giving the parameter R. If R is lower than a user-specified critical value R<sub>Crit</sub>, this is interpreted as excessive local heterogeneity in fiber orientation and the tract is terminated. The value of R is written into <i>stop_list<i>.
   
The variable <i>stop_list</i> is useful to diagnose the reasons of tract propagation failure and to optimize the stop criteria. 

Step 5 – Add the point and continue tracking: If all criteria are successfully met, the next point is calculated and added to <i>fiber_all</i>. The fiber counter is incremented and the tract is recorded as a successful tracking result in the variable <i>roi_flag</i>. Steps 4 and 5 occur within the <i>while</i> loop and continue until a stop criterion is met.  At that point, the <i>while</i> loop breaks and the programs advances to the new [row column] coordinate on the aponeurosis mesh

## 4. Syntax

[fiber_all, roi_flag, stop_list, fiber_len, fa_all, md_all] = fiber_track(tensor_m, mask, roi_mesh, ft_options, plot_options, anat_image);

The input arguments are:

* <i>tensor_m</i>: A 5D matrix with the first-third dimensions matching the [rows columns slice] size of the DTI images and the fourth and fifth dimensions holding the 33 diffusion tensor at each voxel

* <i>mask</i>: The mask delimiting the muscle to be fiber-tracked.    

* <i>roi_mesh</i>: The mesh reconstruction of the aponeurosis of muscle fiber insertion, output from <i>define_roi</i>.  

* <i>ft_options</i>: A structure containing the following fields:

  <i>ref_frame</i>: The frame of reference in which the diffusion directions are specified. For example, set ft_options.ref_frame='LPS'; if the left, posterior, and superior anatomical positions are (+).

  <i>image_orient</i>: The orientation of the images. Specify the anatomical directions at the top and right edges of the image as A (anterior) or P (posterior) and right (R) or left (L).  Input the result as a 2-element string variable (ft_options.image_orient='RA', 'AL', etc.).

  <i>mesh_dist</i>: The number of pixels to shift the mesh into the muscle, prior to fiber tracking. This can be a (+) or (-) number, depending on the desired direction of the shift.  

  <i>depth_ratio</i>: The ratio of slice thickness/in-plane resolution. Note that the function assumes equal in-plane voxel dimensions.

  <i>prop_alg</i>: A string variable that specifies the method for determining the direction of fiber tract propagation. The available options include 'euler', 'rk4', and 'fact'.  The available options include 1) <i>euler</i>: Diagonalization of the observed diffusion tensor at the current fiber tracking point, followed by magnitude-sorting and Euler integration of the first eigenvector. The user must specify the step size in the field ft_options.step_size. 2) <i>rk4</i>: Diagonalization of the observed diffusion tensor, diagonalization to find the firsteigenvector, and integration according to a 4th-order Runge-Kutta method. The user must specify the step size in the field ft_options.step_size. 3) <i>fact</i>: an implementation of the FACT algorithm.

  <i>step_size</i>: The Euler and 4th-order Runge-Kutta methods require the user to set the fiber-tracking step size, in pixels. A step size of 1 reflects the voxel width.

  <i>term_mthd</i>: A string variable that specifies the method for determining whether to terminate a fiber tract or not. Any fiber tracking point that falls outside of the image mask will terminate the tract. Other criteria using the inter-segment angle and the FA, processed according to either of two algorithms. 1) bin1: Inter- segment angle and FA are used as binary criteria to decide whether to terminate the tract. The angle used is the angle formed by two fiber tracking steps. The user can decide whether to calculate this angle between the current step and its immediate predecessor (1-back) or between the current step and a step that looks back by M points. Using the lookback option allows a tract to correct its direction following an initially aberrant result. When bin1 is used the tract terminates if either the angle or FA value is disallowed for a single point. 2) bin2: At each point, the inter-segment angle and FA data are treated as for bin1, but two consecutive points must fail in order to terminate the tract.
  
  The FA and inter-segment angle criteria are set in fields called <i>ft_options.angle_thrsh</i> and <i>ft_options.fa_thrsh</i>. <i>angle_thrsh</i> is a two-element vector containing the angle threshold in degrees and the number of look-back steps. <i>fa_thrsh</i> is a two-element vector containing the lower and upper bounds of allowable FA values

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

2.	DTI images having matrix size 192x192x44, field of view 192x192 mm, and slice thickness 7 mm and with the diffusion tensor data stored in a matrix called tensor_m;

3.	The muscle mask stored in a variable called mask and the aponeurosis mesh stored in a variable called roi_mesh;

4.	All imaging data acquired with a laboratory frame of reference having left/anterior/superior as the positive X, Y, and Z directions and an image orientation with the top and right edges being the anatomical right and anterior directions; 

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

[fiber_all, roi_flag, stop_list, fiber_len, fa_all, md_all] = fiber_track(tensor_m, mask, roi_mesh, ft_options, plot_options, anat_image);

 
###	Example 2

1.	As for Example 1, except using FACT with RCrit equal to 0.8 and examined over 20 voxels.

the following code will allow the user to 

1.	Generate the fiber tracts and

2.	Visualize the fiber tracts and the aponeurosis mesh.

% Tracking options:

ft_options.ref_frame = ‘LAS'; %left-anterior-superior directions are +X, +Y, +Z

ft_options.image_orient = ‘RA'; %image north is right side, image east is anterior 

ft_options.mesh_dist = 0; %don’t shift the mesh

ft_options.prop_algo = 'fact'; %use FACT

ft_options.num_fact_voxels = 5; %calculate R over 5 neighboring voxels

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

[fiber_all, roi_flag, stop_list, fiber_len, fa_all, md_all] = fiber_track(tensor_m, mask, roi_mesh, ft_options, plot_options, anat_image);

## 6. Acknowledgements

 People: Zhaohua Ding, Adam Anderson, Amanda Buck, Anneriet Heemskerk, and Justin Montenegro
 
 Grant support: NIH/NIAMS R01 AR050101, NIH/NIAMS R01 AR073831
