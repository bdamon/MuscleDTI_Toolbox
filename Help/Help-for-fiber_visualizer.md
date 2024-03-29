# Help for the function [<i>fiber_visualizer</i>](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Tractography-Functions/fiber_visualizer.m), v. 1.0.0

## Introduction

This help file contains information about
1) [Purpose of the Program](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_visualizer.md#1-Purpose)
2) [Usage of the Program](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_visualizer.md#2-Usage)
3) [Syntax](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_visualizer.md#3-Syntax)
5) [Example Code](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_visualizer.md#4-Example-Code)
5) [Acknowledgements](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_visualizer.md#5-Acknowledgements)

## 1. Purpose
<i>fiber_ visualizer</i> is used to visualize anatomical images and other structures, including the muscle mask, aponeurosis mesh, and/or the fiber tracts.

[Back to the Top](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_visualizer.md)

## 2. Usage
The user can call <i>fiber_visualizer</i> from the command line. In addition, [<i>define_muscle</i>](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-define_muscle.md), [<i>define_roi</i>](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-define_muscle.md), and [<i>fiber_track</i>](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-define_muscle.md) can be configured to call <i>fiber_visualizer</i> from within the functions, so that the mask, mesh, and fiber tracts can be automatically plotted. The user must supply the anatomical images, a structure with some plotting options, and the other variables to be plotted as input arguments. Fields of view, matrix sizes, slice thickness, etc. are appropriately considered so that all structures are plotted using a consistent measurement scale.

[Back to the Top](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_visualizer.md)

## 3. Syntax 
fiber_figure = fiber_visualizer(anat_image, fv_options, roi_mesh, mask, fiber_all)

The input arguments are:

<i>anat_image</i>: An anatomical image having matrix size = Number of Anatomical Image Rows (N<sub>R,A</sub>) x Number of Anatomical Image Columns (N<sub>C,A</sub>) x Number of Anatomical Image Slices (N<sub>S,A</sub>)

<i>fv_options</i>: A structure containing the following required fields:

  * <i>anat_dims</i>: A two-element vector containing the field of view (FOV) and the slice thickness (∆Z) of the anatomical images, in mm.

  * <i>anat_slices</i>: A vector containing the slice numbers of the anatomical images to be plotted.

   * <i>plot_mesh</i>: If set to 1, the aponeurosis mesh will be plotted.  The default setting is 0.

   * <i>plot_mask</i>: If set to 1, the mask will be plotted.  The default setting is 0.

   * <i>plot_fibers</i>: If set to 1, a single set of fiber tracts will be plotted. If set to 2, two sets of fiber tracts will be plotted.  The default setting is 0.

Depending on the plot options selected, additional fields may be required in <i>fv_options</i>. If <i>plot_mesh</i> equals 1, the user must also specify:
   * <i>mesh_size</i>: This two-element vector specifies the in-plane matrix size of the images used to generate the mesh.

   * <i>mesh_dims</i>: This two-element vector specifies the FOV and ∆Z of the images used to create the mesh, in mm.

   * <i>mesh_color</i>: The user creates a color scale using either of the following two options.
      * If mesh_color is a 3 element vector of values ranging from 0-1, the vector is interpreted as RGB levels.
      * If mesh_color is a matrix with size of (#mesh rows) x (#mesh columns) x 3, and if these values range from 0-1, the matrix will be interpreted as RGB levels specific to each tract. This could be used to represent the distribution of architectural parameters across the aponeurosis
   * <i>mesh_dist</i>: If the mesh was shifted for fiber tracking, the user should set this to the value used during fiber tracking.

If <i>plot_mask</i> equals 1, the user must also specify:
   * <i>mask_size</i>: This two-element vector specifies the in-plane matrix size of the images used to generate the mask.

   * <i>mask_dims</i>: This two-element vector specifies the FOV and ST of the images used to create the mask.

   * <i>mask_color</i>: This three-element vector contains the RGB levels that determine the color of the mask.
     
If <i>plot_fibers</i> equals 1 or 2, you must also specify:
   * <i>dti_size</i>: This two-element vector specifies the matrix size of the images used for fiber tracking.

   * <i>dti_dim</i>s: This two-element vector specifies the FOV and ST of the DT images. The FOV is assumed to be square.

   * <i>fiber_color</i>: Defines the color of the tracts. Several options are available:
      * If plot_fibers equals 1 and fiber_color is a 3 element vector of values ranging from 0-1, the vector is interpreted as RGB levels.
      * If plot_fibers equals 2 and fiber_color is a 2x3 matrix of values ranging from 0-1, each row of the matrix is interpreted as RGB
       levels for the respective sets of tracts.
      * If plot_fibers equals 1 and fiber_color is a matrix with size of (#mesh rows) x (#mesh columns) x 3, and if these values range from 0-1, the third dimension of the matrix will be interpreted as RGB levels specific to each tract. This could be used to represent fiber-specific anatomical properties using a color scale.
      * If plot_fibers equals 2 and fiber_color is a matrix with size of (#mesh rows) x (#mesh columns) x 3 x 2, and if these values range from 0-1, the third dimension of the matrix will be interpreted as RGB levels specific to each tract, separately for sets 1 and 2

   * <i>fiber_skip</i>: Setting fiber_skip to integer values > 1 will skip over fiber tracts when plotting. This may improve visualization and will decrease time for rendering. If not specified, all fibers will be plotted.

 <i>roi_mesh</i>: The output of [<i>define_roi</i>](https://github.com/bdamon/MuscleDTI_Toolbox/edit/master/Help/Help-for-define_roi.md). It is only needed if fv_options.plot_mesh is set to 1.

 <i>mask</i>: A binary mask around the muscle of interest. It could be the output of [<i>define_muscle</i>](https://github.com/bdamon/MuscleDTI_Toolbox/edit/master/Help/Help-for-define_muscle.md) or it could have been defined in another program. It is only needed if fv_options.plot_mask is set to 1.

 <i>fiber_all</i>: The output of [<i>fiber_track</i>](https://github.com/bdamon/MuscleDTI_Toolbox/edit/master/Help/Help-for-fiber_track.md) (original fiber tracts) or [<i>fiber_smoother</i>](https://github.com/bdamon/MuscleDTI_Toolbox/edit/master/Help/Help-for-fiber_smoother.md) (smoothed fiber tracts), or [<i>fiber_goodness</i>]((https://github.com/bdamon/MuscleDTI_Toolbox/edit/master/Help/Help-for-fiber_goodness.md)) (quality-selected fiber tracts). It is only needed if fv_options.plot_fibers is set to 1. If plot_fibers equals 1, the size should be (#mesh rows) x (#mesh columns) x (#fiber tract points) x 3. If plot_fibers equals 1, the size should be (#mesh rows) x (#mesh columns) x (#fiber tract points) x 3 x 2.

The output arguments are:
 <i>fiber_figure</i>: A Matlab figure structure

[Back to the Top](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_visualizer.md)

## 4. Example Code

### Example 1

Given:

1.	An anatomical image with variable name anat_image and having matrix size 192x192x44, field of view 192x192 mm, and slice thickness 7 mm;

2.	The muscle mask stored in a variable called <i>mask</i>;  

the following code will allow the user to 

1.	Visualize the mask; and

2.	Return a MATLAB figure structure called <i>mask_figure</i>.
 
 
% Set visualization options

fv_options.anat_dims = [192 7]; %FOV and slice thickness of the images to be displayed, in mm

fv_options.anat_slices = 14:10:44; %display slices 14, 24, 34, and 44 

fv_options.plot_mesh = 0; %don’t plot an aponeurosis mesh

fv_options.plot_mask = 1; %do plot the mask

fv_options.plot_fibers = 0; %don’t plot any fiber tracts

fv_options.mask_size = [192 192]; %rows x columns of the images used to generate the mask

fv_options.mask_dims = [192 7]; %FOV and ST of the images used to create the mask, in mm

fv_options.mask_color = [1 0 0]; %make the mask a red, semi-transparent overlay

% Call the function:

mask_figure = fiber_visualizer(anat_image, fv_options, [], mask, []);
 
### Example 2
Given:

1.	An anatomical image with variable name anat_image and having matrix size 192x192x44, field of view 192x192 mm, and slice thickness 7 mm;

2.	The aponeurosis mesh stored in a variable called <i>roi_mesh</i>;  

the following code will allow the user to:

1.	Visualize the mesh; and

2.	Return a MATLAB figure structure called <i>mesh_figure</i>.
 
 
% Set visualization options
fv_options.anat_dims = [192 7]; %FOV and slice thickness of the images to be displayed, in mm

fv_options.anat_slices = 14:10:44; %display slices 14, 24, 34, and 44 

fv_options.anat_dims = [192 7]; %FOV and slice thickness of the images to be displayed, in mm

fv_options.anat_slices = 14:10:44; %display slices 14, 24, 34, and 44 

fv_options.plot_mesh = 1; %do plot an aponeurosis mesh

fv_options.plot_mask = 0; %don’t plot the mask

fv_options.plot_fibers = 0; %don’t plot any fiber tracts

fv_options.mesh_size = [192 192]; %rows x columns of the images used to generate the mesh

fv_options.mesh_dims = [192 7]; %FOV and ST of the images used to create the mesh, in mm

fv_options.mesh_color = [0.75 0.75 0.75]; %make the mesh light gray

% Call the function:

mesh_figure = fiber_visualizer(anat_image, fv_options, roi_mesh, [], []);
 
 
###	Example 3
Given:
1.	An anatomical image with variable name anat_image and having matrix size 192x192x44, field of view 192x192 mm, and slice thickness 7 mm;

2.	The aponeurosis mesh stored in a variable called roi_mesh;  

the following code will allow the user to:

1.	Visualize the mesh; 

2.	Adjust the color scale to represent the pennation angle; and

3.	Return a MATLAB figure structure called mesh_figure.

% Set visualization options

fv_options.anat_dims = [192 7]; %FOV and slice thickness of the images to be displayed, in mm


fv_options.anat_slices = 14:10:44; %display slices 14, 24, 34, and 44 

fv_options.anat_dims = [192 7]; %FOV and slice thickness of the images to be displayed, in mm

fv_options.anat_slices = 14:10:44; %display slices 14, 24, 34, and 44 

fv_options.plot_mesh = 1; %do plot an aponeurosis mesh

fv_options.plot_mask = 0; %don’t plot the mask

fv_options.plot_fibers = 0; %don’t plot any fiber tracts

fv_options.mesh_size = [192 192]; %rows x columns of the images used to generate the mesh

fv_options.mesh_dims = [192 7]; %FOV and dZ of the images used to create the mesh, in mm

fv_options.mesh_color = [0.75 0.75 0.75]; %make the mesh light gray

fv_options.plot_fibers = 0; %don’t plot any fiber tracts

fv_options.mesh_size = [192 192]; %rows x columns of the images used to generate the mesh

fv_options.mesh_dims = [192 7]; %FOV and dZ of the images used to create the mesh, in mm

%show increasing pennation angles as decreasing redness/increasing blueness

fv_options.mesh_color = zeros(size(roi_mesh, 1), size(roi_mesh, 2), 3); 

fv_options.mesh_color(:,:,1) = 1-mean_fiber_props(:,:,2)/max(max(mean_fiber_props(:,:,2)));

fv_options.mesh_color(:,:,3) = mean_fiber_props(:,:,2)/max(max(mean_fiber_props(:,:,2))); 

%get rid of NaN and Inf values:

fv_options.mesh_color(isnan(fv_options.mesh_color))=0;

fv_options.mesh_color(isinf(fv_options.mesh_color))=0;

% Call the function:

mesh_figure = fiber_visualizer(anat_image, fv_options, roi_mesh, [], []);


### Example 4
Given:

1.	An anatomical image with variable name anat_image and having matrix size 192x192x44, field of view 192x192 mm, and slice thickness 7 mm;

2.	The aponeurosis mesh stored in a variable called <i>roi_mesh</i>; and 

3.	The fiber tracts stored in a variable called <i>fiber_all</i>.

the following code will allow the user to:

1.	Visualize the mesh; 

2.	Visualize the fibers; and

3.	Return a MATLAB figure structure called <i>fiber_mesh_figure</i>.

 
% Set visualization options

fv_options.anat_dims = [192 7]; %FOV and slice thickness of the images to be displayed, in mm

fv_options.anat_slices = 14:10:44; %display slices 14, 24, 34, and 44 

fv_options.plot_mesh = 1; %do plot an aponeurosis mesh

fv_options.plot_mask = 0; %do plot the mask

fv_options.plot_fibers = 1; %do plot the fiber tracts

fv_options.mesh_size = [192 192]; %rows x columns of the images used to generate the mesh

fv_options.mesh_dims = [192 7]; %FOV and ST of the images used to create the mesh

fv_options.mesh_color = [0.75 0.75 0.75]; %make the mesh light gray

fv_options.fiber_color = [.8 .2 .2];

fv_options.dti_size = [192 192];

fv_options.dti_dims = [192 7];

% Call the function:

fiber_mesh_figure = fiber_visualizer(anat_image, fv_options, roi_mesh, [], fiber_all);

[Back to the Top](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_visualizer.md)

## 5. Acknowledgements
People: Zhaohua Ding, Hannah Kilpatrick

Grant support: NIH/NIAMS R01 AR073831

[Back to the Top](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_visualizer.md)
 
