# Help for the function <i>fiber_visualizer</i>, v. 0.1.x

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
The user can call <i>fiber_visualizer</i> from the command line. The user must supply the anatomical images, a structure with some plotting options, and the other variables to be plotted as input arguments. In addition, <i>define_muscle</i>, <i>define_roi</i>, and <i>fiber_track</i> can be configured to call <i>fiber_visualizer</i> from within the functions, so that the mask, mesh, and fiber tracts can be automatically plotted.  Fields of view, matrix sizes, slice thickness, etc. are appropriately considered so that all structures are plotted using a consistent measurement scale.

[Back to the Top](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_visualizer.md)

## 3. Syntax 
fiber_figure = fiber_visualizer(anat_image, plot_options, roi_mesh, mask, fiber_all)

The input arguments are:

<i>anat_image</i>: An anatomical image having matrix size = Number of Image Rows (N<sub>R,I</sub>) x Number of Image Columns (N<sub>C,I</sub>) x Number of Image Slices (N<sub>S,I</sub>)

<i>plot_options</i>: A structure containing the following required fields:

  * <i>anat_dims</i>: A two-element vector containing the field of view (FOV) and the slice thickness (ST) of the anatomical images, in mm.

  * <i>anat_slices</i>: A vector containing the slice numbers of the anatomical images to be plotted.

   * <i>plot_mesh</i>: If set to 1, this field will allow plotting of the aponeurosis mesh. Otherwise, the mesh will not be plotted.

   * <i>plot_mask</i>: If set to 1, this field will allow plotting of the mask.  Otherwise, the mask will not be plotted.

   * <i>plot_fibers</i>:  If set to 1, this field will allow plotting of a single set of fiber tracts. If set to 2, this will allow plotting of two sets of fiber tracts. Otherwise, the fiber tracts will not be plotted.

Depending on the plot options selected, additional fields may be required in <i>plot_options</i>. If <i>plot_mesh</i> equals 1, the user must also specify:
   * <i>mesh_size</i>: This two-element vector specifies the in-plane matrix size of the images used to generate the mesh.

   * <i>mesh_dims</i>: This two-element vector specifies the FOV and ST of the images used to create the mesh, in mm.

   * <i>mesh_color</i>: The user creates a color scale using either of the following two options.
      * If mesh_color is a 3 element vector of values ranging from 0-1, the vector is interpreted as RGB levels.
      * If mesh_color is a matrix with size of (#mesh rows) x (#mesh columns) x 3, and if these values range from 0-1, the matrix will be interpreted as RGB levels specific to each tract. This could be used to represent the distribution of architectural parameters across the aponeurosis
   * <i>mesh_dist</i>: If the mesh was shifted for fiber tracking, the user should set this to the value used during fiber tracking.

If <i>plot_mask</i> equals 1, the user must also specify:
   * <i>mask_size</i>: This two-element vector specifies the in-plane matrix size of the images used to generate the mask.

   * <i>mask_dims</i>: This two-element vector specifies the FOV and slice thickness of the images used to create the mask.

   * <i>mask_color</i>: This three-element vector contains the RGB levels that determine the color of the mask.
     
If <i>plot_fibers</i> equals 1 or 2, you must also specify:
   * <i>dti_size</i>: This two-element vector specifies the matrix size of the images used for fiber tracking.

   * <i>dti_dim</i>s: This two-element vector specifies the FOV and slice thickness of the DT images. The FOV is assumed to be square.

   * <i>fiber_color</i>: Defines the color of the tracts. Several options are available:
      * If plot_fibers equals 1 and fiber_color is a 3 element vector of values ranging from 0-1, the vector is interpreted as RGB levels.
      * If plot_fibers equals 2 and fiber_color is a 2x3 matrix of values ranging from 0-1, each row of the matrix is interpreted as RGB
       levels for the respective sets of tracts.
      * If plot_fibers equals 1 and fiber_color is a matrix with size of (#mesh rows) x (#mesh columns) x 3, and if these values range from 0-1, the third dimension of the matrix will be interpreted as RGB levels specific to each tract
      * If plot_fibers equals 2 and fiber_color is a matrix with size of (#mesh rows) x (#mesh columns) x 3 x 2, and if these values range from 0-1, the third dimension of the matrix will be interpreted as RGB levels specific to each tract, separately for sets 1 and 2

   * <i>fiber_skip</i>: Setting fiber_skip to integer values > 1 will skip over fiber tracts when plotting. This may improve visualization and will decrease time for rendering. If not specified, all fibers will be plotted.

 <i>roi_mesh</i>: The output of <i>define_roi</i>. It is only needed if plot_options.plot_mesh is set to 1.

 <i>mask</i>: A binary mask around the muscle of interest. It could be the output of <i>define_muscle</i> or it could have been defined in another program. It is only needed if plot_options.plot_mask is set to 1.

 <i>fiber_all</i>: The output of <i>fiber_track</i> (original fiber tracts) or <i>fiber_smoother</i> (smoothed fiber tracts), or <i>fiber_goodness</i> (quality-selected fiber tracts). It is only needed if plot_options.plot_fibers is set to 1. If plot_fibers equals 1, the size should be (#mesh rows) x (#mesh columns) x (#fiber tract points) x 3. If plot_fibers equals 1, the size should be (#mesh rows) x (#mesh columns) x (#fiber tract points) x 3 x 2.

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

plot_options.anat_dims = [192 7]; %FOV and slice thickness of the images to be displayed, in mm

plot_options.anat_slices = 14:10:44; %display slices 14, 24, 34, and 44 

plot_options.plot_mesh = 0; %don’t plot an aponeurosis mesh

plot_options.plot_mask = 1; %do plot the mask

plot_options.plot_fibers = 0; %don’t plot any fiber tracts

plot_options.mask_size = [192 192]; %rows x columns of the images used to generate the mask

plot_options.mask_dims = [192 7]; %FOV and ST of the images used to create the mask, in mm

plot_options.mask_color = [1 0 0]; %make the mask a red, semi-transparent overlay

% Call the function:

mask_figure = fiber_visualizer(anat_image, plot_options, [], mask, []);
 
### Example 2
Given:

1.	An anatomical image with variable name anat_image and having matrix size 192x192x44, field of view 192x192 mm, and slice thickness 7 mm;

2.	The aponeurosis mesh stored in a variable called <i>roi_mesh</i>;  

the following code will allow the user to:

1.	Visualize the mesh; and

2.	Return a MATLAB figure structure called <i>mesh_figure</i>.
 
 
% Set visualization options
plot_options.anat_dims = [192 7]; %FOV and slice thickness of the images to be displayed, in mm

plot_options.anat_slices = 14:10:44; %display slices 14, 24, 34, and 44 

plot_options.anat_dims = [192 7]; %FOV and slice thickness of the images to be displayed, in mm

plot_options.anat_slices = 14:10:44; %display slices 14, 24, 34, and 44 

plot_options.plot_mesh = 1; %do plot an aponeurosis mesh

plot_options.plot_mask = 0; %don’t plot the mask

plot_options.plot_fibers = 0; %don’t plot any fiber tracts

plot_options.mesh_size = [192 192]; %rows x columns of the images used to generate the mesh

plot_options.mesh_dims = [192 7]; %FOV and ST of the images used to create the mesh, in mm

plot_options.mesh_color = [0.75 0.75 0.75]; %make the mesh light gray

% Call the function:

mesh_figure = fiber_visualizer(anat_image, plot_options, roi_mesh, [], []);
 
 
###	Example 3

Given:

1.	An anatomical image with variable name anat_image and having matrix size 192x192x44, field of view 192x192 mm, and slice thickness 7 mm;

2.	The aponeurosis mesh stored in a variable called <i>roi_mesh</i>; and 

3.	The fiber tracts stored in a variable called <i>fiber_all</i>.

the following code will allow the user to:

1.	Visualize the mesh; 

2.	Visualize the fibers; and

3.	Return a MATLAB figure structure called <i>fiber_mesh_figure</i>.

 
% Set visualization options

plot_options.anat_dims = [192 7]; %FOV and slice thickness of the images to be displayed, in mm

plot_options.anat_slices = 14:10:44; %display slices 14, 24, 34, and 44 

plot_options.plot_mesh = 1; %don’t plot an aponeurosis mesh

plot_options.plot_mask = 0; %do plot the mask

plot_options.plot_fibers = 1; %do plot the fiber tracts

plot_options.mesh_size = [192 192]; %rows x columns of the images used to generate the mesh

plot_options.mesh_dims = [192 7]; %FOV and ST of the images used to create the mesh

plot_options.mesh_color = [0.75 0.75 0.75]; %make the mesh light gray

plot_options.fiber_color = [.8 .2 .2];

plot_options.dti_size = [192 192];

plot_options.dti_dims = [192 7];

% Call the function:

fiber_mesh_figure = fiber_visualizer(anat_image, plot_options, roi_mesh, [], fiber_all);

[Back to the Top](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_visualizer.md)

## 5. Acknowledgements
People: Zhaohua Ding, Hannah Kilpatrick

Grant support: NIH/NIAMS R01 AR073831

[Back to the Top](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_visualizer.md)

