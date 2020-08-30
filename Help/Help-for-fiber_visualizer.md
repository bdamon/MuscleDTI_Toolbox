# Help for the function <i>fiber_visualizer</i>, v. 0.1.x

## Introduction

This help file contains information about
1) [Purpose of the Program](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_visualizer.md#1-Purpose)
2) [Usage of the program](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_visualizer.md#2-Usage)
3) [Syntax](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_visualizer.md#3-Syntax)
5) [Example Code](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_visualizer.md#4-Example-Code)
5) [Acknowledgements](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_visualizer.md#5-Acknowledgements)

## 1. Purpose
<i>fiber_ visualizer</i> is used to visualize anatomical images and other structures, including the muscle mask, aponeurosis mesh, and/or the fiber tracts.

## 2. Usage
The user can call <i>fiber_ visualizer</i> from the command line using the [syntax below](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_visualizer.md#3-Syntax). In addition, <i>define_muscle</i>, <i>define_roi</i>, and <i>fiber_track</i> can be configured to call <i>fiber_visualizer</i> from within the functions, so that the mask, mesh, and fiber tracts can be automatically plotted.  The user must supply the anatomical images, a structure with some plotting options, and the other variables to be plotted as input arguments. Fields of view, matrix sizes, slice thickness, etc. are appropriately considered so that all structures are plotted using a consistent measurement scale.

## 3. Syntax 
fiber_figure = fiber_visualizer(anat_image, plot_options, roi_mesh, mask, fiber_all)

The input arguments are:
<i>anat_image</i>: An anatomical image having matrix size = Number of Image Rows (N<sub>R,I</sub>)xNumber of Image Columns (N<sub>C,I</sub>)xNumber of Image Slices (N<sub>S,I</sub>)

<i>plot_options</i>: A structure containing the following required fields:

  * <i>anat_dims</i>: A two element vector containing the FOV and the slice thickness of the anatomical images.

  * <i>anat_slices</i>: A vector containing the slice numbers of the anatomical images to be plotted.

   * <i>plot_mesh</i>: If set to 1, this field will allow plotting of the roi mesh. Otherwise, set to 0.

   * <i>plot_mask</i>: If set to 1, this field will allow plotting of the mask.  Otherwise, set to 0.

   * <i>plot_fibers</i>:  If set to 1, this field will allow plotting of a single set of fiber tracts. If set to 2, this will allow plotting of two sets of fiber tracts. Otherwise, set to 0.

Depending on the plot options selected, additional fields may be required in <i>plot_options</i>. If <i>plot_mesh</i> equals 1, you must also specify:
   * <i>mesh_size</i>: This specifies the in-plane matrix size of the images used to generate the mesh.

   * <i>mesh_dims</i>: This specifies the FOV and slice thickness of the images used to create the mesh.

   * <i>mesh_color</i>: This is a three element vector containing the RGB levels to be used when plotting the mesh

   * <i>mesh_dist</i>: If you shifted the mesh for fiber tracking and want to show add this field and set to the value used during fiber tracking.

If <i>plot_mask</i> equals 1, you must also specify:
   * <i>mask_size</i>: This specifies the in-plane matrix size of the images used to generate the mask.

   * <i>mask_dims</i>: This two-element vector specifies the FOV and slice thickness of the images used to create the mask.

   * <i>mask_color</i>: If the mask is to be plotted, create color scale using either of the following two options:
      * If mesh_color is a 3 element vector of values ranging from 0-1, the vector is interpreted as RGB levels.
      * If mesh_color is a matrix with size of (#mesh rows) x (#mesh columns) x 3, and if these values range from 0-1, the matrix will be interpreted as RGB levels specific to each tract. This could be used to represent the distribution of architectural parameters across the aponeurosis

If <i>plot_fibers</i> equals 1 or 2, you must also specify:
   * <i>dti_size</i>: A2-element vector that specifies the matrix size of the images used for fiber tracking.

   * <i>dti_dim</i>s: This two-element vector specifies the FOV and slice thickness of the DT images. The FOV is assumed to be square.

   * <i>fiber_color</i>: Defines the color of the tracts. Several options are available:
      * If plot_fibers==1 and fiber_color is a 3 element vector of values ranging from 0-1, the vector is interpreted as RGB levels.
      * If plot_fibers==2 and fiber_color is a 2x3 matrix of values ranging from 0-1, each row of the matrix is interpreted as RGB
       levels for the respective sets of tracts.
      * If plot_fibers==1 and fiber_color is a matrix with size of (#mesh rows) x (#mesh columns) x 3, and if these values range from 0-1, the third dimension of the matrix will be interpreted as RGB levels specific to each tract
      * If plot_fibers==2 and fiber_color is a matrix with size of (#mesh rows) x (#mesh columns) x 3 x 2, and if these values range from 0-1, the third dimension of the matrix will be interpreted as RGB levels specific to each tract, separately for sets 1 and 2

   * <i>fiber_skip</i>: Setting fiber_skip to integer values > 1 will skip over fiber tracts when plotting. This may improve visualization and will decrease time for rendering. If not specified, all fibers will be plotted.

   * <i>contrast_multiplier (optional)</i>: Used to adjust the brightness of the images being displayed. Set to 1 for default brightness, <1 to darken the images, or >1 to brighten them.

 <i>roi_mesh</i>: The roi mesh, defined in define_roi. It is only needed if plot_options.plot_mesh is set to 1.

 <i>mask</i>: A binary mask around the muscle of interest. It could be the output of define_muscle or it could have been defined in another program. It is only needed if plot_options.plot_mask is set to 1.

 <i>fiber_all</i>: The output of fiber_track (original fiber tracts) or fiber_smoother (smoothed fiber tracts). It is only needed if
   plot_options.plot_fibers is set to 1. If plot_fibers equals 1, the size should be (#mesh rows) x (#mesh columns) x (#fiber tract points) x 3. If plot_fibers equals 1, the size should be (#mesh rows) x (#mesh columns) x (#fiber tract points) x 3 x 2.

## 4. Output Arguments
 <i>fiber_figure</i>: A Matlab figure structure

## 5. Acknowledgements
People: Zhaohua Ding, Hannah Kilpatrick

Grant support: NIH/NIAMS R01 AR073831

## 6. Example Code
