# Help for the function <i>fiber_visualizer</i>, v. 0.1.x

## Introduction

This help file contains information about
1) [Usage of the program](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-dfiber_visualizer.md#1-usage)
2) [Syntax](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_visualizer.md#2-Syntax)
3) [Input Arguments](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_visualizer.md#3-Input-Arguments)
4) [Output Arguments](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_visualizer.md#4-Output-Arguments)
5) [Acknowledgements](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_visualizer.md#5-Acknowledgements)
6) [Example Code](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_visualizer.md#6-Example-Code)


## 1. Usage
The function fiber_visualizer is used to visualize images and the muscle mask, roi mesh, and/or fiber tracts formed using the MuscleDTI_Toolbox. Input the images and a structure containing the plotting options. Depending on the data you wish to view, the roi mesh, mask, and/or fiber tracts must also be input and additional plotting options will be required.


## 2. Syntax 
fiber_figure = fiber_visualizer(anat_image, plot_options, roi_mesh, mask, fiber_all)

## 3. Input Arguments
<i>anat_image</i>: The stack of images to be plotted for anatomical reference

<i>plot_options</i>: A structure containing the following required fields:

  * <i>anat_dims</i>: A two element vector containing the FOV and the slice thickness of the anatomical images.

  * <i>anat_slices</i>: A vector containing the slice numbers of the anatomical images to be plotted.

   * <i>plot_mesh</i>: If set to 1, this field will allow plotting of the roi mesh. Otherwise, set to 0.

   * <i>plot_mask</i>: If set to 1, this field will allow plotting of the mask.  Otherwise, set to 0.

   * <i>plot_fibers</i>:  If set to 1, this field will allow plotting of a single set of fiber tracts. If set to 2, this will allow plotting of two sets of fiber tracts. Otherwise, set to 0.

 Depending on the plot options selected, the following other fields may be required:

 If plot_mesh equals 1, you must also specify:
   * <i>mesh_size</i>: This specifies the in-plane matrix size of the images used to generate the mesh.

   * <i>mesh_dims</i>: This specifies the FOV and slice thickness of the images used to create the mesh.

   * <i>mesh_color</i>: This is a three element vector containing the RGB levels to be used when plotting the mesh

   * <i>mesh_dist</i>: If you shifted the mesh for fiber tracking and want to show add this field and set to the value used during fiber tracking.

 If plot_mask equals 1, you must also specify:
   * <i>mask_size</i>: This specifies the in-plane matrix size of the images used to generate the mask.

   * <i>mask_dims</i>: This two-element vector specifies the FOV and slice thickness of the images used to create the mask.

   * <i>mask_color</i>: If the mask is to be plotted, create color scale using either of the following two options:
      -If mesh_color is a 3 element vector of values ranging from 0-1, the vector is interpreted as RGB levels.
      -If mesh_color is a matrix with size of (#mesh rows) x (#mesh columns) x 3, and if these values range from 0-1, the matrix will be interpreted as RGB levels specific to each tract. This could be used to represent the distribution of architectural parameters across the aponeurosis

 If plot_fibers equals 1 or 2, you must also specify:
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
