# Help for the function define_roi, v. 0.1

## Introduction

This help file contains information about
1) [Usage of the program](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help%20for%20define_roi.md#1-usage)
2) [Syntax](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help%20for%20define_roi.md#2-Syntax)
3) [Input Arguments](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help%20for%20define_roi.md#3-Input-Arguments)
4) [Output Arguments](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help%20for%20define_roi.md#4-Output-Arguments)
5) [Acknowledgements](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help%20for%20define_roi.md#5-Acknowledgements)
6) [Example Code](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help%20for%20define_roi.md#6-Example-Code)

## 1. Usage

The function define_roi is used to digitize the aponeurosis of muscle fiber insertion in the MuscleDTI_Toolbox.  The digitized points are used to reconstruct a mesh; the mesh is used as the seed surface for fiber tracking.

There are two options for defining the aponeurosis:
1) Manual: The user is prompted initially to select two points, which define the level of zoom to be used throughout the entire process. Then the user advances through the slices to select the aponeurosis. The selected points can form a line or close to form a polygon. At each slice, the user is given the option of repeating the procedure in case of error.  For each figure window, an interactive tool is opened that allows the user to adjust the image's window and level settings.  Eventually, the manual option will be removed.
2) Automatic: The aponeurosis is automatically segmented from within the region of the image represented by the muscle mask. Two automated segmentation methods (edge detection and k-means clustering) are applied. In addition, the region fromthe preceding slice, if available, is used.  The consensus region from these three methods is presented to the user; the user is allowed to correct misassignments by adding points (right mouse button) or deleting points (left mouse button). The boundaries of the segmented region are smoothed using a Savitsky-Golay filter and used to form the mesh. These points are previwed, and the user can further correct points before advancing to the next slice.

The mesh is initially formed with resolution n_row x n_col.  To smooth the mesh, it is then downsampled by a size factor of four. Finally, the smoothed mesh is used to create a high resolution mesh at the desired size. A file called roi_mesh_file.mat is automatically saved in the working directory. 

If the input argument plot_options is included, the mesh and mask are plotted using the function fiber_visualizer.

## 2. Syntax

roi_mesh=define_roi(anat_image, mask, defroi_options, plot_options);

## 3. Input Arguments
anat_image: The imaging data. If input as a structure, then the imaging data are assumed to exist in a field called anat_image.Data.  If specified as a matrix, the data are used directly.

mask: The mask, as defined by the function define_mask or other method.

defroi_options: A structure containing the following fields:

  -slices: A two-element vector containing the first and last slices that the user wishes to digitize.
  
  -dti_size: The size of the DTI image dataset (rows x columns x slices), input as a three element vector.
  
  -mesh_size: A two-element vector containing the numbers of rows (n_row) and columns (n_col) desired in the output mesh.
  
  -method: a string variable set either to 'manual' or 'auto'. The manual and automatic options are described above.

plot_options: Optional. If specified, this calls the fiber_visualizer function to plot the mask and roi mesh.

## 4. Output Arguments
roi_mesh: a 3D matrix containing the reconstructed mesh with size rows x columns x 6. In the 3rd dimension, levels 1-3 hold the row-column-slice data and levels 4-6 hold the normal vector to the roi surface at the point {row column slice}.
   
   
## 5. Acknowledgements

People: Zhaohua Ding

Grant support: NIH/NIAMS R01 AR050101, NIH/NIAMS R01 AR073831

## 6. Example Code


### Example 1:

Given 1) an anatomical image with variable name anat_image and having matrix size 256 x 256 x 50, field of view 192 x 192 mm, and slice thickness = 6 mm and 2) DTI images with variable name dti_images and having matrix size 192 x 192 x 50, field of view 192 x 192 mm, and slice thickness = 6 mm, the code below will allow the user to 
  1) Use automated selection of the aponeurosis;
  2) Define the muscle aponeurosis in slices 15-30; and
  3) Create a mesh of size 175 rows x 50 columns

% Set mesh options:

defroi_options.slices = [15 30];

defroi_options.dti_size = [192 192 50];

defroi_options.mesh_size = [175 50];

defroi_options.method='auto';

% call the function:

roi_mesh=define_roi(anat_image, mask, defroi_options, []);

### Example 2

Example 2 matches Example 1, except that the mesh is automatically called from within the define_roi function:

% Set mesh options:

defroi_options.slices = [15 30];

defroi_options.dti_size = [192 192 50];

defroi_options.mesh_size = [175 50];

defroi_options.method='auto';

% Set options for plotting:

plot_options.plot_fibers=0;                  %don't plot any fiber tracts

plot_options.plot_mesh=1;                    %don't plot an aponeurosis mesh

plot_options.plot_mask=0;                    %do plot the mask

plot_options.anat_dims=[192 6];              %FOV and slice thickness of the images to be displayed, in mm

plot_options.anat_slices=15:10:45;           %slices to display for anatomical reference in fiber_visualizer

plot_options.mesh_size=[256 256];            %in-plane matrix size of the images used to generate the mask

plot_options.mesh_dims=[192 6];              %FOV and slice thickness of the images used to generate the mask, in mm

plot_options.mesh_color=[.75 .75 .75];       %make the mesh gray

% call the function:

roi_mesh=define_roi(anat_image, mask, defroi_options, []);
