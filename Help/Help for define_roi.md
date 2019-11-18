# Help for the function <i>define_roi</i>, v. 0.1.x

## Introduction

This help file contains information about
1) [Usage of the program](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help%20for%20define_roi.md#1-usage)
2) [Syntax](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help%20for%20define_roi.md#2-Syntax)
3) [Input Arguments](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help%20for%20define_roi.md#3-Input-Arguments)
4) [Output Arguments](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help%20for%20define_roi.md#4-Output-Arguments)
5) [Acknowledgements](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help%20for%20define_roi.md#5-Acknowledgements)
6) [Example Code](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help%20for%20define_roi.md#6-Example-Code)

## 1. Usage

The function <i>define_roi</i> is used to digitize the aponeurosis of muscle fiber insertion in the MuscleDTI_Toolbox.  The digitized points are used to reconstruct a mesh; the {row, column, slice} coordinates on the mesh are used as the seed points for fiber tracking using <i>fiber_track</i>.

There are two options for defining the aponeurosis. 
1) <i>Manual</i>: Three windows are opened: the middle (main) window displaying the current slice, the left window displaying the preceding slice, and the right window displaying the upcoming slice. In the first slice, the user is prompted to zoom the window; the user should select Enter on their keyboard when the zoom level is set.  The user can also adjust the contrast and brightness of the image in the main figure window.  The user then uses the left mouse button to select points along the aponeurosis. The selected points can form a line, or they can close to form a polygon. The right mouse button completes the selection. At each slice, the user is given the option of repeating the procedure in case of error (righr mouse button) or advancing to the next slice (left mouse button).  Selection continues in this manner until all slices of interest have been define. The final set of points is used to form the mesh (see below). The manual option is the current best approach for digitizing only a single side of the aponeurosis (such as for unipennate muscles). An [instructional video](https://youtu.be/Q8A8h6JEKBM) is available.

2) <i>Automatic</i>: The aponeurosis is automatically segmented from within the region of the image represented by the muscle mask. Two automated segmentation methods (edge detection and k-means clustering) are applied. In addition, the region from the preceding slice, if available, is used.  The consensus region from these three methods is presented to the user in the middle of three windows; the user is allowed to correct misassignments by adding points (right mouse button) or deleting points (left mouse button). The lefthand window shows the selected pixels as points; this is a helpful way to preview the image intensity data and ensure consistent selection criteria for the apoenurosis. The boundaries of the segmented region are smoothed using a Savitsky-Golay filter and previewed in the righthand window; these are the points that will be used to form the mesh (see next paragraph).  The user can correct points before advancing to the next slice.

Initially, the mesh is formed at the resolution specified by the user in the defroi_options structure.  To smooth the mesh, it is then downsampled by a size factor of four. Finally, the smoothed mesh is used to create a high resolution mesh at the desired size. A file called roi_mesh_file.mat is automatically saved in the working directory. The user is advised to rename the file.

If the input argument plot_options is included, the mesh is automatically plotted using the function <i>fiber_visualizer</i>.

## 2. Syntax

roi_mesh = define_roi(anat_image, mask, defroi_options, plot_options);

## 3. Input Arguments
<i>anat_image</i>: The imaging data. If input as a structure, then the imaging data are assumed to exist in a field called anat_image.Data.  If specified as a matrix, the data are used directly.

* <i>mask</i>: The mask, as defined by the function define_mask or other method.

* <i>defroi_options</i>: A structure containing the following fields:

    <i>.slices</i>: A two-element vector containing the first and last slices that the user wishes to digitize.
  
    <i>.dti_size</i>: The size of the DTI image dataset (rows x columns x slices), input as a three element vector.
  
    <i>.mesh_size</i>: A two-element vector containing the numbers of rows (n_row) and columns (n_col) desired in the output mesh.
  
    <i>.method</i>: a string variable set either to 'manual' or 'auto'. The manual and automatic options are described above.

* <i>plot_options</i>: Optional. If specified, this calls the <i>fiber_visualizer</i> function to plot the mask and roi mesh.

## 4. Output Arguments
* <i>roi_mesh</i>: A 3D matrix containing the reconstructed mesh with size n_row x n_col x 6. In the 3rd dimension, levels 1-3 hold {row column slice} coordinates and levels 4-6 hold the {row column slice} components of the normal vector to the mesh surface at the point {row column slice}.
   
   
## 5. Acknowledgements

People: Zhaohua Ding

Grant support: NIH/NIAMS R01 AR050101, NIH/NIAMS R01 AR073831

## 6. Example Code


### Example 1:

Given 1) an anatomical image with variable name anat_image and having matrix size 192 x 192 x 44, field of view 192 x 192 mm, and slice thickness 7 mm and 2) DTI images having matrix size 192 x 192 x 44, field of view 192 x 192 mm, and slice thickness 7 mm, the code below will allow the user to:
  1) Manually select aponeurosis in slices 6-30;
  2) Create a mesh of size 150 rows x 30 column; and
  3) Visualize the outcome, using slices 14, 24, 34, and 44 of teh anatomical image stack for reference.

%% Example 1 in the define_roi help page:

% Set plotting options:

plot_options.plot_fibers=0;                         %don't plot any fiber tracts

plot_options.plot_mesh=1;                           %do plot an aponeurosis mesh

plot_options.plot_mask=0;                           %don't plot the mask

plot_options.anat_dims=[192 7];                     %FOV and slice thickness of the images to be displayed, in mm

plot_options.anat_slices=14:10:44;                  %display slices 14, 24, 34, and 44 when viewingthe mesh in fiber_visualizer

plot_options.mesh_size=[192 192];                   %in-plane matrix size of the images used to generate the mask

plot_options.mesh_dims=[192 7];                     %FOV and slice thickness of the images used to generate the mask, in mm

plot_options.mesh_color=[.75 .75 .75];              %make the mesh light gray

% Set mesh options:

defroi_options.slices = [4 31];                    %analyze slices 4-31

defroi_options.dti_size = [192 192 44];             %matrix size and # of slices in DTI images

defroi_options.mesh_size = [150 30];                %mesh will have 150 rows and 30 columns

defroi_options.method='manual';                     %digitize it manually

% call the function:

roi_mesh=define_roi(anat_image, mask, defroi_options, plot_options);


### Example 2

Example 2 matches Example 1, except that the mesh is automatically segmented:

%% Example 1 in the define_roi help page:

% Set plotting options:

plot_options.plot_fibers=0;                         %don't plot any fiber tracts

plot_options.plot_mesh=1;                           %do plot an aponeurosis mesh

plot_options.plot_mask=0;                           %don't plot the mask

plot_options.anat_dims=[192 7];                     %FOV and slice thickness of the images to be displayed, in mm

plot_options.anat_slices=14:10:44;                  %display slices 14, 24, 34, and 44 when viewingthe mesh in fiber_visualizer

plot_options.mesh_size=[192 192];                   %in-plane matrix size of the images used to generate the mask

plot_options.mesh_dims=[192 7];                     %FOV and slice thickness of the images used to generate the mask, in mm

plot_options.mesh_color=[.75 .75 .75];              %make the mesh light gray

% Set mesh options:

defroi_options.slices = [4 31];                    %analyze slices 4-31

defroi_options.dti_size = [192 192 44];             %matrix size and # of slices in DTI images

defroi_options.mesh_size = [150 30];                %mesh will have 150 rows and 30 columns

defroi_options.method='auto';                     %digitize it manually

% call the function:

roi_mesh=define_roi(anat_image, mask, defroi_options, plot_options);
