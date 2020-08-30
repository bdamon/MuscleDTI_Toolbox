# Help for the function <i>define_roi</i>, v. 0.1.x

## Introduction

This help file contains information about
1) [Purpose of the program](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-define_roi.md#1-purpose)
2) [Usage of the program](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-define_roi.md#2-usage)
3) [Automated Segmentation Algorithm](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-define_roi.md#2-Automated-Segmentation-Algorithm)
4) [Syntax](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-define_roi.md#4-Syntax)
5) [Example Code](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-define_roi.md#5-Example-Code)
6) [Acknowledgements](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-define_roi.md#6-Acknowledgements)


## 1. Purpose

The function <i>define_roi</i> is used to digitize the aponeurosis of muscle fiber insertion. The digitized points are used to reconstruct a mesh; the mesh is used as the seed surface for fiber tracking.  It is a required input to <i>fiber_track</i> and <i>fiber_quantifier</i>; it may be visualized using <i>fiber_visualizer</i>.

## 2. Usage
There are two options for defining the aponeurosis. 
1) <i>Manual selection</i>: Three figure windows are displayed. The center figure shows the current slice, the left-hand figure shows the preceding slice, and the right-hand figure shows the upcoming slice. An interactive tool is opened that allows the user to adjust the center figure's window and level settings. In the center figure, the edge locations of the mask are indicated. Text prompts in the command window guide the user through the following steps.  First, the user may zoom the image to the area of interest, selecting Enter when finished. Then the user defines the aponeurosis with a series of left mouse clicks. The selected points can form a line or close to form a polygon. A right mouse click is used to complete the definition. At each slice, the user is given the option of repeating the procedure in case of error. This is the current method to use for unipennate muscles. An [instructional video](https://youtu.be/HfQeS_bruQM) is available.

2) <i>Automatic selection</i>: Three figure windows are displayed. The center figure shows the current slice, the left-hand figure shows the preceding slice, and the right-hand figure shows the upcoming slice. An interactive tool is opened that allows the user to adjust the center figure's window and level settings. In the center figure, the edge locations of the mask are indicated. For each slice to be analyzed, the algorithm described in §3.3 is used to present an initial estimate of the aponeurosis’s location and its boundary pixels. The user can correct erroneous assignments in the initial estimate by using the left mouse button to select voxels for removal from the initial assignment and the right mouse button to select voxels to add to the initial assignment. The user selects Enter to proceed. In subsequent slices, the finalized assignment from the preceding slice and the results of an edge detection algorithm are incorporated into the initial segmentation estimate and the process is repeated.

For both selection processes, after all slices are analyzed, the user is given the option of repeating erroneous slices.  The mesh is initially formed with dimensions of N<sub>R,A</sub> x N<sub>C,A</sub>.  To smooth the mesh, it is then down-sampled by a factor of four. Finally, the smoothed mesh is interpolated at the desired resolution. A file called <i>roi_mesh_file.mat</i> is automatically saved in the working directory. The user is advised to rename this file promptly. 

The mesh may be viewed using <i>fiber_visualizer</i>, either as part of the function call to <i>define_roi</i> or directly from the command line

[Back to the top](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-define_roi.md)

## 3. Automated Segmentation Algorithm
In the first slice analyzed, the muscle mask is eroded (boundary pixels are removed) and then multiplied by the anat_image slice to remove image regions outside of the muscle of interest. Then, a k-means clustering algorithm is used to segment the remaining image data into three clusters.  The cluster corresponding to the highest signal intensities is assumed to represent the muscle.  The other clusters are presented to the user as the initial estimate of the aponeurosis’s location. The edge pixels are found and a Savitsky-Golay filter is used to form a smoothed curve that indicates the location of the roi_mesh points for that slice.  The user can correct pixel locations, as described in §3.2; the roi_mesh points are automatically updated.

In subsequent slices, additional information is incorporated into the initial estimate, including the preceding aponeurosis segmentation and the results of an edge detection within the muscle mask.  These images are combined using the weights [1, 1, 2] for the k-means, edge, and previous regions, respectively; voxels with sums greater than 2 are included in the initial estimate. The supervision and manual correction steps occur as described above

[Back to the top](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-define_roi.md)

## 4. Syntax

roi_mesh = define_roi(anat_image, mask, defroi_options, plot_options);

<i>anat_image</i>: The imaging data. 

* <i>mask</i>: The mask, as defined by the function <i>define_mask</i> or another source.

* <i>defroi_options</i>: A structure containing the following fields:

    <i>slices</i>: A two-element vector containing the first and last slices that the user wishes to digitize.
  
    <i>dti_size</i>: The size of the DTI image dataset (rows x columns x slices), input as a three element vector.
  
    <i>mesh_size</i>: A two-element vector containing the numbers of rows (N<sub>R,A</sub>) and columns (N<sub>C,A</sub>) desired in the output mesh.
  
    <i>method</i>: a string variable set either to 'manual' or 'auto'. 

* <i>plot_options</i>: Optional. If specified, this calls the <i>fiber_visualizer</i> function to plot the mask and roi mesh.

The output argument is:
* <i>roi_mesh</i>: A 3D matrix containing the reconstructed mesh with size n_row x n_col x 6. In the 3rd dimension, levels 1-3 hold {row column slice} coordinates and levels 4-6 hold the {row column slice} components of the normal vector to the mesh surface at the point {row column slice}.
   
[Back to the top](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-define_roi.md)

## 6. Example Code

### Example 1
Given 

1.	An anatomical image with variable name anat_image and having matrix size 19219244, field of view 192192 mm, and slice thickness 7 mm;

2.	The muscle mask, stored in a variable called mask; and

3.	DTI images, having matrix size 19219244, field of view 192192 mm, and slice thickness 7 mm

the code below will allow the user to:

1.	Manually select aponeurosis in slices 4-31;

2.	Create a mesh of size 150 rows30 columns; and

3.	Visualize the outcome, using slices 14, 24, 34, and 44 of the anatomical image stack for reference.

% Set mesh options:

defroi_options.slices = [4 31]; %analyze slices 4-31

defroi_options.dti_size = [192 192 44]; %matrix size and # of slices in DTI images

defroi_options.mesh_size = [150 30]; %mesh will have 150 rows and 30 columns

defroi_options.method = ‘manual'; %digitize it manually

% Set visualization options

plot_options.anat_dims = [192 7]; %FOV and slice thickness of the images to be displayed, in mm

plot_options.anat_slices = 14:10:44; %display slices 14, 24, 34, and 44 

plot_options.plot_mesh = 1; %do plot the aponeurosis mesh

plot_options.plot_mask = 0; %don’t plot the mask

plot_options.plot_fibers = 0; %don’t plot any fiber tracts

plot_options.mesh_size = [192 192]; %rows x columns of the images used to generate the mesh

plot_options.mesh_dims = [192 7]; %FOV and ST of the images used to create the mesh

plot_options.mesh_color = [0.75 0.75 0.75]; %make the mesh light gray

% call the function:

roi_mesh = define_roi(anat_image, mask, defroi_options, plot_options);
 
### Example 2
Example 2 matches Example 1, except that the mesh is automatically segmented:

% Set mesh options:

defroi_options.slices = [4 31]; 	%analyze slices 4-31

defroi_options.dti_size = [192 192 44]; %matrix size and # of slices in DTI images

defroi_options.mesh_size = [150 30]; %mesh will have 150 rows and 30 columns

defroi_options.method = ‘auto'; %automatic definition

% Set visualization options

plot_options.anat_dims = [192 7]; %FOV and slice thickness of the images to be displayed, in mm

plot_options.anat_slices = 14:10:44; %display slices 14, 24, 34, and 44 

plot_options.plot_mesh = 1; %do plot the aponeurosis mesh

plot_options.plot_mask = 0; %don’t plot the mask

plot_options.plot_fibers = 0; %don’t plot any fiber tracts

plot_options.mesh_size = [192 192]; %rows x columns of the images used to generate the mesh

plot_options.mesh_dims = [192 7]; %FOV and ST of the images used to create the mesh

plot_options.mesh_color = [0.75 0.75 0.75]; %make the mesh light gray

% call the function:

roi_mesh = define_roi(anat_image, mask, defroi_options, plot_options);

[Back to the top](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-define_roi.md)

## 6. Acknowledgements
People: Zhaohua Ding

Grants: NIH/NIAMS R01 AR050101, NIH/NIAMS R01 AR073831

[Back to the top](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-define_roi.md)
