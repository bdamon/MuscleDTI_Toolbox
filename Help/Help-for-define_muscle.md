# Help for the function [<i>define_muscle</i>](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Tractography-Functions/define_muscle.m), v. 1.0.0

## Introduction

This help file contains information about
1) [Purpose of the Program](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-define_muscle.md#1-Purpose)
2) [Usage of the Program](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-define_muscle.md#2-Usage)
3) [Syntax](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-define_muscle.md#3-Syntax)
4) [Example Code](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-define_muscle.md#4-Example-Code)
5) [Acknowledgements](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-define_muscle.md#5-Acknowledgements)


## 1. Purpose

The function <i>define_muscle</i> is used to define the boundary of a muscle and return its binary image mask. This mask is needed by the functions [<i>define_roi</i>](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-define_roi.md) and [<i>fiber_track</i>](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_track.md). It may be visualized using [<i>fiber_visualizer</i>](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_visualizer.md). 

[Back to the Top](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-define_muscle.md)

## 2. Usage
The user provides the anatomical images to be segmented and defines the slice numbers of interest.  After calling the function, a single figure window is opened. The initial slice of interest is displayed in the middle panel; the preceding two slices (if present) are displayed in the left-most column; and the next two slices (if present) are displayed in the column at the immediate left. For the center panel, the <i>zoom</i> tool is enabled; the user clicks and drags the left mouse button to zoom to the muscle of interest.  To close the zoom tool, the user selects Enter on their keyboard. All images are then zoomed to this level. 

Using the <i>roipoly</i> tool, the user uses a series of left mouse clicks to define the muscle as a region of interest (ROI) in the middle panel. After defining the vertices and closing the ROI, the user can adjust the vertex positions.  To complete ROI selection, the user right-clicks the mouse; this brings up bring up a menu, from which the user selects Create Mask to complete ROI selection.

Then the program advances to the next slice. In this slice and all subsequent slices, the level of zoom is automatically set as ±20 pixels beyond the previous ROI’s row and column limits. In the lower left panel, the preceding slice and its ROI are shown. Also shown are gold and red lines depicting the center row and column, respectively, of this ROI. In the column to the immediate right of the main panel, the projections of the image stack and the ROI along the red line are shown.  In the far-right column, the projections of the image stack and the ROI along the gold line are shown. Examining these windows can help the user maintain consistency in ROI selection. ROI selection continues in this manner until all slices of interest have been defined. These steps are illustrated in the figure below.

<img src="https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Images/FigureS1.png" 
alt="Figure S1" width="960" height="316" border="5" />
<b><i>Example window from <i>define_muscle</i></b>, showing (from left) the two previous slices; the two upcoming slices; the current slice; a sagittal reconstruction of the segmentation; and a coronal reconstruction of the segmentation.</i>

By default, the mask has the same dimensions as the input image. If the DTI images and structural images have different dimensions from each other, an alternatively sized mask may also be calculated.  A MATLAB data file named <i>mask_file.mat</i>, containing the mask and the alternatively sized mask (if present), is automatically saved in the working directory. The user is advised to rename this file promptly.

The mask may be viewed using [<i>fiber_visualizer</i>](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_visualizer.md), either as part of the function call to <i>define_muscle</i> or directly from the command line.

[Back to the Top](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-define_muscle.md)

## 3. Syntax

The function define_muscle is called using the following syntax:

[mask, alt_mask] = define_muscle(anat_image, slices, alt_mask_size, fv_options);

The input arguments are:

* <i>anat_image</i>: A N<sub>R,A</sub> x N<sub>C,A</sub> x N<sub>S,A</sub> stack of images, which the user will use to segment the muscle of interest;

* <i>slices</i>: A two-element vector containing the first and last slices to be analyzed, entered as [first last];

* <i>alt_mask_size</i>: If specified, this is a three-element vector containing the row x column x slice size of a second mask; the same center position, FOV, and foot-head distance of the image stack are assumed. This would be used if the numbers of rows, columns, and slices in the anatomical image dataset differed from the numbers of rows, columns, and slices in the DTI dataset.

* <i>fv_options</i>: If specified, this calls the <i>fiber_visualizer</i> function to plot the mask.

The output arguments are:

* <i>mask</i>: the binary image mask, with dimensions matching that of the original image; and

* <i>alt_mask</i>: If <i>alt_mask_size</i> was specified as in input argument, a second binary image mask is returned.

[Back to the Top](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-define_muscle.md)
 

## 4. Example Code

### Example 1:

Given 

1.	An anatomical image with variable name anat_image and having matrix size 192 x 192 x 44, field of view 192 x 192 mm, and slice thickness 7 mm;

the code below will allow the user to:

1.	Define the muscle mask in slices 4-41 and 

2.	Return a mask of size 192 x 192 x 44:

% define options

slices = [4 41];

% call the function:

[mask, ~] = define_muscle(anat_image, slices, [], []);
 
### Example 2
Given 

1.	An anatomical image with variable name anat_image and having matrix size 192 x 192 x 44, field of view 192 x 192 mm, and slice thickness 7 mm;

the code below will allow the user to 

1.	Define the muscle mask in slices 4-41; 

2.	Return a mask of size 192 x 192 x 44; and 

3.	Return a second mask of size 128 x 128 x 44

% define options:

slices = [4 41];

alt_mask_size = [128 128 44];

% call the function: 

[mask, alt_mask] = define_muscle(anat_image, slices, alt_mask_size, []);

### Example 3

Given 

1.	An anatomical image with variable name anat_image and having matrix size 19219244, field of view 192192 mm, and slice thickness 7 mm;

the code below will allow the user to 

1.	Define the muscle mask in slices 4-41; 

2.	Return a mask of size 192 x 192 x 44; 

3.	Return a second mask of size 128 x 128 x 44; and 

4.	Visualize the result using <i>fiber_visualizer</i>, using slices 14, 24, 34, and 44 for anatomical reference:

% Set mask definition options

alt_mask_size = [128 128 44];

slices = [4 41];

% Set visualization options

fv_options.anat_dims = [192 7]; %FOV and slice thickness of the images to be displayed, in mm

fv_options.anat_slices = 14:10:44; %display slices 14, 24, 34, and 44 

fv_options.plot_mesh = 0; %don’t plot an aponeurosis mesh

fv_options.plot_mask = 1; %do plot the mask

fv_options.plot_fibers = 0; %don’t plot any fiber tracts

fv_options.mask_size = [192 192]; %rows x columns of the images used to generate the mask

fv_options.mask_dims = [192 7]; %FOV and ST of the images used to create the mask

fv_options.mask_color = [1 0 0]; %make the mask a red, semi-transparent overlay

% call the function

 [mask, alt_mask] = define_muscle(anat_image, slices, alt_mask_size, fv_options); 

[Back to the Top](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-define_muscle.md)

## 5. Acknowledgements

Grant support: NIH/NIAMS R01 AR050101, NIH/NIAMS R01 AR073831

[Back to the Top](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-define_muscle.md)
 
