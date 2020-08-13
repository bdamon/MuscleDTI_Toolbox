# Help for the function <i>define_muscle</i>, v. 0.1.x

## Introduction

This help file contains information about
1) [Overview and Usage of the Program](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-define_muscle.md#1-Overview-and-Usage)
2) [Syntax](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-define_muscle.md#2-Syntax)
3) [Example Code](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-define_muscle.md#3-Example-Code)
4) [Help](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-define_muscle.md#4-Help)
5) [Acknowledgements](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-define_muscle.md#5-Acknowledgements)


## 1. Overview and Usage

The function <i>define_muscle</i> is used to define the boundary of a muscle and return its binary image mask. An image mask is a matrix of zeros and ones having the same matrix size as the original image. The ones indicate the location of the muscle in the image stack.This mask is needed by the functions [<i>define_roi</i>](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help%20for%20define_roi.md) and [<i>fiber_track</i>](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help%20for%20fiber_track.md). 

Upon calling <i>define_muscle</i>, three images are displayed: the current slice in the main (center) window; the preceding slice (left window); and the next slice (right window). For the main figure window, an interactive tool is opened that allows the user to adjust the image's contrast and brightness levels. Also, the zoom tool is enabled; the user clicks and drags the left mouse button to zoom in on the muscle of interest. To close the zoom tool, the user selects Enter on their keyboard. 

Using the roipoly tool, the user uses a series of left mouse clicks to define the muscle as a region of interest (ROI) in the middle window. After defining the vertices and closing the ROI, the user can adjust vertex positions.  To complete ROI selection, the user right-clicks the mouse; this brings up bring up a menu, from which the user selects Create Mask to complete ROI selection. Then the program advances to the next slice. In this slice and all subsequent slices, the level of zoom is automatically set as ±5 pixels beyond the previous ROI’s row and column limits. In the left figure, the preceding slice and its ROI are shown.  The righthand figure, the next slice is shown. Examining these windows can help to maintain consistency in ROI selection. ROI selection continues in this manner until all slices of interest have been defined. This procedure continues until ROIs in all slices of interest have been defined.

By default, the program returns a mask of the same dimensions as the image used to select the ROIs. If desired, the user can also create an alternatively sized mask.  This is necessary if the structural images and the DTI images have different matrix sizes.  
   
A MATLAB data file named mask_file.mat, containing the mask and (if present) the alternatively sized mask, is automatically saved in the working directory. The user is advised to rename this file to replace the generic name with something more informative.

The mask may be viewed using [<i>fiber_visualizer</i>](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help%20for%20fiber_visualizer.md).

## 2. Syntax

The function define_muscle is called using the following syntax:

[mask, alt_mask] = define_muscle(anat_image, slices, alt_mask_size, plot_options);

The input arguments are:

* <i>anat_image</i>: A row x column x slices stack of images, which the user will use to segment the muscle of interest;

* <i>slices</i>: A two-element vector containing the first and last slices to be analyzed;

* <i>alt_mask_size</i>: If specified, this is a three-element vector containing the row x column x slice size of a second mask; the same center position of the image stack is assumed; and

* <i>plot_options</i>: If specified, this calls the <i>fiber_visualizer</i> function to plot the mask.

The output arguments are:

* <i>mask</i>: the binary image mask, with dimensions matching that of the original image; and

* <i>alt_mask</i>: If <i>alt_mask_size</i> was specified as in input argument, a second binary image mask is returned.
 

## 3. Example Code

### Example 1:

Given an image with variable name anat_image and having matrix size 192 x 192 x 44, field of view 192 x 192 mm, and slice thickness = 7 mm, the code below will allow the user to:
  1) Define the muscle mask in slices 4-41 and
  2) Return a mask of size 192 x 192 x 44.

% define input options

slices = [4 41];

% call the function:

[mask, ~] = define_muscle(anat_image, slices, [], []);



### Example 2:

Given an image with variable name anat_image and having matrix size 192 x 192 x 44, field of view 192 x 192 mm, and slice thickness = 7 mm, the code below will allow the user to 
  1) Define the muscle mask in slices 4-41; 
  2) Return a mask of size 192 x 192 x 44; and
  3) Return a second mask of size 128 x 128 x 44

% define the size of the alternative mask:

alt_mask_size = [128 128 44];

% define the slices of interest and call the function

slices = [4 41];

[mask, alt_mask] = define_muscle(anat_image, slices, alt_mask_size, []);


### Example 3: 

Given an image with variable name anat_image and having matrix size 192 x 192 x 44, field of view 192 x 192 mm, and slice thickness = 7 mm, the code below will allow the user to 
  1) Define the muscle mask in slices 4-41;
  2) Return a mask of size 192 x 192 x 44; 
  3) Return a second mask of size 128 x 128 x 44; and
  4) Automatically visualize the result using fiber_visualizer, using slices 14, 24, 34, and 44 for anatomical reference

% Set options for plotting:

plot_options.plot_fibers=0;                  %don't plot any fiber tracts

plot_options.plot_mesh=0;                    %don't plot an aponeurosis mesh

plot_options.plot_mask=1;                    %do plot the mask

plot_options.anat_dims=[192 7];              %FOV and slice thickness of the images to be displayed, in mm

plot_options.anat_slices=14:10:44;           %display slices 14, 24, 34, and 44 for anatomical reference in fiber_visualizer

plot_options.mask_size=[192 192];            %in-plane matrix size of the images used to generate the mask

plot_options.mask_dims=[192 7];              %FOV and slice thickness of the images used to generate the mask, in mm

plot_options.mask_color=[1 0 0];             %make the mask a red, semi-transparent overlay


% define the size of the alternative mask:

alt_mask_size = [128 128 44];

% call the function

slices = [4 41];

[mask, alt_mask] = define_muscle(anat_image, slices, alt_mask_size, plot_options);

## 4. Help

<i>define_muscle</i> contains help comments that can be viewed within MATLAB by typing:

help define_muscle

on the command line. An [instructional video](https://youtu.be/TWTZvgVWoB4) detailing the use of this function is available.

## 5. Acknowledgements

Grant support: NIH/NIAMS R01 AR050101, NIH/NIAMS R01 AR073831
