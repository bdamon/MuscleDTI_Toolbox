# Help for the function <i>define_muscle</i>, v. 0.1

## Introduction

This help file contains information about
1) [Usage of the program](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help%20for%20define_muscle.md#1-usage)
2) [Syntax](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help%20for%20define_muscle.md#2-Syntax)
3) [Input Arguments](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help%20for%20define_muscle.md#3-Input-Arguments)
4) [Output Arguments](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help%20for%20define_muscle.md#4-Output-Arguments)
5) [Acknowledgements](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help%20for%20define_muscle.md#5-Acknowledgements)
6) [Example Code](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help%20for%20define_muscle.md#6-Example-Code)


## 1. Usage

The function <i>define_muscle</i> is used to define the boundary of a muscle and return its binary image mask. An image mask is a matrix of zeros and ones having the same matrix size as the original image. The ones indicate the location of the muscle in the image stack.

This mask is needed as an input to the functions <i>define_roi</i> and <i>fiber_track</i>. Upon calling <i>define_muscle</i>, three images are displayed: the current slice in the main (center) window; the preceding slice (left window); and the next slice (right window). For the main figure window, an interactive tool is opened that allows the user to adjust the image's contrast and brightness levels. Also, the zoom tool is enabled; the user should zoom to the muscle of interest. To close the zoom tool, the user selects Enter on their keyboard. 

In the main figure window, the user uses the roipoly tool to define a region of interest (ROI). After closing the ROI, the user may adjust the position of the vertices. When satisfied with the ROI, the user should right-click the mouse to bring up a menu, selecting Create Mask to complete the ROI selection. Then the program displays the next slice in the main window; the level of zoom is automatically set. In the left figure, the preceding slice and its ROI are shown.  The righthand figure, the next slice is shown. Examining htese windows can help to maintain consistency in ROI selection. ROI selection continues in this manner until all slices of interest have been defined.

By default, the program returns a mask of the same dimensions as the image used to select the ROIs. If desired, the user can also create an alternatively sized mask.  This is necessary if the structural images and the DTI images had different matrix sizes.  
   
A MATLAB data file named mask_file.mat, containing the mask and (if present) the alternatively sized mask, is automatically saved in the working directory. The user is advised to rename this file to replace the generic name with more informative name.

An [instructional video](https://youtu.be/Ot-cvL3oRso) is available.


## 2. Syntax

[mask, alt_mask] = define_muscle(anat_image, slices, alt_mask_size, plot_options);
 
## 3. Input Arguments

* <i>anat_image</i>: A row x column x slices stack of images, to be used for selecting the ROIs

* <i>slices</i>: A two element vector containing the first and last slices to be analyzed

* <i>alt_mask_size</i>: If specified, this is a two element vector containing the row x column size of a second mask; the same number of slices is assumed.

* <i>plot_options</i>: If specified, this calls the <i>fiber_visualizer</i> function to plot the mask.


## 4. Output Arguments

* <i>mask</i>: the binary image mask, with size matching that of the original image

* <i>alt_mask</i>: a second binary image mask, with size matching that of the vector alt_mask_size
 
 
## 5. Acknowledgements

Grant support: NIH/NIAMS R01 AR050101, NIH/NIAMS R01 AR073831

 

## 6. Example code

### Example 1:

Given an image with variable name anat_image and having matrix size 192 x 192 x 44, field of view 192 x 192 mm, and slice thickness = 7 mm, the code below will allow the user to:
  1) Define the muscle mask in slices 4-41

% define the slices of interest and call the function

slices = [4 41];

[mask, ~] = define_muscle(anat_image, slices, [], []);



### Example 2:

Given an image with variable name anat_image and having matrix size 192 x 192 x 44, field of view 192 x 192 mm, and slice thickness = 7 mm, the code below will allow the user to 
  1) Define the muscle mask in slices 4-41; and
  2) Return a second mask of size 128 x 128 x 44

% define the size of the alternative mask:

alt_mask_size = [128 128];

% define the slices of interest and call the function

slices = [4 41];

[mask, alt_mask] = define_muscle(anat_image, slices, alt_mask_size, []);


### Example 3: 

Given an image with variable name anat_image and having matrix size 192 x 192 x 44, field of view 192 x 192 mm, and slice thickness = 7 mm, the code below will allow the user to 
  1) Define the muscle mask in slices 4-41;
  2) Return a mask of size 192 x 192 x 44; 
  3) Return a second mask of size 128 x 128 x 44; and
  3) Automatically visualize the result using fiber_visualizer, using slices 14, 24, 34, and 44 for anatomical reference

% Set options for plotting:

plot_options.plot_fibers=0;                  %don't plot any fiber tracts

plot_options.plot_mesh=0;                    %don't plot an aponeurosis mesh

plot_options.plot_mask=1;                    %do plot the mask

plot_options.anat_dims=[192 7];              %FOV and slice thickness of the images to be displayed, in mm

plot_options.anat_slices=14:10:44;           %display slices 15, 25, 35, and 45 for anatomical reference in fiber_visualizer

plot_options.mask_size=[192 192];            %in-plane matrix size of the images used to generate the mask

plot_options.mask_dims=[192 7];              %FOV and slice thickness of the images used to generate the mask, in mm

plot_options.mask_color=[1 0 0];             %make the mask a red, semi-transparent overlay


% define the size of the alternative mask:

alt_mask_size = [128 128];

% call the function

slices = [4 41];

[mask, alt_mask] = define_muscle(anat_image, slices, alt_mask_size, plot_options);
