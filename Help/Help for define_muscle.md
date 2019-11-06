# Help for the function define_muscle, v. 0.1

## Introduction

This help file contains information about
1) [Usage of the program](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help%20for%20define_muscle.md#1-usage)
2) [Syntax](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help%20for%20define_muscle.md#2-Syntax)
3) [Input Arguments](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help%20for%20define_muscle.md#3-Input-Arguments)
4) [Output Arguments](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help%20for%20define_muscle.md#4-Output-Arguments)
5) [Acknowledgements](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help%20for%20define_muscle.md#5-Acknowledgements)
6) [Example Code](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help%20for%20define_muscle.md#6-Example-Code)


## 1. Usage

The function <i>define_muscle</i> is used to define the boundary of a muscle and return its binary image mask. This mask is needed as an input to the functions <i>define_roi</i> and <i>fiber_track</i>. Upon calling <i>define_muscle</i>, three images are displayed: the current slice (main window); the preceding slice with its region of interest (ROI) (if present); and the next slice. For the main figure window, an interactive tool is opened that allows the user to adjust the image's window and level settings. 

In the main figure window, the user uses the roipoly tool to define an ROI. After closing the ROI, the user may adjust the position of the vertices. After completing the ROI selection, the program advances to the next slice. ROI selection continues in this manner until all slices have been defined.

By default, the program returns a mask of the same dimensions as the image used to select the ROIs. If desired, the user can also create an alternatively sized mask.  This would be useful if the structural images and the DTI images had different matrix sizes.  
   
A file named mask_file, containing the mask and (if present) the alternatively sized mask, is automatically saved in the working directory.


## 2. Syntax

[mask, alt_mask] = define_muscle(anat_image, slices, alt_mask_size, plot_options);
 
## 3. Input Arguments

<i>anat_image</i>: A row x column x slices stack of images, to be used for selecting the ROIs

<i>slices</i>: A two element vector containing the first and last slices to be analyzed

<i>alt_mask_size</i>: If specified, this is a two element vector containing the row x column size of a second mask; the same number of slices is assumed.

<i>plot_options</i>: If specified, this calls the <i>fiber_visualizer</i> function to plot the mask.


## 4. Output Arguments

<i>mask</i>: the binary image mask, with size matching that of the original image

<i>alt_mask</i>: a second binary image mask, with size matching that of the vector alt_mask_size
 
 
## 5. Acknowledgements

Grant support: NIH/NIAMS R01 AR050101, NIH/NIAMS R01 AR073831

 

## 6. Example code

### Example 1:

Given an image with variable name anat_image and having matrix size 256 x 256 x 50, field of view 192 x 192 mm, and slice thickness = 6 mm, the code below will allow the user to:
  1) Define the muscle mask in slices 15-40

% call the function

slices = [15 40];

[mask, ~] = define_muscle(anat_image, slices, [], []);



### Example 2:

Given an image with variable name anat_image and having matrix size 256 x 256 x 50, field of view 192 x 192 mm, and slice thickness = 6 mm, the code below will allow the user to 
  1) Define the muscle mask in slices 15-40; and
  2) Return a second mask of size 192 x 192 x 50

% define the size of the alternative mask:

alt_mask_size = [192 192];

% call the function

slices = [15 40];

[mask, alt_mask] = define_muscle(anat_image, slices, alt_mask_size, []);


### Example 3: 

Given an image with variable name anat_image and having matrix size 256 x 256 x 50, field of view 192 x 192 mm, and slice thickness = 6 mm, the code below will allow the user to 
  1) Define the muscle mask in slices 15-40;
  2) Return a mask of size 192 x 192 x 50; and
  3) Automatically visualize the result using fiber_visualizer, using slices 15, 25, 35, and 45 for anatomical reference

% Set options for plotting:

plot_options.plot_fibers=0;                  %don't plot any fiber tracts

plot_options.plot_mesh=0;                    %don't plot an aponeurosis mesh

plot_options.plot_mask=1;                    %do plot the mask

plot_options.anat_dims=[192 6];              %FOV and slice thickness of the images to be displayed, in mm

plot_options.anat_slices=15:10:45;           %display slices 15, 25, 35, and 45 for anatomical reference in fiber_visualizer

plot_options.mask_size=[256 256];            %in-plane matrix size of the images used to generate the mask

plot_options.mask_dims=[192 6];              %FOV and slice thickness of the images used to generate the mask, in mm

plot_options.mask_color=[1 0 0];             %make the mask a red, semi-transparent overlay


% define the size of the alternative mask:

alt_mask_size = [192 192];

% call the function

slices = [15 40];

[mask, alt_mask] = define_muscle(anat_image, slices, alt_mask_size, plot_options);
