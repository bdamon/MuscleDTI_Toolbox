USAGE

The function define_muscle is used to define the boundary of a muscle and return its binary image mask. This mask is needed for muscle fiber tracking using the MuscleDTI_Toolbox. Three images are displayed: the current slice (main panel), the preceding slice with its region of interest (ROI; if present), and the next slice. For the main figure window, an interactive tool is opened that allows the user to adjust the image's window and level settings. 

In the main figure window, the user uses the roipoly tool to define an ROI. After closing the ROI, the user may adjust the position of the vertices. After completing the ROI selection, the program advances to the next slice. ROI selection continues in this manner until all slices have been defined.

By default, the program returns a mask of the same dimensions as the image used to select the ROIs.  If desired, the user can also create an alternatively sized mask.  This would be useful if the structural images and the DTI images had different matrix sizes.  
   
A file named mask_file, containing the mask and (if present) the alternatively sized mask, is automatically saved in the working directory.

SYNTAX

[mask, alt_mask] = define_muscle(anat_image, slices, alt_mask_size, plot_options);

INPUT ARGUMENTS

anat_image: A row x column x slices stack of images, to be used for selecting the ROIs

slices: A two element vector containing the first and last slices to be analyzed

alt_mask_size: If specified, this is a two element vector containing the row x column size of a second mask; the same number of slices is assumed.

plot_options: If specified, this calls the fiber_visualizer function to plot the mask.

OUTPUT ARGUMENTS

mask: the binary image mask, with size matching that of the original image

alt_mask: a second binary image mask, with size matching that of the vector alt_mask_size

EXAMPLE USAGE

Given an image with matrix size 256 x 256 x 50, field of view 192 x 192 mm, and slice thickness = 6 mm, the code below will allow the user to 
  1) Define the muscle mask in slices 15-40;
  2) Return a mask of size 192 x 192 x 50; and
  3) Automatically visualize the result using fiber_visualizer, using slices 15, 25, 35, and 45 for anatomical reference and displaying the mask as a semi-transparent red structure (RGB scale = [1 0 0])

plot_options.plot_fibes=0;

plot_options.plot_mesh=0;

plot_options.plot_mask=1;

plot_options.anat_dims=[192 6];

plot_options.anat_slices=15:10:45;

plot_options.mask_size=[256 256];

plot_options.mask_dims=[192 6];

plot_options.mask_color=[1 0 0];

[mask, alt_mask] = define_muscle(anat_image, slices, alt_mask_size, plot_options);

VERSION INFORMATION

v. 0.5

ACKNOWLEDGMENTS

Grant support: NIH/NIAMS R01 AR050101, NIH/NIAMS R01 AR073831
