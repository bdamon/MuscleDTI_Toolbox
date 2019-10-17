# Help File for define_muscle

## Usage
The function define_muscle is used to define the of a muscle and return the corresponding binary image mask. This mask is used in the fiber_track function of the MuscleDTI_Toolbox. It may be visualized using the fiber_visualizer function.

The user must input a matrix containing the structural images and a two-element vector indicating the first and last slice of this stack the the user wishes to analyze. Three images are displayed: the current slice (middle panel), the preceding slice with its ROI (if present), and the next slice. Using the roipoly tool, the user defines the ROI in the middle slice. For the main figure window, an interactive tool is opened that allows the user to adjust the image's window and level settings. After closing the ROI, the program advances to the next slice until all slices have been defined. A mask of the same size as the anatomical image is returned.

If the user has included the optional input argument mask_size, then an alternativelh sized mask will also be returned.  If the user has inluded the optional input argument plot_options, the fiber_visualizer function will be called and the completed mask will be displayed.

A file named mask_file, containing the mask and (if present) the alternatively sized mask, is automatically saved in the working directory.

## Syntax
INPUT ARGUMENTS
anat_image: A row x column x slices stack of images
slices: A two element vector containing the first and last slices to be
analyzed
mask_size: If specified, this is a two element vector containing the row  
x column size of a second mask
plot_options: If specified, this calls the fiber_visualizer function to
plot the mask.

OUTPUT ARGUMENTS
mask: the binary image mask, with size matching that of the original
image
alt_mask: a second binary image mask, with size matching that of the
vector alt_mask_size


