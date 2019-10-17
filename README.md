# MuscleDTI_Toolbox
## A Matlab Toolbox for Skeletal Muscle Diffusion Tensor MRI Fiber Tractography 

The MuscleDTI_Toolbox consists of a series of custom-written Matlab functions for performing diffusion-tensor MRI fiber tractography in skeletal muscle. The toolbox consists of:
  1) Conventions assumed regarding data acquisition;
  2) A sample script for pre-processing the data;
  3) Functions for generating and quantifying the fiber tracts and a sample script that calls the functions.
Further information about each of the scripts/functions, including video tutorials, appears below.

This work was supported by NIH grants NIH/NIAMS R01 AR050101 and NIH/NIAMS R01 AR073831. By using this software, users agree to acknowledge the active grant (NIH/NIAMS R01 AR073831) in presentations and publications and to adhere to NIH policies regarding open access. Most of the functions were written by Bruce Damon, but reflect the collective contributions of many individuals, including: Adam Anderson, Amanda Buck, Crystal Coolbaugh, Zhaohua Ding, Hannah Kilpatrick, Anneriet Heemskerk, Melissa Hooijmans, and Justin Montenegro. 
  
## Conventions Assumed Regarding Data Acquisition
  1) Slices cover the entire muscle of interest and are numbered from distal to proximal;
  2) The in-plane field of view is square; and
  3) The in-plane reconstructed matrix is square.

## Functions and Sample Script for Pre-processing the Data


## Fiber-tracking Functions

### define_muscle
The function define_muscle is used to define the boundary of a muscle and return the corresponding binary image mask. This mask is used in the fiber_track function of the MuscleDTI_Toolbox.

Three figure windows are open, each displaying a different image: the current slice (middle panel), the preceding slice with its ROI (if present), and the next slice. Using the roipoly tool, the user defines the ROI in the middle slice. For the main figure window, an interactive tool is opened that allows the user to adjust the image's window and level settings. After closing the ROI, the program advances to the next slice until all slices have been defined.

A file named mask_file, containing the mask and (if present) the alternatively sized mask, is automatically saved in the working directory.

### define_roi

### fiber_track

### retrieve_tensor

### fiber_smoother
The function fiber_smoother is used to smooth fiber tracts and increase the spatial resolution of fiber tracts generated using the MuscleDTI_Toolbox. It is recommended for use following fiber_track and prior to calling the fiber_quantifier function.

The x, y, and z positions are fitted to Nth order polynomials as functions of distance along the tract and are uniformly solved at interpolation distances of interpolation_step. The user selects the polynomial order.  This procedure was originally described in Damon et al, Magn Reson Imaging, 2012.  In the MuscleDTI_Toolbox, it has been modified to: 
  1) Fit the tract positions as functions of distance rather than point number and 
  2) Allow selection of the polynomial order.  

The former option is required for tracking algorithms that use variable step sizes,such as FACT.

### fiber_quantifier

### fiber_selector

### fiber_visualizer

