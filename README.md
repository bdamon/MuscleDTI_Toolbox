# MuscleDTI_Toolbox
A Matlab Toolbox for Skeletal Muscle Diffusion Tensor MRI Fiber Tractography 

The MuscleDTI_Toolbox consists of a series of custom-written Matlab functions for performing diffusion-tensor MRI fiber tractography in skeletal muscle. This README file contains
  1) A list of MATLAB requirements;
  2) A list of the conventions assumed regarding data acquisition;
  3) An overview of a typical workflow using the toolbox;
  4) Links to other resources in the toolbox and elsewhere online.

## Acknowledgements
This work was supported by NIH grants NIH/NIAMS R01 AR050101 and NIH/NIAMS R01 AR073831. By using this software, users agree to acknowledge the active grant (NIH/NIAMS R01 AR073831) in presentations and publications; to adhere to NIH policies regarding open access; and to abide by the terms of the GNU General Public License. This work reflects the collective contributions of many individuals, including: Adam Anderson, Amanda Buck, Crystal Coolbaugh, Bruce Damon, Zhaohua Ding, Hannah Kilpatrick, Anneriet Heemskerk, Melissa Hooijmans, and Justin Montenegro. Details regarding authorship and specific contributions are noted in each function.
  
## MATLAB Requirements

## Data Acquisition Conventions Assumed 
  1) Slices cover the entire muscle of interest and are numbered in ascending order from distal to proximal;
  2) The in-plane field of view is square; and
  3) The in-plane reconstructed matrix is square.

## Overview of a Typical Workflow

### define_muscle
The function define_muscle is used to define the boundary of a muscle and return the corresponding binary image mask. This mask is used in the fiber_track function of the MuscleDTI_Toolbox.  Follow [this link](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help%20for%20define_muscle) for detailed help on this function, including an instructional video.

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

