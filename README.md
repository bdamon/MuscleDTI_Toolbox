# MuscleDTI_Toolbox
A Matlab Toolbox for Skeletal Muscle Diffusion Tensor MRI Fiber Tractography 

The MuscleDTI_Toolbox consists of a series of custom-written Matlab functions for performing diffusion-tensor MRI fiber tractography in skeletal muscle. This README file contains
  1) A list of MATLAB requirements;
  2) A list of the conventions assumed regarding data acquisition;
  3) An overview of a typical workflow using the toolbox;
  4) Links to other resources in the toolbox and elsewhere online.

## Acknowledgements
This work was supported by NIH grants NIH/NIAMS R01 AR050101 and NIH/NIAMS R01 AR073831. By using this software, users agree to acknowledge the active grant (NIH/NIAMS R01 AR073831) in presentations and publications and to adhere to NIH policies regarding open access. This work reflects the collective contributions of many individuals, including: Adam Anderson, Amanda Buck, Crystal Coolbaugh, Bruce Damon, Zhaohua Ding, Hannah Kilpatrick, Anneriet Heemskerk, Melissa Hooijmans, and Justin Montenegro. Details regarding authorship and individual contributions are noted in each function.

## License Information
This work is covered under a [GNU General Public License](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/LICENSE.md), v. 3 or later.

## MATLAB Requirements


## Data Acquisition Conventions Assumed 
  1) Slices cover the entire muscle of interest and are numbered in ascending order from distal to proximal;
  2) The in-plane field of view is square; and
  3) The in-plane reconstructed matrix is square.

## Overview of a Typical Workflow

### Pre-processing


### Define muscle boundaries using the function define_muscle
Real muscle fibers, and consequently the DTMRI fiber tracts, are assumed to be contained entirely within a single muscle of interest. The fiber_track function therefore requires the user to input a binary image mask demarcating the muscle boundaries. The function define_muscle is used to define this mask. Follow this link for detailed help on this function, including an instructional video.

### define_roi
Fiber tracts are propagated from a set of points, commonly called "seed points." In the MuscleDTI_Toolbox, the tendinous structure into which the muscle fibers insert is used to define these points. The function define_roi is used to digitize the {row, column, slice} coordinates of the tendon; these points are used to define the seed surface. Follow this link for detailed help on this function, including an instructional video.

### fiber_track

### retrieve_tensor

### fiber_smoother
The function fiber_smoother is used to smooth fiber tracts and increase the spatial resolution of fiber tracts generated using the MuscleDTI_Toolbox. It is recommended for use following fiber_track and prior to calling the fiber_quantifier function.

### fiber_quantifier

### fiber_selector

### fiber_visualizer

