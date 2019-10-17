# MuscleDTI_Toolbox
A Matlab Toolbox for Skeletal Muscle Diffusion Tensor MRI Fiber Tractography 

The MuscleDTI_Toolbox consists of a series of custom-written Matlab functions for performing diffusion-tensor MRI fiber tractography in skeletal muscle. This README file contains
  1) Acknowledgements;
  2) A link to the license
  3) A list of MATLAB requirements;
  4) A list of the conventions assumed regarding data acquisition;
  5) An overview of a typical workflow using the toolbox;
  6) Links to other resources in the toolbox and elsewhere online.

## 1. Acknowledgements
This work was supported by NIH grants NIH/NIAMS R01 AR050101 and NIH/NIAMS R01 AR073831. By using this software, users agree to acknowledge the active grant (NIH/NIAMS R01 AR073831) in presentations and publications and to adhere to NIH policies regarding open access. This work reflects the collective contributions of many individuals, including: Adam Anderson, Amanda Buck, Crystal Coolbaugh, Bruce Damon, Zhaohua Ding, Hannah Kilpatrick, Anneriet Heemskerk, Melissa Hooijmans, and Justin Montenegro. Details regarding authorship and individual contributions are noted in each function.

## 2. License Information
This work is covered under a [GNU General Public License](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/LICENSE.md), v. 3 or later.

## 3. MATLAB Requirements
The functions have been tested using MATLAB v. 2019b.  The toolbox consists primarily of custom-written functions, but also calls functions contained within MATLAB'S Image Processing Toolbox.

## 4. Data Acquisition Conventions Assumed 
  1) Slices cover the entire muscle of interest and are numbered in ascending order from distal to proximal;
  2) The in-plane field of view is square; and
  3) The in-plane reconstructed matrix is square.

## 5. Overview of a Typical Workflow

### A. Pre-processing


### B. Define muscle boundaries using the function define_muscle
Real muscle fibers are assumed to be contained entirely within a single muscle of interest. The fiber_track function therefore requires the user to input a binary image mask demarcating the muscle boundaries; this mask is used to prevent fiber tracts from exiting the muscle of interest. The function define_muscle is used to define this mask. Follow this link for detailed help on this function, including an instructional video.

### C. Define the initial using the function define_roi
Fiber tracts are propagated from a set of points, commonly called "seed points." In the MuscleDTI_Toolbox, the tendinous structure into which the muscle fibers insert is used to define these points. The function define_roi is used to digitize the {row, column, slice} coordinates of the tendon; these points are used to define the seed surface. Follow this link for detailed help on this function, including an instructional video.

### D. Generate the fiber tracts using the function fiber_track
Fiber tracts are propagated from the seed points by following the direction indicated by the first eigenvector of the diffusion tensor. The function fiber_track is used to perform this integration. The user can select from several propagation algorithms and several methods for determining when to stop propagating a tract. The major output of this function is a matrix containing the {row, column, slice} coordinates of each fiber tract.  fiber_track calls the function retrieve_tensor, which finds the diffusion tensor in each voxel of the image.

### E. Smooth the fiber tracts using the function fiber_smoother
Fiber tract points are subject to small errors in position because of the presence of noise and artifacts in the images. To mitigate these effects, the function fiber_smoother performs a polynomial fit to each fiber tract. This also allows the interpolation of the fiber tract positions at a resolution higher than the original tracts.  This step is not required, but is strongly recommended prior to calling the fiber_quantifier function.

### F. Quantify the tracts' structural properties using the function fiber_quantifier
After the fiber tracts have been polynomial-fitted, their structural properties are quantified using the function fiber_quantifier.  The properties quantified include the pennation angle, curvature, and length.  

### fiber_selector

### fiber_visualizer

## 6. Other Resources
