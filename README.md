# MuscleDTI_Toolbox
## A Matlab Toolbox for Skeletal Muscle Diffusion Tensor MRI Fiber Tractography 

The MuscleDTI_Toolbox consists of a series of custom-written Matlab functions for performing diffusion-tensor MRI fiber tractography in skeletal muscle. This README file contains
  1) [Acknowledgements](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/README.md#1-acknowledgements)
  2) [Information about the license](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/README.md#2-license-information)
  3) [A list of MATLAB requirements](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/README.md#3-matlab-requirements)
  4) [A list of the conventions assumed regarding data acquisition](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/README.md#4-data-acquisition-conventions-assumed)
  5) [An overview of a typical workflow using the toolbox](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/README.md#5-overview-of-a-typical-workflow)
  6) [Links to other resources in the toolbox and online](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/README.md#6-other-resources)

## 1. Acknowledgements
This work reflects the collective contributions of many individuals, including: Adam Anderson, Amanda Buck, Crystal Coolbaugh, Bruce Damon, Zhaohua Ding, Hannah Kilpatrick, Anneriet Heemskerk, Melissa Hooijmans, and Justin Montenegro. Details regarding authorship and individual contributions are noted in each function.

This work was supported by NIH grants NIH/NIAMS R01 AR050101 and NIH/NIAMS R01 AR073831. By using this software, users agree to acknowledge the active grant (NIH/NIAMS R01 AR073831) in presentations and publications and to adhere to NIH policies regarding open access to their publications. 

## 2. License Information
This work is covered under a [GNU General Public License](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/LICENSE.md), v. 3 or later.

## 3. MATLAB Requirements
The functions have been tested using MATLAB v. 2019b.  The toolbox consists primarily of custom-written functions, but also calls functions contained within MATLAB's Image Processing Toolbox.

## 4. Data Acquisition Conventions Assumed 
  1) Slices cover the entire muscle of interest and are numbered in ascending order from distal to proximal;
  2) The in-plane field of view is square; and
  3) The in-plane reconstructed matrix is square.

## 5. Overview of a Typical Workflow

### A. Pre-processing
Before performing fiber tractography, several pre-processing steps muscle be performed.  These may include:
  1) File input;
  2) Correction of eddy current-induced distortions in the images;
  3) Concatenation of multiple image acqusitions into a single dataset;
  4) Image registration;
  5) De-noising; and
  6) Calculation of the diffusion tensor throughout the muscle of interest.

Follow this link for help on these steps, including a MATLAB script that performs many of these tasks.

### B. Define muscle boundaries using the function <i>define_muscle</i>
Real muscle fibers are assumed to be contained entirely within a single muscle of interest. The fiber_track function therefore requires the user to input a binary image mask demarcating the muscle boundaries; this mask is used to prevent fiber tracts from exiting the muscle of interest. The function <i>define_muscle</i> is used to define this mask. Follow [this link](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help%20for%20define_muscle.md) for detailed help on this function, including an instructional video.

### C. Define the seed points using the function <i>define_roi</i>
Fiber tracts are propagated from a set of points, commonly called "seed points." In the MuscleDTI_Toolbox, the structure into which the muscle fibers insert (a flattened tendinous structure called an aponeurosis) is used to define these points. The function <i>define_roi</i> is used to digitize the {row, column, slice} coordinates of the aponeurosis; these points are used to define the seed surface. Follow [this link](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help%20for%20define_roi.md) for detailed help on this function, including an instructional video.

### D. Generate the fiber tracts using the function <i>fiber_track</i>
Fiber tracts are propagated from the seed points by following the direction indicated by the first eigenvector of the diffusion tensor. The function <i>fiber_track</i> is used to perform this integration. The user can select from several propagation algorithms and several methods for determining when to stop propagating a tract. The major output of this function is a matrix containing the {row, column, slice} coordinates of each point along each fiber tract.  <i>fiber_track</i> calls the function <i>retrieve_tensor</>, which finds the diffusion tensor in each voxel of the image. Follow [this link](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help%20for%20fiber_track.md) for detailed help on this function.

### E. Smooth the fiber tracts using the function <i>fiber_smoother</i>
Fiber tract points are subject to errors in position because of the presence of noise and artifacts in the images. To mitigate these effects, the function <i>fiber_smoother</i> performs a polynomial fit to each fiber tract. This also allows the interpolation of the fiber tract positions at a resolution higher than the original tracts.  This step is not required, but is strongly recommended prior to calling the <i>fiber_quantifier</i> function. Follow this link for detailed help on this function.

### F. Quantify the tracts' structural properties using the function <i>fiber_quantifier</i>
After the fiber tracts have been polynomial-fitted, their structural properties are quantified using the function <i>fiber_quantifier</i>.  The properties quantified include the pennation angle, curvature, and length. These properties are calculated in a pointwise manner along the fiber tracts. Follow this link for detailed help on this function.

### G. Eliminate erroneuous results using the function <i>fiber_selector</i>
Finally, the quantitative results are examined and obviously wrong results are eliminated from the dataset. The function <i>fiber_selector</i> eliminates tracts that ascend and descend (an error due to overfitting); that have architectural properties that exceed certain limits; or that vary greatly from their neighbors. A final dataset is calulated, including the mean properties for each tract and for the entire muscle. Follow this link for detailed help on this function.

### H. Visualize the results using the function <i>fiber_visualizer</i>
At any stage, the results can be visualized using the function <i>fiber_visualizer</i>. The user can select the mask, seed surface, and/or fiber tracts for display.  The user can also select which image slices to display for anatomical reference.

## 6. Other Resources
Several recent reviews on muscle DTMRI include:

[Skeletal muscle DT-MRI fiber tracking: Rationale, data acquisition and analysis methods, applications, and future directions](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5136336/)

[Techniques and applications of skeletal muscle diffusion tensor imaging: A review](https://www.ncbi.nlm.nih.gov/pubmed/26221741)
