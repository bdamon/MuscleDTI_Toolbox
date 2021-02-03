# MuscleDTI_Toolbox
## A MATLAB Toolbox for Skeletal Muscle Diffusion Tensor MRI Fiber Tractography 

The MuscleDTI_Toolbox consists of a series of custom-written MATLAB functions for performing diffusion-tensor MRI fiber tractography in skeletal muscle. This README file contains
  1) [Acknowledgements](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/README.md#1-acknowledgements)
  2) [License information](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/README.md#2-license-information)
  3) [A list of MATLAB requirements](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/README.md#3-matlab-requirements)
  4) [A list of the conventions assumed regarding data acquisition](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/README.md#4-data-acquisition-conventions-assumed)
  5) [An overview of a typical workflow using the toolbox](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/README.md#5-overview-of-a-typical-workflow)
  6) [Links to other resources in the toolbox and online](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/README.md#6-other-resources)

## 1. Acknowledgements
The functions in this toolbox reflect the collective contributions of many individuals over many years, including: Adam Anderson, Amanda Buck, Crystal Coolbaugh, Bruce Damon, Zhaohua Ding, Hannah Kilpatrick, Anneriet Heemskerk, Melissa Hooijmans, Drew Lansdown, and Justin Montenegro. Details regarding authorship and individual contributions are noted in each function.

This work was supported by NIH grants NIH/NIAMS R01 AR050101 and NIH/NIAMS R01 AR073831. By using this software, users agree to acknowledge the active grant (NIH/NIAMS R01 AR073831) in presentations and publications and to adhere to NIH policies regarding open access to their publications. 

## 2. License Information
This work is covered under a [GNU General Public License](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/LICENSE.md), v. 3 or later.

## 3. MATLAB Requirements
The functions have been tested using MATLAB v. 2019b.  The toolbox consists primarily of custom-written functions, but also calls some built-in MATLAB functions.

## 4. Data Acquisition Conventions Assumed 
In-plane Geometry
* In-plane field of view is square
* In-plane reconstructed imaging matrix size is square

Slice Geometry
* Slices were acquired in the axial anatomical plane
* Slices numbers ascend in the +Z direction 
* Slice geometries of structural and DT images match

Whole-image Properties
* Signal patterns of structural and DT images match
* The slices cover the entire muscle of interest
* Acquisition volumes of structural and DT images match

## 5. Overview of a Typical Workflow
Muscle DTI tractography includes pre-processing and fiber-tract processing steps, as elucidated below.

### A. Pre-processing
Before performing fiber tractography, several pre-processing steps must be performed:
* <i>File input</i>: Depending on the image format, this could be accomplished using the built-in MATLAB function <i>dicomread</i>, the built-in MATLAB function <i>niftiread</i>, or custom-written functions for proprietary image formats such as PAR/REC or XML/REC (Philips Medical Systems). 
* <i>Concatenation of multiple image acqusitions into a single dataset</i>: The need for this depends on the details of the user's image acquisition scheme. 
* <i>Image registration</i>: Switching of the diffusion-encoding gradients induces distortions in the images. These distortions are corrected using image registration; the example in this toolbox uses the Demons registration technique, called using the <i>imregdemons</i> function in MATLAB.
* <i>De-noising</i>: Some level of noise in the data is inevitable. This noise adds variability to the estimation of the diffusion tensor; at a sufficiently low signal-to-noise ratio, it can also add bias. The denoising method used here uses the custom-written function anisotropic smoothing method, [<i>aniso4d_smooth</i>](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Preprocessing-Functions/aniso4D_smoothing.m). Help on this function is available [here](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-signal2tensor2.md).
* <i>Estimation of the diffusion tensor throughout the muscle of interest</i>: The example given here uses a weighted least squares method to estimate the diffusion tensor that best matches the observed signals, given the diffusion-encoding matrix and diffusion encoding (b-) value. This is performed in the function [<i>signal2tensor2</i>](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Preprocessing-Functions/signal2tensor2.m), help for which is found [here](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-signal2tensor2.md).

Follow [this link](https://github.com/bdamon/MuscleDTI_Toolbox/tree/master/Sample-Scripts) for a MATLAB script and [this link](https://github.com/bdamon/MuscleDTI_Toolbox/tree/master/Preprocessing-Functions) for the custom-written MATLAB functions that perform these tasks. Other required functions are part of MATLAB's proprietary toolboxes and cannot be distributed here.

Switching of the phase-encoding gradients can induce distortions in the images. These can be corrected.  If the user wishes to perform this correction, they must have obtained a second set of non-diffusion weighted images with the phase-encoding direction reversed. The eddy-current correction scheme ilustrated here uses freeware called [FSL](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki). 

### B. Fiber-tract Processing
* Define muscle boundaries using the function <i>define_muscle</i>: Real muscle fibers are assumed to be contained entirely within a single muscle of interest. The <i>fiber_track</i> function therefore requires the user to input a binary image mask demarcating the muscle boundaries; this mask is used to prevent fiber tracts from exiting the muscle of interest. The function [<i>define_muscle</i>](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Tractography-Functions/define_muscle.m) is used to define this mask. Follow [this link](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-define_muscle.md) for detailed help on this function. Other programs, such as [ITK-SNAP](http://www.itksnap.org/pmwiki/pmwiki.php), can also be used.
* Define the seed points using the function <i>define_roi</i>: In DTI fiber-tracking, the tracts are propagated from a set of points, commonly called "seed points." In the MuscleDTI_Toolbox, the anatomical structure into which the muscle fibers insert (a flattened tendinous structure called an aponeurosis) is used to define these points. The function [<i>define_roi</i>](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Tractography-Functions/define_roi.m) is used to digitize the {row, column, slice} coordinates of the aponeurosis; these points are used to define the seed surface. Follow [this link](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-define_roi.md) for detailed help on this function.
* Generate the fiber tracts using the function <i>fiber_track</i>: Fiber tracts are propagated from the seed points by integrating the first eigenvector of the diffusion tensor. The function <i>fiber_track</i> is used for this purpose. The user can select from several propagation algorithms and several methods for determining when to stop propagating a tract. The major output of this function is a matrix containing the {row, column, slice} coordinates of each point along each fiber tract.  The function [<i>fiber_track</i>](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Tractography-Functions/fiber_track.m) calls the function [<i>retrieve_tensor</i>](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Tractography-Functions/retrieve_tensor.m), which finds the diffusion tensor in each voxel of the image. Follow [this link](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_track.md) for detailed help on this function.
* Smooth the fiber tracts using the function <i>fiber_smoother</i>: Fiber tract points are subject to errors in position because of the presence of noise and artifacts in the images. To mitigate these effects, the function [<i>fiber_smoother</i>](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Tractography-Functions/fiber_smoother.m) performs a polynomial fit to each fiber tract. This also allows the interpolation of the fiber tract positions at a resolution higher than the original tracts.  This step is not required, but is strongly recommended prior to calling the [<i>fiber_quantifier</i>](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Tractography-Functions/fiber_quantifier.m) function. Follow [this link](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_smoother.md) for detailed help on this function.
* Quantify the tracts' structural properties using the function <i>fiber_quantifier</i>: After the fiber tracts have been polynomial-fitted, their structural properties are quantified using the function [<i>fiber_quantifier</i>](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Tractography-Functions/fiber_quantifier.m).  The properties quantified include the pennation angle, curvature, and length. These properties are calculated in a pointwise manner along the fiber tracts. Follow [this link](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_quantifier.md) for detailed help on this function.
* Eliminate erroneous results using the function <i>fiber_goodness</i>: The quantitative results are examined and obviously wrong results are eliminated from the dataset. The function [<i>fiber_goodness</i>](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Tractography-Functions/fiber_goodness.m) eliminates tracts that ascend and descend (an error due to overfitting); that have architectural properties that exceed certain limits; and/or that vary greatly from their neighbors. A final dataset is calulated, including the mean properties for each tract and for the entire muscle. Follow [this link](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_goodness.md) for detailed help on this function.
* Visualize the results using the function <i>fiber_visualizer</i>: At several of these stages, the results can be visualized using the function [<i>fiber_visualizer</i>](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Tractography-Functions/fiber_visualizer.m). The user can select the mask, seed surface, and/or fiber tracts for display.  The user can also select which image slices to display for anatomical reference. Follow [this link](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-fiber_visualizer.md) for detailed help on this function.

## 6. Other Resources
### A. Within the toolbox:
* [Here's a link to the sample scripts](https://github.com/bdamon/MuscleDTI_Toolbox/tree/master/Sample-Scripts)
* [Here's a link to all of the preprocessing functions](https://github.com/bdamon/MuscleDTI_Toolbox/tree/master/Preprocessing-Functions)
* [Here's a link to all of the tractography functions](https://github.com/bdamon/MuscleDTI_Toolbox/tree/master/Tractography-Functions)
* [Here's a link to all of the help files](https://github.com/bdamon/MuscleDTI_Toolbox/tree/master/Help)
* [Here's a link to templates for submitting feature requests and bug reports](https://github.com/bdamon/MuscleDTI_Toolbox/tree/master/Issues)

### B. External to the toolbox:
Several recent reviews on muscle DTMRI include:
* [Skeletal muscle DT-MRI fiber tracking: Rationale, data acquisition and analysis methods, applications, and future directions](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5136336/)

* [Techniques and applications of skeletal muscle diffusion tensor imaging: A review](https://www.ncbi.nlm.nih.gov/pubmed/26221741)
