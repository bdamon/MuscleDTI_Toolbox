# Help for the function [<i>aniso4d_smoothing</i>](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Tractography-Functions/aniso4d_smoothing.m), v. 1.0.1

## Introduction

This help file contains information about
1) [Purpose of the Program](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-aniso4d_smoothing.md#1-Purpose)
2) [Usage of the Program](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-aniso4d_smoothing.md#2-Usage)
3) [Syntax](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-aniso4d_smoothing.md#3-Syntax)
4) [Example Code](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-aniso4d_smoothing.md#4-Example-Code)
5) [Acknowledgements](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-aniso4d_smoothing.md#5-Acknowledgements)


## 1. Purpose

The function <i>aniso4d_smoothing</i> is used to smooth diffusion-weighted images.

[Back to the Top](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-aniso4d_smoothing.md)

## 2. Usage
The user provides the diffusion weighted images and options that control the degree of smoothing. The images are smoothed using the method described by Ding et al. and Xu et al., which is based on the following partial differential equation:

  δI/δt=div(T∙∇I)

where I is the image intensity, <b>T</b> is a structure tensor that provides smoothing anisotropy, and δt is the iteration time parameter. <b>T</b> is the normalized inverse of the gradient tensor, <b>G</b>. <b>G</b> is the convolution of the outer product of ∇I with a Gaussian kernel having a standard deviation ρ, where the gradient ∇I is estimated from the diffusion-weighted images convolved with a Gaussian kernel with a variance σ.  The anisotropic smoothing method uses a common definition of <b>G</b> for all diffusion-weighting directions, allowing boundary information that is missing in one weighting direction to be captured from another weighting direction. Thus, the algorithm provides isotropic smoothing inside of structures and anisotropic smoothing at their boundaries. The smoothed images are returned.

[Back to the Top](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-aniso4d_smoothing.md)

## 3. Syntax

The function <i>aniso4d_smoothing</i> is called using the following syntax:

[image_smoothed = aniso4D_smoothing(image_unsmoothed, sigma, rho, delta_t, sr, type, isnormg, isfasteig);

The input arguments are:

* <i>image_unsmoothed</i>: The images to be smoothed;

* <i>sigma, rho</i>: The variances of the Gaussian kernel used to smooth the original images and the structure tensor, respectively;

* <i>delta_t</i>: The time step;

* <i>sr</i>: The spatial resolution of the images, input as a two-element vector containing the in-plane resolution and the slice thickness;

* <i>type</i>: The type of smoothing that is used, input as a string variable. The options include ‘explicit’, ‘implicit_multisplitting’, ‘implicit_AOS’, ‘implicit_ADI’, and ‘implicit_AOS_improved’;

* <i>isnormg</i>: Determines whether the gradient in the partial differential equation is normalized; set to either true or false; and

* <i>isfasteig</i>: Determines whether the fast eigenvector solver is used; set to either true or false.

The output argument is:

* <i>image_smoothed</i>: The smoothed images.

[Back to the Top](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-aniso4d_smoothing.md)
 

## 4. Example Code

### Example 1:

Given 

1.	Diffusion-weighted images with a spatial resolution of 1 x 1 x 7 

the code below will allow the user to:

1.	Smooth the images using the parameters specified in Buck et al. (PLoS, 2015, 10(5):e0126953):

%set smoothing options

noise = 5;

sigma = noise/100;

rho = 2*sigma;

delta_t = noise*3/44;

schemetype = 'implicit_multisplitting';

isfasteig  = true;

isnormg = false;

dti_res = [1 7];

dti_all_reg_smooth = aniso4D_smoothing(dti_all_reg, sigma, rho, delta_t, dti_res, schemetype, isnormg, isfasteig);
 
[Back to the Top](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-aniso4d_smoothing.md)

## 5. Acknowledgements

Grant support: NIH/NIBIB R01 EB000461, NIH/NIBIB R01 EB02777

[Back to the Top](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-aniso4d_smoothing.md)
 
