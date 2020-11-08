# Help for <i>signal2tensor2</i>, v. 1.0

## Introduction

This help file contains information about
1) [Purpose of the Program](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-signal2tensor2.md#1-purpose)
2) [Usage of the Program](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-signal2tensor2.md#2-usage)
3) [Syntax](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-signal2tensor2.md#3-Syntax)
4) [Example Code](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-signal2tensor2.md#4-Example-Code)
5) [Acknowledgements](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-signal2tensor2.md#5-Acknowledgements)

## 1. Purpose
<i>signal2tensor2</i> finds the tensor that best fits the observed signal.

[Back to the top](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-signal2tensor2.md)

## 2. Usage
The user inputs a column vector containing unweighted and diffusion-weighted signals, a matrix describing the diffusion-sensitizing directions, and the b-value (assumed to be the same for all directions). The function outputs the calculated tensor, with units of 1/[b].

[Back to the top](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-signal2tensor2.md)

## 3. Syntax

d_m = signal2tensor2(signal_v, dir_m, b)

### The input arguments are:
* <i>signal_v</i>: A column vector of image signals, the first element of which is the unweighted signal and the other elements of which correspond to the directions specified by the rows of <i>dir_m </i>

* <i>dir_m</i>: A matrix containing the X, Y, and Z components of unit-length vectors describing the diffusion-sensitizing directions; it has dimensions of (Number of Directions x 3)

b: The diffusion-weighting (b-) value

### The output argument is:
* <i>d_m</i>: The diffusion tensor

[Back to the top](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-signal2tensor2.md)

## 4. Example code

Given 

1.	A b-value of 450 s/mm<sup>2</sup>; 

2.	A 24-direction diffusion-encoding scheme; and

3.	A 4-dimensional set of DTI images (rows x columns x slices x diffusion-weighting directions, of which the first element is the non-diffusion weighted signal)

the code below will allow the user to:

1.	Define the diffusion-encoding direction scheme;

2.	Obtain the image signals at row 96, column 90, and slice 12; and

3.	Calculate the diffusion-tensor.

% Define diffusion-encoding variables:

b = 450;          %s/mm^2

dir_m = [
        0.1890    0.6940    0.6940

        0.6940    0.6940    0.1890
        
        0.1890    0.6940   -0.6940
        
        -0.6940   0.6940    0.1890
        
        -0.6940   0.1890   -0.6940
        
        0.6940    0.1890   -0.6940
        
        -0.6340   0.4420    0.6340
        
        0.6340    0.4420    0.6340
        
        0.6340    0.6340   -0.4420
        
        -0.4420   0.6340   -0.6340
        
        0.4420    0.6340    0.6340
        
        0.4420    0.6340   -0.6340
        
        -0.6940   0.6940   -0.1890
        
        0.6340    0.6340    0.4420
        
        -0.6940   0.1890    0.6940
        
        -0.6340   0.4420   -0.6340
        
        0.6940    0.1890    0.6940
        
        0.6340    0.4420   -0.6340
        
        -0.1890   0.6940   -0.6940
        
        -0.4420   0.6340    0.6340
        
        0.6940    0.6940   -0.1890
        
        -0.1890   0.6940    0.6940
        
        -0.6340   0.6340   -0.4420
        
        -0.6340   0.6340    0.4420];

% Get signals and call the function

signal_v = squeeze(dti_images(96,90,12,:);

d_m = signal2tensor2(signal_v, dir_m, b)

## 5. Acknowledgements
People: Adam Anderson

[Back to the top](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-signal2tensor2.md)
