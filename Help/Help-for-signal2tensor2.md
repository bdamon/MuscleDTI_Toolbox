# Help for the function <i>signal2tensor2</i>, v. 1.0.1

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

* b: The diffusion-weighting (b-) value

### The output argument is:
* <i>d_m</i>: The diffusion tensor

[Back to the top](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-signal2tensor2.md)

## 4. Example code

Given 

1. Diffusion-weighted images stored in the variable dti_images (A 4-dimensional set of images (rows x columns x slices x diffusion-weighting directions, of which the first element of dimension 4 is the non-diffusion weighted signal and the remaining elements correspond to the diffusion-weighting directions specified in <i>dir_m</i>); 
2. The row, column, and slice coordinates <i>dw_row</i>, <i>dw_col</i>, and <i>dw_sl</i>;
3. The diffusion-encoding vectors stored in the variable <i>dir_m</i>; and
4. A b-value of 450 s/mm<sup>2</sup>, specified in the variable <i>b</i>; 

the code below will allow the user to:

1. Calculate the diffusion-tensor.

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

signal_v = squeeze(dti_images(dw_row,dw_col,dw_slc,:));

d_m = signal2tensor2(signal_v, dir_m, b);

## 5. Acknowledgements
People: Adam Anderson

[Back to the top](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-signal2tensor2.md)
