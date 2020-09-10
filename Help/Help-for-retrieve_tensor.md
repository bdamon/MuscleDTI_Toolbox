# Help for the function <i>retrieve_tensor</i>, v. 0.1.x

## Introduction

This help file contains information about
1) [Purpose of the Program](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-retrieve_tensor.md#1-purpose)
2) [Usage of the Program](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-retrieve_tensor.md#2-usage)
3) [Syntax](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-retrieve_tensor.md#3-Syntax)
4) [Example Code](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-retrieve_tensor.md#4-Example-Code)
5) [Acknowledgements](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-retrieve_tensor.md#5-Acknowledgements)

## 1. Purpose

The function <i>retrieve_tensor</i> is used to retrieve the diffusion tensor from a matrix of diffusion tensors.

[Back to the top](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-retrieve_tensor.md)

## 2. Usage
The user does not interact with this function; it is called from within the <i>fiber_track</i> function.

[Back to the top](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-retrieve_tensor.md)

## 3. Syntax
d_tensor = retrieve_tensor(tensor_m, img_indx);

The input arguments are:

* <i>tensor_m</i>: A 5D matrix with the first-third dimensions matching the [rows columns slice] size of the DTI images and the fourth and fifth dimensions holding the 33 diffusion tensor at each voxel.

* <i>img_indx</i> The row, column, and slice indices into the tensor_m matrix.

The output argument is:

* <i>d_tensor</i>: The diffusion tensor.
   
[Back to the top](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-retrieve_tensor.md)

## 4. Example Code
An example is not given because the function is called from within <i>fiber_track</i>. The user does not interact with this function.

[Back to the top](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-retrieve_tensor.md)

## 5. Acknowledgements
Grants: NIH/NIAMS R01 AR073831

[Back to the top](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Help-for-retrieve_tensor.md)
