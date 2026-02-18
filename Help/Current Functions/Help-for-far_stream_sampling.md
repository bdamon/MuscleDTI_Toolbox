# Help for the function [<i>far_stream_sampling</i>](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Tractography-Functions/far_stream_sampling.m), v. 1.0.0

## Introduction

This help file contains information about
1) [Purpose of the Program](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Current%20Functions/Help-for-fiber_visualizer.md#1-Purpose)
2) [Usage of the Program](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Current%20Functions/Help-for-fiber_visualizer.md#2-Usage)
3) [Syntax](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Current%20Functions/Help-for-fiber_visualizer.md#3-Syntax)
5) [Example Code](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Current%20Functions/Help-for-fiber_visualizer.md#4-Example-Code)
5) [Acknowledgements](https://github.com/bdamon/MuscleDTI_Toolbox/blob/master/Help/Current%20Functions/Help-for-fiber_visualizer.md#5-Acknowledgements)

## 1. Purpose
The function <i>far_stream_sampling</i> is used to uniformly sample fiber-tracts using the farthest streamline sampling method proposed by Li et al., IEEE Trans Biomed Eng, 2024.

## 2. Usage
The user calls <i>far_stream_sampling</i> from the command line. The user must supply a fiber tract matrix having units of mm and specifiy the final number of tracts they wish to preserve (K) and the number of points to use when resampling the tracts (N) during the sampling process. The function preserves the set of K tracts that have the maximum total inter-tract distance, thereby maximizing the sampling uniformity across the muscle.

## 3. Syntax 
[sampled_fiber_all_mm, sampled_fiber_all, sampled_fiber_all_idx] = far_stream_sampling(fiber_all_mm, fiber_all, sample_size, resampled_points)

The input arguments are:

<i>fiber_all_mm</i>: The fiber tracts. The rows and columns correspond to locations on the roi_mesh or the different seeding locations if dimension 2 is equal to 1. Dimension 3 gives point numbers on the tract, and the fourth dimension has row, column, and slice coordinates, with units of mm.

<i>fiber_all</i>: The fiber tracts. The rows and columns correspond to locations on the roi_mesh or the different seeding locations if dimension 2 is equal to 1. Dimension 3 gives point numbers on the tract, and the fourth dimension has row, column, and slice coordinates, with units of voxels.

<i>sample_size</i>: The number of fiber-tracts contained in the output fiber tract matrix
     
 <i>resampled_points</i>: Number of points to which the fiber-tracts will be resampled prior to the computation of the resampling distance; the default is 12

The output arguments are:

 <i>sampled_fiber_all_mm</i>: The matrix including the K sampled fiber tracts, with units of mm.

  <i>sampled_fiber_all</i>: The matrix including the K sampled fiber tracts, with units of voxels.

 <i>sampled_fiber_all_idx</i>: The indices of the sampled fiber-tracts into the original fiber tract matrix.  This is a  1D array for voxel-seeded tracts or a 2D array for roi mesh-seeded tracts.

## 4. Example Code

Given:

1.	The fiber tracking matrices, with point coordinates specified in mm and voxels

the following code will allow the user to 

1.	Resample each fiber tract to 15 points; and

2.	Select the K=500 fiber tracts that most uniformly sample the muscle's volume; and
  
3.	Return the K resampled tracts and their indies into the original fiber tract matrix.
  
% Call the function:

sample_size = 500;

resampled_points = 15;

[sampled_fiber_all_mm, sampled_fiber_all, sampled_fiber_all_idx] = far_stream_sampling(fiber_all_mm, fiber_all, sample_size, resampled_points)
 
## 5. Acknowledgements
People: Roberto Pineda Guzman

Grant support: NIH/NIAMS R01 AR073831
