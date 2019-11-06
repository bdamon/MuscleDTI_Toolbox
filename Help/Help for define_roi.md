# Help for define_roi, v. 0.5

## Usage

The function define_roi is used to digitize the aponeurosis of muscle fiber insertion in the MuscleDTI_Toolbox.  The digitized points are used to reconstruct a mesh; the mesh is used as the seed surface for fiber tracking.

There are two options for defining the aponeurosis:
1) Manual: The user is prompted initially to select two points, which define the level of zoom to be used throughout the entire process. Then the user advances through the slices to select the aponeurosis. The mask is provided as a guide; the aponeurosis must fall within its boundaries. The selected points can form a line or close to form a polygon. At each slice, the user is given the option of repeating the procedure in case of error.  For each figure window, an interactive tool is opened that allows the user to adjust the image's window and level settings.
2) Automatic: The aponeurosis is automatically segmented from within the region of the image represented by the muscle mask. TWo segmentation methods (edge detection and k-means clustering) are used, and the segmented region is defined as the consensus of the two results. The region is displayed and the user is allowed to correct misassignments. The boundaries of the segmented region are smoothed using a Savitsky-Golay filter and used to form the mesh. 

The mesh is initially formed with resolution n_row x n_col.  To smooth the mesh, it is then downsampled by a size factor of four. Finally, the smoothed mesh is used to create a high resolution mesh at the desired size. A file called roi_mesh_file.mat is automatically saved in the working directory. If the input argument plot_options is included, the mesh and mask are plotted using the function fiber_visualizer.

## Input Arguments
anat_image: The imaging data.  If input as a structure, then the imaging data are assumed to exist in a field called anat_image.Data.  If specified as a matrix, the data are used directly.

mask: The mask, as defined by the function define_mask

defroi_options: A structure containing the following fields:

  -slices: A two-element vector containing the first and last slices that the user wishes to digitize.
  
  -dti_size: The size of the DTI image dataset (rows x columns x slices), input as a three element vector.
  
  -mesh_size: A two-element vector containing the numbers of rows (n_row) and columns (n_col) desired in the output mesh.
  
  -method: a string variable set either to 'manual' or 'auto'. The manual and automatic options are described above.

plot_options: Optional. If specified, this calls the fiber_visualizer function to plot the mask and roi mesh.

## Output Arguments
roi_mesh: a 3D matrix containing the reconstructed mesh with size rows x columns x 6. In the 3rd dimension, levels 1-3 hold the row-column-slice data and levels 4-6 hold the normal vector to the roi surface at the point {row column slice}.
   
   
## Acknowledgements
 People: Zhaohua Ding
 Grant support: NIH/NIAMS R01 AR050101, NIH/NIAMS R01 AR073831
