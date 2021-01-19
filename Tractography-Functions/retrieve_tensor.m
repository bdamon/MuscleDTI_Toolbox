function d_tensor = retrieve_tensor(tensor_m, img_indx)
%
%FUNCTION retrieve_tensor
%  d_tensor = retrieve_tensor(tensor_m, img_indx);
%
%USAGE
%  The function retrieve_tensor is used to retrieve the diffusion tensor
%  from a 5-D (#Rows x #Columns x #Slices x 3 x 3) matrix of diffusion tensors.
%  It is called by the fiber_track function in the MuscleDTI_Toolbox.
%    
%  Future development of retrieve_tensor will include the ability to
%  describe the matrix as a continuous field of tensors, allowing improved
%  interpolation.
%
%INPUT ARGUMENTS
%  tensor_m: a 5D matrix containing rows, columns, slices, and the
%    diffusion tensor.
%
%  img_indx: the row, column, and slice indices into the tensor_m matrix.
%
%OUTPUT ARGUMENTS
%  d_tensor: the diffusion tensor
%
%OTHER FUNCTIONS IN THE MUSCLE DTI FIBER-TRACKING TOOLBOX
%  For help with anisotropic smoothing, see <a href="matlab: help aniso4D_smoothing">aniso4D_smoothing</a>.
%  For help defining the muscle mask, see <a href="matlab: help define_muscle">define_muscle</a>.
%  For help defining the aponeurosis ROI, see <a href="matlab: help define_roi">define_roi</a>.
%  For help with fiber tracking, see <a href="matlab: help fiber_track">fiber_track</a>.
%  For help smoothing fiber tracts, see <a href="matlab: help fiber_smoother">fiber_smoother</a>.
%  For help quantifying fiber tracts, see <a href="matlab: help fiber_quantifier">fiber_quantifier</a>.
%  For help selecting fiber tracts following their quantification, see <a href="matlab: help fiber_goodness">fiber_goodness</a>.
%  For help visualizing fiber tracts and other structures, see <a href="matlab: help fiber_visualizer">fiber_visualizer</a>.
%
%VERSION INFORMATION
%  v. 1.0.0 (initial release), 15 Jan 2021, Bruce Damon
%
%ACKNOWLEDGEMENTS
%  Grant support: NIH/NIAMS R01 AR073831


%% get the data


img_r=img_indx(1);
img_c=img_indx(2);
img_s=img_indx(3);
d_tensor = squeeze(tensor_m(img_r, img_c, img_s, :, :));                %get the tensor from teh indicated voxel

%% end the function
return;