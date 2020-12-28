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
%VERSION INFORMATION
%  v 1.0 (initial release), 28 Dec 2020, Bruce Damon
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
