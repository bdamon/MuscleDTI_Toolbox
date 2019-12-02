function d_tensor = retrieve_tensor(tensor_m, img_indx, weight_matrix)
%
%FUNCTION retrieve_tensor
%  d_tensor = retrieve_tensor(tensor_m, img_indx, weight_matrix)
%
%USAGE
%    The function retrieve_tensor is used to retrive the diffusion tensor
%  from a 5-D (#Rows x #Columns x #Slices x 3 x 3) matrix of diffusion tensors.   
%  It is called by the fiber_track function in the MuscleDTI_Toolbox.
%    Future development of retrieve_tensor will include the ability to
%  describe the matrix as a continuous field of tensors, allowing improved
%  interpolation.
%
%INPUT ARGUMENTS
%  tensor_m: a 5D matrix containing rows, columns, slices, and the
%    diffusion tensor.
%
%  img_indx: the row, column, and slice indices into the tensor_m matrix.
%    This can either be a 1x3 vector or a 3x3x3 matrix of row, column, and
%    slice indices.
%
%  weight_matrix: if img_indx is 3x3x3, a matrix of weights must also be
%    specified.
%
%OUTPUT ARGUMENTS
%  d_tensor: the diffusion tensor
%
%VERSION INFORMATION
%  In beta-testing mode
%
%ACKNOWLEDGEMENTS
%  Grant support: NIH/NIAMS R01 AR073831


%% get the data

if nargin==2                                                                %if only tensor_m and img_indx are provided,
    
    img_r=img_indx(1);
    img_c=img_indx(2);
    img_s=img_indx(3);
    d_tensor = squeeze(tensor_m(img_r, img_c, img_s, :, :));                %get the tensor from teh indicated voxel
    
elseif nargin==3
    
    local_tensors = zeros(3,3,3,3,3);                                       %initialize a matrix to hold the local tensors
    
    for r = 1:3
        for c = 1:3
            for s = 1:3
                local_tensors(r, c, s, :, :) = ...                %get out each tensor, multiplying by its weight
                    squeeze(tensor_m(img_indx(r), img_indx(c), img_indx(s), :, :)) .* weight_matrix(r, c, s);
            end
        end
    end
    
    d_tensor = squeeze(sum(sum(sum(local_tensors))))./sum(sum(sum(weight_matrix)));
end


