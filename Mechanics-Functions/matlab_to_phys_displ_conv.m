function dcm_displ_ortho = matlab_to_phys_displ_conv(matlab_dcm_displ,spacing)
%
%FUNCTION matlab_to_phys_displ_conv
% The function converts the MATLAB's imwarp displacement field input to a
% displacement field in physical units
%
%INPUT ARGUMENTS
% matlab_dcm_displ: displacement field that can be used to warp an image
%  opened with dicomread() using the imwarp() function
% spacing: a 3 x 1 array that includes the physical spacing between the
%  elements of the array in the x,y, and z directions. If the spacing is
%  not specified, an isotropic spacing of 1x1x1 units is assumed
%
%OUTPUT ARGUMENTS
% dcm_displ_ortho: displacement field defined in an orthonormal
%  basis, similar to matlab_dcm_displ but with displacement values in 
%  scaled units and multiplied by -1 (forward mapping)
%
%VERSION INFORMATION
% v1.0.0, Roberto Pineda Guzman, October 7 2024.
%
% ACKNOWLEDGMENTS
% Grant support: NIH/NIAMS R01 AR073831

%% 

if nargin<2
    spacing=[1 1 1]; %Assume isotropic voxels unless previously specified
end

displacement_x_dicom=squeeze(matlab_dcm_displ(:,:,:,1));
displacement_y_dicom=squeeze(matlab_dcm_displ(:,:,:,2));
displacement_z_dicom=squeeze(matlab_dcm_displ(:,:,:,3));

dcm_displ_ortho=-1*cat(4,displacement_x_dicom*spacing(1),...
    displacement_y_dicom*spacing(2),...
    displacement_z_dicom*spacing(3));

