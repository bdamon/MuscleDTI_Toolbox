function det_J = jacobian(ux,uy,uz,spacing)
%
%FUNCTION jacobian
% The function jacobian is used to find the determinat of the jacobian of
% a displacement field. This function takes into account MATLAB's
% cootdinate system orthonormal basis used to analyze mechanical data
% in 3D.
%
%INPUT ARGUMENTS
% ux: a 3D array containing the x-component of a displacement field, 
%  defined by positive x in the increasing column direction
% uy: a 3D array containing the y-component of a displacement field, 
%  defined by positive y in the increasing row direction 
% uz: a 3D array containing the z-component of a displacement field, 
%  defined by positive z in the increasing slice direction (increasing 
%  page direction)
% spacing: a 3 x 1 array that includes the physical spacing between the
%  elements of the array in the x,y, and z directions. If the spacing is
%  not specified, an isotropic spacing of 1x1x1 units is assumed
%
%VERSION INFORMATION
% v1.0.0, Roberto Pineda Guzman, October 7 2024.
%
% ACKNOWLEDGMENTS
% Most of the code was taken from Herve Lombaert's Diffeomorphic Log Demons 
% Image Registration 'jacobian' function 
% https://www.mathworks.com/matlabcentral/fileexchange/39194-diffeomorphic-log-demons-image-registration
%
%Grant support: NIH/NIAMS R01 AR073831

%%    
if nargin<4
    spacing=[1 1 1]; %Assume isotropic voxels unless previously specified
end

% Gradients
[gx_x,gx_y,gx_z] = gradient(ux,spacing(1),spacing(2),spacing(3));
[gy_x,gy_y,gy_z] = gradient(uy,spacing(1),spacing(2),spacing(3));
[gz_x,gz_y,gz_z] = gradient(uz,spacing(1),spacing(2),spacing(3));

% Add identity
gx_x = gx_x + 1;
gy_y = gy_y + 1;
gz_z = gz_z + 1;

% Determinant
det_J = gx_x.*gy_y.*gz_z + ...
        gy_x.*gz_y.*gx_z + ...
        gz_x.*gx_y.*gy_z - ...
        gz_x.*gy_y.*gx_z - ...
        gy_x.*gx_y.*gz_z - ...
        gx_x.*gz_y.*gy_z;
end
