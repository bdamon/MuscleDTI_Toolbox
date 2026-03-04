function [F_deform_grad,C_cauchy_strain,E_green_strain] = ...
    strain_tensors(ux,uy,uz,spacing)
%%
%FUNCTION strain_tensors
% The function strain_tensors is used to compute the deformation gradient
% (F), right Cauchy-Green strain tensor (C), and Green Strain tensor (E) from a
% 3D displacement field u = u(x,y,z). 
% The deformation gradient F is defined as: F = grad(u) + I, where I is the
% identity tensor.
% The Cauchy-Green strain tensor (C) is defined as C = F'F.
% The Green strain tensor is defined as E = 0.5*(C-I).
% This function uses MATLAB's image coordinate system.
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
%OUTPUT ARGUMENTS
% F_deform_grad: Deformation gradient a 3D displacement field
% C_cauchy_strain: right Cauchy-Green strain tensor of a 3D displacement field 
% E_green_strain: Green strain tensor of a 3D displacement field 
%
%VERSION INFORMATION
% v1.0.0, Roberto Pineda Guzman, October 7 2024.
%
% ACKNOWLEDGMENTS
% Grant support: NIH/NIAMS R01 AR073831

%% 

if nargin<4
    spacing=[1 1 1]; %Assume isotropic voxels unless previously specified
end

%Compute displacement gradients
[dux_dx,dux_dy,dux_dz] = gradient(ux,spacing(1),spacing(2),spacing(3)); 
[duy_dx,duy_dy,duy_dz] = gradient(uy,spacing(1),spacing(2),spacing(3)); 
[duz_dx,duz_dy,duz_dz] = gradient(uz,spacing(1),spacing(2),spacing(3));

F_deform_grad=zeros(size(dux_dx,1),size(dux_dx,2),size(dux_dx,3),3,3);
C_cauchy_strain=zeros(size(dux_dx,1),size(dux_dx,2),size(dux_dx,3),3,3);
E_green_strain=zeros(size(dux_dx,1),size(dux_dx,2),size(dux_dx,3),3,3);

%Compute deformation tensor F = grad(u) + I
F_deform_grad(:,:,:,1,1)=dux_dx+1; % Add 1 to diagonal elements (Identity tensor addition)
F_deform_grad(:,:,:,1,2)=dux_dy;
F_deform_grad(:,:,:,1,3)=dux_dz;
F_deform_grad(:,:,:,2,1)=duy_dx;
F_deform_grad(:,:,:,2,2)=duy_dy+1;
F_deform_grad(:,:,:,2,3)=duy_dz;
F_deform_grad(:,:,:,3,1)=duz_dx;
F_deform_grad(:,:,:,3,2)=duz_dy;
F_deform_grad(:,:,:,3,3)=duz_dz+1;

%Compute Cauchy C = F'F and Green-Strain Tensor E = 1/2(C-I)
for r=1:size(dux_dx,1)
    for c=1:size(dux_dx,2)
        for s=1:size(dux_dx,3)
                F=squeeze(F_deform_grad(r,c,s,:,:));
                C=F'*F; 
                E=(1/2)*(C-eye(3));
                C_cauchy_strain(r,c,s,:,:)=C; 
                E_green_strain(r,c,s,:,:)=E; 
        end
    end
end