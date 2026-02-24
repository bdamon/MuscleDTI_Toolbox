function [mesh_transform, mesh_disp] = transformmesh(displacementfield, roi_mesh, mesh_img_size, method)
%%
%FUNCTION transformmesh
% The function transformmesh is used to transform a 3D aponeurosis mesh from
% the undeformed to the deformed state using a 3D displacement field.
%
%INPUT ARGUMENTS
% displacementfield: a 4D array containing a displacement field that
%  describes the deformation of the whole muscle from the undeformed to
%  the deformed state (i.e. contraction). This displacement field must be 
%  compatible with MATLAB's imwarp function. 
% roi_mesh: an aponeurosis mesh defined using the define_roi function of
%  the muscleDTI_Toolbox 
% mesh_img_size: size of the image used to define the aponeurosis mesh,
%  i.e. the first 3 dimensions of the DTI images,
%  mesh_img_size = size(dti_img,[1 2 3])
% method: interpolation method used to compute the displacementfield in
%  each point of the mesh (linear is the default), could be any method 
%  compatible with MATLAB's interp3 function 
%
%OUTPUT ARGUMENTS
% mesh_transform: the transformed aponeurosis mesh
% mesh_disp: the displacements applied to each point of the mesh to
%  transform the mesh from the undeformed to the deformed state
%
%VERSION INFORMATION
% v1.0.0, Melissa Hooijmans - Info added by Roberto Pineda Guzman, 
%   October 7 2024.
%
% ACKNOWLEDGMENTS
% Grant support: NIH/NIAMS R01 AR073831

%% 

if nargin<4
    method='linear'; %Assume linear interpolation if not specified 
end

D_field=zeros(size(displacementfield));

% Resize displacement field to match roi_mesh image size
for i=1:3
    D_field(:,:,:,i)=imresize3(displacementfield(:,:,:,i), mesh_img_size,'linear','Antialiasing',false);
end

XD=D_field(:,:,:,1);
YD=D_field(:,:,:,2);
ZD=D_field(:,:,:,3);

[dy,dx,slices,~]=size(D_field);

X=1:1:dy;
Y=1:1:dx;
Z=1:1:slices;

[n_row,n_col,~]=size(roi_mesh); 
mesh_transform=zeros(n_row,n_col,3);
mesh_disp=zeros(n_row,n_col,3);

for i=1:n_row
    for j=1:n_col
        %Interpolate to find displacements of each mesh point
        mesh_disp(i,j,2)=interp3(X,Y,Z,XD,roi_mesh(i,j,2),roi_mesh(i,j,1),roi_mesh(i,j,3),method);
        mesh_disp(i,j,1)=interp3(X,Y,Z,YD,roi_mesh(i,j,2),roi_mesh(i,j,1),roi_mesh(i,j,3),method);
        mesh_disp(i,j,3)=interp3(X,Y,Z,ZD,roi_mesh(i,j,2),roi_mesh(i,j,1),roi_mesh(i,j,3),method);
    end
end

mesh_transform(:,:,2)=roi_mesh(:,:,2)-mesh_disp(:,:,2);
mesh_transform(:,:,1)=roi_mesh(:,:,1)-mesh_disp(:,:,1);
mesh_transform(:,:,3)=roi_mesh(:,:,3)-mesh_disp(:,:,3);


end