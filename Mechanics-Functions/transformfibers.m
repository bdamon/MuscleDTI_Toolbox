function [fiber_all_transform, fiber_disp] = transformfibers(displacementfield, fiber_all, dti_size, method)

%%
%FUNCTION transformfibers
% The function transformfibers is used to transform the 3D fiber-tracts
% representing a muscle's architecture from the undeformed to the deformed 
% state. This function allows for the estimation of 3D muscle architecture
% during contraction (deformed state).
%
%INPUT ARGUMENTS
% displacementfield: a 4D array containing a displacement field that
%  describes the deformation of the whole muscle from the undeformed to
%  the deformed state (i.e. contraction). This displacement field must be 
%  compatible with MATLAB's imwarp function. 
% fiber_all: the 3D fiber-tracts obtained from the fiber_track function
% dti_size: the first three dimensions of the DTI images used to obtain the
%  fiber-tract data, i.e. dti_size = size(dti_img,[1 2 3])
% method: interpolation method used to compute the displacementfield in
%  each point of the fiber-tracts (linear is the default), could be any method 
%  compatible with MATLAB's interp3 function 
%
%OUTPUT ARGUMENTS
% fiber_all_transform: the transformed fiber-tracts
% fiber_disp: the displacements applied to each point of the fiber-tracts to
%  transform them from the undeformed to the deformed state (i.e. contraction)
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
% Resize displacement field to match DTI image size
for i=1:3
    D_field(:,:,:,i)=imresize3(displacementfield(:,:,:,i), dti_size,'linear','Antialiasing',false);
end

XD=D_field(:,:,:,1);
YD=D_field(:,:,:,2);
ZD=D_field(:,:,:,3);

[dy,dx,slices,~]=size(D_field);

X=1:1:dy;
Y=1:1:dx;
Z=1:1:slices;

[n_row,n_col,n_points,~] = size(fiber_all); 
fiber_all_transform=zeros(n_row,n_col,n_points,3);
fiber_disp=zeros(n_row,n_col,n_points,3);

t_total = tic; %Function timer - can be slow if many fiber-tracts are being transformed
for i=1:n_row
    for j=1:n_col
        for k=1:n_points
            if fiber_all(i,j,k,1)>0
                %Interpolate to find displacements of each fiber-tract point
                fiber_disp(i,j,k,2)=interp3(X,Y,Z,XD,fiber_all(i,j,k,2),fiber_all(i,j,k,1),fiber_all(i,j,k,3),method);
                fiber_disp(i,j,k,1)=interp3(X,Y,Z,YD,fiber_all(i,j,k,2),fiber_all(i,j,k,1),fiber_all(i,j,k,3),method);
                fiber_disp(i,j,k,3)=interp3(X,Y,Z,ZD,fiber_all(i,j,k,2),fiber_all(i,j,k,1),fiber_all(i,j,k,3),method);
            end
        end
    end
end
toc(t_total)

% Apply displacements to fiber-tracts
fiber_all_transform(:,:,:,2)=fiber_all(:,:,:,2)-fiber_disp(:,:,:,2);
fiber_all_transform(:,:,:,1)=fiber_all(:,:,:,1)-fiber_disp(:,:,:,1);
fiber_all_transform(:,:,:,3)=fiber_all(:,:,:,3)-fiber_disp(:,:,:,3);


end
