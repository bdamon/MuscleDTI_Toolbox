%% compute_muscle_strain.m 
% This version of compute_tracts_strain.m produces the fiber-tract strain measurements
% shown in Pineda Guzman et al., J Biomech, 2026 

%% Clean slate
clear
close all
clc

%% Load structural images obtained during rest and contraction 

load struct_img_rest % Both struct_imgs must share the same FOV, in-plane resolution, and slice thickness
load struct_img_contr

struct_dims = [struct_fov(1) struct_slcthick];    %FOV and slice thickness of the structural images, in mm
struct_matrix_size = size(struct_img_rest);     %Matrix size of structural images

%% Load muscle fiber direction maps and muscle mask

load muscle_mask
load smoothed_fiber_all

%% Register structural images using the additive Demons registration method 

[displ_field, ~] = imregdemons(struct_img_rest,...
            struct_img_contr,2000,'PyramidLevels',4,'AccumulatedFieldSmoothing',1.5);

%% Process displacement field data

% Convert displacement field to physical units
struct_resol = [struct_dims(1)/struct_matrix_size(1), struct_dims(2)/struct_matrix_size(2) struct_slcthick];
displ_field = matlab_to_phys_displ_conv(displ_field,struct_resol);

displ_field_x = squeeze(displ_field(:,:,:,1));
displ_field_y = squeeze(displ_field(:,:,:,2));
displ_field_z = squeeze(displ_field(:,:,:,3));

% Resize displacement field to isotropic 1 mm spacing

struct_matrix_size_iso = [struct_dims(1) struct_dims(1) struct_dims(2)*struct_matrix_size(3)];    %Image size that will result in isotropic 1 mm voxels

displ_field_x_iso = imresize3(displ_field_x,struct_matrix_size_iso,"linear","Antialiasing",false);
displ_field_y_iso = imresize3(displ_field_y,struct_matrix_size_iso,"linear","Antialiasing",false);
displ_field_z_iso = imresize3(displ_field_z,struct_matrix_size_iso,"linear","Antialiasing",false);

% Smooth displacement field 

for i=1:10
    displ_field_x_iso=smooth3(displ_field_x_iso,"gaussian",3);
    displ_field_y_iso=smooth3(displ_field_y_iso,"gaussian",3);
    displ_field_z_iso=smooth3(displ_field_z_iso,"gaussian",3);
end

% Resize displacement field to DTI data resolution

dti_dims = [dti_fov(1) dti_slcthick];                       %FOV and slice thickness of the DTI images, in mm, FOV must be the same as struct_fov(1)
dti_matrix_size = [dti_numrows dti_numcols dti_numslcs];    %matrix size of DTI images
dti_resol = [dti_dims(1)/dti_matrix_size(1), dti_dims(2)/dti_matrix_size(2) dti_slcthick]; % DTI data resolution

displ_field_x = imresize3(displ_field_x_iso,dti_matrix_size,"linear","Antialiasing",false);
displ_field_y = imresize3(displ_field_y_iso,dti_matrix_size,"linear","Antialiasing",false);
displ_field_z = imresize3(displ_field_z_iso,dti_matrix_size,"linear","Antialiasing",false);

clear displ_field_x_iso %save memory space
clear displ_field_y_iso
clear displ_field_z_iso

displ_field = zeros([size(displ_field_x) 3]);     %Create displ_field variable with smooth displ_field
displ_field(:,:,:,1) = displ_field_x;
displ_field(:,:,:,2) = displ_field_y;
displ_field(:,:,:,3) = displ_field_z;

clear displ_field_x % save memory space
clear displ_field_y
clear displ_field_z

%% Compute Cauchy-Green strain field

%Compute Cauchy-Strain and Green-Strain tensor
[~,C_cauchy_strain,~] = strain_tensors(...
    displ_field_x,displ_field_y,...
    displ_field_z,dti_resol);

%% Compute strain invariants in each fiber-tract

[along_fiber_strain,along_fiber_shear,cross_fiber_shear] = ...
    strain_invariants_fiber_tract(smoothed_fiber_all,C_cauchy_strain,dti_resol);