%% transform_tracts.m 
% This version of transform_tracts.m produces the transformed fiber-tracts
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

%% Load muscle fiber-tracts, mask and aponeurosis information

load muscle_mask
load roi_mesh
load fiber_all

%% Register structural images using the additive Demons registration method 

[displ_field, ~] = imregdemons(struct_img_rest,...
            struct_img_contr,2000,'PyramidLevels',4,'AccumulatedFieldSmoothing',1.5);

%% Smooth displacement field for fiber transformation

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

%% Transform mesh and fiber-tracts

roi_mesh_transformed= transformmesh(displ_field,roi_mesh,dti_matrix_size,'linear');  %Transform mesh

% find normal to mesh at each point:
mesh_row_vec = circshift(roi_mesh_transformed(:, :, 1:3), [0 -1 0]) - roi_mesh_transformed(:, :, 1:3);
mesh_col_vec = circshift(roi_mesh_transformed(:, :, 1:3), [-1 0 0]) - roi_mesh_transformed(:, :, 1:3);
mesh_row_vec(:, n_col, :) = mesh_row_vec(:, n_col-1, :);
mesh_col_vec(n_row, :, :) = mesh_col_vec(n_row-1, :, :);

roi_mesh_transformed(:, :, 4:6) = cross(mesh_row_vec, mesh_col_vec, 3);
roi_mesh_transformed(:, :, 4:6) = smooth3(roi_mesh_transformed(:, :, 4:6));
roi_norm = (roi_mesh_transformed(:, :, 4).^2 + roi_mesh_transformed(:, :, 5).^2 + roi_mesh_transformed(:, :,6).^2).^0.5;
roi_mesh_transformed(:, :, 4) = roi_mesh_transformed(:, :, 4)./roi_norm;
roi_mesh_transformed(:, :, 5) = roi_mesh_transformed(:, :, 5)./roi_norm;
roi_mesh_transformed(:, :, 6) = roi_mesh_transformed(:, :, 6)./roi_norm;
roi_mesh_transformed(:, :, 4:6) = smooth3(roi_mesh_transformed(:, :, 4:6));
roi_norm = (roi_mesh_transformed(:, :, 4).^2 + roi_mesh_transformed(:, :, 5).^2 + roi_mesh_transformed(:, :,6).^2).^0.5;
roi_mesh_transformed(:, :, 4) = roi_mesh_transformed(:, :, 4)./roi_norm;
roi_mesh_transformed(:, :, 5) = roi_mesh_transformed(:, :, 5)./roi_norm;
roi_mesh_transformed(:, :, 6) = roi_mesh_transformed(:, :, 6)./roi_norm;

fiber_all_transformed = transformfibers(displ_field,fiber_all,dti_img_size,'linear');  %Transform fiber-tracts

%% Smooth transformed fiber-tracts

%set fiber_smoother options
fs_options.interpolation_step = 1;                                          %interpolate at 1 pixel widths (=1 mm)
fs_options.p_order = [3 3 3];                                               %fit row, column, and slice positions to 4th, 4th,and 3rd order polynomials
fs_options.tract_units = 'vx';                                              %units - 'vx' or 'mm'
fs_options.dwi_res = [dti_fov(1) dti_numrows dti_slcthick];                 %DTI FOV, matrix size, slice thickness
                
% call the function:
[smoothed_fiber_all_transformed, ~, smoothed_fiber_all_transformed_mm, pcoeff_r, pcoeff_c, pcoeff_s,...
    n_points_smoothed, residuals, residuals_mm] = fiber_smoother(fiber_all_transformed, fs_options);

load anat_image                                                             %Load anatomical image with same FOV, resolution, and matrix size as DTI images

% set fiber visualizer options for the mask
fv_options.anat_dims = [anat_fov(1) anat_slcthick];                         %FOV and slice thickness of the images to be displayed, in mm
fv_options.anat_slices = 5:9:44;                                            %display slices 14, 24, 34, and 44
fv_options.plot_mesh = 0;                                                   %don't plot an aponeurosis mesh
fv_options.plot_mask = 1;                                                   %do plot the mask
fv_options.plot_fibers = 0;                                                 %don't plot any fiber tracts
fv_options.mask_size = [anat_numrows anat_numcols];                         %rows x columns of the images used to generate the mask
fv_options.mask_dims = [anat_fov(1) anat_slcthick];                         %FOV and slice thickness of the images used to create the mask, in mm
fv_options.mask_color = [1 0 0];                                            %make the mask a red, semi-transparent overlay

%visualize transformed fiber tracts
sampled_smoothed_fiber_fig = fiber_visualizer(anat_image, fv_options, roi_mesh_transformed, [], smoothed_fiber_all_transformed);

%% Quantify architectural properties

fq_options.dwi_res = [dti_fov(1) dti_numrows dti_slcthick];                 %DTI FOV, matrix size, slice thickness
fq_options.filt_kernel = 3;                                                 %size of smoothing kernel for determining aponeurosis normal vectors
fq_options.mesh_units = 'vx';                                               %tracts are in units of voxels
fq_options.tract_units = 'vx';

[angle_list, distance_list, curvature_list, fiber_all_mm, n_points, apo_area] = ...
     fiber_quantifier(fq_options, sampled_smoothed_fiber_all_vx, roi_mesh, mask);