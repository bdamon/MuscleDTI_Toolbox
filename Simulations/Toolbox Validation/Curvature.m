%% Curvature
% v 1.0.0 (1/9/2021), Bruce Damon
% 
% This code simulates a set of circular fibers that emerge from an 
% aponeurosis on the anatomical left side. 
%
% Images are not formed
%
% The apparent curvature can be varied by changing the apparent image
% resolution (1x1x1 mm^3 and 2x2x2 mm^3)


%% clean slate
clear
close all
clc

%% form synthetic image; isotropic voxels
anat_image = zeros(35, 35, 120);
anat_image(2:34, 2:34, 2:125) = 1;

%% form roi mesh

roi_mesh=zeros(20,33,6);
for r=1:20
    roi_mesh(r,1:33,1)=2:34;
    roi_mesh(r,1:33,3)= r+13;
end
roi_mesh(:,1:33,2)=2;
roi_mesh(:,:,1:3) = ...              %need to add a little randomness to prevent /0 in fiber_quantifier
    roi_mesh(:,:,1:3)+randn(size(squeeze(roi_mesh(:,:,1:3))))*0.000001;

% find normal to mesh at each point:
mesh_row_vec = circshift(roi_mesh(:, :, 1:3), [0 -1 0]) - roi_mesh(:, :, 1:3);
mesh_col_vec = circshift(roi_mesh(:, :, 1:3), [-1 0 0]) - roi_mesh(:, :, 1:3);
mesh_row_vec(:, end, :) = mesh_row_vec(:, end-1, :);
mesh_col_vec(end, :, :) = mesh_col_vec(end-1, :, :);

roi_mesh(:, :, 4:6) = cross(mesh_row_vec, mesh_col_vec, 3);
roi_mesh(:, :, 4:6) = smooth3(roi_mesh(:, :, 4:6));
roi_norm = (roi_mesh(:, :, 4).^2 + roi_mesh(:, :, 5).^2 + roi_mesh(:, :,6).^2).^0.5;
roi_mesh(:, :, 4) = roi_mesh(:, :, 4)./roi_norm;
roi_mesh(:, :, 5) = roi_mesh(:, :, 5)./roi_norm;
roi_mesh(:, :, 6) = roi_mesh(:, :, 6)./roi_norm;
roi_mesh(:, :, 4:6) = smooth3(roi_mesh(:, :, 4:6));
roi_norm = (roi_mesh(:, :, 4).^2 + roi_mesh(:, :, 5).^2 + roi_mesh(:, :,6).^2).^0.5;
roi_mesh(:, :, 4) = roi_mesh(:, :, 4)./roi_norm;
roi_mesh(:, :, 5) = roi_mesh(:, :, 5)./roi_norm;
roi_mesh(:, :, 6) = roi_mesh(:, :, 6)./roi_norm;

% set up to visualize the meshes, starting with AR1:
plot_options.anat_dims = [35 1]; %FOV and slice thickness of the images to be displayed, in mm
plot_options.anat_slices = 7:21:126;  
plot_options.plot_mesh = 1; %do plot an aponeurosis mesh
plot_options.plot_mask = 0; %don’t plot the mask
plot_options.plot_fibers = 0; %don’t plot any fiber tracts
plot_options.mesh_size = [35 35]; %rows x columns of the images used to generate the mesh
plot_options.mesh_dims = [35 1]; %FOV and ST of the images used to create the mesh, in mm
plot_options.mesh_color = [0.75 0.75 0.75]; %make the mesh light gray
plot_options.mesh_size = [35 35]; %rows x columns of the images used to generate the mesh
plot_options.mesh_dims = [35 1]; %FOV and ST of the images used to create the mesh, in mm
plot_options.mesh_color = [0.75 0.75 0.75]; %make the mesh light gray

% view the mesh
mesh_figure = fiber_visualizer(anat_image, plot_options, roi_mesh, [], []);
xlabel('Column Position (mm)')
ylabel('Row Position (mm)')
zlabel('Slice Position (mm)')
set(gca, 'CameraPosition', [188 -255 330], 'CameraTarget', [18 18 75], 'CameraViewAngle', 10, ...
    'CameraUpVector', [0 0 1], 'DataAspectRatio', [1 1 3], ...
    'xlim', [0 36], 'xtick', 0:9:36, 'ylim', [0 36], 'ytick', 0:9:36, 'zlim', [0 127], 'ztick', 0:21:126)
set(gcf, 'Position', [1000 100 700 700])
title('Aponeurosis Mesh')


%% form fiber_all

% note that equation of a circle is r^2 = (x-a)^2 + (y-b)^2 
column_points = 0:0.5:33;
row_points = zeros(size(column_points));
slice_points = 33 - (33^2 - column_points.^2).^(1/2);
figure
plot(column_points, slice_points)

%put points into fiber_all matrix
fiber_all = zeros(length(roi_mesh(:,1,1)), length(roi_mesh(1,:,1)), length(column_points), 3);
for r=1:length(roi_mesh(:,1,1))
    for c=1:length(roi_mesh(1,:,1))
        seed_point = squeeze(roi_mesh(r,c,1:3));
        
        loop_row_points = row_points + seed_point(1);
        loop_column_points = column_points + seed_point(2);
        loop_slice_points = slice_points + seed_point(3);
        
        fiber_all(r,c,:,1) = loop_row_points;
        fiber_all(r,c,:,2) = loop_column_points;
        fiber_all(r,c,:,3) = loop_slice_points;
    end
end

% Set visualization options
plot_options_fibers.anat_dims = [35 1]; %FOV and slice thickness of the images to be displayed, in mm
plot_options_fibers.anat_slices = 7:21:126; 
plot_options_fibers.plot_mesh = 1; %do plot an aponeurosis mesh
plot_options_fibers.plot_mask = 0; %don’t plot the mask
plot_options_fibers.plot_fibers = 1; %do plot any fiber tracts
plot_options_fibers.mesh_size = [35 35]; %rows x columns of the images used to generate the mesh
plot_options_fibers.mesh_dims = [35 1]; %FOV and ST of the images used to create the mesh
plot_options_fibers.mesh_color = [0.75 0.75 0.75]; %make the mesh light gray
plot_options_fibers.mask_size = [192 192]; %rows x columns of the images used to generate the mask
plot_options_fibers.mask_dims = [192 1]; %FOV and ST of the images used to create the mask
plot_options_fibers.mask_color = [1 0 0]; %make the mask a red, semi-transparent overlay
plot_options_fibers.fiber_color = [.8 .2 .2]; %make the fibers red
plot_options_fibers.dti_size = [35 35]; %rows x columns of the DTI data
plot_options_fibers.dti_dims = [35 1]; %FOV and ST of the DTI data

% view the fiber tracts
fiber_figure1 = fiber_visualizer(anat_image, plot_options_fibers, roi_mesh, [], fiber_all);
xlabel('Column Position (mm)')
ylabel('Row Position (mm)')
zlabel('Slice Position (mm)')
set(gca, 'CameraPosition', [188 -255 330], 'CameraTarget', [18 18 75], 'CameraViewAngle', 10, ...
    'CameraUpVector', [0 0 1], 'DataAspectRatio', [1 1 3], ...
    'xlim', [0 36], 'xtick', 0:9:36, 'ylim', [0 36], 'ytick', 0:9:36, 'zlim', [0 127], 'ztick', 0:21:126)
set(gcf, 'Position', [1000 100 700 700])
title('Aponeurosis Mesh and Fiber Tracts (AL)')

        
%% fiber_quantifier

%first, imagine that the image resolution is 1x1x1:
fq_options.dwi_res=[35 35 1];
fq_options.filt_kernel=3;
fq_options.mesh_units='vx';
fq_options.tract_units='vx';

[angle_list_111, distance_list_111, curvature_list_111, fiber_all_mm_111, n_points_111, apo_area_111] = ...
     fiber_quantifier(fiber_all, roi_mesh, fq_options);
full_fiber_length_111 = max(distance_list_111, [],3);


%first, imagine that the image resolution is 1x1x1:
fq_options.dwi_res=[70 35 2];
fq_options.filt_kernel=3;

[angle_list_222, distance_list_222, curvature_list_222, fiber_all_mm_222, n_points_222, apo_area_222] = ...
     fiber_quantifier(fiber_all, roi_mesh, fq_options);
full_fiber_length_222 = max(distance_list_222, [],3);
        

%% summary data
clc

summary_111 = [mean(curvature_list_111(curvature_list_111~=0)) std(curvature_list_111(curvature_list_111~=0)) 1/.033 ...
mean(full_fiber_length_111(full_fiber_length_111~=0)) std(full_fiber_length_111(full_fiber_length_111~=0)) 2*pi*33/4]


summary_222 = [mean(curvature_list_222(curvature_list_222~=0)) std(curvature_list_222(curvature_list_222~=0)) 1/.066 ...
mean(full_fiber_length_222(full_fiber_length_222~=0)) std(full_fiber_length_222(full_fiber_length_222~=0)) 2*pi*66/4]

close all
save Validation_Curvature



