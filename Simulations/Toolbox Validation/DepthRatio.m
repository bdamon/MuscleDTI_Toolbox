%% DepthRatio
% v 1.0.0 (1/9/2021), Bruce Damon
% 
% This code simulates a set of straight fibers that emerge from an 
% aponeurosis on the anatomical left side. 
%
% The known values for pennation angle, curvature, and fiber length are 20
% degrees, 0 m^-1, and 96.5 mm, respectively.
%
% In all sets of simulations, the frame of reference for the diffusion 
% directions is LAS.
% 
% The image orientation is AL.
%
% Voxels are simulated at isotropic resolution (1x1x1 mm^3) and with a
% depth ratio of 7 (1x1x7 mm^3)
%
% Euler integration is used
%
% The step size is 1 voxel widths

%% start with a clean slate
clear
close all
clc

%% form diffusion tensors for muscles having different fiber orientations and using different frames of reference

% all muscles: principal diffusivities
L = [2.1 0 0; 0 1.6 0; 0 0 1.4] * 1E-3;

% fibers go right and upward 
E1 = [-sind(20) 0 cosd(20)]';
E2 = [-cosd(20) 0 -sind(20)]';
E3 = [0 1 0]';
E = [E1 E2 E3];
D = E*L*E';


%% set DTI encoding parameters

DTI_dir = [0    0   0
    0.1890    0.6940    0.6940      % philips 24 direction DTI encoding scheme
    0.6940    0.6940    0.1890
    0.1890    0.6940   -0.6940
    -0.6940    0.6940    0.1890
    -0.6940    0.1890   -0.6940
    0.6940    0.1890   -0.6940
    -0.6340    0.4420    0.6340
    0.6340    0.4420    0.6340
    0.6340    0.6340   -0.4420
    -0.4420    0.6340   -0.6340
    0.4420    0.6340    0.6340
    0.4420    0.6340   -0.6340
    -0.6940    0.6940   -0.1890
    0.6340    0.6340    0.4420
    -0.6940    0.1890    0.6940
    -0.6340    0.4420   -0.6340
    0.6940    0.1890    0.6940
    0.6340    0.4420   -0.6340
    -0.1890    0.6940   -0.6940
    -0.4420    0.6340    0.6340
    0.6940    0.6940   -0.1890
    -0.1890    0.6940    0.6940
    -0.6340    0.6340   -0.4420
    -0.6340    0.6340    0.4420];
b_value = 450;                        % b-value, s/mm^2

%% form isotropic diffusion-weighted images

% form mask
mask_DR1 = zeros(35, 35, 126);                                               %RU
mask_DR1(2:34,2:34,2:125) = 1;                                               %main part of the muscle is square

%in slice 1, create markers to track flip/rotation
mask_DR1(1:35,2,1) = 2;                                                      %seed point from left side
mask_DR1(34,1:35,1) = 2;                                                     %seed point from posterior

% signal generation loop
dwi = zeros(35, 35, 126, 25);
for d=1:25
    
    % get encoding direction
    loop_dir = squeeze(DTI_dir(d,:));
    
    % form images
    loop_signal = 1*exp(-b_value*loop_dir*D*loop_dir');
    dwi(:,:,:,d) = loop_signal;
    loop_signal = 1*exp(-b_value*loop_dir*D*loop_dir');
    dwi(:,:,:,d) = loop_signal;
    
    % mask the image
    dwi(:,:,:,d) = dwi(:,:,:,d).*mask_DR1;
    dwi(:,:,:,d) = dwi(:,:,:,d).*mask_DR1;
    
end

% form tensor matrix 
tensor = zeros(35, 35, 126, 3, 3);
for r=1:35
    for c=1:35
        for s=1:126
            if mask_DR1(r, c, s)==1
                
                loop_signal = squeeze(dwi(r, c, s, :));
                loop_D = signal2tensor2(loop_signal, DTI_dir(2:end,:), b_value);
                tensor(r, c, s, :, :) = loop_D;
                
            end
        end
    end
end

% flip images to form AL image orientation
tensor_AL_DR1 = zeros(35, 35, 126, 3, 3);
for s=1:126
    
    mask_DR1(:, :, s) = fliplr(mask_DR1(:, :, s));
    
    for r=1:3
        for c=1:3
            
            tensor_AL_DR1(:, :, s, r, c) = fliplr(tensor(:, :, s, r, c));
            
        end
    end
end

%compress in slice direction
mask_DR7 = imresize3(mask_DR1, [35 35 126/7]);
mask_DR7(mask_DR7>0.5)=1;

tensor_AL_DR7 = zeros(35, 35, 126/7, 3,3);
for r=1:3
    for c=1:3
        loop_images = squeeze(tensor_AL_DR1(:,:,:,r,c));
        loop_images_DR7 = imresize3(loop_images, [35 35 126/7]);
        tensor_AL_DR7(:,:,:,r,c) = loop_images_DR7;
    end
end


%% form roi_meshes - depth ratio of 1

% set up to visualize the meshes, starting with DR1:
plot_options_DR1.anat_dims = [35 1];                                            %FOV and slice thickness of the images to be displayed, in mm
plot_options_DR1.anat_slices = 14:14:126;                                       %display slices [14 28 42...126]
plot_options_DR1.plot_mesh = 1;                                                 %do plot an aponeurosis mesh
plot_options_DR1.plot_mask = 0;                                                 %don’t plot the mask
plot_options_DR1.plot_fibers = 0;                                               %don’t plot any fiber tracts
plot_options_DR1.mesh_size = [35 35];                                           %rows x columns of the images used to generate the mesh
plot_options_DR1.mesh_dims = [35 1];                                            %FOV and ST of the images used to create the mesh, in mm
plot_options_DR1.mesh_color = [0.75 0.75 0.75];                                 %make the mesh light gray
plot_options_DR1.mesh_size = [35 35];                                           %rows x columns of the images used to generate the mesh
plot_options_DR1.mesh_dims = [35 1];                                            %FOV and ST of the images used to create the mesh, in mm
plot_options_DR1.mesh_color = [0.75 0.75 0.75];                                 %make the mesh light gray

% set up to visualize the meshes, D71:
plot_options_DR7.anat_dims = [35 7];                                            %FOV and slice thickness of the images to be displayed, in mm
plot_options_DR7.anat_slices = 1:3:18;                                          %display slices [1 4 7...16]
plot_options_DR7.plot_mesh = 1;                                                 %do plot an aponeurosis mesh
plot_options_DR7.plot_mask = 0;                                                 %don’t plot the mask
plot_options_DR7.plot_fibers = 0;                                               %don’t plot any fiber tracts
plot_options_DR7.mesh_size = [35 35];                                           %rows x columns of the images used to generate the mesh
plot_options_DR7.mesh_dims = [35 7];                                            %FOV and ST of the images used to create the mesh, in mm
plot_options_DR7.mesh_color = [0.75 0.75 0.75];                                 %make the mesh light gray
plot_options_DR7.mesh_size = [35 35];                                           %rows x columns of the images used to generate the mesh
plot_options_DR7.mesh_dims = [35 7];                                            %FOV and ST of the images used to create the mesh, in mm
plot_options_DR7.mesh_color = [0.75 0.75 0.75];                                 %make the mesh light gray


% roi_mesh
roi_mesh_DR1=zeros(15,33,6);
for r=1:15
    roi_mesh_DR1(r,1:33,1)=2:34;
    roi_mesh_DR1(r,1:33,3)= r+13;
end
roi_mesh_DR1(:,1:33,2)=34;
roi_mesh_DR1(:,:,1:3) = ...              %need to add a little randomness to prevent /0 in fiber_quantifier
    roi_mesh_DR1(:,:,1:3)+randn(size(squeeze(roi_mesh_DR1(:,:,1:3))))*0.000001;

% find normal to mesh at each point:
mesh_row_vec = circshift(roi_mesh_DR1(:, :, 1:3), [0 -1 0]) - roi_mesh_DR1(:, :, 1:3);
mesh_col_vec = circshift(roi_mesh_DR1(:, :, 1:3), [-1 0 0]) - roi_mesh_DR1(:, :, 1:3);
mesh_row_vec(:, end, :) = mesh_row_vec(:, end-1, :);
mesh_col_vec(end, :, :) = mesh_col_vec(end-1, :, :);

roi_mesh_DR1(:, :, 4:6) = cross(mesh_row_vec, mesh_col_vec, 3);
roi_mesh_DR1(:, :, 4:6) = smooth3(roi_mesh_DR1(:, :, 4:6));
roi_norm = (roi_mesh_DR1(:, :, 4).^2 + roi_mesh_DR1(:, :, 5).^2 + roi_mesh_DR1(:, :,6).^2).^0.5;
roi_mesh_DR1(:, :, 4) = roi_mesh_DR1(:, :, 4)./roi_norm;
roi_mesh_DR1(:, :, 5) = roi_mesh_DR1(:, :, 5)./roi_norm;
roi_mesh_DR1(:, :, 6) = roi_mesh_DR1(:, :, 6)./roi_norm;
roi_mesh_DR1(:, :, 4:6) = smooth3(roi_mesh_DR1(:, :, 4:6));
roi_norm = (roi_mesh_DR1(:, :, 4).^2 + roi_mesh_DR1(:, :, 5).^2 + roi_mesh_DR1(:, :,6).^2).^0.5;
roi_mesh_DR1(:, :, 4) = roi_mesh_DR1(:, :, 4)./roi_norm;
roi_mesh_DR1(:, :, 5) = roi_mesh_DR1(:, :, 5)./roi_norm;
roi_mesh_DR1(:, :, 6) = roi_mesh_DR1(:, :, 6)./roi_norm;

% view the mesh
mesh_figure = fiber_visualizer(mask_DR1, plot_options_DR1, roi_mesh_DR1, [], []);
xlabel('Column Position (mm)')
ylabel('Row Position (mm)')
zlabel('Slice Position (mm)')
set(gca, 'CameraPosition', [188 -255 330], 'CameraTarget', [18 18 75], 'CameraViewAngle', 10, ...
    'CameraUpVector', [0 0 1], 'DataAspectRatio', [1 1 3], ...
    'xlim', [0 36], 'xtick', 0:9:36, 'ylim', [0 36], 'ytick', 0:9:36, 'zlim', [0 127], 'ztick', 0:21:126)
set(gcf, 'Position', [1000 100 700 700])
title('Aponeurosis Mesh')


% roi_mesh
roi_mesh_DR7 = roi_mesh_DR1;
roi_mesh_DR7(:,:,3) = roi_mesh_DR1(:,:,3)/7;

% find normal to mesh at each point:
mesh_row_vec_DR7 = circshift(roi_mesh_DR1(:, :, 1:3), [0 -1 0]) - roi_mesh_DR1(:, :, 1:3);
mesh_col_vec_DR7 = circshift(roi_mesh_DR1(:, :, 1:3), [-1 0 0]) - roi_mesh_DR1(:, :, 1:3);
mesh_row_vec_DR7(:, end, :) = mesh_row_vec_DR7(:, end-1, :);
mesh_col_vec_DR7(end, :, :) = mesh_col_vec_DR7(end-1, :, :);

roi_mesh_DR7(:, :, 4:6) = cross(mesh_row_vec_DR7, mesh_col_vec_DR7, 3);
roi_mesh_DR7(:, :, 4:6) = smooth3(roi_mesh_DR7(:, :, 4:6));
roi_norm = (roi_mesh_DR7(:, :, 4).^2 + roi_mesh_DR7(:, :, 5).^2 + roi_mesh_DR7(:, :,6).^2).^0.5;
roi_mesh_DR7(:, :, 4) = roi_mesh_DR7(:, :, 4)./roi_norm;
roi_mesh_DR7(:, :, 5) = roi_mesh_DR7(:, :, 5)./roi_norm;
roi_mesh_DR7(:, :, 6) = roi_mesh_DR7(:, :, 6)./roi_norm;
roi_mesh_DR7(:, :, 4:6) = smooth3(roi_mesh_DR7(:, :, 4:6));
roi_norm = (roi_mesh_DR7(:, :, 4).^2 + roi_mesh_DR7(:, :, 5).^2 + roi_mesh_DR7(:, :,6).^2).^0.5;
roi_mesh_DR7(:, :, 4) = roi_mesh_DR7(:, :, 4)./roi_norm;
roi_mesh_DR7(:, :, 5) = roi_mesh_DR7(:, :, 5)./roi_norm;
roi_mesh_DR7(:, :, 6) = roi_mesh_DR7(:, :, 6)./roi_norm;

% view the mesh
mesh_figure_DR7 = fiber_visualizer(mask_DR7, plot_options_DR7, roi_mesh_DR7, [], []);
xlabel('Column Position (mm)')
ylabel('Row Position (mm)')
zlabel('Slice Position (mm)')
set(gca, 'CameraPosition', [188 -255 330], 'CameraTarget', [18 18 75], 'CameraViewAngle', 10, ...
    'CameraUpVector', [0 0 1], 'DataAspectRatio', [1 1 3], ...
    'xlim', [0 36], 'xtick', 0:9:36, 'ylim', [0 36], 'ytick', 0:9:36, 'zlim', [0 127], 'ztick', 0:21:126)
set(gcf, 'Position', [1000 100 700 700])
title('Aponeurosis Mesh')


%% fiber tracking: DR1, DR7

% Set fiber tracking options - DR of 1
ft_options_DR1.ref_frame = 'LAS';                                            %left-anterior-superior directions are +X, +Y, +Z
ft_options_DR1.image_orient = 'AL';                                          %image top is anterior, image right is left
ft_options_DR1.mesh_dist = 0;                                                %don’t shift the mesh
ft_options_DR1.prop_algo = 'euler';                                          %Euler integration
ft_options_DR1.step_size = 1;                                                %1 pixel width step
ft_options_DR1.term_mthd = 'bin2';                                           %BIN2 stop algorithm
ft_options_DR1.angle_thrsh = [25 3];                                         %>=25 degree inter-segment angle disallowed; look back three points
ft_options_DR1.fa_thrsh = [.1 .4];                                           %0.1<FA<0.4 range of FA values allowed
ft_options_DR1.depth_ratio = 1;                                              %ratio of ST/in-plane resolution of reconstructed images
 
% Set fiber tracking options - DR of 1
ft_options_DR7.ref_frame = 'LAS';                                            %left-anterior-superior directions are +X, +Y, +Z
ft_options_DR7.image_orient = 'AL';                                          %image top is anterior, image right is left
ft_options_DR7.mesh_dist = 0;                                                %don’t shift the mesh
ft_options_DR7.prop_algo = 'euler';                                          %Euler integration
ft_options_DR7.step_size = 1;                                                %1 pixel width step
ft_options_DR7.term_mthd = 'bin2';                                           %BIN2 stop algorithm
ft_options_DR7.angle_thrsh = [25 3];                                         %>=25 degree inter-segment angle disallowed; look back three points
ft_options_DR7.fa_thrsh = [.1 .4];                                           %0.1<FA<0.4 range of FA values allowed
ft_options_DR7.depth_ratio = 7;                                              %ratio of ST/in-plane resolution of reconstructed images

% Set visualization options
plot_options_DR1.plot_fibers = 1;                                               %do plot any fiber tracts
plot_options_DR1.fiber_color = [.8 .2 .2];                                      %make the fibers red
plot_options_DR1.dti_size = [35 35];                                            %rows x columns of the DTI data
plot_options_DR1.dti_dims = [35 1];                                             %FOV and ST of the DTI data

plot_options_DR7.plot_fibers = 1;                                               %do plot any fiber tracts
plot_options_DR7.fiber_color = [.8 .2 .2];                                      %make the fibers red
plot_options_DR7.dti_size = [35 35];                                            %rows x columns of the DTI data
plot_options_DR7.dti_dims = [35 7];                                             %FOV and ST of the DTI data

% Fiber track - depth ratio of 1
[fiber_all_DR1, roi_flag_DR1, stop_list_DR1, fiber_length_DR1, FA_all_DR1, MD_DR1] = fiber_track...
    (tensor_AL_DR1, mask_DR1, roi_mesh_DR1, ft_options_DR1);

% Fiber track - depth ratio of 7
[fiber_all_DR7, roi_flag_DR7, stop_list_DR7, fiber_length_DR7, FA_all_DR7, MD_DR7] = fiber_track...
    (tensor_AL_DR7, mask_DR7, roi_mesh_DR7, ft_options_DR7);

% view the fiber tracts
fiber_figure_DR1 = fiber_visualizer(mask_DR1, plot_options_DR1, roi_mesh_DR1, [], fiber_all_DR1);
set(gca, 'CameraPosition', [339 -553 262], 'CameraTarget', [18 18 64], 'CameraViewAngle', 10, ...
    'CameraUpVector', [0 0 1], 'DataAspectRatio', [1 1 1], ...
    'xlim', [-.01 36.01], 'xtick', 0:9:36, 'ylim', [0 36.01], 'ytick', 0:9:36, 'zlim', [0 127], 'ztick', 0:14:126, ...
    'fontsize', 12)
xlabel('Column Position (mm)', 'fontweight', 'bold', 'fontsize', 14)
ylabel('Row Position (mm)', 'fontweight', 'bold', 'fontsize', 14)
zlabel('Slice Position (mm)', 'fontweight', 'bold', 'fontsize', 14)
set(gcf, 'Position', [1000 100 700 850])
title('Aponeurosis Mesh and Fiber Tracts (DR=1)')

% view the fiber tracts
fiber_figure_DR7 = fiber_visualizer(mask_DR7, plot_options_DR7, roi_mesh_DR7, [], fiber_all_DR7);
set(gca, 'CameraPosition', [339 -553 262], 'CameraTarget', [18 18 64], 'CameraViewAngle', 10, ...
    'CameraUpVector', [0 0 1], 'DataAspectRatio', [1 1 1], ...
    'xlim', [-.01 36.01], 'xtick', 0:9:36, 'ylim', [0 36.01], 'ytick', 0:9:36, 'zlim', [0 127], 'ztick', 0:14:126, ...
    'fontsize', 12)
xlabel('Column Position (mm)', 'fontweight', 'bold', 'fontsize', 14)
ylabel('Row Position (mm)', 'fontweight', 'bold', 'fontsize', 14)
zlabel('Slice Position (mm)', 'fontweight', 'bold', 'fontsize', 14)
set(gcf, 'Position', [1000 100 700 850])
title('Aponeurosis Mesh and Fiber Tracts (DR=7)')


%% fiber quantifier

% set options
fq_options_DR7.dwi_res=[35 35 7];
fq_options_DR7.filt_kernel=3;
fq_options_DR7.mesh_units='vx';
fq_options_DR7.tract_units='vx';

fq_options_DR1.dwi_res=[35 35 1];
fq_options_DR1.filt_kernel=3;
fq_options_DR1.mesh_units='vx';
fq_options_DR1.tract_units='vx';

% Depth ratio of 7
[angle_list_DR7, distance_list_DR7, curvature_list_DR7, fiber_all_mm_DR7, n_points_DR7, apo_area_DR7] = ...
     fiber_quantifier(fiber_all_DR7, roi_mesh_DR7, fq_options_DR7);
full_fiber_length_DR7 = max(distance_list_DR7, [],3);

%Depth ratio of 1
[angle_list_DR1, distance_list_DR1, curvature_list_DR1, fiber_all_mm_DR1, n_points_DR1, apo_area_DR1] = ...
     fiber_quantifier(fiber_all_DR1, roi_mesh_DR1, fq_options_DR1);
full_fiber_length_DR1 = max(distance_list_DR1, [],3);


%% output data

clc

summary_DR1 = [mean(angle_list_DR1(angle_list_DR1~=0)) std(angle_list_DR1(angle_list_DR1~=0)) 20 ...
mean(curvature_list_DR1(curvature_list_DR1~=0)) std(curvature_list_DR1(curvature_list_DR1~=0)) 0 ...
mean(full_fiber_length_DR1(full_fiber_length_DR1~=0)) std(full_fiber_length_DR1(full_fiber_length_DR1~=0)) 33/sind(20)]

summary_DR7 = [mean(angle_list_DR7(angle_list_DR7~=0)) std(angle_list_DR7(angle_list_DR7~=0)) 20 ...
mean(curvature_list_DR7(curvature_list_DR7~=0)) std(curvature_list_DR7(curvature_list_DR7~=0)) 0 ...
mean(full_fiber_length_DR7(full_fiber_length_DR7~=0)) std(full_fiber_length_DR7(full_fiber_length_DR7~=0)) 33/sind(20)]

close all
save Validation_DepthRatio
