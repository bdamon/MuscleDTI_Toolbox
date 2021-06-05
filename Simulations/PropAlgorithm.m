%% PropAlgorithm
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
% Voxels are simulated at isotropic resolution (1x1x1 mm^3).
%
% Three integration methods are used: Euler, 4th order Runge Kutta, and
% FACT
%
% The step size is 1 voxel widths

%% start with a clean slate
clear
close all
clc

%% form diffusion tensors for muscles having different fiber orientations and using different frames of reference

% all muscles: principal diffusivities
L = [2.1 0 0; 0 1.6 0; 0 0 1.4] * 1E-3;

% fibers go right and upward at 20 degrees. Adjust vectors for trhe frame
% of reference
E1 = [-sind(20) 0 cosd(20)]';
E2 = [-cosd(20) 0 -sind(20)]';
E3 = [0 1 0]';
E = [E1 E2 E3];
D = E*L*E';

%% set DTI encoding parameters

DTI_dir = [0    0   0
    0.1890    0.6940    0.6940                                              % philips 24 direction DTI encoding scheme: LAS
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
mask = zeros(35, 35, 126);
mask(2:34,2:34,2:125) = 1;      %main part of the muscle is square

%in slice 1, create markers to track flip/rotation
mask(1:35,2,1) = 2;             %seed point from A-P on left side

% signal generation loop
dwi = zeros(35, 35, 126, 25);
for d=1:25
    
    % get encoding direction
    loop_dir = squeeze(DTI_dir(d,:));
    
    % form images
    loop_signal = 1*exp(-b_value*loop_dir*D*loop_dir');
    dwi(:,:,:,d) = loop_signal;
    
    % mask the image
    dwi(:,:,:,d) = dwi(:,:,:,d).*mask;
    
end

% form tensor matrix
tensor = zeros(35, 35, 126, 3, 3);
for r=1:35
    for c=1:35
        for s=1:126
            if mask(r, c, s)==1
                
                loop_signal = squeeze(dwi(r, c, s, :));
                loop_D = signal2tensor2(loop_signal, DTI_dir(2:end,:), b_value);
                tensor(r, c, s, :, :) = loop_D;
                
            end
        end
    end
end

% flip images to form AL image orientation
tensor_AL = zeros(35, 35, 126, 3, 3);
tensor_AL_LPS = zeros(35, 35, 126, 3, 3);
tensor_AL_RAS = zeros(35, 35, 126, 3, 3);
for s=1:126
    mask(:, :, s) = fliplr(mask(:, :, s));
    for r=1:3
        for c=1:3
            
            tensor_AL(:, :, s, r, c) = fliplr(tensor(:, :, s, r, c));
            
        end
    end
end


%% form roi_mesh

% set up to visualize the meshes, starting with AR1:
plot_options.anat_dims = [35 1];                                            %FOV and slice thickness of the images to be displayed, in mm
plot_options.anat_slices = 14:14:126;                                       %display slices [14 28 42...126]
plot_options.plot_mesh = 1;                                                 %do plot an aponeurosis mesh
plot_options.plot_mask = 0;                                                 %don’t plot the mask
plot_options.plot_fibers = 0;                                               %don’t plot any fiber tracts
plot_options.mesh_size = [35 35];                                           %rows x columns of the images used to generate the mesh
plot_options.mesh_dims = [35 1];                                            %FOV and ST of the images used to create the mesh, in mm
plot_options.mesh_color = [0.75 0.75 0.75];                                 %make the mesh light gray
plot_options.mesh_size = [35 35];                                           %rows x columns of the images used to generate the mesh
plot_options.mesh_dims = [35 1];                                            %FOV and ST of the images used to create the mesh, in mm
plot_options.mesh_color = [0.75 0.75 0.75];                                 %make the mesh light gray

% roi_mesh
roi_mesh=zeros(15,33,6);
for r=1:15
    roi_mesh(r,1:33,1)=2:34;
    roi_mesh(r,1:33,3)= r+13;
end
roi_mesh(:,1:33,2)=34;
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

% view the mesh
mesh_figure = fiber_visualizer(mask, plot_options, roi_mesh, [], []);
set(gca, 'CameraPosition', [339 -553 262], 'CameraTarget', [18 18 64], 'CameraViewAngle', 10, ...
    'CameraUpVector', [0 0 1], 'DataAspectRatio', [1 1 1], ...
    'xlim', [-.01 36.01], 'xtick', 0:9:36, 'ylim', [0 36.01], 'ytick', 0:9:36, 'zlim', [0 127], 'ztick', 0:14:126, ...
    'fontsize', 12)
xlabel('Column Position (mm)', 'fontweight', 'bold', 'fontsize', 14)
ylabel('Row Position (mm)', 'fontweight', 'bold', 'fontsize', 14)
zlabel('Slice Position (mm)', 'fontweight', 'bold', 'fontsize', 14)
set(gcf, 'Position', [1000 100 700 850])
title('Aponeurosis Mesh')


%% fiber tracking

% Set fiber tracking options - Eul
ft_options_Eul.ref_frame = 'LAS';                                           %left-anterior-superior directions are +X, +Y, +Z
ft_options_Eul.image_orient = 'AL';                                         %image top is anterior, image right is left
ft_options_Eul.mesh_dist = 0;                                               %don’t shift the mesh
ft_options_Eul.prop_algo = 'euler';                                         %Euler integration
ft_options_Eul.step_size = 1;                                               %1 pixel width step
ft_options_Eul.term_mthd = 'bin2';                                          %BIN2 stop algorithm
ft_options_Eul.angle_thrsh = [25 3];                                        %>=25 degree inter-segment angle disallowed; look back three points
ft_options_Eul.fa_thrsh = [.1 .4];                                          %0.1<FA<0.4 range of FA values allowed
ft_options_Eul.depth_ratio = 1;                                             %ratio of ST/in-plane resolution of reconstructed images

% Set fiber tracking options - RK4
ft_options_RK4.ref_frame = 'LAS';                                           %left-anterior-superior directions are +X, +Y, +Z
ft_options_RK4.image_orient = 'AL';                                         %image top is anterior, image right is left
ft_options_RK4.mesh_dist = 0;                                               %don’t shift the mesh
ft_options_RK4.prop_algo = 'rk4';                                           %4th order Runge Kutta integration
ft_options_RK4.step_size = 1;                                               %1 pixel width step
ft_options_RK4.term_mthd = 'bin2';                                          %BIN2 stop algorithm
ft_options_RK4.angle_thrsh = [25 3];                                        %>=25 degree inter-segment angle disallowed; look back three points
ft_options_RK4.fa_thrsh = [.1 .4];                                          %0.1<FA<0.4 range of FA values allowed
ft_options_RK4.depth_ratio = 1;                                             %ratio of ST/in-plane resolution of reconstructed images

% Set fiber tracking options - FACT
ft_options_FACT.ref_frame = 'LAS';                                          %left-anterior-superior directions are +X, +Y, +Z
ft_options_FACT.image_orient = 'AL';                                        %image top is anterior, image right is left
ft_options_FACT.mesh_dist = 0;                                              %don’t shift the mesh
ft_options_FACT.prop_algo = 'fact';                                         %Use FACT
ft_options_FACT.r_crit = 0.95;                                              %critical R value
ft_options_FACT.num_fact_voxels = 5;                                        %number of voxels over which to calculate R
ft_options_FACT.depth_ratio = 1;                                            %ratio of ST/in-plane resolution of reconstructed images

% Set visualization options
plot_options.plot_fibers = 1;                                               %do plot any fiber tracts
plot_options.fiber_color = [.8 .2 .2];                                      %make the fibers red
plot_options.dti_size = [35 35];                                            %rows x columns of the DTI data
plot_options.dti_dims = [35 1];                                             %FOV and ST of the DTI data

% Fiber track - Euler:
[fiber_all_Eul, roi_flag_Eul, stop_list_Eul, fiber_length_Eul, FA_all_Eul, MD_Eul] = fiber_track...
    (tensor, mask, roi_mesh, ft_options_Eul);

% Fiber track - RK4:
[fiber_all_RK4, roi_flag_RK4, stop_list_RK4, fiber_length_RK4, FA_all_RK4, MD_RK4] = fiber_track...
    (tensor, mask, roi_mesh, ft_options_RK4);

% Fiber track - FACT:
[fiber_all_FACT, roi_flag_FACT, stop_list_FACT, fiber_length_FACT, FA_all_FACT, MD_FACT] = fiber_track...
    (tensor, mask, roi_mesh, ft_options_FACT);


% view the fiber tracts - Euler
fiber_figure_Eul = fiber_visualizer(mask, plot_options, roi_mesh, [], fiber_all_Eul);
set(gca, 'CameraPosition', [339 -553 262], 'CameraTarget', [18 18 64], 'CameraViewAngle', 10, ...
    'CameraUpVector', [0 0 1], 'DataAspectRatio', [1 1 1], ...
    'xlim', [-.01 36.01], 'xtick', 0:9:36, 'ylim', [0 36.01], 'ytick', 0:9:36, 'zlim', [0 127], 'ztick', 0:14:126, ...
    'fontsize', 12)
xlabel('Column Position (mm)', 'fontweight', 'bold', 'fontsize', 14)
ylabel('Row Position (mm)', 'fontweight', 'bold', 'fontsize', 14)
zlabel('Slice Position (mm)', 'fontweight', 'bold', 'fontsize', 14)
set(gcf, 'Position', [1000 100 700 850])
title('Aponeurosis Mesh and Fiber Tracts (Euler)')

% view the fiber tracts - RK4
fiber_figure_RK4 = fiber_visualizer(mask, plot_options, roi_mesh, [], fiber_all_RK4);
set(gca, 'CameraPosition', [339 -553 262], 'CameraTarget', [18 18 64], 'CameraViewAngle', 10, ...
    'CameraUpVector', [0 0 1], 'DataAspectRatio', [1 1 1], ...
    'xlim', [-.01 36.01], 'xtick', 0:9:36, 'ylim', [0 36.01], 'ytick', 0:9:36, 'zlim', [0 127], 'ztick', 0:14:126, ...
    'fontsize', 12)
xlabel('Column Position (mm)', 'fontweight', 'bold', 'fontsize', 14)
ylabel('Row Position (mm)', 'fontweight', 'bold', 'fontsize', 14)
zlabel('Slice Position (mm)', 'fontweight', 'bold', 'fontsize', 14)
set(gcf, 'Position', [1000 100 700 850])
title('Aponeurosis Mesh and Fiber Tracts (RK4)')

% view the fiber tracts - FACT
fiber_figure_FACT = fiber_visualizer(mask, plot_options, roi_mesh, [], fiber_all_FACT);
set(gca, 'CameraPosition', [339 -553 262], 'CameraTarget', [18 18 64], 'CameraViewAngle', 10, ...
    'CameraUpVector', [0 0 1], 'DataAspectRatio', [1 1 1], ...
    'xlim', [-.01 36.01], 'xtick', 0:9:36, 'ylim', [0 36.01], 'ytick', 0:9:36, 'zlim', [0 127], 'ztick', 0:14:126, ...
    'fontsize', 12)
xlabel('Column Position (mm)', 'fontweight', 'bold', 'fontsize', 14)
ylabel('Row Position (mm)', 'fontweight', 'bold', 'fontsize', 14)
zlabel('Slice Position (mm)', 'fontweight', 'bold', 'fontsize', 14)
set(gcf, 'Position', [1000 100 700 850])
title('Aponeurosis Mesh and Fiber Tracts (FACT)')


%% fiber quantifier

% set options
fq_options.dwi_res=[35 35 1];
fq_options.filt_kernel=3;
fq_options.mesh_units='vx';
fq_options.tract_units='vx';

% Euler
[angle_list_Eul, distance_list_Eul, curvature_list_Eul, fiber_all_mm_Eul, n_points_Eul, apo_area_Eul] = ...
     fiber_quantifier(fiber_all_Eul, roi_mesh, fq_options);
full_fiber_length_Eul = max(distance_list_Eul, [],3);

% 4th order RK
[angle_list_RK4, distance_list_RK4, curvature_list_RK4, fiber_all_mm_RK4, n_points_RK4, apo_area_RK4] = ...
     fiber_quantifier(fiber_all_RK4, roi_mesh, fq_options);
full_fiber_length_RK4 = max(distance_list_RK4, [],3);

% FACT
[angle_list_FACT, distance_list_FACT, curvature_list_FACT, fiber_all_mm_FACT, n_points_FACT, apo_area_FACT] = ...
     fiber_quantifier(fiber_all_FACT, roi_mesh, fq_options);
full_fiber_length_FACT = max(distance_list_FACT, [],3);


%% output data

clc

summary_Eul = [mean(angle_list_Eul(angle_list_Eul~=0)) std(angle_list_Eul(angle_list_Eul~=0)) 20 ...
mean(curvature_list_Eul(curvature_list_Eul~=0)) std(curvature_list_Eul(curvature_list_Eul~=0)) 0 ...
mean(full_fiber_length_Eul(full_fiber_length_Eul~=0)) std(full_fiber_length_Eul(full_fiber_length_Eul~=0)) 33/sind(20)]

summary_RK4 = [mean(angle_list_RK4(angle_list_RK4~=0)) std(angle_list_RK4(angle_list_RK4~=0)) 20 ...
mean(curvature_list_RK4(curvature_list_RK4~=0)) std(curvature_list_RK4(curvature_list_RK4~=0)) 0 ...
mean(full_fiber_length_RK4(full_fiber_length_RK4~=0)) std(full_fiber_length_RK4(full_fiber_length_RK4~=0)) 33/sind(20)]

summary_FACT = [mean(angle_list_FACT(angle_list_FACT~=0)) std(angle_list_FACT(angle_list_FACT~=0)) 20 ...
mean(curvature_list_FACT(curvature_list_FACT~=0)) std(curvature_list_FACT(curvature_list_FACT~=0)) 0 ...
mean(full_fiber_length_FACT(full_fiber_length_FACT~=0)) std(full_fiber_length_FACT(full_fiber_length_FACT~=0)) 33/sind(20)]

close all
save Validation_PropAlgorithm

