%% FrameOfReference
% v 1.0.0 (1/9/2021), Bruce Damon
% 
% This code simulates a set of straight fibers that emerge from an 
% aponeurosis on the anatomical left side. 
%
% The known values for pennation angle, curvature, and fiber length are 20
% degrees, 0 m^-1, and 96.5 mm, respectively.
%
% In the three sets of simulations, the frames of reference for the diffusion 
% directions are LAS, RAS, LPS.  The image orientation is AL.
% 
% The image orientation is AL.
%
% Voxels are simulated at isotropic resolution (1x1x1 mm^3).
%
% Euler integration is used.
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
E1_LAS = [-sind(20) 0 cosd(20)]';
E2_LAS = [-cosd(20) 0 -sind(20)]';
E3_LAS = [0 1 0]';
E_LAS = [E1_LAS E2_LAS E3_LAS];
D_LAS = E_LAS*L*E_LAS';

E1_LPS = [-sind(20) -0 cosd(20)]';
E2_LPS = [-cosd(20) -0 -sind(20)]';
E3_LPS = [0 -1 0]';
E_LPS = [E1_LPS E2_LPS E3_LPS];
D_LPS = E_LPS*L*E_LPS';

E1_RAS = [sind(20) 0 cosd(20)]';
E2_RAS = [-cosd(20) 0 sind(20)]';
E3_RAS = [0 1 0]';
E_RAS = [E1_RAS E2_RAS E3_RAS];
D_RAS = E_RAS*L*E_RAS';

%% set DTI encoding parameters

DTI_dir_LAS = [0    0   0
    0.1890    0.6940    0.6940      % philips 24 direction DTI encoding scheme: XYZ
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
DTI_dir_LPS = DTI_dir_LAS;
DTI_dir_LPS(:,2) = -DTI_dir_LPS(:,2);
DTI_dir_RAS = DTI_dir_LAS;
DTI_dir_RAS(:,1) = -DTI_dir_RAS(:,1);

b_value = 450;                        % b-value, s/mm^2

%% form isotropic diffusion-weighted images

% form mask
mask = zeros(35, 35, 126);
mask(2:34,2:34,2:125) = 1;      %main part of the muscle is square

%in slice 1, create markers to track flip/rotation
mask(1:35,2,1) = 2;             %seed point from A-P on left side

% signal generation loop
dwi_LAS = zeros(35, 35, 126, 25);
dwi_LPS = zeros(35, 35, 126, 25);
dwi_RAS = zeros(35, 35, 126, 25);
for d=1:25
    
    % get encoding direction
    loop_dir_LAS = squeeze(DTI_dir_LAS(d,:));
    loop_dir_LPS = squeeze(DTI_dir_LPS(d,:));
    loop_dir_RAS = squeeze(DTI_dir_RAS(d,:));
    
    % form images
    loop_signal_LAS = 1*exp(-b_value*loop_dir_LAS*D_LAS*loop_dir_LAS');
    dwi_LAS(:,:,:,d) = loop_signal_LAS;
    loop_signal_LPS = 1*exp(-b_value*loop_dir_LPS*D_LPS*loop_dir_LPS');
    dwi_LPS(:,:,:,d) = loop_signal_LPS;
    loop_signal_RAS = 1*exp(-b_value*loop_dir_RAS*D_RAS*loop_dir_RAS');
    dwi_RAS(:,:,:,d) = loop_signal_RAS;
    
    % mask the image
    dwi_LAS(:,:,:,d) = dwi_LAS(:,:,:,d).*mask;
    dwi_LPS(:,:,:,d) = dwi_LPS(:,:,:,d).*mask;
    dwi_RAS(:,:,:,d) = dwi_RAS(:,:,:,d).*mask;
    
end

% form tensor matrix
tensor_LAS = zeros(35, 35, 126, 3, 3);
tensor_LPS = zeros(35, 35, 126, 3, 3);
tensor_RAS = zeros(35, 35, 126, 3, 3);
for r=1:35
    for c=1:35
        for s=1:126
            if mask(r, c, s)==1
                
                loop_signal_LAS = squeeze(dwi_LAS(r, c, s, :));
                loop_D_LAS = signal2tensor2(loop_signal_LAS, DTI_dir_LAS(2:end,:), b_value);
                tensor_LAS(r, c, s, :, :) = loop_D_LAS;
                
                loop_signal_LPS = squeeze(dwi_LPS(r, c, s, :));
                loop_D_LPS = signal2tensor2(loop_signal_LPS, DTI_dir_LPS(2:end,:), b_value);
                tensor_LPS(r, c, s, :, :) = loop_D_LPS;
                
                loop_signal_RAS = squeeze(dwi_RAS(r, c, s, :));
                loop_D_RAS = signal2tensor2(loop_signal_RAS, DTI_dir_RAS(2:end,:), b_value);
                tensor_RAS(r, c, s, :, :) = loop_D_RAS;
                
            end
        end
    end
end

% flip images to form AL image orientation
tensor_AL_LAS = zeros(35, 35, 126, 3, 3);
tensor_AL_LPS = zeros(35, 35, 126, 3, 3);
tensor_AL_RAS = zeros(35, 35, 126, 3, 3);
for s=1:126
    mask(:, :, s) = fliplr(mask(:, :, s));
    for r=1:3
        for c=1:3
            
            tensor_AL_LAS(:, :, s, r, c) = fliplr(tensor_LAS(:, :, s, r, c));
            tensor_AL_LPS(:, :, s, r, c) = fliplr(tensor_LPS(:, :, s, r, c));
            tensor_AL_RAS(:, :, s, r, c) = fliplr(tensor_RAS(:, :, s, r, c));
            
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
xlabel('Column Position (mm)')
ylabel('Row Position (mm)')
zlabel('Slice Position (mm)')
set(gca, 'CameraPosition', [188 -255 330], 'CameraTarget', [18 18 75], 'CameraViewAngle', 10, ...
    'CameraUpVector', [0 0 1], 'DataAspectRatio', [1 1 3], ...
    'xlim', [0 36], 'xtick', 0:9:36, 'ylim', [0 36], 'ytick', 0:9:36, 'zlim', [0 127], 'ztick', 0:21:126)
set(gcf, 'Position', [1000 100 700 700])
title('Aponeurosis Mesh')


%% fiber tracking

% Set fiber tracking options - LAS
ft_options_LAS.ref_frame = 'LAS';                                           %left-anterior-superior directions are +X, +Y, +Z
ft_options_LAS.image_orient = 'AL';                                         %image top is anterior, image right is left
ft_options_LAS.mesh_dist = 0;                                               %don’t shift the mesh
ft_options_LAS.prop_algo = 'euler';                                         %Euler integration
ft_options_LAS.step_size = 1;                                               %1 pixel width step
ft_options_LAS.term_mthd = 'bin2';                                          %BIN2 stop algorithm
ft_options_LAS.angle_thrsh = [25 3];                                        %>=25 degree inter-segment angle disallowed; look back three points
ft_options_LAS.fa_thrsh = [.1 .4];                                          %0.1<FA<0.4 range of FA values allowed
ft_options_LAS.depth_ratio = 1;                                             %ratio of ST/in-plane resolution of reconstructed images

% Set fiber tracking options - LPS
ft_options_LPS.ref_frame = 'LPS';                                           %left-posterior-superior directions are +X, +Y, +Z
ft_options_LPS.image_orient = 'AL';                                         %image top is anterior, image right is left
ft_options_LPS.mesh_dist = 0;                                               %don’t shift the mesh
ft_options_LPS.prop_algo = 'euler';                                         %Euler integration
ft_options_LPS.step_size = 1;                                               %1 pixel width step
ft_options_LPS.term_mthd = 'bin2';                                          %BIN2 stop algorithm
ft_options_LPS.angle_thrsh = [25 3];                                        %>=25 degree inter-segment angle disallowed; look back three points
ft_options_LPS.fa_thrsh = [.1 .4];                                          %0.1<FA<0.4 range of FA values allowed
ft_options_LPS.depth_ratio = 1;                                             %ratio of ST/in-plane resolution of reconstructed images

% Set fiber tracking options - RAS
ft_options_RAS.ref_frame = 'RAS';                                           %right-anterior-superior directions are +X, +Y, +Z
ft_options_RAS.image_orient = 'AL';                                         %image top is anterior, image right is left
ft_options_RAS.mesh_dist = 0;                                               %don’t shift the mesh
ft_options_RAS.prop_algo = 'euler';                                         %Euler integration
ft_options_RAS.step_size = 1;                                               %1 pixel width step
ft_options_RAS.term_mthd = 'bin2';                                          %BIN2 stop algorithm
ft_options_RAS.angle_thrsh = [25 3];                                        %>=25 degree inter-segment angle disallowed; look back three points
ft_options_RAS.fa_thrsh = [.1 .4];                                          %0.1<FA<0.4 range of FA values allowed
ft_options_RAS.depth_ratio = 1;                                             %ratio of ST/in-plane resolution of reconstructed images

% Set visualization options
plot_options.plot_fibers = 1;                                               %do plot any fiber tracts
plot_options.fiber_color = [.8 .2 .2];                                      %make the fibers red
plot_options.dti_size = [35 35];                                            %rows x columns of the DTI data
plot_options.dti_dims = [35 1];                                             %FOV and ST of the DTI data

% Fiber track - LAS:
[fiber_all_LAS, roi_flag_LAS, stop_list_LAS, fiber_length_LAS, FA_all_LAS, MD_LAS] = fiber_track...
    (tensor_LAS, mask, roi_mesh, ft_options_LAS);

% Fiber track - LPS:
[fiber_all_LPS, roi_flag_LPS, stop_list_LPS, fiber_length_LPS, FA_all_LPS, MD_LPS] = fiber_track...
    (tensor_LPS, mask, roi_mesh, ft_options_LPS);

% Fiber track - RAS:
[fiber_all_RAS, roi_flag_RAS, stop_list_RAS, fiber_length_RAS, FA_all_RAS, MD_RAS] = fiber_track...
    (tensor_RAS, mask, roi_mesh, ft_options_RAS);


% view the fiber tracts - LAS
fiber_figure_LAS = fiber_visualizer(mask, plot_options, roi_mesh, [], fiber_all_LAS);
xlabel('Column Position (mm)')
ylabel('Row Position (mm)')
zlabel('Slice Position (mm)')
set(gca, 'CameraPosition', [188 -255 330], 'CameraTarget', [18 18 75], 'CameraViewAngle', 10, ...
    'CameraUpVector', [0 0 1], 'DataAspectRatio', [1 1 3], ...
    'xlim', [0 36], 'xtick', 0:9:36, 'ylim', [0 36], 'ytick', 0:9:36, 'zlim', [0 127], 'ztick', 0:21:126)
set(gcf, 'Position', [1000 100 700 700])
title('Aponeurosis Mesh and Fiber Tracts (LAS)')

% view the fiber tracts - LPS
fiber_figure_LPS = fiber_visualizer(mask, plot_options, roi_mesh, [], fiber_all_LPS);
xlabel('Column Position (mm)')
ylabel('Row Position (mm)')
zlabel('Slice Position (mm)')
set(gca, 'CameraPosition', [188 -255 330], 'CameraTarget', [18 18 75], 'CameraViewAngle', 10, ...
    'CameraUpVector', [0 0 1], 'DataAspectRatio', [1 1 3], ...
    'xlim', [0 36], 'xtick', 0:9:36, 'ylim', [0 36], 'ytick', 0:9:36, 'zlim', [0 127], 'ztick', 0:21:126)
set(gcf, 'Position', [1000 100 700 700])
title('Aponeurosis Mesh and Fiber Tracts (LPS)')

% view the fiber tracts - RAS
fiber_figure_RAS = fiber_visualizer(mask, plot_options, roi_mesh, [], fiber_all_RAS);
xlabel('Column Position (mm)')
ylabel('Row Position (mm)')
zlabel('Slice Position (mm)')
set(gca, 'CameraPosition', [188 -255 330], 'CameraTarget', [18 18 75], 'CameraViewAngle', 10, ...
    'CameraUpVector', [0 0 1], 'DataAspectRatio', [1 1 3], ...
    'xlim', [0 36], 'xtick', 0:9:36, 'ylim', [0 36], 'ytick', 0:9:36, 'zlim', [0 127], 'ztick', 0:21:126)
set(gcf, 'Position', [1000 100 700 700])
title('Aponeurosis Mesh and Fiber Tracts (RAS)')


%% fiber quantifier

% set options
fq_options.dwi_res=[35 35 1];
fq_options.filt_kernel=3;
fq_options.mesh_units='vx';
fq_options.tract_units='vx';

% LAS
[angle_list_LAS, distance_list_LAS, curvature_list_LAS, fiber_all_mm_LAS, n_points_LAS, apo_area_LAS] = ...
     fiber_quantifier(fiber_all_LAS, roi_mesh, fq_options);
full_fiber_length_LAS = max(distance_list_LAS, [],3);

% LPS
[angle_list_LPS, distance_list_LPS, curvature_list_LPS, fiber_all_mm_LPS, n_points_LPS, apo_area_LPS] = ...
     fiber_quantifier(fiber_all_LPS, roi_mesh, fq_options);
full_fiber_length_LPS = max(distance_list_LPS, [],3);

% RAS
[angle_list_RAS, distance_list_RAS, curvature_list_RAS, fiber_all_mm_RAS, n_points_RAS, apo_area_RAS] = ...
     fiber_quantifier(fiber_all_RAS, roi_mesh, fq_options);
full_fiber_length_RAS = max(distance_list_RAS, [],3);


%% output data

clc

summary_LAS = [mean(angle_list_LAS(angle_list_LAS~=0)) std(angle_list_LAS(angle_list_LAS~=0)) 20 ...
mean(curvature_list_LAS(curvature_list_LAS~=0)) std(curvature_list_LAS(curvature_list_LAS~=0)) 0 ...
mean(full_fiber_length_LAS(full_fiber_length_LAS~=0)) std(full_fiber_length_LAS(full_fiber_length_LAS~=0)) 33/sind(20)]

summary_LPS = [mean(angle_list_LPS(angle_list_LPS~=0)) std(angle_list_LPS(angle_list_LPS~=0)) 20 ...
mean(curvature_list_LPS(curvature_list_LPS~=0)) std(curvature_list_LPS(curvature_list_LPS~=0)) 0 ...
mean(full_fiber_length_LPS(full_fiber_length_LPS~=0)) std(full_fiber_length_LPS(full_fiber_length_LPS~=0)) 33/sind(20)]

summary_RAS = [mean(angle_list_RAS(angle_list_RAS~=0)) std(angle_list_RAS(angle_list_RAS~=0)) 20 ...
mean(curvature_list_RAS(curvature_list_RAS~=0)) std(curvature_list_RAS(curvature_list_RAS~=0)) 0 ...
mean(full_fiber_length_RAS(full_fiber_length_RAS~=0)) std(full_fiber_length_RAS(full_fiber_length_RAS~=0)) 33/sind(20)]

close all
save Validation_FrameOfReference