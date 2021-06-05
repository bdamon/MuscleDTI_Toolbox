%% Image Orientation
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
% The image orientations are RA and AL.
%
% Voxels are simulated at isotropic resolution (1x1x1 mm^3)
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

% fibers go right and upward at 20 degrees in an LAS frame of reference
E1 = [0 sind(20) cosd(20)]';
E2 = [0 -cosd(20) sind(20)]';
E3 = [1 0 0]';
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
mask(34,1:35,1) = 2;             %seed point from A-P on left side

% signal generation loop
dwi = zeros(35, 35, 126, 25);
for d=1:25
    
    % get encoding direction
    loop_DTI_dir = squeeze(DTI_dir(d,:));
    
    % form images
    loop_signal = 1*exp(-b_value*loop_DTI_dir*D*loop_DTI_dir');
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
mask_AL = zeros(35, 35, 126);
for s=1:126
    mask_AL(:, :, s) = fliplr(mask(:, :, s));
    for r=1:3
        for c=1:3
            tensor_AL(:, :, s, r, c) = fliplr(tensor(:, :, s, r, c));
        end
    end
end

% rotate images to form RA image orientation
tensor_RA = zeros(35, 35, 126, 3, 3);
mask_RA = zeros(35, 35, 126);
for s=1:126
    mask_RA(:, :, s) = imrotate(mask_AL(:, :, s), -90);
    for r=1:3
        for c=1:3
            tensor_RA(:,:,s, r, c) = imrotate(tensor_AL(:, :, s, r, c), -90);
        end
    end
end


%% form roi_meshes

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

% roi_mesh - AL
roi_mesh_AL=zeros(20,33,6);
for r=1:20
    roi_mesh_AL(r,1:33,2)=2:34;
    roi_mesh_AL(r,1:33,3)= r+13;
end
roi_mesh_AL(:,1:33,1)=34;
roi_mesh_AL(:,:,1:3) = ...              %need to add a little randomness to prevent /0 in fiber_quantifier
    roi_mesh_AL(:,:,1:3)+randn(size(squeeze(roi_mesh_AL(:,:,1:3))))*0.000001;

% find normal to mesh at each point:
mesh_row_vec = circshift(roi_mesh_AL(:, :, 1:3), [0 -1 0]) - roi_mesh_AL(:, :, 1:3);
mesh_col_vec = circshift(roi_mesh_AL(:, :, 1:3), [-1 0 0]) - roi_mesh_AL(:, :, 1:3);
mesh_row_vec(:, end, :) = mesh_row_vec(:, end-1, :);
mesh_col_vec(end, :, :) = mesh_col_vec(end-1, :, :);

roi_mesh_AL(:, :, 4:6) = cross(mesh_row_vec, mesh_col_vec, 3);
roi_mesh_AL(:, :, 4:6) = smooth3(roi_mesh_AL(:, :, 4:6));
roi_norm = (roi_mesh_AL(:, :, 4).^2 + roi_mesh_AL(:, :, 5).^2 + roi_mesh_AL(:, :,6).^2).^0.5;
roi_mesh_AL(:, :, 4) = roi_mesh_AL(:, :, 4)./roi_norm;
roi_mesh_AL(:, :, 5) = roi_mesh_AL(:, :, 5)./roi_norm;
roi_mesh_AL(:, :, 6) = roi_mesh_AL(:, :, 6)./roi_norm;
roi_mesh_AL(:, :, 4:6) = smooth3(roi_mesh_AL(:, :, 4:6));
roi_norm = (roi_mesh_AL(:, :, 4).^2 + roi_mesh_AL(:, :, 5).^2 + roi_mesh_AL(:, :,6).^2).^0.5;
roi_mesh_AL(:, :, 4) = roi_mesh_AL(:, :, 4)./roi_norm;
roi_mesh_AL(:, :, 5) = roi_mesh_AL(:, :, 5)./roi_norm;
roi_mesh_AL(:, :, 6) = roi_mesh_AL(:, :, 6)./roi_norm;

% view the mesh
mesh_figure_AL = fiber_visualizer(mask, plot_options, roi_mesh_AL, [], []);
xlabel('Column Position (mm)')
ylabel('Row Position (mm)')
zlabel('Slice Position (mm)')
set(gca, 'CameraPosition', [188 -255 330], 'CameraTarget', [18 18 75], 'CameraViewAngle', 10, ...
    'CameraUpVector', [0 0 1], 'DataAspectRatio', [1 1 3], ...
    'xlim', [0 36], 'xtick', 0:9:36, 'ylim', [0 36], 'ytick', 0:9:36, 'zlim', [0 127], 'ztick', 0:21:126)
set(gcf, 'Position', [1000 100 700 700])
title('Aponeurosis Mesh (AL)')

% roi_mes - RA
roi_mesh_RA = zeros(20,33,6);
for r=1:20
    roi_mesh_RA(r,1:33,1)=2:34;
    roi_mesh_RA(r,1:33,3)= r+13;
end
roi_mesh_RA(:,1:33,2)=2;
roi_mesh_RA(:,:,1:3) = ...              %need to add a little randomness to prevent /0 in fiber_quantifier
    roi_mesh_RA(:,:,1:3)+randn(size(squeeze(roi_mesh_RA(:,:,1:3))))*0.000001;

% find normal to mesh at each point:
mesh_row_vec = circshift(roi_mesh_RA(:, :, 1:3), [0 -1 0]) - roi_mesh_RA(:, :, 1:3);
mesh_col_vec = circshift(roi_mesh_RA(:, :, 1:3), [-1 0 0]) - roi_mesh_RA(:, :, 1:3);
mesh_row_vec(:, end, :) = mesh_row_vec(:, end-1, :);
mesh_col_vec(end, :, :) = mesh_col_vec(end-1, :, :);

roi_mesh_RA(:, :, 4:6) = cross(mesh_row_vec, mesh_col_vec, 3);
roi_mesh_RA(:, :, 4:6) = smooth3(roi_mesh_RA(:, :, 4:6));
roi_norm = (roi_mesh_RA(:, :, 4).^2 + roi_mesh_RA(:, :, 5).^2 + roi_mesh_RA(:, :,6).^2).^0.5;
roi_mesh_RA(:, :, 4) = roi_mesh_RA(:, :, 4)./roi_norm;
roi_mesh_RA(:, :, 5) = roi_mesh_RA(:, :, 5)./roi_norm;
roi_mesh_RA(:, :, 6) = roi_mesh_RA(:, :, 6)./roi_norm;
roi_mesh_RA(:, :, 4:6) = smooth3(roi_mesh_RA(:, :, 4:6));
roi_norm = (roi_mesh_RA(:, :, 4).^2 + roi_mesh_RA(:, :, 5).^2 + roi_mesh_RA(:, :,6).^2).^0.5;
roi_mesh_RA(:, :, 4) = roi_mesh_RA(:, :, 4)./roi_norm;
roi_mesh_RA(:, :, 5) = roi_mesh_RA(:, :, 5)./roi_norm;
roi_mesh_RA(:, :, 6) = roi_mesh_RA(:, :, 6)./roi_norm;

% view the mesh
mesh_figure_RA = fiber_visualizer(mask, plot_options, roi_mesh_RA, [], []);
xlabel('Column Position (mm)')
ylabel('Row Position (mm)')
zlabel('Slice Position (mm)')
set(gca, 'CameraPosition', [188 -255 330], 'CameraTarget', [18 18 75], 'CameraViewAngle', 10, ...
    'CameraUpVector', [0 0 1], 'DataAspectRatio', [1 1 3], ...
    'xlim', [0 36], 'xtick', 0:9:36, 'ylim', [0 36], 'ytick', 0:9:36, 'zlim', [0 127], 'ztick', 0:21:126)
set(gcf, 'Position', [1000 100 700 700])
title('Aponeurosis Mesh (RA)')


%% fiber tracking 

% Set fiber tracking options
ft_options_AL.ref_frame = 'LAS';                                           %left-anterior-superior directions are +X, +Y, +Z
ft_options_AL.image_orient = 'AL';                                         %image top is anterior, image right is left
ft_options_AL.mesh_dist = 0;                                               %don’t shift the mesh
ft_options_AL.prop_algo = 'euler';                                         %Euler integration
ft_options_AL.step_size = 1;                                               %1 pixel width step
ft_options_AL.term_mthd = 'bin2';                                          %BIN2 stop algorithm
ft_options_AL.angle_thrsh = [25 3];                                        %>=25 degree inter-segment angle disallowed; look back three points
ft_options_AL.fa_thrsh = [.1 .4];                                          %0.1<FA<0.4 range of FA values allowed
ft_options_AL.depth_ratio = 1;                                             %ratio of ST/in-plane resolution of reconstructed images

% Set fiber tracking options - RA
ft_options_RA.ref_frame = 'LAS';                                           %left-anterior-superior directions are +X, +Y, +Z
ft_options_RA.image_orient = 'RA';                                         %image top is anterior, image right is left
ft_options_RA.mesh_dist = 0;                                               %don’t shift the mesh
ft_options_RA.prop_algo = 'euler';                                         %Euler integration
ft_options_RA.step_size = 1;                                               %1 pixel width step
ft_options_RA.term_mthd = 'bin2';                                          %BIN2 stop algorithm
ft_options_RA.angle_thrsh = [25 3];                                        %>=25 degree inter-segment angle disallowed; look back three points
ft_options_RA.fa_thrsh = [.1 .4];                                          %0.1<FA<0.4 range of FA values allowed
ft_options_RA.depth_ratio = 1;                                             %ratio of ST/in-plane resolution of reconstructed images

% Set visualization options
plot_options.plot_fibers = 1;                                               %do plot any fiber tracts
plot_options.fiber_color = [.8 .2 .2];                                      %make the fibers red
plot_options.dti_size = [35 35];                                            %rows x columns of the DTI data
plot_options.dti_dims = [35 1];                                             %FOV and ST of the DTI data

% Fiber track - AL
[fiber_all_AL, roi_flag_AL_RK4_AR1, stop_list_AL, fiber_length_AL, FA_all_AL, MD_AL] = fiber_track...
    (tensor_AL, mask, roi_mesh_AL, ft_options_AL);

% Fiber track - RA
[fiber_all_RA, roi_flag_RA, stop_list_RA, fiber_length_RA, FA_all_RA, MD_RA] = fiber_track...
    (tensor_RA, mask, roi_mesh_RA, ft_options_RA);

% view the fiber tracts
fiber_figure_AL = fiber_visualizer(mask, plot_options, roi_mesh_AL, [], fiber_all_AL);
xlabel('Column Position (mm)')
ylabel('Row Position (mm)')
zlabel('Slice Position (mm)')
set(gca, 'CameraPosition', [188 -255 330], 'CameraTarget', [18 18 75], 'CameraViewAngle', 10, ...
    'CameraUpVector', [0 0 1], 'DataAspectRatio', [1 1 3], ...
    'xlim', [0 36], 'xtick', 0:9:36, 'ylim', [0 36], 'ytick', 0:9:36, 'zlim', [0 127], 'ztick', 0:21:126)
set(gcf, 'Position', [1000 100 700 700])
title('Aponeurosis Mesh and Fiber Tracts (LAS, AU, AL)')


% view the fiber tracts
fiber_figure_RA = fiber_visualizer(mask, plot_options, roi_mesh_RA, [], fiber_all_RA);
xlabel('Column Position (mm)')
ylabel('Row Position (mm)')
zlabel('Slice Position (mm)')
set(gca, 'CameraPosition', [188 -255 330], 'CameraTarget', [18 18 75], 'CameraViewAngle', 10, ...
    'CameraUpVector', [0 0 1], 'DataAspectRatio', [1 1 3], ...
    'xlim', [0 36], 'xtick', 0:9:36, 'ylim', [0 36], 'ytick', 0:9:36, 'zlim', [0 127], 'ztick', 0:21:126)
set(gcf, 'Position', [1000 100 700 700])
title('Aponeurosis Mesh and Fiber Tracts (LAS, AU, RA)')


%% fiber quantifier

% set options
fq_options.dwi_res=[35 35 1];
fq_options.filt_kernel=3;
fq_options.mesh_units='vx';
fq_options.tract_units='vx';

% AL
[angle_list_AL, distance_list_AL, curvature_list_AL, fiber_all_mm_AL, n_points_AL, apo_area_AL] = ...
     fiber_quantifier(fiber_all_AL, roi_mesh_AL, fq_options);
full_fiber_length_AL = max(distance_list_AL, [],3);

% RA
[angle_list_RA, distance_list_RA, curvature_list_RA, fiber_RA_mm_RA, n_points_RA, apo_area_RA] = ...
     fiber_quantifier(fiber_all_RA, roi_mesh_RA, fq_options);
full_fiber_length_RA = max(distance_list_RA, [],3);



%% output data

clc

summary_AL = [mean(angle_list_AL(angle_list_AL~=0)) std(angle_list_AL(angle_list_AL~=0)) 20 ...
mean(curvature_list_AL(curvature_list_AL~=0)) std(curvature_list_AL(curvature_list_AL~=0)) 0 ...
mean(full_fiber_length_AL(full_fiber_length_AL~=0)) std(full_fiber_length_AL(full_fiber_length_AL~=0)) 33/cosd(70)]

summary_RA = [mean(angle_list_RA(angle_list_RA~=0)) std(angle_list_RA(angle_list_RA~=0)) 20 ...
mean(curvature_list_RA(curvature_list_RA~=0)) std(curvature_list_RA(curvature_list_RA~=0)) 0 ... 
mean(full_fiber_length_RA(full_fiber_length_RA~=0)) std(full_fiber_length_RA(full_fiber_length_RA~=0)) 33/cosd(70)]

save Validation