%% PennationLength
% v 1.0.0 (1/9/2021), Bruce Damon
% 
% This code simulates two sets of straight fibers that emerge from an 
% aponeurosis on the anatomical left side. One set of fibers projects right
% and upward at an angle of 30 degrees.  The other set projects right and
% upward at an angle of 20 degrees.  
% 
% The known values for pennation angle, curvature, and fiber length are 
% 30/0/66 and 20/0/96.49 for the two sets fo fibers.  The units are degrees,
% m^-1, and mm, respectively.
%
% In both sets of simulations, the frame of reference for the diffusion 
% directions is LAS.
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

% fibers go right and upward at 30 degrees in an LAS frame of reference
E1_30 = [-sind(30) 0 cosd(30)]';
E2_30 = [-cosd(30) 0 -sind(30)]';
E3_30 = [0 1 0]';
E_30 = [E1_30 E2_30 E3_30];
D_30 = E_30*L*E_30';

% fibers go right and upward at 20 degrees
E1_20 = [-sind(20) 0 cosd(20)]';
E2_20 = [-cosd(20) 0 -sind(20)]';
E3_20 = [0 1 0]';
E_20 = [E1_20 E2_20 E3_20];
D_20 = E_20*L*E_20';


%% set DTI encoding parameters

DTI_dir = [0    0   0
    0.1890    0.6940    0.6940                                              % philips 24 direction DTI encoding scheme:LAS
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
b_value = 450;                                                              % b-value, s/mm^2

%% form isotropic diffusion-weighted images

% form mask
mask = zeros(35, 35, 126);
mask(2:34,2:34,2:125) = 1;      %main part of the muscle is square

%in slice 1, create markers to track flip/rotation
mask(1:35,2,1) = 2;             %seed point from A-P on left side

% signal generation loop
dwi_30 = zeros(35, 35, 126, 25);                                            %DT images for 30 degree pennation angle
dwi_20 = zeros(35, 35, 126, 25);                                            %DT images for 20 degree pennation angle
for d=1:25
    
    % get encoding direction
    Loop_DTI_Dir = squeeze(DTI_dir(d,:));
    
    % form images
    loop_signal_30 = 1*exp(-b_value*Loop_DTI_Dir*D_30*Loop_DTI_Dir');       %generate signals - 30 degrees
    dwi_30(:,:,:,d) = loop_signal_30;                                       %into the DWI matrix

    loop_signal_20 = 1*exp(-b_value*Loop_DTI_Dir*D_20*Loop_DTI_Dir');       %generate signals - 20 degrees
    dwi_20(:,:,:,d) = loop_signal_20;                                       %into the DWI matrix

    % mask the images
    dwi_30(:,:,:,d) = dwi_30(:,:,:,d).*mask;
    dwi_20(:,:,:,d) = dwi_20(:,:,:,d).*mask;
    
end

% form tensor matrix
tensor_30 = zeros(35, 35, 126, 3, 3);
tensor_20 = zeros(35, 35, 126, 3, 3);
for r=1:35
    for c=1:35
        for s=1:126
            if mask(r, c, s)==1
                
                loop_signal_30 = squeeze(dwi_30(r, c, s, :));
                loop_D = signal2tensor2(loop_signal_30, DTI_dir(2:end,:), b_value);
                tensor_30(r, c, s, :, :) = loop_D;
                
                loop_signal_20 = squeeze(dwi_20(r, c, s, :));
                loop_D = signal2tensor2(loop_signal_20, DTI_dir(2:end,:), b_value);
                tensor_20(r, c, s, :, :) = loop_D;
                
            end
        end
    end
end

% flip images to form AL image orientation
tensor_AL_30 = zeros(35, 35, 126, 3, 3);
tensor_AL_20 = zeros(35, 35, 126, 3, 3);
for s=1:126
    
    mask(:, :, s) = fliplr(mask(:, :, s));
    
    for r=1:3
        for c=1:3
            
            tensor_AL_30(:, :, s, r, c) = fliplr(tensor_30(:, :, s, r, c));
            tensor_AL_20(:, :, s, r, c) = fliplr(tensor_20(:, :, s, r, c));
            
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

% Set fiber tracking options
ft_options.ref_frame = 'LAS';                                           %left-anterior-superior directions are +X, +Y, +Z
ft_options.image_orient = 'AL';                                         %image top is anterior, image right is left
ft_options.mesh_dist = 0;                                               %don’t shift the mesh
ft_options.prop_algo = 'euler';                                         %Euler integration
ft_options.step_size = 1;                                               %1 pixel width step
ft_options.term_mthd = 'bin2';                                          %BIN2 stop algorithm
ft_options.angle_thrsh = [25 3];                                        %>=25 degree inter-segment angle disallowed; look back three points
ft_options.fa_thrsh = [.1 .4];                                          %0.1<FA<0.4 range of FA values allowed
ft_options.depth_ratio = 1;                                             %ratio of ST/in-plane resolution of reconstructed images

% Set visualization options
plot_options.plot_fibers = 1;                                               %do plot any fiber tracts
plot_options.fiber_color = [.8 .2 .2];                                      %make the fibers red
plot_options.dti_size = [35 35];                                            %rows x columns of the DTI data
plot_options.dti_dims = [35 1];                                             %FOV and ST of the DTI data

% Fiber track - 30 degrees
[fiberall_30, roiflag_30, stoplist_30, fiber_length_30, FA_30, MD_30] = fiber_track...
    (tensor_AL_30, mask, roi_mesh, ft_options);

% Fiber track - 20 degrees
[fiberall_20, roiflag_20, stoplist_20, fiberlength_20, FA_20, MD_20] = fiber_track...
    (tensor_AL_20, mask, roi_mesh, ft_options);

% view the fiber tracts - 30 degrees
fiber_figure_30 = fiber_visualizer(mask, plot_options, roi_mesh, [], fiberall_30);
set(gca, 'CameraPosition', [339 -553 262], 'CameraTarget', [18 18 64], 'CameraViewAngle', 10, ...
    'CameraUpVector', [0 0 1], 'DataAspectRatio', [1 1 1], ...
    'xlim', [-.01 36.01], 'xtick', 0:9:36, 'ylim', [-.01 36.01], 'ytick', 0:9:36, 'zlim', [-.01 126.01], 'ztick', 0:14:126, ...
    'fontsize', 12)
xlabel('Column Position (mm)', 'fontweight', 'bold', 'fontsize', 14)
ylabel('Row Position (mm)', 'fontweight', 'bold', 'fontsize', 14)
zlabel('Slice Position (mm)', 'fontweight', 'bold', 'fontsize', 14)
set(gcf, 'Position', [1000 100 700 850])
title('Aponeurosis Mesh and Fiber Tracts (30^o)')
axis image

% view the fiber tracts
fiber_figure_20 = fiber_visualizer(mask, plot_options, roi_mesh, [], fiberall_20);
set(gca, 'CameraPosition', [339 -553 262], 'CameraTarget', [18 18 64], 'CameraViewAngle', 10, ...
    'CameraUpVector', [0 0 1], 'DataAspectRatio', [1 1 1], ...
    'xlim', [-.01 36.01], 'xtick', 0:9:36, 'ylim', [0 36.01], 'ytick', 0:9:36, 'zlim', [0 127], 'ztick', 0:14:126, ...
    'fontsize', 12)
xlabel('Column Position (mm)', 'fontweight', 'bold', 'fontsize', 14)
ylabel('Row Position (mm)', 'fontweight', 'bold', 'fontsize', 14)
zlabel('Slice Position (mm)', 'fontweight', 'bold', 'fontsize', 14)
set(gcf, 'Position', [1000 100 700 850])
title('Aponeurosis Mesh and Fiber Tracts (20^o)')
axis image


%% Quantify fiber tract properties

% set options
fq_options.dwi_res=[35 35 1];
fq_options.filt_kernel=3;
fq_options.mesh_units='vx';
fq_options.tract_units='vx';

% 30 degrees
[angle_list_30, distance_list_30, curvature_list_30, fiberall_mm_30, n_points_30, apo_area_30] = ...
     fiber_quantifier(fiberall_30, roi_mesh, fq_options);
full_fiber_length_30 = max(distance_list_30, [],3);

% 20 degrees
[angle_list_20, distance_list_20, curvature_list_20, fiber_all_mm_20, n_points_20, apo_area_20] = ...
     fiber_quantifier(fiberall_20, roi_mesh, fq_options);
full_fiber_length_20 = max(distance_list_20, [],3);

%% output data

clc

summary_20 = [mean(angle_list_20(angle_list_20~=0)) std(angle_list_20(angle_list_20~=0)) 30 ...
mean(curvature_list_20(curvature_list_20~=0)) std(curvature_list_20(curvature_list_20~=0)) 0 ...
mean(full_fiber_length_20(full_fiber_length_20~=0)) std(full_fiber_length_20(full_fiber_length_20~=0)) 33/sind(20)]

summary_30 = [mean(angle_list_30(angle_list_30~=0)) std(angle_list_30(angle_list_30~=0)) 30 ...
mean(curvature_list_30(curvature_list_30~=0)) std(curvature_list_30(curvature_list_30~=0)) 0 ...
mean(full_fiber_length_30(full_fiber_length_30~=0)) std(full_fiber_length_30(full_fiber_length_30~=0)) 33/sind(30)]

close all
save Validation_PennationLength
