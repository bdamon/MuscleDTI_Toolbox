%% ProjectionDirection
% v 1.0.0 (1/9/2021), Bruce Damon
% 
% This code simulates two sets of straight fibers that project upwards from an 
% aponeurosis. One set projects rightwards from an aponeurosis on the left 
% side and the other set projects anteriorly from a posteiorly placed aponeurosis 
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

% fibers go right and upward 
E1_RU = [-sind(20) 0 cosd(20)]';
E2_RU = [-cosd(20) 0 -sind(20)]';
E3_RU = [0 1 0]';
E_RU = [E1_RU E2_RU E3_RU];
D_RU = E_RU*L*E_RU';

% fibers go anterior and upward
E1_AU = [0 sind(20) cosd(20)]';
E2_AU = [0 -cosd(20) sind(20)]';
E3_AU = [1 0 0]';
E_AU = [E1_AU E2_AU E3_AU];
D_AU = E_AU*L*E_AU';


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
mask_RU = zeros(35, 35, 126);                                               %RU
mask_RU(2:34,2:34,2:125) = 1;                                               %main part of the muscle is square
mask_AU = zeros(35, 35, 126);                                               %AU
mask_AU(2:34,2:34,2:125) = 1;                                               

%in slice 1, create markers to track flip/rotation
mask_RU(1:35,2,1) = 2;                                                      %seed point from left side
mask_AU(34,1:35,1) = 2;                                                     %seed point from posterior

% signal generation loop
DWI_RU = zeros(35, 35, 126, 25);
DWI_AU = zeros(35, 35, 126, 25);
for d=1:25
    
    % get encoding direction
    loop_dir = squeeze(DTI_dir(d,:));
    
    % form images
    loop_signal_RU = 1*exp(-b_value*loop_dir*D_RU*loop_dir');
    DWI_RU(:,:,:,d) = loop_signal_RU;
    loop_signal_AU = 1*exp(-b_value*loop_dir*D_AU*loop_dir');
    DWI_AU(:,:,:,d) = loop_signal_AU;
    
    % mask the image
    DWI_RU(:,:,:,d) = DWI_RU(:,:,:,d).*mask_RU;
    DWI_AU(:,:,:,d) = DWI_AU(:,:,:,d).*mask_AU;
    
end

% form tensor matrix 
tensor_RU = zeros(35, 35, 126, 3, 3);
tensor_AU = zeros(35, 35, 126, 3, 3);
for r=1:35
    for c=1:35
        for s=1:126
            if mask_RU(r, c, s)==1
                
                loop_signal_RU = squeeze(DWI_RU(r, c, s, :));
                loop_D_RU = signal2tensor2(loop_signal_RU, DTI_dir(2:end,:), b_value);
                tensor_RU(r, c, s, :, :) = loop_D_RU;
                
                loop_signal_AU = squeeze(DWI_AU(r, c, s, :));
                loop_D_AU = signal2tensor2(loop_signal_AU, DTI_dir(2:end,:), b_value);
                tensor_AU(r, c, s, :, :) = loop_D_AU;
                
            end
        end
    end
end

% flip images to form AL image orientation
tensor_AL_RU = zeros(35, 35, 126, 3, 3);
tensor_AL_AU = zeros(35, 35, 126, 3, 3);
for s=1:126
    
    mask_RU(:, :, s) = fliplr(mask_RU(:, :, s));
    mask_AU(:, :, s) = fliplr(mask_AU(:, :, s));
    
    for r=1:3
        for c=1:3
            
            tensor_AL_RU(:, :, s, r, c) = fliplr(tensor_RU(:, :, s, r, c));
            tensor_AL_AU(:, :, s, r, c) = fliplr(tensor_AU(:, :, s, r, c));
            
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

% roi_mesh
roi_mesh_RU=zeros(15,33,6);
for r=1:15
    roi_mesh_RU(r,1:33,1)=2:34;
    roi_mesh_RU(r,1:33,3)= r+13;
end
roi_mesh_RU(:,1:33,2)=34;
roi_mesh_RU(:,:,1:3) = ...              %need to add a little randomness to prevent /0 in fiber_quantifier
    roi_mesh_RU(:,:,1:3)+randn(size(squeeze(roi_mesh_RU(:,:,1:3))))*0.000001;

% find normal to mesh at each point:
mesh_row_vec = circshift(roi_mesh_RU(:, :, 1:3), [0 -1 0]) - roi_mesh_RU(:, :, 1:3);
mesh_col_vec = circshift(roi_mesh_RU(:, :, 1:3), [-1 0 0]) - roi_mesh_RU(:, :, 1:3);
mesh_row_vec(:, end, :) = mesh_row_vec(:, end-1, :);
mesh_col_vec(end, :, :) = mesh_col_vec(end-1, :, :);

roi_mesh_RU(:, :, 4:6) = cross(mesh_row_vec, mesh_col_vec, 3);
roi_mesh_RU(:, :, 4:6) = smooth3(roi_mesh_RU(:, :, 4:6));
roi_norm = (roi_mesh_RU(:, :, 4).^2 + roi_mesh_RU(:, :, 5).^2 + roi_mesh_RU(:, :,6).^2).^0.5;
roi_mesh_RU(:, :, 4) = roi_mesh_RU(:, :, 4)./roi_norm;
roi_mesh_RU(:, :, 5) = roi_mesh_RU(:, :, 5)./roi_norm;
roi_mesh_RU(:, :, 6) = roi_mesh_RU(:, :, 6)./roi_norm;
roi_mesh_RU(:, :, 4:6) = smooth3(roi_mesh_RU(:, :, 4:6));
roi_norm = (roi_mesh_RU(:, :, 4).^2 + roi_mesh_RU(:, :, 5).^2 + roi_mesh_RU(:, :,6).^2).^0.5;
roi_mesh_RU(:, :, 4) = roi_mesh_RU(:, :, 4)./roi_norm;
roi_mesh_RU(:, :, 5) = roi_mesh_RU(:, :, 5)./roi_norm;
roi_mesh_RU(:, :, 6) = roi_mesh_RU(:, :, 6)./roi_norm;

% view the mesh
mesh_figure_RU = fiber_visualizer(mask_RU, plot_options, roi_mesh_RU, [], []);
xlabel('Column Position (mm)')
ylabel('Row Position (mm)')
zlabel('Slice Position (mm)')
set(gca, 'CameraPosition', [188 -255 330], 'CameraTarget', [18 18 75], 'CameraViewAngle', 10, ...
    'CameraUpVector', [0 0 1], 'DataAspectRatio', [1 1 3], ...
    'xlim', [0 36], 'xtick', 0:9:36, 'ylim', [0 36], 'ytick', 0:9:36, 'zlim', [0 127], 'ztick', 0:21:126)
set(gcf, 'Position', [1000 100 700 700])
title('Aponeurosis Mesh')


% roi_mesh
roi_mesh_AU=zeros(15,33,6);
for r=1:15
    roi_mesh_AU(r,1:33,2)=2:34;
    roi_mesh_AU(r,1:33,3)= r+13;
end
roi_mesh_AU(:,1:33,1)=34;
roi_mesh_AU(:,:,1:3) = ...              %need to add a little randomness to prevent /0 in fiber_quantifier
    roi_mesh_AU(:,:,1:3)+randn(size(squeeze(roi_mesh_AU(:,:,1:3))))*0.000001;

% find normal to mesh at each point:
mesh_row_vec = circshift(roi_mesh_AU(:, :, 1:3), [0 -1 0]) - roi_mesh_AU(:, :, 1:3);
mesh_col_vec = circshift(roi_mesh_AU(:, :, 1:3), [-1 0 0]) - roi_mesh_AU(:, :, 1:3);
mesh_row_vec(:, end, :) = mesh_row_vec(:, end-1, :);
mesh_col_vec(end, :, :) = mesh_col_vec(end-1, :, :);

roi_mesh_AU(:, :, 4:6) = cross(mesh_row_vec, mesh_col_vec, 3);
roi_mesh_AU(:, :, 4:6) = smooth3(roi_mesh_AU(:, :, 4:6));
roi_norm = (roi_mesh_AU(:, :, 4).^2 + roi_mesh_AU(:, :, 5).^2 + roi_mesh_AU(:, :,6).^2).^0.5;
roi_mesh_AU(:, :, 4) = roi_mesh_AU(:, :, 4)./roi_norm;
roi_mesh_AU(:, :, 5) = roi_mesh_AU(:, :, 5)./roi_norm;
roi_mesh_AU(:, :, 6) = roi_mesh_AU(:, :, 6)./roi_norm;
roi_mesh_AU(:, :, 4:6) = smooth3(roi_mesh_AU(:, :, 4:6));
roi_norm = (roi_mesh_AU(:, :, 4).^2 + roi_mesh_AU(:, :, 5).^2 + roi_mesh_AU(:, :,6).^2).^0.5;
roi_mesh_AU(:, :, 4) = roi_mesh_AU(:, :, 4)./roi_norm;
roi_mesh_AU(:, :, 5) = roi_mesh_AU(:, :, 5)./roi_norm;
roi_mesh_AU(:, :, 6) = roi_mesh_AU(:, :, 6)./roi_norm;

% view the mesh
mesh_figure_AU = fiber_visualizer(mask_AU, plot_options, roi_mesh_AU, [], []);
set(gca, 'CameraPosition', [339 -553 262], 'CameraTarget', [18 18 64], 'CameraViewAngle', 10, ...
    'CameraUpVector', [0 0 1], 'DataAspectRatio', [1 1 1], ...
    'xlim', [-.01 36.01], 'xtick', 0:9:36, 'ylim', [0 36.01], 'ytick', 0:9:36, 'zlim', [0 127], 'ztick', 0:14:126, ...
    'fontsize', 12)
xlabel('Column Position (mm)', 'fontweight', 'bold', 'fontsize', 14)
ylabel('Row Position (mm)', 'fontweight', 'bold', 'fontsize', 14)
zlabel('Slice Position (mm)', 'fontweight', 'bold', 'fontsize', 14)
set(gcf, 'Position', [1000 100 700 850])
title('Aponeurosis Mesh')


%% fiber tracking: AR1, AL, RK4

% Set fiber tracking options - LAS
ft_options.ref_frame = 'LAS';                                            %left-anterior-superior directions are +X, +Y, +Z
ft_options.image_orient = 'AL';                                          %image top is anterior, image right is left
ft_options.mesh_dist = 0;                                                %don’t shift the mesh
ft_options.prop_algo = 'euler';                                          %Euler integration
ft_options.step_size = 1;                                                %1 pixel width step
ft_options.term_mthd = 'bin2';                                           %BIN2 stop algorithm
ft_options.angle_thrsh = [25 3];                                         %>=25 degree inter-segment angle disallowed; look back three points
ft_options.fa_thrsh = [.1 .4];                                           %0.1<FA<0.4 range of FA values allowed
ft_options.depth_ratio = 1;                                              %ratio of ST/in-plane resolution of reconstructed images
 
% Set visualization options
plot_options.plot_fibers = 1;                                               %do plot any fiber tracts
plot_options.fiber_color = [.8 .2 .2];                                      %make the fibers red
plot_options.dti_size = [35 35];                                            %rows x columns of the DTI data
plot_options.dti_dims = [35 1];                                             %FOV and ST of the DTI data

% Fiber track - RU
[fiber_all_RU, roi_flag_RU, stop_list_RU, fiber_length_RU, FA_all_RU, MD_RU] = fiber_track...
    (tensor_AL_RU, mask_RU, roi_mesh_RU, ft_options);

% Fiber track - AU
[fiber_all_AU, roi_flag_AU, stop_list_AU, fiber_length_AU, FA_all_AU, MD_AU] = fiber_track...
    (tensor_AL_AU, mask_AU, roi_mesh_AU, ft_options);

% view the fiber tracts
fiber_figure_RU = fiber_visualizer(mask_RU, plot_options, roi_mesh_RU, [], fiber_all_RU);
set(gca, 'CameraPosition', [339 -553 262], 'CameraTarget', [18 18 64], 'CameraViewAngle', 10, ...
    'CameraUpVector', [0 0 1], 'DataAspectRatio', [1 1 1], ...
    'xlim', [-.01 36.01], 'xtick', 0:9:36, 'ylim', [0 36.01], 'ytick', 0:9:36, 'zlim', [0 127], 'ztick', 0:14:126, ...
    'fontsize', 12)
xlabel('Column Position (mm)', 'fontweight', 'bold', 'fontsize', 14)
ylabel('Row Position (mm)', 'fontweight', 'bold', 'fontsize', 14)
zlabel('Slice Position (mm)', 'fontweight', 'bold', 'fontsize', 14)
set(gcf, 'Position', [1000 100 700 850])
title('Aponeurosis Mesh and Fiber Tracts (RU)')

% view the fiber tracts
fiber_figure_AU = fiber_visualizer(mask_AU, plot_options, roi_mesh_AU, [], fiber_all_AU);
set(gca, 'CameraPosition', [339 -553 262], 'CameraTarget', [18 18 64], 'CameraViewAngle', 10, ...
    'CameraUpVector', [0 0 1], 'DataAspectRatio', [1 1 1], ...
    'xlim', [-.01 36.01], 'xtick', 0:9:36, 'ylim', [0 36.01], 'ytick', 0:9:36, 'zlim', [0 127], 'ztick', 0:14:126, ...
    'fontsize', 12)
xlabel('Column Position (mm)', 'fontweight', 'bold', 'fontsize', 14)
ylabel('Row Position (mm)', 'fontweight', 'bold', 'fontsize', 14)
zlabel('Slice Position (mm)', 'fontweight', 'bold', 'fontsize', 14)
set(gcf, 'Position', [1000 100 700 850])
title('Aponeurosis Mesh and Fiber Tracts (AU)')


%% fiber quantifier: AR1, RK4

% set options
fq_options.dwi_res=[35 35 1];
fq_options.filt_kernel=3;
fq_options.mesh_units='vx';
fq_options.tract_units='vx';

% RU
[angle_list_RU, distance_list_RU, curvature_list_RU, fiber_all_mm_RU, n_points_RU, apo_area_RU] = ...
     fiber_quantifier(fiber_all_RU, roi_mesh_RU, fq_options);
full_fiber_length_RU = max(distance_list_RU, [],3);

% AU
[angle_list_AU, distance_list_AU, curvature_list_AU, fiber_all_mm_AU, n_points_AU, apo_area_AU] = ...
     fiber_quantifier(fiber_all_AU, roi_mesh_AU, fq_options);
full_fiber_length_AU = max(distance_list_AU, [],3);


%% output data

clc

summary_RU = [mean(angle_list_RU(angle_list_RU~=0)) std(angle_list_RU(angle_list_RU~=0)) 20 ...
mean(curvature_list_RU(curvature_list_RU~=0)) std(curvature_list_RU(curvature_list_RU~=0)) 0 ...
mean(full_fiber_length_RU(full_fiber_length_RU~=0)) std(full_fiber_length_RU(full_fiber_length_RU~=0)) 33/sind(20)]

summary_AU = [mean(angle_list_AU(angle_list_AU~=0)) std(angle_list_AU(angle_list_AU~=0)) 20 ...
mean(curvature_list_AU(curvature_list_AU~=0)) std(curvature_list_AU(curvature_list_AU~=0)) 0 ...
mean(full_fiber_length_AU(full_fiber_length_AU~=0)) std(full_fiber_length_AU(full_fiber_length_AU~=0)) 33/sind(20)]

close all
save Validation_ProjectionDirection