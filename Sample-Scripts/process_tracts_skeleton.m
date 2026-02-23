%% process_tracts.m 
% This version of process_tracts is intended to be used with version 2 of the fiber_track function.
% Completed on 23 Oct 2025, Bruce Damon

%% Clean slate:
clear
close all
clc

%% File Input

% prompt user to select pre-processed data file/load
[pre_fiber_tracking_file, pre_fiber_tracking_path] = uigetfile('*.mat');
cd(pre_fiber_tracking_path)
load(pre_fiber_tracking_file)

% Load muscle mask and aponeurosis ROI
load ta_mask mask
load ta_mesh roi_*

%% view the mask and mesh

% set fiber visualizer options for the mask
fv_options.anat_dims = [anat_fov(1) anat_slcthick];                         %FOV and slice thickness of the images to be displayed, in mm
fv_options.anat_slices = 5:9:44;                                            %display slices 14, 24, 34, and 44
fv_options.plot_mesh = 0;                                                   %don’t plot an aponeurosis mesh
fv_options.plot_mask = 1;                                                   %do plot the mask
fv_options.plot_fibers = 0;                                                 %don’t plot any fiber tracts
fv_options.mask_size = [anat_numrows anat_numcols];                         %rows x columns of the images used to generate the mask
fv_options.mask_dims = [anat_fov(1) anat_slcthick];                         %FOV and slice thickness of the images used to create the mask, in mm
fv_options.mask_color = [1 0 0];                                            %make the mask a red, semi-transparent overlay

% visualize the mask:
mask_figure = fiber_visualizer(anat_image, fv_options, [], mask, []);
xlabel('Left-Right Direction (mm)')
ylabel('Anterior-Posterior Direction (mm)')
zlabel('Head-Foot Direction (mm)')

% set fiber visualizer options to view only the mesh
fv_options.plot_mask =                                                      %don't plot the mask
fv_options.plot_mesh =                                                      %do plot the mesh
fv_options.mesh_color =                                                     % set the mesh color to be a favorite color
fv_options.mesh_size = [anat_numrows anat_numcols];                         %rows x columns of the images used to generate the mesh
fv_options.mesh_dims = [anat_fov(1) anat_slcthick];                         %FOV and slice thickness of the images used to create the mesh, in mm

% visualize the aponeurosis mesh:


% now, visualize both the mask and the mesh simultaneously


%% diagonalize D, find E1 and FA in each voxel

e1fa = zeros(length(tensor_m(:,1,1,1,1)), length(tensor_m(1,:,1,1,1)), length(tensor_m(1,1,:,1,1)), 4);

for r=1:length(tensor_m(:,1,1,1,1))
    for c=1:length(tensor_m(1,:,1,1,1))
        for s=1:length(tensor_m(1,1,:,1,1))
            if mask(r,c,s)

                % get the diffusion tnesor from teh matrix provided
                D = squeeze(tensor_m(r,c,s,:,:));

                % find the eigenvalues and eigenvectors using either svd or eigs
                [E, L] = 

                % sort the eigenvalues in order of descending magnitude

                % find the mean diffusivity and fractional anisotropy
                MD = mean(diag(L));
                L1 = 
                L2 = 
                L3 = 
                FA = 

                % extract the first eigenvector from the E matrix
                E1 = 

                %correct for antipodal symmetry (whazzat?)
                if E1(3)<0
                    
                end

                % write to e1fa matrix
                e1fa(r,c,s,1:3) = E1;
                e1fa(r,c,s,4) = FA;

            end
        end
    end
end

% display slice 30 as an rgb image
figure


% smooth first eigenvector field
E1map_smooth = E1_smoothing(squeeze(e1fa(:,:,:,1:3)));

figure

% create a matrix for the smoothed E1's
e1fa_smooth = e1fa;
e1fa_smooth(:,:,:,1:3) = E1map_smooth;


%% Fiber track

% set fiber-tracking options
ft_options.ref_frame = 'LPS';                                               %LPS frame of reference
ft_options.image_orient = 'AL';                                             %anterior/left image orientation
ft_options.seed_method = 'apo';                                             %track from aponeurosis
ft_options.mesh_dist = 0;                                                   %no shift in mesh position
ft_options.prop_algo = 'rk4';                                               %integrate E1 using 4th order Runge Kutta
ft_options.step_size = 1;                                                   %step size of 1 voxel width
ft_options.term_mthd = 'bin2';                                              %BIN2 requires two points to meet stop criteria
ft_options.angle_thrsh = [30 2];                                            %30 degree angle between current step and the step 2 points earlier
ft_options.fa_thrsh = [.1 .4];                                              %FA limits
ft_options.depth_ratio = dti_slcthick/dti_pixelspace(1);                    %depth ratio of 7

% call the function:
[fiber_all, roi_flag, stop_list, fiber_len] = ...
     fiber_track(ft_options, e1fa, mask, roi_mesh, fv_options);

%update fiber_visualizer options
fv_options.plot_mask = 0;                                                   %don't plot the mask
fv_options.plot_mesh = 1;                                                   %do plot the mesh
fv_options.mesh_dist = ft_options.mesh_dist;                                %match to fiber tracking option/noshift
fv_options.mesh_color=[.7 .7 .7];                                           %gray
fv_options.plot_fibers = 1;                                                 %do plot the mesh
fv_options.fiber_color = [.8 .2 .2];                                        %fibers will be red
fv_options.fiber_width = 0.5;                                               %linewidth for fiber tracts
fv_options.dti_size = [dti_numrows dti_numcols];                            %matrix size of DT images
fv_options.dti_dims = [dti_fov(1) dti_slcthick];                            %FOV and slice thickness of DT images

%visualize fiber tracts
fiber_figure = fiber_visualizer(anat_image, fv_options, roi_mesh, mask, fiber_all);
xlabel('Left-Right Direction (mm)')
ylabel('Anterior-Posterior Direction (mm)')
zlabel('Head-Foot Direction (mm)')


% play with seeding options, tracking options, etc.

%% Smooth tracts

%set fiber_smoother options
fs_options.interpolation_step = 1;                                          %interpolate at 1 pixel widths (=1 mm)
fs_options.p_order = [3 3 3];                                               %fit row, column, and slice positions to 4th, 4th,and 3rd order polynomials
fs_options.tract_units = 'vx';                                              %units - 'vx' or 'mm'
fs_options.dwi_res = [dti_fov(1) dti_numrows dti_slcthick];                 %DTI FOV, matrix size, slice thickness
                
% call the function:
[smoothed_fiber_all, ~, smoothed_fiber_all_mm, pcoeff_r, pcoeff_c, pcoeff_s,...
    n_points_smoothed, residuals, residuals_mm] = fiber_smoother(fiber_all, fs_options);

%visualize fiber tracts


%% select 500 tracts

[sampled_fiber_all_mm, sampled_fiber_all_idx] = far_stream_sampling( , );

% convert to voxels for display
sampled_smoothed_fiber_all_vx = zeros(size(sampled_fiber_all_mm));
for k=1:length(sampled_fiber_all_idx)

    loop_sampled_fiber_all_vx = squeeze(sampled_fiber_all_mm(sampled_fiber_all_idx(k,1), sampled_fiber_all_idx(k,2), :, :));
    loop_sampled_fiber_all_vx(:,1) = loop_sampled_fiber_all_vx(:,1)/(dti_fov(1)/dti_numrows);
    loop_sampled_fiber_all_vx(:,2) = loop_sampled_fiber_all_vx(:,2)/(dti_fov(1)/dti_numrows);
    loop_sampled_fiber_all_vx(:,3) = loop_sampled_fiber_all_vx(:,3)/dti_slcthick;

    sampled_smoothed_fiber_all_vx(sampled_fiber_all_idx(k,1), sampled_fiber_all_idx(k,2), :, :) = loop_sampled_fiber_all_vx;

end

%visualize fiber tracts
sampled_smoothed_fiber_fig = fiber_visualizer(anat_image, fv_options, roi_mesh, [], sampled_smoothed_fiber_all_vx);


%% Quantify architectural properties

fq_options.dwi_res = [dti_fov(1) dti_numrows dti_slcthick];                 %DTI FOV, matrix size, slice thickness
fq_options.filt_kernel = 3;                                                 %size of smoothing kernel for determining aponeurosis normal vectors
fq_options.mesh_units = 'vx';                                               %tracts are in units of voxels
fq_options.tract_units = 'vx';

[angle_list, distance_list, curvature_list, fiber_all_mm, n_points, apo_area] = ...
     fiber_quantifier(fq_options, sampled_smoothed_fiber_all_vx, roi_mesh, mask);


% find average properties for each tract, display a heat map with colorscale



