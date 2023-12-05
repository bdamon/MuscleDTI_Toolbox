function track_simulate(muscle_structure, simulation_vars, isdisp, raw_dir)
%TRACK_SIMULATE Track through diffusion tensors and generate raw tracts
%   muscle_structure: Returned from FORM_MODEL_MUSCLE
%   simulation_vars: Structure to control simulation variables
%       -fit_order: Polynomial order for fitted fiber
%       -poly_params: A 3-by-N vector
%           poly_params(1, : ) - Polynomial coefficients for X positions
%           poly_params(2, : ) - Polynomial coefficients for Y positions
%           poly_params(3, : ) - Polynomial coefficients for Z positions
%           as functions of *NORMALIZED* tracking distance for true fiber
%       -curvature: Scalar, mean curvature along fiber
%       -full_len: Full length of the whole fiber, ground truth
%       -angle_list: Pennation angle along the fiber, ground truth
%       -concave_dir: Concavity direction; 'U' for up and 'D' for down
%       -voxel_dim: A 1-by-3 vector
%           voxel_dim(1) - In-plane resolution (mm)
%           voxel_dim(2) - In-plane resolution (mm)
%           voxel_dim(3) - Slice thickness (mm)
%           assume *SQUARE* in-plane pixels
%       -SNR: Scalar, noise level
%       -step_size: Scalar, tracking step size in *VOXEL WIDTH*
%       -prop_algo: Integration algorithm; 'rk4' or 'euler'
%       -term_mthd: Termination criteria; 'bin1' or 'bin2'
%
%   isdisp: Scalar, whether to display demo figures
%       0 - No
%       1 - Yes
%
%   raw_dir: Directory to read/write raw tractography data

%% 00 - Initialization
% poly_params = simulation_vars.poly_params;
voxel_dim   = simulation_vars.voxel_dim;
SNR         = simulation_vars.SNR;
step_size   = simulation_vars.step_size;
prop_algo   = simulation_vars.prop_algo;
term_mthd   = simulation_vars.term_mthd;

FOV_core_xy_vx = ceil(30 / voxel_dim(1));
FOV_core_xy_mm = FOV_core_xy_vx * voxel_dim(1);
FOV_full_xy_vx = FOV_core_xy_vx + 4;
FOV_full_xy_mm = FOV_full_xy_vx * voxel_dim(1);

fname_core = sprintf('CR%d%c%c_VX%d_%d_SN%d_SS%d.mat', ...
    round(simulation_vars.curvature * 100), ...
    simulation_vars.concave_dir, ...
    simulation_vars.dim_tag, ...
    round(voxel_dim(1) * 100), ...
    round(voxel_dim(3) * 100), ...
    round(SNR), ...
    round(step_size * 100));
fraw_save = sprintf('%s/%s', raw_dir, fname_core);

if isdisp == 1
    anat_image = squeeze(muscle_structure.diff_img( : , : , : , 1));
end

%% 01 - ROI Mesh Generation
% create roi_mesh matrix
mesh_H_vx = ceil(FOV_core_xy_mm / voxel_dim(3));
mesh_W_vx = FOV_core_xy_vx;

% Number of voxels to skip when placing seed points horizontally
if strcmp(simulation_vars.dim_tag, 'A')
    mesh_step_y_vx = 3;
    mesh_step_z_vx = 1;
elseif strcmp(simulation_vars.dim_tag, 'B')
    mesh_step_y_vx = 4;
    mesh_step_z_vx = 1;
end

seed_pt_zv = 1 : mesh_step_z_vx : mesh_H_vx;
seed_pt_yv = 1 : mesh_step_y_vx : mesh_W_vx;

roi_mesh = zeros(length(seed_pt_zv), length(seed_pt_yv), 6);
roi_mesh( : , : , 1) = 3;  %x component

[C, R] = meshgrid(1 : mesh_step_y_vx : mesh_W_vx, ...
                  1 : mesh_step_z_vx : mesh_H_vx);
roi_mesh( : , : , 2) = C + 2;
roi_mesh( : , : , 3) = R + 20 / voxel_dim(3);
roi_mesh( : , : , 5) = -1;  %normal to mesh

%% 01.1 - Visualization
if isdisp == 1
    fv_options.anat_dims = [FOV_full_xy_mm, voxel_dim(3)];
    fv_options.anat_slices  = round((5 : 20 : 120) ./ voxel_dim(3)) + 1;
    fv_options.plot_mask = 0;
    fv_options.plot_fibers = 0;
    fv_options.plot_mesh = 1;
    fv_options.mesh_size = round([FOV_full_xy_vx, FOV_full_xy_vx]);
    fv_options.mesh_dims = [FOV_full_xy_mm, voxel_dim(3)];
    fv_options.mesh_color = [0.5, 0.5, 0.5];
    mesh_figure = fiber_visualizer(anat_image, fv_options, roi_mesh, [], []);
    axis image
    zlim([0, 130]);
    set(gcf, 'position', [100 50 560 420])

    mask = muscle_structure.mask_muscle;

    fv_options.plot_mask = 1;
    fv_options.plot_fibers = 0;
    fv_options.plot_mesh = 1;
    fv_options.mask_size = round([FOV_full_xy_vx, FOV_full_xy_vx]);
    fv_options.mask_dims = [FOV_full_xy_mm, voxel_dim(3)];
    fv_options.mask_color = [1, 0.25, 0.25];
    mask_figure = fiber_visualizer(anat_image, fv_options, roi_mesh, mask, []);
    axis image
    zlim([0, 130]);
    set(gcf, 'position', [100 50 560 420])
end

%% 02 - Fiber Tracking
% Use published parameters
ft_options.ref_frame = 'LPS';                                               %LPS frame of reference
ft_options.image_orient = 'AL';                                             %anterior/left image orientation
ft_options.mesh_dist = 0;                                                   %no shift in mesh position
ft_options.prop_algo = prop_algo;                                           %integrate E1 using 4th order Runge Kutta
ft_options.step_size = step_size;                                           %step size of *VOXEL WIDTH*
ft_options.term_mthd = term_mthd;                                           %BIN2 requires two points to meet stop criteria
ft_options.angle_thrsh = [30 2];                                            %30 degree angle between current step and the step 2 points earlier
ft_options.fa_thrsh = [0.1, 0.4];                                           %FA limits
ft_options.depth_ratio = voxel_dim(3) / voxel_dim(1);                       %depth ratio of

% ft_options.prop_algo = 'fact';
% ft_options.r_crit = 0.95;
% ft_options.num_fact_voxels = 5;

if isdisp == 1
    fv_options.plot_mask = 0;                                               %don't plot the mask
    fv_options.plot_mesh = 1;                                               %do plot the mesh
    fv_options.mesh_dist = ft_options.mesh_dist;                            %match to fiber tracking option/noshift
    fv_options.mesh_color=[0.7, 0.7, 0.7];                                  %gray
    fv_options.plot_fibers = 1;                                             %do plot the mesh
    fv_options.fiber_color = [0.9, 0.2, 0.2];                               %fibers will be red
    fv_options.fiber_width = 0.5;                                           %linewidth for fiber tracts
    fv_options.dti_size = round([FOV_full_xy_vx, FOV_full_xy_vx]);          %matrix size of DT images
    fv_options.dti_dims = [FOV_full_xy_mm, voxel_dim(3)];                   %FOV and slice thickness of DT images

    % Customized
    % fv_options.fiber_skip = 2;

    [fiber_all, roi_flag, stop_list, fiber_len, fa_all, md_all] = ...
        fiber_track(...
            double(muscle_structure.tensor_m), ...
            muscle_structure.mask_muscle, ...
            roi_mesh, ft_options, fv_options, anat_image);
else
    [fiber_all, roi_flag, stop_list, fiber_len, fa_all, md_all] = ...
        fiber_track(...
            double(muscle_structure.tensor_m), ...
            muscle_structure.mask_muscle, ...
            roi_mesh, ft_options);
end



if isdisp == 1
    axis image
    zlim([0, 130]);
end

%% 06 - Save Data
fiber_all_mm = zeros(size(fiber_all));
fiber_all_mm( : , : , : , 1) = ...
    squeeze(fiber_all( : , : , : , 1)) * voxel_dim(1);
fiber_all_mm( : , : , : , 2) = ...
    squeeze(fiber_all( : , : , : , 2)) * voxel_dim(2);
fiber_all_mm( : , : , : , 3) = ...
    squeeze(fiber_all( : , : , : , 3)) * voxel_dim(3);

roi_mesh     = single(roi_mesh);
fiber_all    = single(fiber_all);
fiber_all_mm = single(fiber_all_mm);
roi_flag     = single(roi_flag);
stop_list    = single(stop_list);
fiber_len    = single(fiber_len);
save(fraw_save, ...
    'simulation_vars', 'roi_mesh', ...
    'fiber_all', 'fiber_all_mm', 'roi_flag', 'stop_list', 'fiber_len');
% fprintf('Raw fiber tracts saved to:\n');
% fprintf('  %s\n', fraw_save);
end