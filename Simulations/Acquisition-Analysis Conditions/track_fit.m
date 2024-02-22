function track_fit(muscle_structure, simulation_vars, isdisp, track_base_dir)
%TRACK_FIT Smooth fiber tracts and perform fiber goodness check
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
%   track_base_dir: Directory to read/write data files

%% 00 - Initialization
% poly_params = simulation_vars.poly_params;
voxel_dim   = simulation_vars.voxel_dim;
SNR         = simulation_vars.SNR;
step_size   = simulation_vars.step_size;

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
fraw_read = sprintf('%s/track_raw/%s', track_base_dir, fname_core);
ffit_save = sprintf('%s/fit_p%d/%s', ...
    track_base_dir, simulation_vars.fit_order, fname_core);
if exist(ffit_save, 'file')
    fprintf('Tract fitting results existed:\n  %s\n', ffit_save);
    fprintf('To re-run tract fitting, please delete the data file.\n');
    return;
end

load(fraw_read, ...
    'roi_mesh', 'fiber_all', 'fiber_all_mm', ...
    'roi_flag', 'stop_list', 'fiber_len');
% fprintf('Raw fiber tracts loaded from:\n');
% fprintf('  %s\n', fraw_read);

if isdisp == 1
    anat_image = squeeze(muscle_structure.diff_img( : , : , : , 1));

    if strcmp(simulation_vars.dim_tag, 'A')
        mesh_step_y_vx = 3;
        mesh_step_z_vx = 1;
    elseif strcmp(simulation_vars.dim_tag, 'B')
        mesh_step_y_vx = 4;
        mesh_step_z_vx = 1;
    end

    mesh_H_vx = ceil(FOV_core_xy_mm / voxel_dim(3));
    mesh_W_vx = FOV_core_xy_vx;

    seed_pt_zv = 1 : mesh_step_z_vx : mesh_H_vx;
    seed_pt_yv = 1 : mesh_step_y_vx : mesh_W_vx;
end

%% 03 - Fiber Tract Smoothing
%set fiber_smoother options
fs_options.interpolation_step = 1;                                          %interpolate at 1 pixel widths (=1 mm)
fs_options.p_order = simulation_vars.fit_order;                             %fit row, column, and slice positions to 4th, 4th,and 3rd order polynomials
fs_options.dwi_res = [FOV_full_xy_mm, FOV_full_xy_vx, voxel_dim(3)];        %DTI FOV, matrix size, slice thickness
fs_options.tract_units = 'vx';

% call the function:
[smoothed_fiber_all, ~ , smoothed_fiber_all_mm, pcoeff_r, pcoeff_c, pcoeff_s, n_points_smoothed, residuals, residuals_mm] = ...
    fiber_smoother(double(fiber_all), fs_options);

if isdisp == 1
    % pick some random fiber tracts to plot
    tracts_row = randperm(length(seed_pt_zv));                   %choose a random set of tracts
    tracts_row = tracts_row(1:length(seed_pt_zv));
    tracts_col = [randperm(length(seed_pt_yv)) randperm(length(seed_pt_yv)) randperm(length(seed_pt_yv)) randperm(length(seed_pt_yv))];         %column direction

%     close all
    for k=1:length(seed_pt_zv)
        figure(k * 100)

        if n_points_smoothed(tracts_row(k), tracts_col(k))>5
            subplot(1,3,1)
            plot(((1:n_points_smoothed(tracts_row(k), tracts_col(k)))-1)/(n_points_smoothed(tracts_row(k), tracts_col(k))-1), nonzeros(squeeze(smoothed_fiber_all(tracts_row(k), tracts_col(k), :, 1))), 'b');
            hold on
            plot(((1:fiber_len(tracts_row(k), tracts_col(k)))-0)/(fiber_len(tracts_row(k), tracts_col(k))-1), nonzeros(squeeze(fiber_all(tracts_row(k), tracts_col(k), :, 1))), 'b.');

            subplot(1,3,2)
            plot(((1:n_points_smoothed(tracts_row(k), tracts_col(k)))-1)/(n_points_smoothed(tracts_row(k), tracts_col(k))-1), nonzeros(squeeze(smoothed_fiber_all(tracts_row(k), tracts_col(k), :, 2))), 'b');
            hold on
            plot(((1:fiber_len(tracts_row(k), tracts_col(k)))-0)/(fiber_len(tracts_row(k), tracts_col(k))-1), nonzeros(squeeze(fiber_all(tracts_row(k), tracts_col(k), :, 2))), 'b.');

            subplot(1,3,3)
            plot(((1:n_points_smoothed(tracts_row(k), tracts_col(k)))-1)/(n_points_smoothed(tracts_row(k), tracts_col(k))-1), nonzeros(squeeze(smoothed_fiber_all(tracts_row(k), tracts_col(k), :, 3))), 'b');
            hold on
            plot(((1:fiber_len(tracts_row(k), tracts_col(k)))-0)/(fiber_len(tracts_row(k), tracts_col(k))-1), nonzeros(squeeze(fiber_all(tracts_row(k), tracts_col(k), :, 3))), 'b.');
        end
    end
%     beep                                                                        %get the user's attention
%     pause                                                                       %pause to inspect results

%     close all

    %view all tracts
    fv_options.plot_fibers = 1;                                                 %do plot the mesh
    fv_options.mesh_color = [.7 .7 .7];

    smoothed_fiber_fig = fiber_visualizer(anat_image, fv_options, roi_mesh, [], smoothed_fiber_all);
    axis image
    zlim([0, 130]);
end

%% 04 - Architectural Property Quantification
% roi_mesh1 = roi_mesh + randn(size(roi_mesh)) .* 1e-4;
roi_mesh1 = roi_mesh;
roi_mesh1( : , : , 1) = roi_mesh1( : , : , 1) + randn(size(roi_mesh( : , : , 1))) .* 1e-4;
fq_options.dwi_res = [FOV_full_xy_mm, FOV_full_xy_vx, voxel_dim(3)];       %DTI FOV, matrix size, slice thickness
fq_options.filt_kernel = 3;                                                 %size of smoothing kernel for determining aponeurosis normal vectors
fq_options.mesh_units = 'vx';                                               %tracts are in units of voxels
fq_options.tract_units = 'vx';

% BE VERY CAUTIOUS
% [angle_list, distance_list, curvature_list, fiber_all_mm, n_points, apo_area] = ...
%     fiber_quantifier(fiber_all, roi_mesh1, fq_options);
[angle_list, distance_list, curvature_list, ~, n_points, apo_area] = ...
    fiber_quantifier(smoothed_fiber_all, roi_mesh1, fq_options);

%% 05 - Outcome Evaluation
% Set fiber_goodness options:
fg_options.dwi_res = [FOV_full_xy_mm, FOV_full_xy_vx, voxel_dim(3)];
fg_options.min_distance  = simulation_vars.full_len * 0.90;                 %tracts must be >90% full length
fg_options.min_pennation = 3;                                               %acceptable pennation angles are from 3 - max angle
fg_options.max_pennation = max(simulation_vars.angle_list);
fg_options.max_curvature = 40;                                              %acceptable curvatures are from 0-40 m^1

[final_fibers, final_curvature, final_angle, final_distance, qual_mask, num_tracked, mean_fiber_props, mean_apo_props] = ...
    fiber_goodness(smoothed_fiber_all, angle_list, distance_list, curvature_list, n_points, roi_flag, apo_area, fg_options);

if isdisp == 1
    % view architectural properties of tracts - example is curvature
    fv_options.fiber_color = zeros(size(squeeze(roi_mesh(:,:,1:3))));
    fv_options.fiber_color(:,:,2) = mean_fiber_props(:,:,1)/40;                 %increasing greenness is increasing curvature, 0 to 40 m^-1
    fv_options.fiber_color(:,:,3) = 1-mean_fiber_props(:,:,1)/40;               %increasing blueness is decreasing curvature, 40 to 0 m^-1

    final_fiber_fig = fiber_visualizer(anat_image, fv_options, roi_mesh1, [], final_fibers);
    axis image
    zlim([0, 130]);
end

%% 06 - Save Data
smoothed_fiber_all_mm = single(smoothed_fiber_all_mm);
residuals_mm          = single(residuals_mm);
final_fibers          = single(final_fibers);
final_curvature       = single(final_curvature);
final_angle           = single(final_angle);
final_distance        = single(final_distance);
mean_fiber_props      = single(mean_fiber_props);
mean_apo_props        = single(mean_apo_props);

save(ffit_save, ...
    'smoothed_fiber_all_mm', 'residuals_mm', ...
    'final_fibers', 'final_curvature', 'final_angle', 'final_distance', ...
        'qual_mask', 'num_tracked', 'mean_fiber_props', 'mean_apo_props');
% fprintf('Smoothed fiber tracts saved to:\n');
% fprintf('  %s\n', ffit_save);
end