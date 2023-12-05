function DTI_Simulation_Driver(simu_tag, fit_order, subset_tag)
% DTI_SIMULATION_DRIVER DTI Tractography Simulation Driver
%   To set up parameter values and invoke the core simulation module

%% 00 - Initialization
clc;
%clear all;
close all;
rng('shuffle');

if nargin  == 2
    subset_tag = [];
end

code_dir = 'my_code_directory';
cd(code_dir);

%% 01 - Parameter Assignment
isdisp = 0;
curv_groups = 1 : 4;
opt_set = 1;
param_init();

simu_rw_dir = sprintf('%s/%s', simu_base_dir, simu_tag);
if ~exist(simu_rw_dir, 'dir')
    mkdir(simu_rw_dir);
end

% Simulated fiber geometry and DWI images
simu_dwi_dir = sprintf('%s/%s', simu_rw_dir, 'dwi_simulate');
if ~exist(simu_dwi_dir, 'dir')
    mkdir(simu_dwi_dir);
end

% Tractography results
track_base_dir = sprintf('%s/%s_%s', ...
    simu_rw_dir, prop_algo, term_mthd);
if ~exist(track_base_dir, 'dir')
    mkdir(track_base_dir);
end

track_raw_dir = sprintf('%s/%s', track_base_dir, 'track_raw');
if ~exist(track_raw_dir, 'dir')
    mkdir(track_raw_dir);
end

track_fit_dir = sprintf('%s/fit_p%d', track_base_dir, fit_order);
if ~exist(track_fit_dir, 'dir')
    mkdir(track_fit_dir);
end

simulation_vars.fit_order = fit_order;

%% 02 - Loop through Parameter Levels
for kcurv = 1 : Ncurv
    if mod(kcurv, 4) == 1
        fprintf('[%s] Simulation in progress: %d/%d\n', ...
            simu_tag, kcurv - 1, Ncurv);
    end
    for kvs = 1 : Nvs
        for ksn = 1 : Nsn
            simulation_vars.poly_params = poly_params( : , : , kcurv);
            simulation_vars.curvature   = mean_curv(kcurv);
            simulation_vars.concave_dir = concave_dir(kcurv);
            simulation_vars.dim_tag     = dim_tag(kcurv);
            simulation_vars.voxel_dim   = voxel_dims(kvs, : );
            simulation_vars.SNR         = SN_list(ksn);

            fname_core = sprintf('CR%d%c%c_VX%d_%d_SN%d', ...
                round(simulation_vars.curvature * 100), ...
                simulation_vars.concave_dir, ...
                simulation_vars.dim_tag, ...
                round(simulation_vars.voxel_dim(1) * 100), ...
                round(simulation_vars.voxel_dim(3) * 100), ...
                round(simulation_vars.SNR));
            fdwi_save = sprintf('%s/%s.mat', simu_dwi_dir, fname_core);

%             fprintf('\n====================================\n');
%             fprintf('Curvature: %-6.3f m-1    ', ...
%                 simulation_vars.curvature);
%             fprintf('In-plane resolution: %.2f mm    ', ...
%                 simulation_vars.voxel_dim(1));
%             fprintf('Slice thickness: %.2f mm    ', ...
%                 simulation_vars.voxel_dim(3));
%             fprintf('SNR: %d\n\n', ...
%                 simulation_vars.SNR);

            if ~exist(fdwi_save, 'file')
                %% Generate model muscle
                muscle_structure = ...
                    form_model_muscle(dti_options, simulation_vars, ...
                    isdisp);

                if isdisp == 1
                    fprintf(['Press any key to ' ...
                        'close all figure windows and continue.\n']);
                    pause();
                    close all;
                end

                save(fdwi_save, 'muscle_structure');
%                 fprintf('Simulated dataset saved to:\n');
%                 fprintf('  %s\n', fdwi_save);
            else
                load(fdwi_save, 'muscle_structure');
%                 fprintf('Simulated dataset loaded from:\n');
%                 fprintf('  %s\n', fdwi_save);
            end

            % Update additional structural infomation
            simulation_vars.full_len = ...
                get_full_len(muscle_structure.model_fiber);
            simulation_vars.angle_list = ...
                get_angle(muscle_structure.model_fiber, [0, -1, 0]);

            for kss = 1 : Nss
                simulation_vars.step_size = step_list(kss);

                fraw_save = sprintf('%s/%s_SS%d.mat', ...
                    track_raw_dir, fname_core, ...
                    round(simulation_vars.step_size * 100));
                if ~exist(fraw_save, 'file')
                    %% Compute raw fiber tracts
%                     fprintf('\nRunning tractography simulation...\n');
%                     fprintf('Tracking step size: %.2f voxel widith ', ...
%                         simulation_vars.step_size);
%                     fprintf('(%.2f mm)\n', ...
%                         simulation_vars.step_size * ...
%                         simulation_vars.voxel_dim(1));
%                     tic
                    track_simulate(muscle_structure, simulation_vars, ...
                        isdisp, track_raw_dir);
%                     toc
                    if isdisp == 1
                        fprintf(['Press any key to ' ...
                            'close all figure windows and continue.\n']);
                        pause();
                        close all;
                    end
                end

                %% Post-processing (smoothing, goodness filtering)
%                 fprintf('\nRunning tractography post-processing...\n');
%                 fprintf('Tracking step size: %.2f voxel widith ', ...
%                     simulation_vars.step_size);
%                 fprintf('(%.2f mm)    ', ...
%                     simulation_vars.step_size * ...
%                     simulation_vars.voxel_dim(1));
%                 fprintf('Polynomial fitting order: %d\n\n', ...
%                     simulation_vars.fit_order);

                track_fit(muscle_structure, simulation_vars, ...
                    isdisp, track_base_dir);
                if isdisp == 1
                    fprintf(['Press any key to ' ...
                        'close all figure windows and continue.\n']);
                    pause();
                    close all;
                end
            end
        end
    end
end
fprintf('[%s] Simulation finished.\n', simu_tag);

%% 03 - Remove the unnecessary FIBER_TRACK_DATA.MAT
extra_track_file = sprintf('%s/%s', code_dir, 'fiber_tract_data.mat');
if exist(extra_track_file, 'file')
    delete(extra_track_file);
end
end
