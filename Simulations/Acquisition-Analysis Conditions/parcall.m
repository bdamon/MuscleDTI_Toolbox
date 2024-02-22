% To call DTI_Simulation_Driver and run parallel, the Parallel Computing
% Toolbox (PCT) must be installed.

%% 00 - Initialization
clc;
clear all;
close all;

code_dir = 'my_code_directory';
cd(code_dir);

tagstr = 'Simulation_New';

%% 01 - Simulation
fit_order_v = [2, 3, 4];
k0 = 1;
Nsimu = 1; % 84
Nthread = 1; % 6
% Variable CURV_GROUPS is initialized in DTI_SIMULATION_DRIVER.M
% Variable OPT_SET is initialized in DTI_SIMULATION_DRIVER.M

subset_tag = sprintf('%03d_%03d', k0, k0 + Nsimu * Nthread - 1);

tic
for fit_order = fit_order_v
    for ks = 1 : Nsimu
        parfor kt = 1 : Nthread
            ksimu = k0 - 1 + (ks - 1) * Nthread + kt;
            simu_tag = sprintf('%s_%d', tagstr, ksimu);
            fprintf('Simulation %d started...\n', ksimu);
            DTI_Simulation_Driver(simu_tag, fit_order, subset_tag);
            fprintf('Simulation %d finished.\n\n', ksimu);
        end
    end
end
toc