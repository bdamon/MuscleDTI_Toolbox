%% process_tracts.m
% This version of process_tracts produces the data published in the manuscript.
% Completed on 17 Jan 2021, Bruce Damon


%% File Input

%start with a clean slate:
clear
close all
clc

% prompt user to select pre-processed data file/load
[pre_fiber_tracking_file, pre_fiber_tracking_path] = uigetfile('*.mat');
cd(pre_fiber_tracking_path)
load(pre_fiber_tracking_file)

%% Form muscle mask

% decide which slices to analyze when forming the mask
s = 1;
while 1
    imagesc(anat_image(:,:,s))
    colormap gray
    axis image
    title(['Slice ' num2str(s) '; Press left mouse button to advance, right button to go back, or middle button to quit'])
    [~, ~, b] = ginput(1);
    if b==1
        s = min([anat_numslcs s+1]);
    elseif b==3
        s = max([1 s-1]);
    else
        break
    end
end
clc
slices(1) = input('What is the first slice that you would like to analyze? ');
slices(2) = input('What is the last slice that you would like to analyze? ');

% set visualization options
fv_options.anat_dims = [anat_fov(1) anat_slcthick];                         %FOV and slice thickness of the images to be displayed, in mm
fv_options.anat_slices = 5:9:44;                                            %display slices 14, 24, 34, and 44
fv_options.plot_mesh = 0;                                                   %don’t plot an aponeurosis mesh
fv_options.plot_mask = 1;                                                   %do plot the mask
fv_options.plot_fibers = 0;                                                 %don’t plot any fiber tracts
fv_options.mask_size = [anat_numrows anat_numcols];                         %rows x columns of the images used to generate the mask
fv_options.mask_dims = [anat_fov(1) anat_slcthick];                         %FOV and slice thickness of the images used to create the mask, in mm
fv_options.mask_color = [1 0 0];                                            %make the mask a red, semi-transparent overlay

% call the define_muscle function:
close all
[mask, ~] = define_muscle(anat_image, slices, [], fv_options);
set(gca, 'color', 'k')

% save file with unique name
save ta_mask mask


%% Form mesh reconstruction of aponeurosis

% decide which slices to analyze when forming the mask and which slice has the largest aponeurosis
close all
s = 1;
while 1
    imagesc(anat_image(:,:,s))
    colormap gray
    axis image
    title(['Slice ' num2str(s) '; Press left mouse button to advance, right button to go back, or middle button to quit'])
    [~, ~, b] = ginput(1);
    if b==1
        s = min([anat_numslcs s+1]);
    elseif b==3
        s = max([1 s-1]);
    else
        break
    end
end
clc
dr_options.slices(1) = input('First slice to analyze: ');                   %begin forming define_roi options structure
dr_options.slices(2) = input('Last slice to analyze: ');

% figure out the mesh dimensions
apo_slice = input('In which slice is the aponeurosis the biggest? ');
close all
imagesc(anat_image(:,:,apo_slice));
colormap gray
axis image
title('Digitize the aponeurosis using left mouse clicks; then press return')
[x_points, y_points] = ginput;                                              %select number of points needed to define the aponeurosis, then hit enter
dx = diff(x_points)*dti_pixelspace(1);                                      %calculate inter-point differences in X/multiply by pixel spacing
dy = diff(y_points)*dti_pixelspace(1);                                      %calculate inter-point differences in Y/multiply by pixel spacing
max_apo_length = sum((dx.^2 + dy.^2).^(1/2));                               %length of aponeurosis in mm

% update fiber visualize options for the mesh
fv_options.plot_mask = 0;                                                   %don't plot the mask
fv_options.plot_mesh = 1;                                                   %do plot the mesh
fv_options.mesh_color = [.5 .5 .5];                                         %medium gray
fv_options.mesh_size = [anat_numrows anat_numcols];                         %rows x columns of the images used to generate the mesh
fv_options.mesh_dims = [anat_fov(1) anat_slcthick];                         %FOV and slice thickness of the images used to create the mesh, in mm

% finish setting define_roi options
dr_options.dti_size = [dti_numrows dti_numcols dti_all_numslcs];            %rows x columns x slices of the images used to generate the mesh
dr_options.mesh_size = [ceil(dti_slcthick*diff(dr_options.slices)) 2*round(max_apo_length)];  %about 2/mm^2 at widest point
dr_options.method = 'auto';                                                 %use automated segmentation
dr_options.dilate_mesh = 'y';                                               %dilate the mesh within the slice plane
dr_options.n_steps = 1;                                                     %by one pixel

% call the define_roi function:
[roi_mesh, roi_mask, roi_mesh_dilated_1] = define_roi(anat_image, mask, dr_options, fv_options);

% save mesh file
save ta_mesh roi_mesh* roi_mask

% save pre-tracking data

save pretracking_data mask roi_* anat_image fv_options dr_options tensor_m anat* dti_* max_apo_length

%free up some memory
clear
load pretracking_data

%% Fiber track

% set fiber-tracking options
ft_options.ref_frame = 'LPS';                                               %LPS frame of reference
ft_options.image_orient = 'AL';                                             %anterior/left image orientation
ft_options.mesh_dist = 0;                                                   %no shift in mesh position
ft_options.prop_algo = 'rk4';                                               %integrate E1 using 4th order Runge Kutta
ft_options.step_size = 1;                                                   %step size of 1 voxel width
ft_options.term_mthd = 'bin2';                                              %BIN2 requires two points to meet stop criteria
ft_options.angle_thrsh = [30 2];                                            %30 degree angle between current step and the step 2 points earlier
ft_options.fa_thrsh = [.1 .4];                                              %FA limits
ft_options.depth_ratio = dti_slcthick/dti_pixelspace(1);                    %depth ratio of 7

%update fiber_visualizer options
fv_options.plot_mask = 0;                                                   %don't plot the mask
fv_options.plot_mesh = 1;                                                   %do plot the mesh
fv_options.mesh_dist = ft_options.mesh_dist;                                %match to fiber tracking option/noshift
fv_options.mesh_color=[.7 .7 .7];                                           %gray
fv_options.plot_fibers = 1;                                                 %do plot the mesh
fv_options.fiber_color = [.9 .2 .2];                                        %fibers will be red
fv_options.fiber_width = 0.5;                                               %linewidth for fiber tracts
fv_options.dti_size = [dti_numrows dti_numcols];                            %matrix size of DT images
fv_options.dti_dims = [dti_fov(1) dti_slcthick];                            %FOV and slice thickness of DT images

% call the function:
[fiber_all, roi_flag, stop_list, fiber_len, fa_all, md_all] = ...
    fiber_track(tensor_m, mask-roi_mask, roi_mesh_dilated_1, ft_options, fv_options, anat_image);


%% Smooth tracts

%set fiber_smoother options
fs_options.interpolation_step = 1;                                          %interpolate at 1 pixel widths (=1 mm)
fs_options.p_order = [4 4 2];                                               %fit row, column, and slice positions to 4th, 4th,and 3rd order polynomials

% call the function:
[smoothed_fiber_all, pcoeff_r, pcoeff_c, pcoeff_s, n_points_smoothed] = ...
    fiber_smoother(fiber_all, fs_options);

% pick some random fiber tracts to plot
tracts_row = randperm(length(roi_mesh_dilated_1(:,1,1)));                   %choose a random set of tracts
tracts_row = tracts_row(1:72);                                              %72 of them - row direction
tracts_col = [randperm(18) randperm(18) randperm(18) randperm(18)];         %column direction

close all
for k=1:72
    figure(k)
    
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
beep                                                                        %get the user's attention
pause                                                                       %pause to inspect results

close all

%view all tracts
fv_options.plot_fibers = 1;                                                 %do plot the mesh
fv_options.mesh_color = [.7 .7 .7];

smoothed_fiber_fig = fiber_visualizer(anat_image, fv_options, roi_mesh_dilated_1, [], smoothed_fiber_all);

%% Quantify architectural properties

fq_options.dwi_res = [dti_fov(1) dti_numrows dti_slcthick];                 %DTI FOV, matrix size, slice thickness
fq_options.filt_kernel = 3;                                                 %size of smoothing kernel for determining aponeurosis normal vectors
fq_options.mesh_units = 'vx';                                               %tracts are in units of voxels
fq_options.tract_units = 'vx';

[angle_list, distance_list, curvature_list, fiber_all_mm, n_points, apo_area] = ...
    fiber_quantifier(fiber_all, roi_mesh_dilated_1, fq_options);

%% Evaluate outcomes and select tracts

% Set fiber_goodness options:
fg_options.dwi_res = [dti_fov(1) dti_numrows dti_slcthick];
fg_options.min_distance  =  10;                                             %tracts must be >10 mm long
fg_options.min_pennation = 0;                                               %acceptable pennation angles are from 0-40 degrees
fg_options.max_pennation = 40;
fg_options.max_curvature = 40;                                              %acceptable curvatures are from 0-40 m^1

[final_fibers, final_curvature, final_angle, final_distance, qual_mask, num_tracked, mean_fiber_props, mean_apo_props] = fiber_goodness ...
    (smoothed_fiber_all, angle_list, distance_list, curvature_list, n_points, roi_flag, apo_area, roi_mesh_dilated_1, fg_options);

% view architectural properties of tracts - example is curvature
fv_options.fiber_color = zeros(size(squeeze(roi_mesh_dilated_1(:,:,1:3))));
fv_options.fiber_color(:,:,2) = mean_fiber_props(:,:,1)/40;                 %increasing greenness is increasing curvature, 0 to 40 m^-1
fv_options.fiber_color(:,:,3) = 1-mean_fiber_props(:,:,1)/40;               %increasing blueness is decreasing curvature, 40 to 0 m^-1

final_fiber_fig = fiber_visualizer(anat_image, fv_options, roi_mesh_dilated_1, [], final_fibers);


%% save output

% save full dataset
close all
save all_tracking_data

% to reduce file size for Github upload, clear out variables stored elsewhere and save as new file
clear tensor_m mask roi_mesh* roi_mask
save final_tracking_data




