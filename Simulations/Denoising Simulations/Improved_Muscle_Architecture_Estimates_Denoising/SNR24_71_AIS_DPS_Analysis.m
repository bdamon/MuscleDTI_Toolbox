%% This code contains the analyses made on the simulated dataset presented
% in the manuscript 'Improved DTI-based skeletal muscle architecture
% estimation via diffusion weighted image denoising and first eigenvector
% map smoothing'

% Written by: Roberto Pineda Guzman, Carle Foundation Hospital

% Supporting functions:
% smoothn: https://www.mathworks.com/matlabcentral/fileexchange/25634-smoothn/

clear
clc
close all
 
% Load simulated data of model muscle

load NoiseFreeImages_RestrictedMesh
load NoiseFree_Euler_BIN2_30_2_AllApoMethods apo_ft_options all_fs_options...
     apo_fq_options composite_mask muscle_mask muscle_mask_filled apo_mask ...
     roi_mesh_restricted roi_mesh_restricted_mm
close all

%% generate and analyze noisy data

% define number of Monte Carle trials
num_trials = 1000;

%create zeros matrices
all_snr24 = zeros(num_trials, 1);
all_snr71 = zeros(num_trials, 1);


%% fiber tract computation on noise-free images

%%%%%%%%%%%%%%%% fiber tracking - aponeurosis seeding %%%%%%%%%%%%%%%%%%
[apo_fiber_all, ~, apo_stop_list] = ...
    fiber_track_v20(apo_ft_options, e1fa, composite_mask, roi_mesh_restricted);
[apo_smoothed_fiber_all, ~, apo_smoothed_fiber_all_mm, ~, ~, ~, ~, ~, ~] = ...
    fiber_smoother_v14a(apo_fiber_all, all_fs_options);

% prep for filtering based on same seed location
apo_seeds_mm_rounded = round(roi_mesh_restricted_mm(:,:,1:3));
apo_smoothed_fiber_all_mm_final = zeros(size(apo_smoothed_fiber_all_mm));   %will hold tracts after seed-location filtering
apo_smoothed_fiber_all_final = zeros(size(apo_smoothed_fiber_all));
apo_seeds_final = zeros(10,10,5);

n=1;
for r=min(min(apo_seeds_mm_rounded(:,:,1))):1:max(max(apo_seeds_mm_rounded(:,:,1)))
    for c=min(min(apo_seeds_mm_rounded(:,:, 2))):max(max(apo_seeds_mm_rounded(:,:, 2)))
        for s=min(min(apo_seeds_mm_rounded(:,:, 3))):max(max(apo_seeds_mm_rounded(:,:, 3)))
            idx = find(apo_seeds_mm_rounded(:,:,1)==r...
                & apo_seeds_mm_rounded(:,:, 2)==c...
                & apo_seeds_mm_rounded(:,:, 3)==s);

            if isempty(idx)
                continue

            else

                [loop_r, loop_c] = ind2sub(size(squeeze(roi_mesh_restricted(:,:,1))), idx);
                randm_v = randperm(length(loop_r));
                loop_r = loop_r(randm_v);
                loop_c = loop_c(randm_v);

                for k=1:length(loop_r)
                    if apo_stop_list(loop_r(k), loop_c(k))==4 && ...
                            length(nonzeros(apo_smoothed_fiber_all_mm(loop_r(k),loop_c(k),:,1)))>5
                        apo_seeds_final(loop_r(k),loop_c(k),:) = [r c s loop_r(k) loop_c(k)];
                        apo_smoothed_fiber_all_mm_final(loop_r(k),loop_c(k),:,1:3) = ...
                            apo_smoothed_fiber_all_mm(loop_r(k),loop_c(k),:,:);
                        apo_smoothed_fiber_all_final(loop_r(k),loop_c(k),1:length(apo_smoothed_fiber_all(loop_r(k),loop_c(k),:,1)),1:3) = ...
                            apo_smoothed_fiber_all(loop_r(k),loop_c(k),:,:);
                        break;
                    end
                end

                n=n+1;

            end  %of if isempty(idx)
        end %of s loop
    end %of c loop
end % of r loop

[apo_angle_list_final, apo_distance_list_final, apo_curvature_list_final, ~, apo_n_points_final, apo_area_final] = ...
    fiber_quantifier_v20a(apo_fq_options, apo_smoothed_fiber_all_final, roi_mesh_restricted, muscle_mask);

apo_seeds_final_mask = apo_seeds_final(:,:,4); % Create mask of apo_seeds_final to conserve same seed points on other datasets
apo_seeds_final_mask (apo_seeds_final_mask > 0) = 1;

% calculate raw summary statistics
apo_total_distance_final = squeeze(max(apo_distance_list_final, [], 3));
apo_total_distance_final(apo_total_distance_final<5)=0;

apo_mean_alpha_final = sum(squeeze(apo_angle_list_final(:,:,:,1)), 3)./squeeze(apo_n_points_final(:,:,2));
apo_mean_alpha_final(apo_total_distance_final<5)=0;

apo_mean_curvature_final = sum(apo_curvature_list_final, 3)./squeeze(apo_n_points_final(:,:,3));
apo_mean_curvature_final(apo_total_distance_final<5)=0;

%get rid of zeros and NaNs
apo_mean_alpha_final(isnan(apo_mean_alpha_final)) = 0;
apo_mean_alpha_final = nonzeros(apo_mean_alpha_final);
apo_mean_curvature_final(isnan(apo_mean_curvature_final)) = 0;
apo_mean_curvature_final = nonzeros(apo_mean_curvature_final);
apo_total_distance_final(isnan(apo_total_distance_final)) = 0;
apo_total_distance_final = nonzeros(apo_total_distance_final);


% save summary data for each time through the loop
all_mean_alpha_deg_noise_free = mean(nonzeros(apo_mean_alpha_final));
all_mean_curvature_m1_noise_free = mean(nonzeros(apo_mean_curvature_final));
all_total_distance_mm_noise_free = mean(nonzeros(apo_total_distance_final));

apo_smoothed_fiber_all_mm_noise_free=apo_smoothed_fiber_all_mm_final;

%% Monte Carlo simulation

angle_noise_24 = zeros([size(noise_free_img,1:3) num_trials]);
angle_noise_aniso_24 = zeros([size(noise_free_img,1:3) num_trials]);
angle_noise_71 = zeros([size(noise_free_img,1:3) num_trials]);
angle_noise_aniso_71 = zeros([size(noise_free_img,1:3) num_trials]);
angle_noise_24_smooth = zeros([size(noise_free_img,1:3) num_trials]);
angle_noise_aniso_24_smooth = zeros([size(noise_free_img,1:3) num_trials]);
angle_noise_71_smooth = zeros([size(noise_free_img,1:3) num_trials]);
angle_noise_aniso_71_smooth = zeros([size(noise_free_img,1:3) num_trials]);
fiber_similarity_24 = zeros([size(apo_smoothed_fiber_all_mm_noise_free,1:2) num_trials]);
fiber_similarity_aniso_24 = zeros([size(apo_smoothed_fiber_all_mm_noise_free,1:2) num_trials]);
fiber_similarity_71 = zeros([size(apo_smoothed_fiber_all_mm_noise_free,1:2) num_trials]);
fiber_similarity_aniso_71 = zeros([size(apo_smoothed_fiber_all_mm_noise_free,1:2) num_trials]);
fiber_similarity_24_smooth = zeros([size(apo_smoothed_fiber_all_mm_noise_free,1:2) num_trials]);
fiber_similarity_aniso_24_smooth = zeros([size(apo_smoothed_fiber_all_mm_noise_free,1:2) num_trials]);
fiber_similarity_71_smooth = zeros([size(apo_smoothed_fiber_all_mm_noise_free,1:2) num_trials]);
fiber_similarity_aniso_71_smooth = zeros([size(apo_smoothed_fiber_all_mm_noise_free,1:2) num_trials]);
euclid_distance_24 = zeros([size(apo_smoothed_fiber_all_mm_noise_free,1:2) num_trials]);
euclid_distance_aniso_24 = zeros([size(apo_smoothed_fiber_all_mm_noise_free,1:2) num_trials]);
euclid_distance_71 = zeros([size(apo_smoothed_fiber_all_mm_noise_free,1:2) num_trials]);
euclid_distance_aniso_71 = zeros([size(apo_smoothed_fiber_all_mm_noise_free,1:2) num_trials]);
euclid_distance_24_smooth = zeros([size(apo_smoothed_fiber_all_mm_noise_free,1:2) num_trials]);
euclid_distance_aniso_24_smooth = zeros([size(apo_smoothed_fiber_all_mm_noise_free,1:2) num_trials]);
euclid_distance_71_smooth = zeros([size(apo_smoothed_fiber_all_mm_noise_free,1:2) num_trials]);
euclid_distance_aniso_71_smooth = zeros([size(apo_smoothed_fiber_all_mm_noise_free,1:2) num_trials]);
cor_segment_24 = zeros([size(apo_smoothed_fiber_all_mm_noise_free,1:2) num_trials]);
cor_segment_aniso_24 = zeros([size(apo_smoothed_fiber_all_mm_noise_free,1:2) num_trials]);
cor_segment_71 = zeros([size(apo_smoothed_fiber_all_mm_noise_free,1:2) num_trials]);
cor_segment_aniso_71 = zeros([size(apo_smoothed_fiber_all_mm_noise_free,1:2) num_trials]);
cor_segment_24_smooth = zeros([size(apo_smoothed_fiber_all_mm_noise_free,1:2) num_trials]);
cor_segment_aniso_24_smooth = zeros([size(apo_smoothed_fiber_all_mm_noise_free,1:2) num_trials]);
cor_segment_71_smooth = zeros([size(apo_smoothed_fiber_all_mm_noise_free,1:2) num_trials]);
cor_segment_aniso_71_smooth = zeros([size(apo_smoothed_fiber_all_mm_noise_free,1:2) num_trials]);
all_mean_alpha_deg_noise_24 = zeros([num_trials 1]);
all_mean_curvature_m1_noise_24 = zeros([num_trials 1]);
all_total_distance_mm_noise_24 = zeros([num_trials 1]);
all_mean_alpha_deg_noise_71 = zeros([num_trials 1]);
all_mean_curvature_m1_noise_71 = zeros([num_trials 1]);
all_total_distance_mm_noise_71 = zeros([num_trials 1]);
all_mean_alpha_deg_noise_aniso_24 = zeros([num_trials 1]);
all_mean_curvature_m1_noise_aniso_24 = zeros([num_trials 1]);
all_total_distance_mm_noise_aniso_24 = zeros([num_trials 1]);
all_mean_alpha_deg_noise_aniso_71 = zeros([num_trials 1]);
all_mean_curvature_m1_noise_aniso_71 = zeros([num_trials 1]);
all_total_distance_mm_noise_aniso_71 = zeros([num_trials 1]);
all_mean_alpha_deg_noise_24_smooth = zeros([num_trials 1]);
all_mean_curvature_m1_noise_24_smooth = zeros([num_trials 1]);
all_total_distance_mm_noise_24_smooth = zeros([num_trials 1]);
all_mean_alpha_deg_noise_71_smooth = zeros([num_trials 1]);
all_mean_curvature_m1_noise_71_smooth = zeros([num_trials 1]);
all_total_distance_mm_noise_71_smooth = zeros([num_trials 1]);
all_mean_alpha_deg_noise_aniso_24_smooth = zeros([num_trials 1]);
all_mean_curvature_m1_noise_aniso_24_smooth = zeros([num_trials 1]);
all_total_distance_mm_noise_aniso_24_smooth = zeros([num_trials 1]);
all_mean_alpha_deg_noise_aniso_71_smooth = zeros([num_trials 1]);
all_mean_curvature_m1_noise_aniso_71_smooth = zeros([num_trials 1]);
all_total_distance_mm_noise_aniso_71_smooth = zeros([num_trials 1]);
fiber_similarity_24_mean = zeros([num_trials 1]);
fiber_similarity_aniso_24_mean = zeros([num_trials 1]);
fiber_similarity_smooth_24_mean = zeros([num_trials 1]);
fiber_similarity_aniso_smooth_24_mean = zeros([num_trials 1]);
fiber_similarity_71_mean = zeros([num_trials 1]);
fiber_similarity_aniso_71_mean = zeros([num_trials 1]);
fiber_similarity_smooth_71_mean = zeros([num_trials 1]);
fiber_similarity_aniso_smooth_71_mean = zeros([num_trials 1]);
euclid_distance_24_mean = zeros([num_trials 1]);
euclid_distance_aniso_24_mean = zeros([num_trials 1]);
euclid_distance_smooth_24_mean = zeros([num_trials 1]);
euclid_distance_aniso_smooth_24_mean = zeros([num_trials 1]);
euclid_distance_71_mean = zeros([num_trials 1]);
euclid_distance_aniso_71_mean = zeros([num_trials 1]);
euclid_distance_smooth_71_mean = zeros([num_trials 1]);
euclid_distance_aniso_smooth_71_mean = zeros([num_trials 1]);
fiber_similarity_24_std = zeros([num_trials 1]);
fiber_similarity_aniso_24_std = zeros([num_trials 1]);
fiber_similarity_smooth_24_std = zeros([num_trials 1]);
fiber_similarity_aniso_smooth_24_std = zeros([num_trials 1]);
fiber_similarity_71_std = zeros([num_trials 1]);
fiber_similarity_aniso_71_std = zeros([num_trials 1]);
fiber_similarity_smooth_71_std = zeros([num_trials 1]);
fiber_similarity_aniso_smooth_71_std = zeros([num_trials 1]);
euclid_distance_24_std = zeros([num_trials 1]);
euclid_distance_aniso_24_std = zeros([num_trials 1]);
euclid_distance_smooth_24_std = zeros([num_trials 1]);
euclid_distance_aniso_smooth_24_std = zeros([num_trials 1]);
euclid_distance_71_std = zeros([num_trials 1]);
euclid_distance_aniso_71_std = zeros([num_trials 1]);
euclid_distance_smooth_71_std = zeros([num_trials 1]);
euclid_distance_aniso_smooth_71_std = zeros([num_trials 1]);
percent_empty_fiber_tracts_24 = zeros([num_trials 1]);
percent_empty_fiber_tracts_aniso_24 = zeros([num_trials 1]);
percent_empty_fiber_tracts_smooth_24 = zeros([num_trials 1]);
percent_empty_fiber_tracts_aniso_smooth_24 = zeros([num_trials 1]);
percent_empty_fiber_tracts_71 = zeros([num_trials 1]);
percent_empty_fiber_tracts_aniso_71 = zeros([num_trials 1]);
percent_empty_fiber_tracts_smooth_71 = zeros([num_trials 1]);
percent_empty_fiber_tracts_aniso_smooth_71 = zeros([num_trials 1]);
smooth_parameter_24 = zeros([num_trials 1]);
smooth_parameter_71 = zeros([num_trials 1]);

tic
for trial_cntr=1:num_trials

    % form a noisy image
    noise=randn(size(noise_free_img));
    noisy_img_24 = abs(noise_free_img + noise/24);
    noisy_img_71 = abs(noise_free_img + noise/71);
    loop_b0_snr24 = squeeze(noisy_img_24(:,:, 25,1));
    loop_b0_snr24 = loop_b0_snr24(muscle_mask(:,:, 25)==1);
    all_snr24(trial_cntr) = 1/std(loop_b0_snr24);
    loop_b0_snr71 = squeeze(noisy_img_71(:,:, 25,1));
    loop_b0_snr71 = loop_b0_snr71(muscle_mask(:,:, 25)==1);
    all_snr71(trial_cntr) = 1/std(loop_b0_snr71);
    alpha_noise24_loop = zeros(size(noise_free_img,1:3));
    alpha_noise24_loop_aniso = zeros(size(noise_free_img,1:3));
    alpha_noise71_loop = zeros(size(noise_free_img,1:3));
    alpha_noise71_loop_aniso = zeros(size(noise_free_img,1:3));
    e1map_24 = zeros([size(noise_free_img,1:3) 3]);
    e1map_71 = zeros([size(noise_free_img,1:3) 3]);
    e1map_aniso24 = zeros([size(noise_free_img,1:3) 3]);
    e1map_aniso71 = zeros([size(noise_free_img,1:3) 3]);
    famap_24 = zeros(size(noise_free_img,1:3));
    famap_71 = zeros(size(noise_free_img,1:3));
    famap_aniso24 = zeros(size(noise_free_img,1:3));
    famap_aniso71 = zeros(size(noise_free_img,1:3));
    alpha_noise24_loop_smooth = zeros(size(noise_free_img,1:3));
    alpha_noise24_loop_aniso_smooth = zeros(size(noise_free_img,1:3));
    alpha_noise71_loop_smooth = zeros(size(noise_free_img,1:3));
    alpha_noise71_loop_aniso_smooth = zeros(size(noise_free_img,1:3));

    % Denoise data

    %set parameters based on Buck et al
    noise = 5;
    sigma = noise/100;
    rho = 2*sigma;
    delta_t = noise*3/44;
    schemetype = 'implicit_multisplitting';
    isfasteig  = true;
    isnormg = false;
    dti_res = [1 1 6];
    
    noisy_img_aniso_24 = aniso4D_smoothing(noisy_img_24, sigma, rho, delta_t, dti_res, schemetype, isnormg, isfasteig);
    noisy_img_aniso_71 = aniso4D_smoothing(noisy_img_71, sigma, rho, delta_t, dti_res, schemetype, isnormg, isfasteig);

    % get noisy diffusion tensors
    for r_img=1:50
        for c_img=1:50
            for s_img=1:40
                if muscle_apo_mask(r_img,c_img,s_img)>0

                    signal_v_24 = squeeze(noisy_img_24(r_img, c_img, s_img, :));
                    signal_v_71 = squeeze(noisy_img_71(r_img, c_img, s_img, :));
                    signal_v_24aniso = squeeze(noisy_img_aniso_24(r_img, c_img, s_img, :));
                    signal_v_71aniso = squeeze(noisy_img_aniso_71(r_img, c_img, s_img, :));

                    D_24 = signal2tensor2(signal_v_24, dir_m, 475);
                    [E, L] = svd(D_24);

                    L_v = diag(L);
                    md = mean(L_v);
                    fa = sqrt(3/2) * sqrt(((L_v(1) - md)^2 + (L_v(2) - md)^2 + (L_v(3) - md)^2) / ...
                        (L_v(1)^2 + L_v(2)^2 + L_v(3)^2));

                    loop_E1 = E(:,L_v==max(L_v));
                    if loop_E1(3)<0
                        loop_E1 = -loop_E1;
                    end
                    
                    e1map_24(r_img,c_img,s_img,:) = loop_E1;
                    famap_24(r_img,c_img,s_img,:) = fa;
                    alpha_noise24_loop(r_img,c_img,s_img) = dot(loop_E1,squeeze(e1_map(r_img,c_img,s_img,:)));

                    D_24aniso = signal2tensor2(signal_v_24aniso, dir_m, 475);
                    [E, L] = svd(D_24aniso);

                    L_v = diag(L);
                    fa = sqrt(3/2) * sqrt(((L_v(1) - md)^2 + (L_v(2) - md)^2 + (L_v(3) - md)^2) / ...
                        (L_v(1)^2 + L_v(2)^2 + L_v(3)^2));

                    loop_E1 = E(:,L_v==max(L_v));
                    if loop_E1(3)<0
                        loop_E1 = -loop_E1;
                    end
                    
                    e1map_aniso24(r_img,c_img,s_img,:) = loop_E1;
                    famap_aniso24(r_img,c_img,s_img,:) = fa;
                    alpha_noise24_loop_aniso(r_img,c_img,s_img) = dot(loop_E1,squeeze(e1_map(r_img,c_img,s_img,:)));

                    D_71 = signal2tensor2(signal_v_71, dir_m, 475);
                    [E, L] = svd(D_71);

                    L_v = diag(L);
                    fa = sqrt(3/2) * sqrt(((L_v(1) - md)^2 + (L_v(2) - md)^2 + (L_v(3) - md)^2) / ...
                        (L_v(1)^2 + L_v(2)^2 + L_v(3)^2));

                    loop_E1 = E(:,L_v==max(L_v));
                    if loop_E1(3)<0
                        loop_E1 = -loop_E1;
                    end
                    
                    e1map_71(r_img,c_img,s_img,:) = loop_E1;
                    famap_71(r_img,c_img,s_img,:) = fa;
                    alpha_noise71_loop(r_img,c_img,s_img) = dot(loop_E1,squeeze(e1_map(r_img,c_img,s_img,:)));

                    D_71aniso = signal2tensor2(signal_v_71aniso, dir_m, 475);
                    [E, L] = svd(D_71aniso);

                    L_v = diag(L);
                    md = mean(L_v);
                    fa = sqrt(3/2) * sqrt(((L_v(1) - md)^2 + (L_v(2) - md)^2 + (L_v(3) - md)^2) / ...
                        (L_v(1)^2 + L_v(2)^2 + L_v(3)^2));

                    loop_E1 = E(:,L_v==max(L_v));
                    if loop_E1(3)<0
                        loop_E1 = -loop_E1;
                    end
                    
                    e1map_aniso71(r_img,c_img,s_img,:) = loop_E1;
                    famap_aniso71(r_img,c_img,s_img,:) = fa;
                    alpha_noise71_loop_aniso(r_img,c_img,s_img) = dot(loop_E1,squeeze(e1_map(r_img,c_img,s_img,:)));

                end
            end
        end
    end

    angle_noise_24(:,:,:,trial_cntr) = acosd(alpha_noise24_loop).*muscle_mask;
    angle_noise_aniso_24(:,:,:,trial_cntr) = acosd(alpha_noise24_loop_aniso).*muscle_mask;
    angle_noise_71(:,:,:,trial_cntr) = acosd(alpha_noise71_loop).*muscle_mask;
    angle_noise_aniso_71(:,:,:,trial_cntr) = acosd(alpha_noise71_loop_aniso).*muscle_mask;

    % Smooth eigenvector field 

    [e1map_24_smooth,smooth_parameter]=smoothn({e1map_24(:,:,:,1),e1map_24(:,:,:,2),e1map_24(:,:,:,3)});
    smooth_parameter_24(trial_cntr) = smooth_parameter;
    [e1map_71_smooth,smooth_parameter]=smoothn({e1map_71(:,:,:,1),e1map_71(:,:,:,2),e1map_24(:,:,:,3)});
    smooth_parameter_71(trial_cntr) = smooth_parameter;
    e1map_24_aniso_smooth=smoothn({e1map_aniso24(:,:,:,1),e1map_aniso24(:,:,:,2),e1map_aniso24(:,:,:,3)});
    e1map_71_aniso_smooth=smoothn({e1map_aniso71(:,:,:,1),e1map_aniso71(:,:,:,2),e1map_aniso71(:,:,:,3)});

    e1map_24_smooth=cat(4,e1map_24_smooth{1},e1map_24_smooth{2},e1map_24_smooth{3});
    e1map_71_smooth=cat(4,e1map_71_smooth{1},e1map_71_smooth{2},e1map_71_smooth{3});
    e1map_24_aniso_smooth=cat(4,e1map_24_aniso_smooth{1},e1map_24_aniso_smooth{2},e1map_24_aniso_smooth{3});
    e1map_71_aniso_smooth=cat(4,e1map_71_aniso_smooth{1},e1map_71_aniso_smooth{2},e1map_71_aniso_smooth{3});

    for r_img=1:50
        for c_img=1:50
            for s_img=1:40
                if muscle_apo_mask(r_img,c_img,s_img)>0

                    E1_24_smooth = squeeze(e1map_24_smooth(r_img,c_img,s_img,:));
                    E1_71_smooth = squeeze(e1map_71_smooth(r_img,c_img,s_img,:));
                    E1_24_aniso_smooth = squeeze(e1map_24_aniso_smooth(r_img,c_img,s_img,:));
                    E1_71_aniso_smooth = squeeze(e1map_71_aniso_smooth(r_img,c_img,s_img,:));

                    E1_24_smooth = E1_24_smooth/norm(E1_24_smooth);
                    E1_71_smooth = E1_71_smooth/norm(E1_71_smooth);
                    E1_24_aniso_smooth = E1_24_aniso_smooth/norm(E1_24_aniso_smooth);
                    E1_71_aniso_smooth = E1_71_aniso_smooth/norm(E1_71_aniso_smooth);

                    e1map_24_smooth(r_img,c_img,s_img,:)=E1_24_smooth;
                    e1map_71_smooth(r_img,c_img,s_img,:)=E1_71_smooth;
                    e1map_24_aniso_smooth(r_img,c_img,s_img,:)=E1_24_aniso_smooth;
                    e1map_71_aniso_smooth(r_img,c_img,s_img,:)=E1_71_aniso_smooth;

                    alpha_noise24_loop_smooth(r_img,c_img,s_img) = dot(E1_24_smooth,squeeze(e1_map(r_img,c_img,s_img,:)));
                    alpha_noise24_loop_aniso_smooth(r_img,c_img,s_img) = dot(E1_24_aniso_smooth,squeeze(e1_map(r_img,c_img,s_img,:)));
                    alpha_noise71_loop_smooth(r_img,c_img,s_img) = dot(E1_71_smooth,squeeze(e1_map(r_img,c_img,s_img,:)));
                    alpha_noise71_loop_aniso_smooth(r_img,c_img,s_img) = dot(E1_71_aniso_smooth,squeeze(e1_map(r_img,c_img,s_img,:)));

                end
            end
        end
    end

    angle_noise_24_smooth(:,:,:,trial_cntr) = acosd(alpha_noise24_loop_smooth).*muscle_mask;
    angle_noise_aniso_24_smooth(:,:,:,trial_cntr) = acosd(alpha_noise24_loop_aniso_smooth).*muscle_mask;
    angle_noise_71_smooth(:,:,:,trial_cntr) = acosd(alpha_noise71_loop_smooth).*muscle_mask;
    angle_noise_aniso_71_smooth(:,:,:,trial_cntr) = acosd(alpha_noise71_loop_aniso_smooth).*muscle_mask;

        %%%%%%%%%%%%%%%% fiber tracking - aponeurosis seeding %%%%%%%%%%%%%%%%%%
    % SNR 24
    [apo_fiber_all, ~, apo_stop_list] = ...
        fiber_track_v20(apo_ft_options, cat(4,e1map_24,famap_24), composite_mask, roi_mesh_restricted);
    [apo_smoothed_fiber_all, ~, apo_smoothed_fiber_all_mm, ~, ~, ~, ~, ~, ~] = ...
        fiber_smoother_v14a(apo_fiber_all, all_fs_options);

    % prep for filtering based on same seed location
    apo_smoothed_fiber_all_mm_final = zeros(size(apo_smoothed_fiber_all_mm));   %will hold tracts after seed-location filtering
    apo_smoothed_fiber_all_final = zeros(size(apo_smoothed_fiber_all));

    for step_cntr = 1:size(apo_smoothed_fiber_all,3)
        for dim_cntr = 1:size(apo_smoothed_fiber_all,4)
            apo_smoothed_fiber_all_final(:,:,step_cntr,dim_cntr) = ...
                squeeze(apo_smoothed_fiber_all(:,:,step_cntr,dim_cntr)).*apo_seeds_final_mask;
            apo_smoothed_fiber_all_mm_final(:,:,step_cntr,dim_cntr) = ...
                squeeze(apo_smoothed_fiber_all_mm(:,:,step_cntr,dim_cntr)).*apo_seeds_final_mask;
        end
    end

    [apo_angle_list_final, apo_distance_list_final, apo_curvature_list_final, ~, apo_n_points_final, apo_area_final] = ...
        fiber_quantifier_v20a(vxl_fq_options, apo_smoothed_fiber_all_final, roi_mesh_restricted, muscle_mask);

    % calculate raw summary statistics
    apo_total_distance_final = squeeze(max(apo_distance_list_final, [], 3));
    apo_total_distance_final(apo_total_distance_final<5)=0;

    apo_mean_alpha_final = sum(squeeze(apo_angle_list_final(:,:,:,1)), 3)./squeeze(apo_n_points_final(:,:,2));
    apo_mean_alpha_final(apo_total_distance_final<5)=0;

    apo_mean_curvature_final = sum(apo_curvature_list_final, 3)./squeeze(apo_n_points_final(:,:,3));
    apo_mean_curvature_final(apo_total_distance_final<5)=0;

    %get rid of zeros and NaNs
    apo_mean_alpha_final(isnan(apo_mean_alpha_final)) = 0;
    apo_mean_alpha_final = nonzeros(apo_mean_alpha_final);
    apo_mean_curvature_final(isnan(apo_mean_curvature_final)) = 0;
    apo_mean_curvature_final = nonzeros(apo_mean_curvature_final);
    apo_total_distance_final(isnan(apo_total_distance_final)) = 0;
    apo_total_distance_final = nonzeros(apo_total_distance_final);

    % save summary data for each time through the loop
    all_mean_alpha_deg_noise_24(trial_cntr,1) = mean(nonzeros(apo_mean_alpha_final));
    all_mean_curvature_m1_noise_24(trial_cntr,1) = mean(nonzeros(apo_mean_curvature_final));
    all_total_distance_mm_noise_24(trial_cntr,1) = mean(nonzeros(apo_total_distance_final));

    % Compute fiber similarity index, euclidean distance, and corresponding
    % segment ratio
    for fiber_cntr = 1:size(apo_smoothed_fiber_all_mm_final,1)*size(apo_smoothed_fiber_all_mm_final,2)
        [row,col] = ind2sub([size(apo_smoothed_fiber_all_mm_final,1) size(apo_smoothed_fiber_all_mm_final,2)],fiber_cntr);
        fiber_noise_free = squeeze(apo_smoothed_fiber_all_mm_noise_free(row,col,:,:));
        fiber_noisy = squeeze(apo_smoothed_fiber_all_mm_final(row,col,:,:));
        if mean(fiber_noise_free,'all') ~=0 
            fiber_noise_free ( all(~fiber_noise_free ,2), : ) = []; % remove rows containing zero elements
            fiber_noisy ( all(~fiber_noisy ,2), : ) = [];
            [S, Rcs, D, Dpos] = fit_similarity(fiber_noise_free,fiber_noisy,1); % 1 mm of in-plane resolution
            fiber_similarity_24(row,col,trial_cntr) = S;
            euclid_distance_24(row,col,trial_cntr) = D;
            cor_segment_24(row,col,trial_cntr) = Rcs;
        end
    end

    % SNR 24 - anisotropic smoothing
    [apo_fiber_all, ~, apo_stop_list] = ...
        fiber_track_v20(apo_ft_options, cat(4,e1map_aniso24,famap_aniso24), composite_mask, roi_mesh_restricted);
    [apo_smoothed_fiber_all, ~, apo_smoothed_fiber_all_mm, ~, ~, ~, ~, ~, ~] = ...
        fiber_smoother_v14a(apo_fiber_all, all_fs_options);

    % prep for filtering based on same seed location
    apo_smoothed_fiber_all_mm_final = zeros(size(apo_smoothed_fiber_all_mm));   %will hold tracts after seed-location filtering
    apo_smoothed_fiber_all_final = zeros(size(apo_smoothed_fiber_all));

    for step_cntr = 1:size(apo_smoothed_fiber_all,3)
        for dim_cntr = 1:size(apo_smoothed_fiber_all,4)
            apo_smoothed_fiber_all_final(:,:,step_cntr,dim_cntr) = ...
                squeeze(apo_smoothed_fiber_all(:,:,step_cntr,dim_cntr)).*apo_seeds_final_mask;
            apo_smoothed_fiber_all_mm_final(:,:,step_cntr,dim_cntr) = ...
                squeeze(apo_smoothed_fiber_all_mm(:,:,step_cntr,dim_cntr)).*apo_seeds_final_mask;
        end
    end

    [apo_angle_list_final, apo_distance_list_final, apo_curvature_list_final, ~, apo_n_points_final, apo_area_final] = ...
        fiber_quantifier_v20a(vxl_fq_options, apo_smoothed_fiber_all_final, roi_mesh_restricted, muscle_mask);

    % calculate raw summary statistics
    apo_total_distance_final = squeeze(max(apo_distance_list_final, [], 3));
    apo_total_distance_final(apo_total_distance_final<5)=0;

    apo_mean_alpha_final = sum(squeeze(apo_angle_list_final(:,:,:,1)), 3)./squeeze(apo_n_points_final(:,:,2));
    apo_mean_alpha_final(apo_total_distance_final<5)=0;

    apo_mean_curvature_final = sum(apo_curvature_list_final, 3)./squeeze(apo_n_points_final(:,:,3));
    apo_mean_curvature_final(apo_total_distance_final<5)=0;

    %get rid of zeros and NaNs
    apo_mean_alpha_final(isnan(apo_mean_alpha_final)) = 0;
    apo_mean_alpha_final = nonzeros(apo_mean_alpha_final);
    apo_mean_curvature_final(isnan(apo_mean_curvature_final)) = 0;
    apo_mean_curvature_final = nonzeros(apo_mean_curvature_final);
    apo_total_distance_final(isnan(apo_total_distance_final)) = 0;
    apo_total_distance_final = nonzeros(apo_total_distance_final);

    % save summary data for each time through the loop
    all_mean_alpha_deg_noise_aniso_24(trial_cntr,1) = mean(nonzeros(apo_mean_alpha_final));
    all_mean_curvature_m1_noise_aniso_24(trial_cntr,1) = mean(nonzeros(apo_mean_curvature_final));
    all_total_distance_mm_noise_aniso_24(trial_cntr,1) = mean(nonzeros(apo_total_distance_final));

    % Compute fiber similarity index, euclidean distance, and corresponding
    % segment ratio
    for fiber_cntr = 1:size(apo_smoothed_fiber_all_mm_final,1)*size(apo_smoothed_fiber_all_mm_final,2)
        [row,col] = ind2sub([size(apo_smoothed_fiber_all_mm_final,1) size(apo_smoothed_fiber_all_mm_final,2)],fiber_cntr);
        fiber_noise_free = squeeze(apo_smoothed_fiber_all_mm_noise_free(row,col,:,:));
        fiber_noisy = squeeze(apo_smoothed_fiber_all_mm_final(row,col,:,:));
        if mean(fiber_noise_free,'all') ~=0 
            fiber_noise_free ( all(~fiber_noise_free ,2), : ) = []; % remove rows containing zero elements
            fiber_noisy ( all(~fiber_noisy ,2), : ) = [];
            [S, Rcs, D, Dpos] = fit_similarity(fiber_noise_free,fiber_noisy,1); % 1 mm of in-plane resolution
            fiber_similarity_aniso_24(row,col,trial_cntr) = S;
            euclid_distance_aniso_24(row,col,trial_cntr) = D;
            cor_segment_aniso_24(row,col,trial_cntr) = Rcs;
        end
    end

    % SNR 24 + smoothn
    [apo_fiber_all, ~, apo_stop_list] = ...
        fiber_track_v20(apo_ft_options, cat(4,e1map_24_smooth,famap_24), composite_mask, roi_mesh_restricted);
    [apo_smoothed_fiber_all, ~, apo_smoothed_fiber_all_mm, ~, ~, ~, ~, ~, ~] = ...
        fiber_smoother_v14a(apo_fiber_all, all_fs_options);

    % prep for filtering based on same seed location
    apo_smoothed_fiber_all_mm_final = zeros(size(apo_smoothed_fiber_all_mm));   %will hold tracts after seed-location filtering
    apo_smoothed_fiber_all_final = zeros(size(apo_smoothed_fiber_all));

    for step_cntr = 1:size(apo_smoothed_fiber_all,3)
        for dim_cntr = 1:size(apo_smoothed_fiber_all,4)
            apo_smoothed_fiber_all_final(:,:,step_cntr,dim_cntr) = ...
                squeeze(apo_smoothed_fiber_all(:,:,step_cntr,dim_cntr)).*apo_seeds_final_mask;
            apo_smoothed_fiber_all_mm_final(:,:,step_cntr,dim_cntr) = ...
                squeeze(apo_smoothed_fiber_all_mm(:,:,step_cntr,dim_cntr)).*apo_seeds_final_mask;
        end
    end

    [apo_angle_list_final, apo_distance_list_final, apo_curvature_list_final, ~, apo_n_points_final, apo_area_final] = ...
        fiber_quantifier_v20a(vxl_fq_options, apo_smoothed_fiber_all_final, roi_mesh_restricted, muscle_mask);

    % calculate raw summary statistics
    apo_total_distance_final = squeeze(max(apo_distance_list_final, [], 3));
    apo_total_distance_final(apo_total_distance_final<5)=0;

    apo_mean_alpha_final = sum(squeeze(apo_angle_list_final(:,:,:,1)), 3)./squeeze(apo_n_points_final(:,:,2));
    apo_mean_alpha_final(apo_total_distance_final<5)=0;

    apo_mean_curvature_final = sum(apo_curvature_list_final, 3)./squeeze(apo_n_points_final(:,:,3));
    apo_mean_curvature_final(apo_total_distance_final<5)=0;

    %get rid of zeros and NaNs
    apo_mean_alpha_final(isnan(apo_mean_alpha_final)) = 0;
    apo_mean_alpha_final = nonzeros(apo_mean_alpha_final);
    apo_mean_curvature_final(isnan(apo_mean_curvature_final)) = 0;
    apo_mean_curvature_final = nonzeros(apo_mean_curvature_final);
    apo_total_distance_final(isnan(apo_total_distance_final)) = 0;
    apo_total_distance_final = nonzeros(apo_total_distance_final);

    % save summary data for each time through the loop
    all_mean_alpha_deg_noise_24_smooth(trial_cntr,1) = mean(nonzeros(apo_mean_alpha_final));
    all_mean_curvature_m1_noise_24_smooth(trial_cntr,1) = mean(nonzeros(apo_mean_curvature_final));
    all_total_distance_mm_noise_24_smooth(trial_cntr,1) = mean(nonzeros(apo_total_distance_final));

    % Compute fiber similarity index, euclidean distance, and corresponding
    % segment ratio
    for fiber_cntr = 1:size(apo_smoothed_fiber_all_mm_final,1)*size(apo_smoothed_fiber_all_mm_final,2)
        [row,col] = ind2sub([size(apo_smoothed_fiber_all_mm_final,1) size(apo_smoothed_fiber_all_mm_final,2)],fiber_cntr);
        fiber_noise_free = squeeze(apo_smoothed_fiber_all_mm_noise_free(row,col,:,:));
        fiber_noisy = squeeze(apo_smoothed_fiber_all_mm_final(row,col,:,:));
        if mean(fiber_noise_free,'all') ~=0 
            fiber_noise_free ( all(~fiber_noise_free ,2), : ) = []; % remove rows containing zero elements
            fiber_noisy ( all(~fiber_noisy ,2), : ) = [];
            [S, Rcs, D, Dpos] = fit_similarity(fiber_noise_free,fiber_noisy,1); % 1 mm of in-plane resolution
            fiber_similarity_24_smooth(row,col,trial_cntr) = S;
            euclid_distance_24_smooth(row,col,trial_cntr) = D;
            cor_segment_24_smooth(row,col,trial_cntr) = Rcs;
        end
    end

    % SNR 24 - anisotropic smoothing + smoothn
    [apo_fiber_all, ~, apo_stop_list] = ...
        fiber_track_v20(apo_ft_options, cat(4,e1map_24_aniso_smooth,famap_aniso24), composite_mask, roi_mesh_restricted);
    [apo_smoothed_fiber_all, ~, apo_smoothed_fiber_all_mm, ~, ~, ~, ~, ~, ~] = ...
        fiber_smoother_v14a(apo_fiber_all, all_fs_options);

    % prep for filtering based on same seed location
    apo_smoothed_fiber_all_mm_final = zeros(size(apo_smoothed_fiber_all_mm));   %will hold tracts after seed-location filtering
    apo_smoothed_fiber_all_final = zeros(size(apo_smoothed_fiber_all));

    for step_cntr = 1:size(apo_smoothed_fiber_all,3)
        for dim_cntr = 1:size(apo_smoothed_fiber_all,4)
            apo_smoothed_fiber_all_final(:,:,step_cntr,dim_cntr) = ...
                squeeze(apo_smoothed_fiber_all(:,:,step_cntr,dim_cntr)).*apo_seeds_final_mask;
            apo_smoothed_fiber_all_mm_final(:,:,step_cntr,dim_cntr) = ...
                squeeze(apo_smoothed_fiber_all_mm(:,:,step_cntr,dim_cntr)).*apo_seeds_final_mask;
        end
    end

    [apo_angle_list_final, apo_distance_list_final, apo_curvature_list_final, ~, apo_n_points_final, apo_area_final] = ...
        fiber_quantifier_v20a(vxl_fq_options, apo_smoothed_fiber_all_final, roi_mesh_restricted, muscle_mask);

    % calculate raw summary statistics
    apo_total_distance_final = squeeze(max(apo_distance_list_final, [], 3));
    apo_total_distance_final(apo_total_distance_final<5)=0;

    apo_mean_alpha_final = sum(squeeze(apo_angle_list_final(:,:,:,1)), 3)./squeeze(apo_n_points_final(:,:,2));
    apo_mean_alpha_final(apo_total_distance_final<5)=0;

    apo_mean_curvature_final = sum(apo_curvature_list_final, 3)./squeeze(apo_n_points_final(:,:,3));
    apo_mean_curvature_final(apo_total_distance_final<5)=0;

    %get rid of zeros and NaNs
    apo_mean_alpha_final(isnan(apo_mean_alpha_final)) = 0;
    apo_mean_alpha_final = nonzeros(apo_mean_alpha_final);
    apo_mean_curvature_final(isnan(apo_mean_curvature_final)) = 0;
    apo_mean_curvature_final = nonzeros(apo_mean_curvature_final);
    apo_total_distance_final(isnan(apo_total_distance_final)) = 0;
    apo_total_distance_final = nonzeros(apo_total_distance_final);

    % save summary data for each time through the loop
    all_mean_alpha_deg_noise_aniso_24_smooth(trial_cntr,1) = mean(nonzeros(apo_mean_alpha_final));
    all_mean_curvature_m1_noise_aniso_24_smooth(trial_cntr,1) = mean(nonzeros(apo_mean_curvature_final));
    all_total_distance_mm_noise_aniso_24_smooth(trial_cntr,1) = mean(nonzeros(apo_total_distance_final));

    % Compute fiber similarity index, euclidean distance, and corresponding
    % segment ratio
    for fiber_cntr = 1:size(apo_smoothed_fiber_all_mm_final,1)*size(apo_smoothed_fiber_all_mm_final,2)
        [row,col] = ind2sub([size(apo_smoothed_fiber_all_mm_final,1) size(apo_smoothed_fiber_all_mm_final,2)],fiber_cntr);
        fiber_noise_free = squeeze(apo_smoothed_fiber_all_mm_noise_free(row,col,:,:));
        fiber_noisy = squeeze(apo_smoothed_fiber_all_mm_final(row,col,:,:));
        if mean(fiber_noise_free,'all') ~=0 
            fiber_noise_free ( all(~fiber_noise_free ,2), : ) = []; % remove rows containing zero elements
            fiber_noisy ( all(~fiber_noisy ,2), : ) = [];
            [S, Rcs, D, Dpos] = fit_similarity(fiber_noise_free,fiber_noisy,1); % 1 mm of in-plane resolution
            fiber_similarity_aniso_24_smooth(row,col,trial_cntr) = S;
            euclid_distance_aniso_24_smooth(row,col,trial_cntr) = D;
            cor_segment_aniso_24_smooth(row,col,trial_cntr) = Rcs;
        end
    end

    % SNR 71
    [apo_fiber_all, ~, apo_stop_list] = ...
        fiber_track_v20(apo_ft_options, cat(4,e1map_71,famap_71), composite_mask, roi_mesh_restricted);
    [apo_smoothed_fiber_all, ~, apo_smoothed_fiber_all_mm, ~, ~, ~, ~, ~, ~] = ...
        fiber_smoother_v14a(apo_fiber_all, all_fs_options);

    % prep for filtering based on same seed location
    apo_smoothed_fiber_all_mm_final = zeros(size(apo_smoothed_fiber_all_mm));   %will hold tracts after seed-location filtering
    apo_smoothed_fiber_all_final = zeros(size(apo_smoothed_fiber_all));

    for step_cntr = 1:size(apo_smoothed_fiber_all,3)
        for dim_cntr = 1:size(apo_smoothed_fiber_all,4)
            apo_smoothed_fiber_all_final(:,:,step_cntr,dim_cntr) = ...
                squeeze(apo_smoothed_fiber_all(:,:,step_cntr,dim_cntr)).*apo_seeds_final_mask;
            apo_smoothed_fiber_all_mm_final(:,:,step_cntr,dim_cntr) = ...
                squeeze(apo_smoothed_fiber_all_mm(:,:,step_cntr,dim_cntr)).*apo_seeds_final_mask;
        end
    end

    [apo_angle_list_final, apo_distance_list_final, apo_curvature_list_final, ~, apo_n_points_final, apo_area_final] = ...
        fiber_quantifier_v20a(vxl_fq_options, apo_smoothed_fiber_all_final, roi_mesh_restricted, muscle_mask);

    % calculate raw summary statistics
    apo_total_distance_final = squeeze(max(apo_distance_list_final, [], 3));
    apo_total_distance_final(apo_total_distance_final<5)=0;

    apo_mean_alpha_final = sum(squeeze(apo_angle_list_final(:,:,:,1)), 3)./squeeze(apo_n_points_final(:,:,2));
    apo_mean_alpha_final(apo_total_distance_final<5)=0;

    apo_mean_curvature_final = sum(apo_curvature_list_final, 3)./squeeze(apo_n_points_final(:,:,3));
    apo_mean_curvature_final(apo_total_distance_final<5)=0;

    %get rid of zeros and NaNs
    apo_mean_alpha_final(isnan(apo_mean_alpha_final)) = 0;
    apo_mean_alpha_final = nonzeros(apo_mean_alpha_final);
    apo_mean_curvature_final(isnan(apo_mean_curvature_final)) = 0;
    apo_mean_curvature_final = nonzeros(apo_mean_curvature_final);
    apo_total_distance_final(isnan(apo_total_distance_final)) = 0;
    apo_total_distance_final = nonzeros(apo_total_distance_final);

    % save summary data for each time through the loop
    all_mean_alpha_deg_noise_71(trial_cntr,1) = mean(nonzeros(apo_mean_alpha_final));
    all_mean_curvature_m1_noise_71(trial_cntr,1) = mean(nonzeros(apo_mean_curvature_final));
    all_total_distance_mm_noise_71(trial_cntr,1) = mean(nonzeros(apo_total_distance_final));

    % Compute fiber similarity index, euclidean distance, and corresponding
    % segment ratio
    for fiber_cntr = 1:size(apo_smoothed_fiber_all_mm_final,1)*size(apo_smoothed_fiber_all_mm_final,2)
        [row,col] = ind2sub([size(apo_smoothed_fiber_all_mm_final,1) size(apo_smoothed_fiber_all_mm_final,2)],fiber_cntr);
        fiber_noise_free = squeeze(apo_smoothed_fiber_all_mm_noise_free(row,col,:,:));
        fiber_noisy = squeeze(apo_smoothed_fiber_all_mm_final(row,col,:,:));
        if mean(fiber_noise_free,'all') ~=0 
            fiber_noise_free ( all(~fiber_noise_free ,2), : ) = []; % remove rows containing zero elements
            fiber_noisy ( all(~fiber_noisy ,2), : ) = [];
            [S, Rcs, D, Dpos] = fit_similarity(fiber_noise_free,fiber_noisy,1); % 1 mm of in-plane resolution
            fiber_similarity_71(row,col,trial_cntr) = S;
            euclid_distance_71(row,col,trial_cntr) = D;
            cor_segment_71(row,col,trial_cntr) = Rcs;
        end
    end

    % SNR 71 - anisotropic smoothing
    [apo_fiber_all, ~, apo_stop_list] = ...
        fiber_track_v20(apo_ft_options, cat(4,e1map_aniso71,famap_aniso71), composite_mask, roi_mesh_restricted);
    [apo_smoothed_fiber_all, ~, apo_smoothed_fiber_all_mm, ~, ~, ~, ~, ~, ~] = ...
        fiber_smoother_v14a(apo_fiber_all, all_fs_options);

    % prep for filtering based on same seed location
    apo_smoothed_fiber_all_mm_final = zeros(size(apo_smoothed_fiber_all_mm));   %will hold tracts after seed-location filtering
    apo_smoothed_fiber_all_final = zeros(size(apo_smoothed_fiber_all));

    for step_cntr = 1:size(apo_smoothed_fiber_all,3)
        for dim_cntr = 1:size(apo_smoothed_fiber_all,4)
            apo_smoothed_fiber_all_final(:,:,step_cntr,dim_cntr) = ...
                squeeze(apo_smoothed_fiber_all(:,:,step_cntr,dim_cntr)).*apo_seeds_final_mask;
            apo_smoothed_fiber_all_mm_final(:,:,step_cntr,dim_cntr) = ...
                squeeze(apo_smoothed_fiber_all_mm(:,:,step_cntr,dim_cntr)).*apo_seeds_final_mask;
        end
    end

    [apo_angle_list_final, apo_distance_list_final, apo_curvature_list_final, ~, apo_n_points_final, apo_area_final] = ...
        fiber_quantifier_v20a(vxl_fq_options, apo_smoothed_fiber_all_final, roi_mesh_restricted, muscle_mask);

    % calculate raw summary statistics
    apo_total_distance_final = squeeze(max(apo_distance_list_final, [], 3));
    apo_total_distance_final(apo_total_distance_final<5)=0;

    apo_mean_alpha_final = sum(squeeze(apo_angle_list_final(:,:,:,1)), 3)./squeeze(apo_n_points_final(:,:,2));
    apo_mean_alpha_final(apo_total_distance_final<5)=0;

    apo_mean_curvature_final = sum(apo_curvature_list_final, 3)./squeeze(apo_n_points_final(:,:,3));
    apo_mean_curvature_final(apo_total_distance_final<5)=0;

    %get rid of zeros and NaNs
    apo_mean_alpha_final(isnan(apo_mean_alpha_final)) = 0;
    apo_mean_alpha_final = nonzeros(apo_mean_alpha_final);
    apo_mean_curvature_final(isnan(apo_mean_curvature_final)) = 0;
    apo_mean_curvature_final = nonzeros(apo_mean_curvature_final);
    apo_total_distance_final(isnan(apo_total_distance_final)) = 0;
    apo_total_distance_final = nonzeros(apo_total_distance_final);

    % save summary data for each time through the loop
    all_mean_alpha_deg_noise_aniso_71(trial_cntr,1) = mean(nonzeros(apo_mean_alpha_final));
    all_mean_curvature_m1_noise_aniso_71(trial_cntr,1) = mean(nonzeros(apo_mean_curvature_final));
    all_total_distance_mm_noise_aniso_71(trial_cntr,1) = mean(nonzeros(apo_total_distance_final));

    % Compute fiber similarity index, euclidean distance, and corresponding
    % segment ratio
    for fiber_cntr = 1:size(apo_smoothed_fiber_all_mm_final,1)*size(apo_smoothed_fiber_all_mm_final,2)
        [row,col] = ind2sub([size(apo_smoothed_fiber_all_mm_final,1) size(apo_smoothed_fiber_all_mm_final,2)],fiber_cntr);
        fiber_noise_free = squeeze(apo_smoothed_fiber_all_mm_noise_free(row,col,:,:));
        fiber_noisy = squeeze(apo_smoothed_fiber_all_mm_final(row,col,:,:));
        if mean(fiber_noise_free,'all') ~=0 
            fiber_noise_free ( all(~fiber_noise_free ,2), : ) = []; % remove rows containing zero elements
            fiber_noisy ( all(~fiber_noisy ,2), : ) = [];
            [S, Rcs, D, Dpos] = fit_similarity(fiber_noise_free,fiber_noisy,1); % 1 mm of in-plane resolution
            fiber_similarity_aniso_71(row,col,trial_cntr) = S;
            euclid_distance_aniso_71(row,col,trial_cntr) = D;
            cor_segment_aniso_71(row,col,trial_cntr) = Rcs;
        end
    end

    % SNR 71 + smoothn
    [apo_fiber_all, ~, apo_stop_list] = ...
        fiber_track_v20(apo_ft_options, cat(4,e1map_71_smooth,famap_71), composite_mask, roi_mesh_restricted);
    [apo_smoothed_fiber_all, ~, apo_smoothed_fiber_all_mm, ~, ~, ~, ~, ~, ~] = ...
        fiber_smoother_v14a(apo_fiber_all, all_fs_options);

     % prep for filtering based on same seed location
    apo_smoothed_fiber_all_mm_final = zeros(size(apo_smoothed_fiber_all_mm));   %will hold tracts after seed-location filtering
    apo_smoothed_fiber_all_final = zeros(size(apo_smoothed_fiber_all));

    for step_cntr = 1:size(apo_smoothed_fiber_all,3)
        for dim_cntr = 1:size(apo_smoothed_fiber_all,4)
            apo_smoothed_fiber_all_final(:,:,step_cntr,dim_cntr) = ...
                squeeze(apo_smoothed_fiber_all(:,:,step_cntr,dim_cntr)).*apo_seeds_final_mask;
            apo_smoothed_fiber_all_mm_final(:,:,step_cntr,dim_cntr) = ...
                squeeze(apo_smoothed_fiber_all_mm(:,:,step_cntr,dim_cntr)).*apo_seeds_final_mask;
        end
    end

    [apo_angle_list_final, apo_distance_list_final, apo_curvature_list_final, ~, apo_n_points_final, apo_area_final] = ...
        fiber_quantifier_v20a(vxl_fq_options, apo_smoothed_fiber_all_final, roi_mesh_restricted, muscle_mask);

    % calculate raw summary statistics
    apo_total_distance_final = squeeze(max(apo_distance_list_final, [], 3));
    apo_total_distance_final(apo_total_distance_final<5)=0;

    apo_mean_alpha_final = sum(squeeze(apo_angle_list_final(:,:,:,1)), 3)./squeeze(apo_n_points_final(:,:,2));
    apo_mean_alpha_final(apo_total_distance_final<5)=0;

    apo_mean_curvature_final = sum(apo_curvature_list_final, 3)./squeeze(apo_n_points_final(:,:,3));
    apo_mean_curvature_final(apo_total_distance_final<5)=0;

    %get rid of zeros and NaNs
    apo_mean_alpha_final(isnan(apo_mean_alpha_final)) = 0;
    apo_mean_alpha_final = nonzeros(apo_mean_alpha_final);
    apo_mean_curvature_final(isnan(apo_mean_curvature_final)) = 0;
    apo_mean_curvature_final = nonzeros(apo_mean_curvature_final);
    apo_total_distance_final(isnan(apo_total_distance_final)) = 0;
    apo_total_distance_final = nonzeros(apo_total_distance_final);

    % save summary data for each time through the loop
    all_mean_alpha_deg_noise_71_smooth(trial_cntr,1) = mean(nonzeros(apo_mean_alpha_final));
    all_mean_curvature_m1_noise_71_smooth(trial_cntr,1) = mean(nonzeros(apo_mean_curvature_final));
    all_total_distance_mm_noise_71_smooth(trial_cntr,1) = mean(nonzeros(apo_total_distance_final));

    % Compute fiber similarity index, euclidean distance, and corresponding
    % segment ratio
    for fiber_cntr = 1:size(apo_smoothed_fiber_all_mm_final,1)*size(apo_smoothed_fiber_all_mm_final,2)
        [row,col] = ind2sub([size(apo_smoothed_fiber_all_mm_final,1) size(apo_smoothed_fiber_all_mm_final,2)],fiber_cntr);
        fiber_noise_free = squeeze(apo_smoothed_fiber_all_mm_noise_free(row,col,:,:));
        fiber_noisy = squeeze(apo_smoothed_fiber_all_mm_final(row,col,:,:));
        if mean(fiber_noise_free,'all') ~=0 
            fiber_noise_free ( all(~fiber_noise_free ,2), : ) = []; % remove rows containing zero elements
            fiber_noisy ( all(~fiber_noisy ,2), : ) = [];
            [S, Rcs, D, Dpos] = fit_similarity(fiber_noise_free,fiber_noisy,1); % 1 mm of in-plane resolution
            fiber_similarity_71_smooth(row,col,trial_cntr) = S;
            euclid_distance_71_smooth(row,col,trial_cntr) = D;
            cor_segment_71_smooth(row,col,trial_cntr) = Rcs;
        end
    end

    % SNR 71 - anisotropic smoothing + smoothn
    [apo_fiber_all, ~, apo_stop_list] = ...
        fiber_track_v20(apo_ft_options, cat(4,e1map_71_aniso_smooth,famap_aniso71), composite_mask, roi_mesh_restricted);
    [apo_smoothed_fiber_all, ~, apo_smoothed_fiber_all_mm, ~, ~, ~, ~, ~, ~] = ...
        fiber_smoother_v14a(apo_fiber_all, all_fs_options);

    % prep for filtering based on same seed location
    apo_smoothed_fiber_all_mm_final = zeros(size(apo_smoothed_fiber_all_mm));   %will hold tracts after seed-location filtering
    apo_smoothed_fiber_all_final = zeros(size(apo_smoothed_fiber_all));

    for step_cntr = 1:size(apo_smoothed_fiber_all,3)
        for dim_cntr = 1:size(apo_smoothed_fiber_all,4)
            apo_smoothed_fiber_all_final(:,:,step_cntr,dim_cntr) = ...
                squeeze(apo_smoothed_fiber_all(:,:,step_cntr,dim_cntr)).*apo_seeds_final_mask;
            apo_smoothed_fiber_all_mm_final(:,:,step_cntr,dim_cntr) = ...
                squeeze(apo_smoothed_fiber_all_mm(:,:,step_cntr,dim_cntr)).*apo_seeds_final_mask;
        end
    end

    [apo_angle_list_final, apo_distance_list_final, apo_curvature_list_final, ~, apo_n_points_final, apo_area_final] = ...
        fiber_quantifier_v20a(vxl_fq_options, apo_smoothed_fiber_all_final, roi_mesh_restricted, muscle_mask);

    % calculate raw summary statistics
    apo_total_distance_final = squeeze(max(apo_distance_list_final, [], 3));
    apo_total_distance_final(apo_total_distance_final<5)=0;

    apo_mean_alpha_final = sum(squeeze(apo_angle_list_final(:,:,:,1)), 3)./squeeze(apo_n_points_final(:,:,2));
    apo_mean_alpha_final(apo_total_distance_final<5)=0;

    apo_mean_curvature_final = sum(apo_curvature_list_final, 3)./squeeze(apo_n_points_final(:,:,3));
    apo_mean_curvature_final(apo_total_distance_final<5)=0;

    %get rid of zeros and NaNs
    apo_mean_alpha_final(isnan(apo_mean_alpha_final)) = 0;
    apo_mean_alpha_final = nonzeros(apo_mean_alpha_final);
    apo_mean_curvature_final(isnan(apo_mean_curvature_final)) = 0;
    apo_mean_curvature_final = nonzeros(apo_mean_curvature_final);
    apo_total_distance_final(isnan(apo_total_distance_final)) = 0;
    apo_total_distance_final = nonzeros(apo_total_distance_final);

    % save summary data for each time through the loop
    all_mean_alpha_deg_noise_aniso_71_smooth(trial_cntr,1) = mean(nonzeros(apo_mean_alpha_final));
    all_mean_curvature_m1_noise_aniso_71_smooth(trial_cntr,1) = mean(nonzeros(apo_mean_curvature_final));
    all_total_distance_mm_noise_aniso_71_smooth(trial_cntr,1) = mean(nonzeros(apo_total_distance_final));

    % Compute fiber similarity index, euclidean distance, and corresponding
    % segment ratio
    for fiber_cntr = 1:size(apo_smoothed_fiber_all_mm_final,1)*size(apo_smoothed_fiber_all_mm_final,2)
        [row,col] = ind2sub([size(apo_smoothed_fiber_all_mm_final,1) size(apo_smoothed_fiber_all_mm_final,2)],fiber_cntr);
        fiber_noise_free = squeeze(apo_smoothed_fiber_all_mm_noise_free(row,col,:,:));
        fiber_noisy = squeeze(apo_smoothed_fiber_all_mm_final(row,col,:,:));
        if mean(fiber_noise_free,'all') ~=0 
            fiber_noise_free ( all(~fiber_noise_free ,2), : ) = []; % remove rows containing zero elements
            fiber_noisy ( all(~fiber_noisy ,2), : ) = [];
            [S, Rcs, D, Dpos] = fit_similarity(fiber_noise_free,fiber_noisy,1); % 1 mm of in-plane resolution
            fiber_similarity_aniso_71_smooth(row,col,trial_cntr) = S;
            euclid_distance_aniso_71_smooth(row,col,trial_cntr) = D;
            cor_segment_aniso_71_smooth(row,col,trial_cntr) = Rcs;
        end
    end

    % Percent of fiber tracts that didn't start due to noisy data violating the
    % termination criteria
    percent_empty_fiber_tracts_24(trial_cntr,1) = sum(isnan(fiber_similarity_24(:,:,trial_cntr)),'all')/sum(apo_seeds_final_mask,'all');
    percent_empty_fiber_tracts_aniso_24(trial_cntr,1) = sum(isnan(fiber_similarity_aniso_24(:,:,trial_cntr)),'all')/sum(apo_seeds_final_mask,'all');
    percent_empty_fiber_tracts_smooth_24(trial_cntr,1) = sum(isnan(fiber_similarity_24_smooth(:,:,trial_cntr)),'all')/sum(apo_seeds_final_mask,'all');
    percent_empty_fiber_tracts_aniso_smooth_24(trial_cntr,1) = sum(isnan(fiber_similarity_aniso_24_smooth(:,:,trial_cntr)),'all')/sum(apo_seeds_final_mask,'all');
    percent_empty_fiber_tracts_71(trial_cntr,1) = sum(isnan(fiber_similarity_71(:,:,trial_cntr)),'all')/sum(apo_seeds_final_mask,'all');
    percent_empty_fiber_tracts_aniso_71(trial_cntr,1) = sum(isnan(fiber_similarity_aniso_71(:,:,trial_cntr)),'all')/sum(apo_seeds_final_mask,'all');
    percent_empty_fiber_tracts_smooth_71(trial_cntr,1) = sum(isnan(fiber_similarity_71_smooth(:,:,trial_cntr)),'all')/sum(apo_seeds_final_mask,'all');
    percent_empty_fiber_tracts_aniso_smooth_71(trial_cntr,1) = sum(isnan(fiber_similarity_aniso_71_smooth(:,:,trial_cntr)),'all')/sum(apo_seeds_final_mask,'all');

    % Mean fiber similarity of all fiber tracts -  noisy tracts that didn't
    % start have fiber similarity = 0
    fiber_similarity_24_mean(trial_cntr,1) = sum(fiber_similarity_24(:,:,trial_cntr),'all','omitnan')/sum(apo_seeds_final_mask,'all');
    fiber_similarity_24_reshape = reshape(fiber_similarity_24(:,:,trial_cntr),[numel(fiber_similarity_24(:,:,trial_cntr)) 1]);
    fiber_similarity_24_reshape(fiber_similarity_24_reshape == 0) = [];
    fiber_similarity_24_reshape(isnan(fiber_similarity_24_reshape)) = 0;
    fiber_similarity_24_std(trial_cntr,1) = std(fiber_similarity_24_reshape);
    fiber_similarity_smooth_24_mean(trial_cntr,1) = sum(fiber_similarity_24_smooth(:,:,trial_cntr),'all','omitnan')/sum(apo_seeds_final_mask,'all');
    fiber_similarity_smooth_24_reshape = reshape(fiber_similarity_24_smooth(:,:,trial_cntr),[numel(fiber_similarity_24_smooth(:,:,trial_cntr)) 1]);
    fiber_similarity_smooth_24_reshape(fiber_similarity_smooth_24_reshape == 0) = [];
    fiber_similarity_smooth_24_reshape(isnan(fiber_similarity_smooth_24_reshape)) = 0;
    fiber_similarity_smooth_24_std(trial_cntr,1) = std(fiber_similarity_smooth_24_reshape);
    fiber_similarity_aniso_24_mean(trial_cntr,1) = sum(fiber_similarity_aniso_24(:,:,trial_cntr),'all','omitnan')/sum(apo_seeds_final_mask,'all');
    fiber_similarity_aniso_24_reshape = reshape(fiber_similarity_aniso_24(:,:,trial_cntr),[numel(fiber_similarity_aniso_24(:,:,trial_cntr)) 1]);
    fiber_similarity_aniso_24_reshape(fiber_similarity_aniso_24_reshape == 0) = [];
    fiber_similarity_aniso_24_reshape(isnan(fiber_similarity_aniso_24_reshape)) = 0;
    fiber_similarity_aniso_24_std(trial_cntr,1) = std(fiber_similarity_aniso_24_reshape);
    fiber_similarity_aniso_smooth_24_mean(trial_cntr,1) = sum(fiber_similarity_aniso_24_smooth(:,:,trial_cntr),'all','omitnan')/sum(apo_seeds_final_mask,'all');
    fiber_similarity_aniso_smooth_24_reshape = reshape(fiber_similarity_aniso_24_smooth(:,:,trial_cntr),[numel(fiber_similarity_aniso_24_smooth(:,:,trial_cntr)) 1]);
    fiber_similarity_aniso_smooth_24_reshape(fiber_similarity_aniso_smooth_24_reshape == 0) = [];
    fiber_similarity_aniso_smooth_24_reshape(isnan(fiber_similarity_aniso_smooth_24_reshape)) = 0;
    fiber_similarity_aniso_smooth_24_std(trial_cntr,1) = std(fiber_similarity_aniso_smooth_24_reshape);
    fiber_similarity_71_mean(trial_cntr,1) = sum(fiber_similarity_71(:,:,trial_cntr),'all','omitnan')/sum(apo_seeds_final_mask,'all');
    fiber_similarity_71_reshape = reshape(fiber_similarity_71(:,:,trial_cntr),[numel(fiber_similarity_71(:,:,trial_cntr)) 1]);
    fiber_similarity_71_reshape(fiber_similarity_71_reshape == 0) = [];
    fiber_similarity_71_reshape(isnan(fiber_similarity_71_reshape)) = 0;
    fiber_similarity_71_std(trial_cntr,1) = std(fiber_similarity_71_reshape);
    fiber_similarity_smooth_71_mean(trial_cntr,1) = sum(fiber_similarity_71_smooth(:,:,trial_cntr),'all','omitnan')/sum(apo_seeds_final_mask,'all');
    fiber_similarity_smooth_71_reshape = reshape(fiber_similarity_71_smooth(:,:,trial_cntr),[numel(fiber_similarity_71_smooth(:,:,trial_cntr)) 1]);
    fiber_similarity_smooth_71_reshape(fiber_similarity_smooth_71_reshape == 0) = [];
    fiber_similarity_smooth_71_reshape(isnan(fiber_similarity_smooth_71_reshape)) = 0;
    fiber_similarity_smooth_71_std(trial_cntr,1) = std(fiber_similarity_smooth_71_reshape);
    fiber_similarity_aniso_71_mean(trial_cntr,1) = sum(fiber_similarity_aniso_71(:,:,trial_cntr),'all','omitnan')/sum(apo_seeds_final_mask,'all');
    fiber_similarity_aniso_71_reshape = reshape(fiber_similarity_aniso_71(:,:,trial_cntr),[numel(fiber_similarity_aniso_71(:,:,trial_cntr)) 1]);
    fiber_similarity_aniso_71_reshape(fiber_similarity_aniso_71_reshape == 0) = [];
    fiber_similarity_aniso_71_reshape(isnan(fiber_similarity_aniso_71_reshape)) = 0;
    fiber_similarity_aniso_71_std(trial_cntr,1) = std(fiber_similarity_aniso_71_reshape);
    fiber_similarity_aniso_smooth_71_mean(trial_cntr,1) = sum(fiber_similarity_aniso_71_smooth(:,:,trial_cntr),'all','omitnan')/sum(apo_seeds_final_mask,'all');
    fiber_similarity_aniso_smooth_71_reshape = reshape(fiber_similarity_aniso_71_smooth(:,:,trial_cntr),[numel(fiber_similarity_aniso_71_smooth(:,:,trial_cntr)) 1]);
    fiber_similarity_aniso_smooth_71_reshape(fiber_similarity_aniso_smooth_71_reshape == 0) = [];
    fiber_similarity_aniso_smooth_71_reshape(isnan(fiber_similarity_aniso_smooth_71_reshape)) = 0;
    fiber_similarity_aniso_smooth_71_std(trial_cntr,1) = std(fiber_similarity_aniso_smooth_71_reshape);

    % Mean euclidean distance of all fiber tracts -  noisy tracts that didn't
    % start have fiber similarity = 0
    euclid_distance_24_mean(trial_cntr,1) = sum(euclid_distance_24(:,:,trial_cntr),'all','omitnan')/sum(apo_seeds_final_mask,'all');
    euclid_distance_24_reshape = reshape(euclid_distance_24(:,:,trial_cntr),[numel(euclid_distance_24(:,:,trial_cntr)) 1]);
    euclid_distance_24_reshape(euclid_distance_24_reshape == 0) = [];
    euclid_distance_24_reshape(isnan(euclid_distance_24_reshape)) = 0;
    euclid_distance_24_std(trial_cntr,1) = std(euclid_distance_24_reshape);
    euclid_distance_smooth_24_mean(trial_cntr,1) = sum(euclid_distance_24_smooth(:,:,trial_cntr),'all','omitnan')/sum(apo_seeds_final_mask,'all');
    euclid_distance_smooth_24_reshape = reshape(euclid_distance_24_smooth(:,:,trial_cntr),[numel(euclid_distance_24_smooth(:,:,trial_cntr)) 1]);
    euclid_distance_smooth_24_reshape(euclid_distance_smooth_24_reshape == 0) = [];
    euclid_distance_smooth_24_reshape(isnan(euclid_distance_smooth_24_reshape)) = 0;
    euclid_distance_smooth_24_std(trial_cntr,1) = std(euclid_distance_smooth_24_reshape);
    euclid_distance_aniso_24_mean(trial_cntr,1) = sum(euclid_distance_aniso_24(:,:,trial_cntr),'all','omitnan')/sum(apo_seeds_final_mask,'all');
    euclid_distance_aniso_24_reshape = reshape(euclid_distance_aniso_24(:,:,trial_cntr),[numel(euclid_distance_aniso_24(:,:,trial_cntr)) 1]);
    euclid_distance_aniso_24_reshape(euclid_distance_aniso_24_reshape == 0) = [];
    euclid_distance_aniso_24_reshape(isnan(euclid_distance_aniso_24_reshape)) = 0;
    euclid_distance_aniso_24_std(trial_cntr,1) = std(euclid_distance_aniso_24_reshape);
    euclid_distance_aniso_smooth_24_mean(trial_cntr,1) = sum(euclid_distance_aniso_24_smooth(:,:,trial_cntr),'all','omitnan')/sum(apo_seeds_final_mask,'all');
    euclid_distance_aniso_smooth_24_reshape = reshape(euclid_distance_aniso_24_smooth(:,:,trial_cntr),[numel(euclid_distance_aniso_24_smooth(:,:,trial_cntr)) 1]);
    euclid_distance_aniso_smooth_24_reshape(euclid_distance_aniso_smooth_24_reshape == 0) = [];
    euclid_distance_aniso_smooth_24_reshape(isnan(euclid_distance_aniso_smooth_24_reshape)) = 0;
    euclid_distance_aniso_smooth_24_std(trial_cntr,1) = std(euclid_distance_aniso_smooth_24_reshape);
    euclid_distance_71_mean(trial_cntr,1) = sum(euclid_distance_71(:,:,trial_cntr),'all','omitnan')/sum(apo_seeds_final_mask,'all');
    euclid_distance_71_reshape = reshape(euclid_distance_71(:,:,trial_cntr),[numel(euclid_distance_71(:,:,trial_cntr)) 1]);
    euclid_distance_71_reshape(euclid_distance_71_reshape == 0) = [];
    euclid_distance_71_reshape(isnan(euclid_distance_71_reshape)) = 0;
    euclid_distance_71_std(trial_cntr,1) = std(euclid_distance_71_reshape);
    euclid_distance_smooth_71_mean(trial_cntr,1) = sum(euclid_distance_71_smooth(:,:,trial_cntr),'all','omitnan')/sum(apo_seeds_final_mask,'all');
    euclid_distance_smooth_71_reshape = reshape(euclid_distance_71_smooth(:,:,trial_cntr),[numel(euclid_distance_71_smooth(:,:,trial_cntr)) 1]);
    euclid_distance_smooth_71_reshape(euclid_distance_smooth_71_reshape == 0) = [];
    euclid_distance_smooth_71_reshape(isnan(euclid_distance_smooth_71_reshape)) = 0;
    euclid_distance_smooth_71_std(trial_cntr,1) = std(euclid_distance_smooth_71_reshape);
    euclid_distance_aniso_71_mean(trial_cntr,1) = sum(euclid_distance_aniso_71(:,:,trial_cntr),'all','omitnan')/sum(apo_seeds_final_mask,'all');
    euclid_distance_aniso_71_reshape = reshape(euclid_distance_aniso_71(:,:,trial_cntr),[numel(euclid_distance_aniso_71(:,:,trial_cntr)) 1]);
    euclid_distance_aniso_71_reshape(euclid_distance_aniso_71_reshape == 0) = [];
    euclid_distance_aniso_71_reshape(isnan(euclid_distance_aniso_71_reshape)) = 0;
    euclid_distance_aniso_71_std(trial_cntr,1) = std(euclid_distance_aniso_71_reshape);
    euclid_distance_aniso_smooth_71_mean(trial_cntr,1) = sum(euclid_distance_aniso_71_smooth(:,:,trial_cntr),'all','omitnan')/sum(apo_seeds_final_mask,'all');
    euclid_distance_aniso_smooth_71_reshape = reshape(euclid_distance_aniso_71_smooth(:,:,trial_cntr),[numel(euclid_distance_aniso_71_smooth(:,:,trial_cntr)) 1]);
    euclid_distance_aniso_smooth_71_reshape(euclid_distance_aniso_smooth_71_reshape == 0) = [];
    euclid_distance_aniso_smooth_71_reshape(isnan(euclid_distance_aniso_smooth_71_reshape)) = 0;
    euclid_distance_aniso_smooth_71_std(trial_cntr,1) = std(euclid_distance_aniso_smooth_71_reshape);
    
end

elapsed_time = toc;

%% Calculate summary data for fiber similarity metrics

percent_empty_fiber_tracts_24_mean = mean(percent_empty_fiber_tracts_24);
percent_empty_fiber_tracts_24_std = std(percent_empty_fiber_tracts_24);
percent_empty_fiber_tracts_aniso_24_mean = mean(percent_empty_fiber_tracts_aniso_24);
percent_empty_fiber_tracts_aniso_24_std = std(percent_empty_fiber_tracts_aniso_24);
percent_empty_fiber_tracts_smooth_24_mean = mean(percent_empty_fiber_tracts_smooth_24);
percent_empty_fiber_tracts_smooth_24_std = std(percent_empty_fiber_tracts_smooth_24);
percent_empty_fiber_tracts_aniso_smooth_24_mean = mean(percent_empty_fiber_tracts_aniso_smooth_24);
percent_empty_fiber_tracts_aniso_smooth_24_std = std(percent_empty_fiber_tracts_aniso_smooth_24);
percent_empty_fiber_tracts_71_mean = mean(percent_empty_fiber_tracts_71);
percent_empty_fiber_tracts_71_std = std(percent_empty_fiber_tracts_71);
percent_empty_fiber_tracts_aniso_71_mean = mean(percent_empty_fiber_tracts_aniso_71);
percent_empty_fiber_tracts_aniso_71_std = std(percent_empty_fiber_tracts_aniso_71);
percent_empty_fiber_tracts_smooth_71_mean = mean(percent_empty_fiber_tracts_smooth_71);
percent_empty_fiber_tracts_smooth_71_std = std(percent_empty_fiber_tracts_smooth_71);
percent_empty_fiber_tracts_aniso_smooth_71_mean = mean(percent_empty_fiber_tracts_aniso_smooth_71);
percent_empty_fiber_tracts_aniso_smooth_71_std = std(percent_empty_fiber_tracts_aniso_smooth_71);

fiber_similarity_24_mean_all=mean(fiber_similarity_24_mean);
fiber_similarity_24_std_all=mean(fiber_similarity_24_std);
fiber_similarity_aniso_24_mean_all=mean(fiber_similarity_aniso_24_mean);
fiber_similarity_aniso_24_std_all=mean(fiber_similarity_aniso_24_std);
fiber_similarity_smooth_24_mean_all=mean(fiber_similarity_smooth_24_mean);
fiber_similarity_smooth_24_std_all=mean(fiber_similarity_smooth_24_std);
fiber_similarity_aniso_smooth_24_mean_all=mean(fiber_similarity_aniso_smooth_24_mean);
fiber_similarity_aniso_smooth_24_std_all=mean(fiber_similarity_aniso_smooth_24_std);
fiber_similarity_71_mean_all=mean(fiber_similarity_71_mean);
fiber_similarity_71_std_all=mean(fiber_similarity_71_std);
fiber_similarity_aniso_71_mean_all=mean(fiber_similarity_aniso_71_mean);
fiber_similarity_aniso_71_std_all=mean(fiber_similarity_aniso_71_std);
fiber_similarity_smooth_71_mean_all=mean(fiber_similarity_smooth_71_mean);
fiber_similarity_smooth_71_std_all=mean(fiber_similarity_smooth_71_std);
fiber_similarity_aniso_smooth_71_mean_all=mean(fiber_similarity_aniso_smooth_71_mean);
fiber_similarity_aniso_smooth_71_std_all=mean(fiber_similarity_aniso_smooth_71_std);

euclid_distance_24_mean_all=mean(euclid_distance_24_mean);
euclid_distance_24_std_all=mean(euclid_distance_24_std);
euclid_distance_aniso_24_mean_all=mean(euclid_distance_aniso_24_mean);
euclid_distance_aniso_24_std_all=mean(euclid_distance_aniso_24_std);
euclid_distance_smooth_24_mean_all=mean(euclid_distance_smooth_24_mean);
euclid_distance_smooth_24_std_all=mean(euclid_distance_smooth_24_std);
euclid_distance_aniso_smooth_24_mean_all=mean(euclid_distance_aniso_smooth_24_mean);
euclid_distance_aniso_smooth_24_std_all=mean(euclid_distance_aniso_smooth_24_std);
euclid_distance_71_mean_all=mean(euclid_distance_71_mean);
euclid_distance_71_std_all=mean(euclid_distance_71_std);
euclid_distance_aniso_71_mean_all=mean(euclid_distance_aniso_71_mean);
euclid_distance_aniso_71_std_all=mean(euclid_distance_aniso_71_std);
euclid_distance_smooth_71_mean_all=mean(euclid_distance_smooth_71_mean);
euclid_distance_smooth_71_std_all=mean(euclid_distance_smooth_71_std);
euclid_distance_aniso_smooth_71_mean_all=mean(euclid_distance_aniso_smooth_71_mean);
euclid_distance_aniso_smooth_71_std_all=mean(euclid_distance_aniso_smooth_71_std);

%% Create table to plot data on R

SNR_Cat = repelem(["SNR_24";"SNR_71"],[4000,4000]);
Data_Cat = repelem(["Raw";"Aniso";"Smoothn";"Aniso_Smoothn"],[1000,1000,1000,1000]);
Data_Cat = [Data_Cat;Data_Cat];
Curvature_m1 = [all_mean_curvature_m1_noise_24;all_mean_curvature_m1_noise_aniso_24;...
    all_mean_curvature_m1_noise_24_smooth;all_mean_curvature_m1_noise_aniso_24_smooth;...
    all_mean_curvature_m1_noise_71;all_mean_curvature_m1_noise_aniso_71;...
    all_mean_curvature_m1_noise_71_smooth;all_mean_curvature_m1_noise_aniso_71_smooth];
Fiber_length_mm = [all_total_distance_mm_noise_24;all_total_distance_mm_noise_aniso_24;...
    all_total_distance_mm_noise_24_smooth;all_total_distance_mm_noise_aniso_24_smooth;...
    all_total_distance_mm_noise_71;all_total_distance_mm_noise_aniso_71;...
    all_total_distance_mm_noise_71_smooth;all_total_distance_mm_noise_aniso_71_smooth];
Pennation_angle_deg = [all_mean_alpha_deg_noise_24;all_mean_alpha_deg_noise_aniso_24;...
    all_mean_alpha_deg_noise_24_smooth;all_mean_alpha_deg_noise_aniso_24_smooth;...
    all_mean_alpha_deg_noise_71;all_mean_alpha_deg_noise_aniso_71;...
    all_mean_alpha_deg_noise_71_smooth;all_mean_alpha_deg_noise_aniso_71_smooth];
Fiber_similarity = [fiber_similarity_24_mean;fiber_similarity_aniso_24_mean;...
    fiber_similarity_smooth_24_mean;fiber_similarity_aniso_smooth_24_mean;
    fiber_similarity_71_mean;fiber_similarity_aniso_71_mean;...
    fiber_similarity_smooth_71_mean;fiber_similarity_aniso_smooth_71_mean];
Simulation_Table = table(SNR_Cat,Data_Cat,Fiber_length_mm,Pennation_angle_deg,...
    Curvature_m1,Fiber_similarity); 
writetable(Simulation_Table,'Sim_Results_APO3_SNR24_71.csv')

%% 
% Results for the manuscript

fiber_length_24_mean=mean(all_total_distance_mm_noise_24)
fiber_length_24_ci=quantile(all_total_distance_mm_noise_24,[0.025 0.975])

fiber_length_aniso_24_mean=mean(all_total_distance_mm_noise_aniso_24)
fiber_length_aniso_24_ci=quantile(all_total_distance_mm_noise_aniso_24,[0.025 0.975])

fiber_length_smooth_24_mean=mean(all_total_distance_mm_noise_24_smooth)
fiber_length_smooth_24_ci=quantile(all_total_distance_mm_noise_24_smooth,[0.025 0.975])

fiber_length_aniso_smooth_24_mean=mean(all_total_distance_mm_noise_aniso_24_smooth)
fiber_length_aniso_smooth_24_ci=quantile(all_total_distance_mm_noise_aniso_24_smooth,[0.025 0.975])

fiber_length_71_mean=mean(all_total_distance_mm_noise_71)
fiber_length_71_ci=quantile(all_total_distance_mm_noise_71,[0.025 0.975])

fiber_length_aniso_71_mean=mean(all_total_distance_mm_noise_aniso_71)
fiber_length_aniso_71_ci=quantile(all_total_distance_mm_noise_aniso_71,[0.025 0.975])

fiber_length_smooth_71_mean=mean(all_total_distance_mm_noise_71_smooth)
fiber_length_smooth_71_ci=quantile(all_total_distance_mm_noise_71_smooth,[0.025 0.975])

fiber_length_aniso_smooth_71_mean=mean(all_total_distance_mm_noise_aniso_71_smooth)
fiber_length_aniso_smooth_71_ci=quantile(all_total_distance_mm_noise_aniso_71_smooth,[0.025 0.975])

pen_angle_24_mean=mean(all_mean_alpha_deg_noise_24)
pen_angle_24_ci=quantile(all_mean_alpha_deg_noise_24,[0.025 0.975])

pen_angle_aniso_24_mean=mean(all_mean_alpha_deg_noise_aniso_24)
pen_angle_aniso_24_ci=quantile(all_mean_alpha_deg_noise_aniso_24,[0.025 0.975])

pen_angle_smooth_24_mean=mean(all_mean_alpha_deg_noise_24_smooth)
pen_angle_smooth_24_ci=quantile(all_mean_alpha_deg_noise_24_smooth,[0.025 0.975])

pen_angle_aniso_smooth_24_mean=mean(all_mean_alpha_deg_noise_aniso_24_smooth)
pen_angle_aniso_smooth_24_ci=quantile(all_mean_alpha_deg_noise_aniso_24_smooth,[0.025 0.975])

pen_angle_71_mean=mean(all_mean_alpha_deg_noise_71)
pen_angle_71_ci=quantile(all_mean_alpha_deg_noise_71,[0.025 0.975])

pen_angle_aniso_71_mean=mean(all_mean_alpha_deg_noise_aniso_71)
pen_angle_aniso_71_ci=quantile(all_mean_alpha_deg_noise_aniso_71,[0.025 0.975])

pen_angle_smooth_71_mean=mean(all_mean_alpha_deg_noise_71_smooth)
pen_angle_smooth_71_ci=quantile(all_mean_alpha_deg_noise_71_smooth,[0.025 0.975])

pen_angle_aniso_smooth_71_mean=mean(all_mean_alpha_deg_noise_aniso_71_smooth)
pen_angle_aniso_smooth_71_ci=quantile(all_mean_alpha_deg_noise_aniso_71_smooth,[0.025 0.975])

curvature_24_mean=mean(all_mean_curvature_m1_noise_24)
curvature_24_ci=quantile(all_mean_curvature_m1_noise_24,[0.025 0.975])

curvature_aniso_24_mean=mean(all_mean_curvature_m1_noise_aniso_24)
curvature_aniso_24_ci=quantile(all_mean_curvature_m1_noise_aniso_24,[0.025 0.975])

curvature_smooth_24_mean=mean(all_mean_curvature_m1_noise_24_smooth)
curvature_smooth_24_ci=quantile(all_mean_curvature_m1_noise_24_smooth,[0.025 0.975])

curvature_aniso_smooth_24_mean=mean(all_mean_curvature_m1_noise_aniso_24_smooth)
curvature_aniso_smooth_24_ci=quantile(all_mean_curvature_m1_noise_aniso_24_smooth,[0.025 0.975])

curvature_71_mean=mean(all_mean_curvature_m1_noise_71)
curvature_71_ci=quantile(all_mean_curvature_m1_noise_71,[0.025 0.975])

curvature_aniso_71_mean=mean(all_mean_curvature_m1_noise_aniso_71)
curvature_aniso_71_ci=quantile(all_mean_curvature_m1_noise_aniso_71,[0.025 0.975])

curvature_smooth_71_mean=mean(all_mean_curvature_m1_noise_71_smooth)
curvature_smooth_71_ci=quantile(all_mean_curvature_m1_noise_71_smooth,[0.025 0.975])

curvature_aniso_smooth_71_mean=mean(all_mean_curvature_m1_noise_aniso_71_smooth)
curvature_aniso_smooth_71_ci=quantile(all_mean_curvature_m1_noise_aniso_71_smooth,[0.025 0.975])

%% Plot 1st eigenvector maps

clim=[-0.5 0.5];

figure
t = tiledlayout(2,6);

nexttile(1,[2 2]);
imagesc(squeeze(squeeze(e1_map(:,:,21,1))).*muscle_mask(:,:,21),clim); pbaspect([1 1 1])
title('Noise-free map')

nexttile(3)
imagesc(squeeze(squeeze(e1map_24(:,:,21,1))).*muscle_mask(:,:,21),clim); pbaspect([1 1 1])
title('SNR = 24')

nexttile(4)
imagesc(squeeze(squeeze(e1map_aniso24(:,:,21,1))).*muscle_mask(:,:,21),clim); pbaspect([1 1 1])
title('SNR = 24, Aniso 4D')

nexttile(5)
imagesc(squeeze(squeeze(e1map_24_smooth(:,:,21,1))).*muscle_mask(:,:,21),clim); pbaspect([1 1 1])
title('SNR = 24, smoothn')

nexttile(6)
imagesc(squeeze(squeeze(e1map_24_aniso_smooth(:,:,21,1))).*muscle_mask(:,:,21),clim); pbaspect([1 1 1])
title('SNR = 24, Aniso 4D + smoothn')

nexttile(9)
imagesc(squeeze(squeeze(e1map_71(:,:,21,1))).*muscle_mask(:,:,21),clim); pbaspect([1 1 1])
title('SNR = 71')

nexttile(10)
imagesc(squeeze(squeeze(e1map_aniso71(:,:,21,1))).*muscle_mask(:,:,21),clim); pbaspect([1 1 1])
title('SNR = 71, Aniso 4D')

nexttile(11)
imagesc(squeeze(squeeze(e1map_71_smooth(:,:,21,1))).*muscle_mask(:,:,21),clim); pbaspect([1 1 1])
title('SNR = 71, smoothn')

nexttile(12)
imagesc(squeeze(squeeze(e1map_71_aniso_smooth(:,:,21,1))).*muscle_mask(:,:,21),clim); pbaspect([1 1 1])
title('SNR = 71, Aniso 4D + smoothn')


cb = colorbar;
cb.Layout.Tile = 'east';
title(t,'1st eigenvector - x component')

clim=[-0.5 0.5];

figure
t = tiledlayout(2,6);

nexttile(1,[2 2]);
imagesc(squeeze(squeeze(e1_map(:,:,21,2))).*muscle_mask(:,:,21),clim); pbaspect([1 1 1])
title('Noise-free map')

nexttile(3)
imagesc(squeeze(squeeze(e1map_24(:,:,21,2))).*muscle_mask(:,:,21),clim); pbaspect([1 1 1])
title('SNR = 24')

nexttile(4)
imagesc(squeeze(squeeze(e1map_aniso24(:,:,21,2))).*muscle_mask(:,:,21),clim); pbaspect([1 1 1])
title('SNR = 24, Aniso 4D')

nexttile(5)
imagesc(squeeze(squeeze(e1map_24_smooth(:,:,21,2))).*muscle_mask(:,:,21),clim); pbaspect([1 1 1])
title('SNR = 24, smoothn')

nexttile(6)
imagesc(squeeze(squeeze(e1map_24_aniso_smooth(:,:,21,2))).*muscle_mask(:,:,21),clim); pbaspect([1 1 1])
title('SNR = 24, Aniso 4D + smoothn')

nexttile(9)
imagesc(squeeze(squeeze(e1map_71(:,:,21,2))).*muscle_mask(:,:,21),clim); pbaspect([1 1 1])
title('SNR = 71')

nexttile(10)
imagesc(squeeze(squeeze(e1map_aniso71(:,:,21,2))).*muscle_mask(:,:,21),clim); pbaspect([1 1 1])
title('SNR = 71, Aniso 4D')

nexttile(11)
imagesc(squeeze(squeeze(e1map_71_smooth(:,:,21,2))).*muscle_mask(:,:,21),clim); pbaspect([1 1 1])
title('SNR = 71, smoothn')

nexttile(12)
imagesc(squeeze(squeeze(e1map_71_aniso_smooth(:,:,21,2))).*muscle_mask(:,:,21),clim); pbaspect([1 1 1])
title('SNR = 71, Aniso 4D + smoothn')


cb = colorbar;
cb.Layout.Tile = 'east';
title(t,'1st eigenvector - y component')

figure
t = tiledlayout(2,6);
clim=[0.7 1];

nexttile(1,[2 2]);
imagesc(squeeze(squeeze(e1_map(:,:,21,3))).*muscle_mask(:,:,21),clim); pbaspect([1 1 1])
title('Noise-free map')

nexttile(3)
imagesc(squeeze(squeeze(e1map_24(:,:,21,3))).*muscle_mask(:,:,21),clim); pbaspect([1 1 1])
title('SNR = 24')

nexttile(4)
imagesc(squeeze(squeeze(e1map_aniso24(:,:,21,3))).*muscle_mask(:,:,21),clim); pbaspect([1 1 1])
title('SNR = 24, Aniso 4D')

nexttile(5)
imagesc(squeeze(squeeze(e1map_24_smooth(:,:,21,3))).*muscle_mask(:,:,21),clim); pbaspect([1 1 1])
title('SNR = 24, smoothn')

nexttile(6)
imagesc(squeeze(squeeze(e1map_24_aniso_smooth(:,:,21,3))).*muscle_mask(:,:,21),clim); pbaspect([1 1 1])
title('SNR = 24, Aniso 4D + smoothn')

nexttile(9)
imagesc(squeeze(squeeze(e1map_71(:,:,21,3))).*muscle_mask(:,:,21),clim); pbaspect([1 1 1])
title('SNR = 71')

nexttile(10)
imagesc(squeeze(squeeze(e1map_aniso71(:,:,21,3))).*muscle_mask(:,:,21),clim); pbaspect([1 1 1])
title('SNR = 71, Aniso 4D')

nexttile(11)
imagesc(squeeze(squeeze(e1map_71_smooth(:,:,21,3))).*muscle_mask(:,:,21),clim); pbaspect([1 1 1])
title('SNR = 71, smoothn')

nexttile(12)
imagesc(squeeze(squeeze(e1map_71_aniso_smooth(:,:,21,3))).*muscle_mask(:,:,21),clim); pbaspect([1 1 1])
title('SNR = 71, Aniso 4D + smoothn')


cb = colorbar;
cb.Layout.Tile = 'east';
title(t,'1st eigenvector - z component')

%% Plot angle of deviation maps

clim=[0 15];

angle_noise_24(angle_noise_24 ==0) = NaN;
angle_noise_aniso_24(angle_noise_aniso_24 ==0) = NaN;
angle_noise_24_smooth(angle_noise_24_smooth ==0) = NaN;
angle_noise_aniso_24_smooth(angle_noise_aniso_24_smooth ==0) = NaN;
angle_noise_71(angle_noise_71 ==0) = NaN;
angle_noise_aniso_71(angle_noise_aniso_71 ==0) = NaN;
angle_noise_71_smooth(angle_noise_71_smooth ==0) = NaN;
angle_noise_aniso_71_smooth(angle_noise_aniso_71_smooth ==0) = NaN;

figure
t = tiledlayout(2,4);
nexttile
imagesc(squeeze(squeeze(angle_noise_24(:,:,21,end))),'AlphaData',~isnan(squeeze(squeeze(angle_noise_24(:,:,21,end)))),clim);set(gca,'color',[0.7 0.7 0.7]); pbaspect([1 1 1])
title('SNR = 24')

nexttile
imagesc(squeeze(squeeze(angle_noise_aniso_24(:,:,21,end))),'AlphaData',~isnan(squeeze(squeeze(angle_noise_24(:,:,21,end)))),clim);set(gca,'color',[0.7 0.7 0.7]); pbaspect([1 1 1])
title('SNR = 24, Aniso 4D')

nexttile
imagesc(squeeze(squeeze(angle_noise_24_smooth(:,:,21,end))),'AlphaData',~isnan(squeeze(squeeze(angle_noise_24(:,:,21,end)))),clim);set(gca,'color',[0.7 0.7 0.7]); pbaspect([1 1 1])
title('SNR = 24, smoothn')

nexttile
imagesc(squeeze(squeeze(angle_noise_aniso_24_smooth(:,:,21,end))),'AlphaData',~isnan(squeeze(squeeze(angle_noise_24(:,:,21,end)))),clim);set(gca,'color',[0.7 0.7 0.7]); pbaspect([1 1 1])
title('SNR = 24, Aniso 4D + smoothn')

nexttile
imagesc(squeeze(squeeze(angle_noise_71(:,:,21,end))),'AlphaData',~isnan(squeeze(squeeze(angle_noise_24(:,:,21,end)))),clim);set(gca,'color',[0.7 0.7 0.7]); pbaspect([1 1 1])
title('SNR = 71')

nexttile
imagesc(squeeze(squeeze(angle_noise_aniso_71(:,:,21,end))),'AlphaData',~isnan(squeeze(squeeze(angle_noise_24(:,:,21,end)))),clim);set(gca,'color',[0.7 0.7 0.7]); pbaspect([1 1 1])
title('SNR = 71, Aniso 4D')

nexttile
imagesc(squeeze(squeeze(angle_noise_71_smooth(:,:,21,end))),'AlphaData',~isnan(squeeze(squeeze(angle_noise_24(:,:,21,end)))),clim);set(gca,'color',[0.7 0.7 0.7]); pbaspect([1 1 1])
title('SNR = 71, smoothn')

nexttile
imagesc(squeeze(squeeze(angle_noise_aniso_71_smooth(:,:,21,end))),'AlphaData',~isnan(squeeze(squeeze(angle_noise_24(:,:,21,end)))),clim);set(gca,'color',[0.7 0.7 0.7]); pbaspect([1 1 1])
title('SNR = 71, Aniso 4D + smoothn')

cb = colorbar;
cb.Layout.Tile = 'east';
title(t,'Angular deviation $\theta wrt 1st eigenvector')

cd('S:\Muscle_DTI\Roberto\DTI_Muscle_Analysis\Processed_Data\Aim1E')
print(gcf,'Angular_Deviation_Images_SNR24_71.png','-dpng','-r1000'); 

figure
imagesc(squeeze(squeeze(angle_noise_24(:,:,21,end))),'AlphaData',~isnan(squeeze(squeeze(angle_noise_24(:,:,21,end)))),clim);set(gca,'color',[0.7 0.7 0.7]); pbaspect([1 1 1])
title('SNR = 24')
print(gcf,'Angular_Deviation_Images_SNR24_Raw.tif','-dtiffn','-r300'); 

figure
imagesc(squeeze(squeeze(angle_noise_aniso_24(:,:,21,end))),'AlphaData',~isnan(squeeze(squeeze(angle_noise_24(:,:,21,end)))),clim);set(gca,'color',[0.7 0.7 0.7]); pbaspect([1 1 1])
title('SNR = 24, Aniso 4D')
print(gcf,'Angular_Deviation_Images_SNR24_Aniso.png','-dpng','-r2000'); 

figure
imagesc(squeeze(squeeze(angle_noise_24_smooth(:,:,21,end))),'AlphaData',~isnan(squeeze(squeeze(angle_noise_24(:,:,21,end)))),clim);set(gca,'color',[0.7 0.7 0.7]); pbaspect([1 1 1])
title('SNR = 24, smoothn')
print(gcf,'Angular_Deviation_Images_SNR24_Smoothn.png','-dpng','-r2000'); 

figure
imagesc(squeeze(squeeze(angle_noise_aniso_24_smooth(:,:,21,end))),'AlphaData',~isnan(squeeze(squeeze(angle_noise_24(:,:,21,end)))),clim);set(gca,'color',[0.7 0.7 0.7]); pbaspect([1 1 1])
title('SNR = 24, Aniso 4D + smoothn')
print(gcf,'Angular_Deviation_Images_SNR24_Aniso_Smoothn.png','-dpng','-r2000'); 

figure
imagesc(squeeze(squeeze(angle_noise_71(:,:,21,end))),'AlphaData',~isnan(squeeze(squeeze(angle_noise_24(:,:,21,end)))),clim);set(gca,'color',[0.7 0.7 0.7]); pbaspect([1 1 1])
title('SNR = 71')
print(gcf,'Angular_Deviation_Images_SNR71_Raw.png','-dpng','-r2000'); 

figure
imagesc(squeeze(squeeze(angle_noise_aniso_71(:,:,21,end))),'AlphaData',~isnan(squeeze(squeeze(angle_noise_24(:,:,21,end)))),clim);set(gca,'color',[0.7 0.7 0.7]); pbaspect([1 1 1])
title('SNR = 71, Aniso 4D')
print(gcf,'Angular_Deviation_Images_SNR71_Aniso.png','-dpng','-r2000'); 

figure
imagesc(squeeze(squeeze(angle_noise_71_smooth(:,:,21,end))),'AlphaData',~isnan(squeeze(squeeze(angle_noise_24(:,:,21,end)))),clim);set(gca,'color',[0.7 0.7 0.7]); pbaspect([1 1 1])
title('SNR = 71, smoothn')
print(gcf,'Angular_Deviation_Images_SNR71_Smoothn.png','-dpng','-r2000'); 

figure
imagesc(squeeze(squeeze(angle_noise_aniso_71_smooth(:,:,21,end))),'AlphaData',~isnan(squeeze(squeeze(angle_noise_24(:,:,21,end)))),clim);set(gca,'color',[0.7 0.7 0.7]); pbaspect([1 1 1])
title('SNR = 71, Aniso 4D + smoothn')
print(gcf,'Angular_Deviation_Images_SNR71_Aniso_Smoothn.png','-dpng','-r2000'); 

%% 
% Results for the manuscript

angle_noise_24_iterations=squeeze(mean(angle_noise_24,[1 2 3],"omitnan"));
angle_noise_24_iterations_mean = mean(angle_noise_24_iterations)
angle_noise_24_ci=quantile(angle_noise_24_iterations,[0.025 0.975])

angle_noise_aniso_24_iterations = squeeze(mean(angle_noise_aniso_24,[1 2 3],"omitnan"));
angle_noise_aniso_24_iterations_mean = mean(angle_noise_aniso_24_iterations)
angle_noise_aniso_24_ci=quantile(angle_noise_aniso_24_iterations,[0.025 0.975])

angle_noise_24_smooth_iterations = squeeze(mean(angle_noise_24_smooth,[1 2 3],"omitnan"));
angle_noise_24_smooth_iterations_mean = mean(angle_noise_24_smooth_iterations)
angle_noise_24_smooth_ci=quantile(angle_noise_24_smooth_iterations,[0.025 0.975])

angle_noise_aniso_24_smooth_iterations = squeeze(mean(angle_noise_aniso_24_smooth,[1 2 3],"omitnan"));
angle_noise_aniso_24_smooth_iterations_mean = mean(angle_noise_aniso_24_smooth_iterations)
angle_noise_aniso_24_smooth_ci=quantile(angle_noise_aniso_24_smooth_iterations,[0.025 0.975])

angle_noise_71_iterations=squeeze(mean(angle_noise_71,[1 2 3],"omitnan"));
angle_noise_71_iterations_mean = mean(angle_noise_71_iterations)
angle_noise_71_ci=quantile(angle_noise_71_iterations,[0.025 0.975])

angle_noise_aniso_71_iterations = squeeze(mean(angle_noise_aniso_71,[1 2 3],"omitnan"));
angle_noise_aniso_71_iterations_mean = mean(angle_noise_aniso_71_iterations)
angle_noise_aniso_71_ci=quantile(angle_noise_aniso_71_iterations,[0.025 0.975])

angle_noise_71_smooth_iterations = squeeze(mean(angle_noise_71_smooth,[1 2 3],"omitnan"));
angle_noise_71_smooth_iterations_mean = mean(angle_noise_71_smooth_iterations)
angle_noise_71_smooth_ci=quantile(angle_noise_71_smooth_iterations,[0.025 0.975])

angle_noise_aniso_71_smooth_iterations = squeeze(mean(angle_noise_aniso_71_smooth,[1 2 3],"omitnan"));
angle_noise_aniso_71_smooth_iterations_mean = mean(angle_noise_aniso_71_smooth_iterations)
angle_noise_aniso_71_smooth_ci=quantile(angle_noise_aniso_71_smooth_iterations,[0.025 0.975])