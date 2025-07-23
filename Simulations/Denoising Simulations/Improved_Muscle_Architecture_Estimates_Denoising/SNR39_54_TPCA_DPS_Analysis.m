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
all_snr39 = zeros(num_trials, 1);
all_snr54 = zeros(num_trials, 1);


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

angle_noise_39 = zeros([size(noise_free_img,1:3) num_trials]);
angle_noise_tpca_39 = zeros([size(noise_free_img,1:3) num_trials]);
angle_noise_54 = zeros([size(noise_free_img,1:3) num_trials]);
angle_noise_tpca_54 = zeros([size(noise_free_img,1:3) num_trials]);
angle_noise_tpca_39_smooth = zeros([size(noise_free_img,1:3) num_trials]);
angle_noise_tpca_54_smooth = zeros([size(noise_free_img,1:3) num_trials]);
fiber_similarity_39 = zeros([size(apo_smoothed_fiber_all_mm_noise_free,1:2) num_trials]);
fiber_similarity_tpca_39 = zeros([size(apo_smoothed_fiber_all_mm_noise_free,1:2) num_trials]);
fiber_similarity_54 = zeros([size(apo_smoothed_fiber_all_mm_noise_free,1:2) num_trials]);
fiber_similarity_tpca_54 = zeros([size(apo_smoothed_fiber_all_mm_noise_free,1:2) num_trials]);
fiber_similarity_tpca_39_smooth = zeros([size(apo_smoothed_fiber_all_mm_noise_free,1:2) num_trials]);
fiber_similarity_tpca_54_smooth = zeros([size(apo_smoothed_fiber_all_mm_noise_free,1:2) num_trials]);
euclid_distance_39 = zeros([size(apo_smoothed_fiber_all_mm_noise_free,1:2) num_trials]);
euclid_distance_tpca_39 = zeros([size(apo_smoothed_fiber_all_mm_noise_free,1:2) num_trials]);
euclid_distance_54 = zeros([size(apo_smoothed_fiber_all_mm_noise_free,1:2) num_trials]);
euclid_distance_tpca_54 = zeros([size(apo_smoothed_fiber_all_mm_noise_free,1:2) num_trials]);
euclid_distance_tpca_39_smooth = zeros([size(apo_smoothed_fiber_all_mm_noise_free,1:2) num_trials]);
euclid_distance_tpca_54_smooth = zeros([size(apo_smoothed_fiber_all_mm_noise_free,1:2) num_trials]);
cor_segment_39 = zeros([size(apo_smoothed_fiber_all_mm_noise_free,1:2) num_trials]);
cor_segment_tpca_39 = zeros([size(apo_smoothed_fiber_all_mm_noise_free,1:2) num_trials]);
cor_segment_54 = zeros([size(apo_smoothed_fiber_all_mm_noise_free,1:2) num_trials]);
cor_segment_tpca_54 = zeros([size(apo_smoothed_fiber_all_mm_noise_free,1:2) num_trials]);
cor_segment_tpca_39_smooth = zeros([size(apo_smoothed_fiber_all_mm_noise_free,1:2) num_trials]);
cor_segment_tpca_54_smooth = zeros([size(apo_smoothed_fiber_all_mm_noise_free,1:2) num_trials]);
all_mean_alpha_deg_noise_39 = zeros([num_trials 1]);
all_mean_curvature_m1_noise_39 = zeros([num_trials 1]);
all_total_distance_mm_noise_39 = zeros([num_trials 1]);
all_mean_alpha_deg_noise_54 = zeros([num_trials 1]);
all_mean_curvature_m1_noise_54 = zeros([num_trials 1]);
all_total_distance_mm_noise_54 = zeros([num_trials 1]);
all_mean_alpha_deg_noise_tpca_39 = zeros([num_trials 1]);
all_mean_curvature_m1_noise_tpca_39 = zeros([num_trials 1]);
all_total_distance_mm_noise_tpca_39 = zeros([num_trials 1]);
all_mean_alpha_deg_noise_tpca_54 = zeros([num_trials 1]);
all_mean_curvature_m1_noise_tpca_54 = zeros([num_trials 1]);
all_total_distance_mm_noise_tpca_54 = zeros([num_trials 1]);
all_mean_alpha_deg_noise_tpca_39_smooth = zeros([num_trials 1]);
all_mean_curvature_m1_noise_tpca_39_smooth = zeros([num_trials 1]);
all_total_distance_mm_noise_tpca_39_smooth = zeros([num_trials 1]);
all_mean_alpha_deg_noise_tpca_54_smooth = zeros([num_trials 1]);
all_mean_curvature_m1_noise_tpca_54_smooth = zeros([num_trials 1]);
all_total_distance_mm_noise_tpca_54_smooth = zeros([num_trials 1]);
fiber_similarity_39_mean = zeros([num_trials 1]);
fiber_similarity_tpca_39_mean = zeros([num_trials 1]);
fiber_similarity_tpca_smooth_39_mean = zeros([num_trials 1]);
fiber_similarity_54_mean = zeros([num_trials 1]);
fiber_similarity_tpca_54_mean = zeros([num_trials 1]);
fiber_similarity_tpca_smooth_54_mean = zeros([num_trials 1]);
euclid_distance_39_mean = zeros([num_trials 1]);
euclid_distance_tpca_39_mean = zeros([num_trials 1]);
euclid_distance_tpca_smooth_39_mean = zeros([num_trials 1]);
euclid_distance_54_mean = zeros([num_trials 1]);
euclid_distance_tpca_54_mean = zeros([num_trials 1]);
euclid_distance_tpca_smooth_54_mean = zeros([num_trials 1]);
fiber_similarity_39_std = zeros([num_trials 1]);
fiber_similarity_tpca_39_std = zeros([num_trials 1]);
fiber_similarity_tpca_smooth_39_std = zeros([num_trials 1]);
fiber_similarity_54_std = zeros([num_trials 1]);
fiber_similarity_tpca_54_std = zeros([num_trials 1]);
fiber_similarity_tpca_smooth_54_std = zeros([num_trials 1]);
euclid_distance_39_std = zeros([num_trials 1]);
euclid_distance_tpca_39_std = zeros([num_trials 1]);
euclid_distance_tpca_smooth_39_std = zeros([num_trials 1]);
euclid_distance_54_std = zeros([num_trials 1]);
euclid_distance_tpca_54_std = zeros([num_trials 1]);
euclid_distance_tpca_smooth_54_std = zeros([num_trials 1]);
percent_empty_fiber_tracts_39 = zeros([num_trials 1]);
percent_empty_fiber_tracts_tpca_39 = zeros([num_trials 1]);
percent_empty_fiber_tracts_tpca_smooth_39 = zeros([num_trials 1]);
percent_empty_fiber_tracts_54 = zeros([num_trials 1]);
percent_empty_fiber_tracts_tpca_54 = zeros([num_trials 1]);
percent_empty_fiber_tracts_tpca_smooth_54 = zeros([num_trials 1]);

tic
for trial_cntr=1:num_trials

    % form a noisy image
    noise=randn(size(noise_free_img));
    noisy_img_39 = abs(noise_free_img + noise/39);
    noisy_img_54 = abs(noise_free_img + noise/54);
    loop_b0_snr39 = squeeze(noisy_img_39(:,:, 25,1));
    loop_b0_snr39 = loop_b0_snr39(muscle_mask(:,:, 25)==1);
    all_snr39(trial_cntr) = 1/std(loop_b0_snr39);
    loop_b0_snr54 = squeeze(noisy_img_54(:,:, 25,1));
    loop_b0_snr54 = loop_b0_snr54(muscle_mask(:,:, 25)==1);
    all_snr54(trial_cntr) = 1/std(loop_b0_snr54);
    alpha_noise39_loop = zeros(size(noise_free_img,1:3));
    alpha_noise39_loop_tpca = zeros(size(noise_free_img,1:3));
    alpha_noise54_loop = zeros(size(noise_free_img,1:3));
    alpha_noise54_loop_tpca = zeros(size(noise_free_img,1:3));
    e1map_39 = zeros([size(noise_free_img,1:3) 3]);
    e1map_54 = zeros([size(noise_free_img,1:3) 3]);
    e1map_tpca39 = zeros([size(noise_free_img,1:3) 3]);
    e1map_tpca54 = zeros([size(noise_free_img,1:3) 3]);
    famap_39 = zeros(size(noise_free_img,1:3));
    famap_54 = zeros(size(noise_free_img,1:3));
    famap_tpca39 = zeros(size(noise_free_img,1:3));
    famap_tpca54 = zeros(size(noise_free_img,1:3));
    alpha_noise39_loop_tpca_smooth = zeros(size(noise_free_img,1:3));
    alpha_noise54_loop_tpca_smooth = zeros(size(noise_free_img,1:3));

    % Denoise data - Use tpca algorithm
    %noisy_img_tpca_39 = TPCA_denoising(noisy_img_39,ones(50,50,40),[5 5 5],'fast',1,ones(50,50,40)/39);
    [noisy_img_tpca_39,n_comps39] = TPCA_denoising(noisy_img_39,ones(50,50,40),[9 9 1],ones(50,50,40)/(39*39));
    %noisy_img_tpca_54 = TPCA_denoising(noisy_img_54,ones(50,50,40),[5 5 5],'fast',1,ones(50,50,40)/54);
    [noisy_img_tpca_54,n_comps54] = TPCA_denoising(noisy_img_54,ones(50,50,40),[9 9 1],ones(50,50,40)/(54*54));

    % get noisy diffusion tensors
    for r_img=1:50
        for c_img=1:50
            for s_img=1:40
                if muscle_apo_mask(r_img,c_img,s_img)>0

                    signal_v_39 = squeeze(noisy_img_39(r_img, c_img, s_img, :));
                    signal_v_54 = squeeze(noisy_img_54(r_img, c_img, s_img, :));
                    signal_v_39tpca = abs(squeeze(noisy_img_tpca_39(r_img, c_img, s_img, :))); %MP-PCA outputs negative signals in the aponeurosis, SNR ~ 2.4
                    signal_v_54tpca = abs(squeeze(noisy_img_tpca_54(r_img, c_img, s_img, :)));

                    D_39 = signal2tensor2(signal_v_39, dir_m, 475);
                    [E, L] = svd(D_39);

                    L_v = diag(L);
                    fa = sqrt(3/2) * sqrt(((L_v(1) - md)^2 + (L_v(2) - md)^2 + (L_v(3) - md)^2) / ...
                        (L_v(1)^2 + L_v(2)^2 + L_v(3)^2));

                    loop_E1 = E(:,L_v==max(L_v));
                    if loop_E1(3)<0
                        loop_E1 = -loop_E1;
                    end
                    
                    e1map_39(r_img,c_img,s_img,:) = loop_E1;
                    famap_39(r_img,c_img,s_img,:) = fa;
                    alpha_noise39_loop(r_img,c_img,s_img) = dot(loop_E1,squeeze(e1_map(r_img,c_img,s_img,:)));
                    
                    D_39tpca = signal2tensor2(signal_v_39tpca, dir_m, 475);
                    [E, L] = svd(D_39tpca);

                    L_v = diag(L);
                    fa = sqrt(3/2) * sqrt(((L_v(1) - md)^2 + (L_v(2) - md)^2 + (L_v(3) - md)^2) / ...
                        (L_v(1)^2 + L_v(2)^2 + L_v(3)^2));

                    loop_E1 = E(:,L_v==max(L_v));
                    if loop_E1(3)<0
                        loop_E1 = -loop_E1;
                    end
                    
                    e1map_tpca39(r_img,c_img,s_img,:) = loop_E1;
                    famap_tpca39(r_img,c_img,s_img,:) = fa;
                    alpha_noise39_loop_tpca(r_img,c_img,s_img) = dot(loop_E1,squeeze(e1_map(r_img,c_img,s_img,:)));

                    D_54 = signal2tensor2(signal_v_54, dir_m, 475);
                    [E, L] = svd(D_54);

                    L_v = diag(L);
                    fa = sqrt(3/2) * sqrt(((L_v(1) - md)^2 + (L_v(2) - md)^2 + (L_v(3) - md)^2) / ...
                        (L_v(1)^2 + L_v(2)^2 + L_v(3)^2));

                    loop_E1 = E(:,L_v==max(L_v));
                    if loop_E1(3)<0
                        loop_E1 = -loop_E1;
                    end
                    
                    e1map_54(r_img,c_img,s_img,:) = loop_E1;
                    famap_54(r_img,c_img,s_img,:) = fa;
                    alpha_noise54_loop(r_img,c_img,s_img) = dot(loop_E1,squeeze(e1_map(r_img,c_img,s_img,:)));
                    
                    D_54tpca = signal2tensor2(signal_v_54tpca, dir_m, 475);
                    [E, L] = svd(D_54tpca);

                    L_v = diag(L);
                    fa = sqrt(3/2) * sqrt(((L_v(1) - md)^2 + (L_v(2) - md)^2 + (L_v(3) - md)^2) / ...
                        (L_v(1)^2 + L_v(2)^2 + L_v(3)^2));

                    loop_E1 = E(:,L_v==max(L_v));
                    if loop_E1(3)<0
                        loop_E1 = -loop_E1;
                    end
                    
                    e1map_tpca54(r_img,c_img,s_img,:) = loop_E1;
                    famap_tpca54(r_img,c_img,s_img,:) = fa;
                    alpha_noise54_loop_tpca(r_img,c_img,s_img) = dot(loop_E1,squeeze(e1_map(r_img,c_img,s_img,:)));

                end
            end
        end
    end

    angle_noise_39(:,:,:,trial_cntr) = acosd(alpha_noise39_loop).*muscle_mask;
    angle_noise_tpca_39(:,:,:,trial_cntr) = acosd(alpha_noise39_loop_tpca).*muscle_mask;
    angle_noise_54(:,:,:,trial_cntr) = acosd(alpha_noise54_loop).*muscle_mask;
    angle_noise_tpca_54(:,:,:,trial_cntr) = acosd(alpha_noise54_loop_tpca).*muscle_mask;

    % Smooth eigenvector field 

    e1map_39_tpca_smooth=smoothn({e1map_tpca39(:,:,:,1),e1map_tpca39(:,:,:,2),e1map_tpca39(:,:,:,3)});
    e1map_54_tpca_smooth=smoothn({e1map_tpca54(:,:,:,1),e1map_tpca54(:,:,:,2),e1map_tpca54(:,:,:,3)});

    e1map_39_tpca_smooth=cat(4,e1map_39_tpca_smooth{1},e1map_39_tpca_smooth{2},e1map_39_tpca_smooth{3});
    e1map_54_tpca_smooth=cat(4,e1map_54_tpca_smooth{1},e1map_54_tpca_smooth{2},e1map_54_tpca_smooth{3});

    for r_img=1:50
        for c_img=1:50
            for s_img=1:40
                if muscle_apo_mask(r_img,c_img,s_img)>0

                    E1_39_tpca_smooth = squeeze(e1map_39_tpca_smooth(r_img,c_img,s_img,:));
                    E1_54_tpca_smooth = squeeze(e1map_54_tpca_smooth(r_img,c_img,s_img,:));

                    E1_39_tpca_smooth = E1_39_tpca_smooth/norm(E1_39_tpca_smooth);
                    E1_54_tpca_smooth = E1_54_tpca_smooth/norm(E1_54_tpca_smooth);

                    e1map_39_tpca_smooth(r_img,c_img,s_img,:)=E1_39_tpca_smooth;
                    e1map_54_tpca_smooth(r_img,c_img,s_img,:)=E1_54_tpca_smooth;

                    alpha_noise39_loop_tpca_smooth(r_img,c_img,s_img) = dot(E1_39_tpca_smooth,squeeze(e1_map(r_img,c_img,s_img,:)));
                    alpha_noise54_loop_tpca_smooth(r_img,c_img,s_img) = dot(E1_54_tpca_smooth,squeeze(e1_map(r_img,c_img,s_img,:)));

                end
            end
        end
    end

    angle_noise_tpca_39_smooth(:,:,:,trial_cntr) = acosd(alpha_noise39_loop_tpca_smooth).*muscle_mask;
    angle_noise_tpca_54_smooth(:,:,:,trial_cntr) = acosd(alpha_noise54_loop_tpca_smooth).*muscle_mask;

        %%%%%%%%%%%%%%%% fiber tracking - aponeurosis seeding %%%%%%%%%%%%%%%%%%
    % SNR 39
    [apo_fiber_all, ~, apo_stop_list] = ...
        fiber_track_v20(apo_ft_options, cat(4,e1map_39,famap_39), composite_mask, roi_mesh_restricted);
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
    all_mean_alpha_deg_noise_39(trial_cntr,1) = mean(nonzeros(apo_mean_alpha_final));
    all_mean_curvature_m1_noise_39(trial_cntr,1) = mean(nonzeros(apo_mean_curvature_final));
    all_total_distance_mm_noise_39(trial_cntr,1) = mean(nonzeros(apo_total_distance_final));

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
            fiber_similarity_39(row,col,trial_cntr) = S;
            euclid_distance_39(row,col,trial_cntr) = D;
            cor_segment_39(row,col,trial_cntr) = Rcs;
        end
    end

    % SNR 39 - TPCA smoothing
    [apo_fiber_all, ~, apo_stop_list] = ...
        fiber_track_v20(apo_ft_options, cat(4,e1map_tpca39,famap_tpca39), composite_mask, roi_mesh_restricted);
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
    all_mean_alpha_deg_noise_tpca_39(trial_cntr,1) = mean(nonzeros(apo_mean_alpha_final));
    all_mean_curvature_m1_noise_tpca_39(trial_cntr,1) = mean(nonzeros(apo_mean_curvature_final));
    all_total_distance_mm_noise_tpca_39(trial_cntr,1) = mean(nonzeros(apo_total_distance_final));

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
            fiber_similarity_tpca_39(row,col,trial_cntr) = S;
            euclid_distance_tpca_39(row,col,trial_cntr) = D;
            cor_segment_tpca_39(row,col,trial_cntr) = Rcs;
        end
    end

    % SNR 39 - TPCA smoothing + smoothn
    [apo_fiber_all, ~, apo_stop_list] = ...
        fiber_track_v20(apo_ft_options, cat(4,e1map_39_tpca_smooth,famap_tpca39), composite_mask, roi_mesh_restricted);
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
    all_mean_alpha_deg_noise_tpca_39_smooth(trial_cntr,1) = mean(nonzeros(apo_mean_alpha_final));
    all_mean_curvature_m1_noise_tpca_39_smooth(trial_cntr,1) = mean(nonzeros(apo_mean_curvature_final));
    all_total_distance_mm_noise_tpca_39_smooth(trial_cntr,1) = mean(nonzeros(apo_total_distance_final));

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
            fiber_similarity_tpca_39_smooth(row,col,trial_cntr) = S;
            euclid_distance_tpca_39_smooth(row,col,trial_cntr) = D;
            cor_segment_tpca_39_smooth(row,col,trial_cntr) = Rcs;
        end
    end

    % SNR 54
    [apo_fiber_all, ~, apo_stop_list] = ...
        fiber_track_v20(apo_ft_options, cat(4,e1map_54,famap_54), composite_mask, roi_mesh_restricted);
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
    all_mean_alpha_deg_noise_54(trial_cntr,1) = mean(nonzeros(apo_mean_alpha_final));
    all_mean_curvature_m1_noise_54(trial_cntr,1) = mean(nonzeros(apo_mean_curvature_final));
    all_total_distance_mm_noise_54(trial_cntr,1) = mean(nonzeros(apo_total_distance_final));

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
            fiber_similarity_54(row,col,trial_cntr) = S;
            euclid_distance_54(row,col,trial_cntr) = D;
            cor_segment_54(row,col,trial_cntr) = Rcs;
        end
    end

    % SNR 54 - TPCA smoothing
    [apo_fiber_all, ~, apo_stop_list] = ...
        fiber_track_v20(apo_ft_options, cat(4,e1map_tpca54,famap_tpca54), composite_mask, roi_mesh_restricted);
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
    all_mean_alpha_deg_noise_tpca_54(trial_cntr,1) = mean(nonzeros(apo_mean_alpha_final));
    all_mean_curvature_m1_noise_tpca_54(trial_cntr,1) = mean(nonzeros(apo_mean_curvature_final));
    all_total_distance_mm_noise_tpca_54(trial_cntr,1) = mean(nonzeros(apo_total_distance_final));

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
            fiber_similarity_tpca_54(row,col,trial_cntr) = S;
            euclid_distance_tpca_54(row,col,trial_cntr) = D;
            cor_segment_tpca_54(row,col,trial_cntr) = Rcs;
        end
    end

    % SNR 54 - TPCA smoothing + smoothn
    [apo_fiber_all, ~, apo_stop_list] = ...
        fiber_track_v20(apo_ft_options, cat(4,e1map_54_tpca_smooth,famap_tpca54), composite_mask, roi_mesh_restricted);
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
    all_mean_alpha_deg_noise_tpca_54_smooth(trial_cntr,1) = mean(nonzeros(apo_mean_alpha_final));
    all_mean_curvature_m1_noise_tpca_54_smooth(trial_cntr,1) = mean(nonzeros(apo_mean_curvature_final));
    all_total_distance_mm_noise_tpca_54_smooth(trial_cntr,1) = mean(nonzeros(apo_total_distance_final));

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
            fiber_similarity_tpca_54_smooth(row,col,trial_cntr) = S;
            euclid_distance_tpca_54_smooth(row,col,trial_cntr) = D;
            cor_segment_tpca_54_smooth(row,col,trial_cntr) = Rcs;
        end
    end

    % Percent of fiber tracts that didn't start due to noisy data violating the
    % termination criteria
    percent_empty_fiber_tracts_39(trial_cntr,1) = sum(isnan(fiber_similarity_39(:,:,trial_cntr)),'all')/sum(apo_seeds_final_mask,'all');
    percent_empty_fiber_tracts_tpca_39(trial_cntr,1) = sum(isnan(fiber_similarity_tpca_39(:,:,trial_cntr)),'all')/sum(apo_seeds_final_mask,'all');
    percent_empty_fiber_tracts_tpca_smooth_39(trial_cntr,1) = sum(isnan(fiber_similarity_tpca_39_smooth(:,:,trial_cntr)),'all')/sum(apo_seeds_final_mask,'all');
    percent_empty_fiber_tracts_54(trial_cntr,1) = sum(isnan(fiber_similarity_54(:,:,trial_cntr)),'all')/sum(apo_seeds_final_mask,'all');
    percent_empty_fiber_tracts_tpca_54(trial_cntr,1) = sum(isnan(fiber_similarity_tpca_54(:,:,trial_cntr)),'all')/sum(apo_seeds_final_mask,'all');
    percent_empty_fiber_tracts_tpca_smooth_54(trial_cntr,1) = sum(isnan(fiber_similarity_tpca_54_smooth(:,:,trial_cntr)),'all')/sum(apo_seeds_final_mask,'all');

    % Mean fiber similarity of all fiber tracts -  noisy tracts that didn't
    % start have fiber similarity = 0
    fiber_similarity_39_mean(trial_cntr,1) = sum(fiber_similarity_39(:,:,trial_cntr),'all','omitnan')/sum(apo_seeds_final_mask,'all');
    fiber_similarity_39_reshape = reshape(fiber_similarity_39(:,:,trial_cntr),[numel(fiber_similarity_39(:,:,trial_cntr)) 1]);
    fiber_similarity_39_reshape(fiber_similarity_39_reshape == 0) = [];
    fiber_similarity_39_reshape(isnan(fiber_similarity_39_reshape)) = 0;
    fiber_similarity_39_std(trial_cntr,1) = std(fiber_similarity_39_reshape);
    fiber_similarity_tpca_39_mean(trial_cntr,1) = sum(fiber_similarity_tpca_39(:,:,trial_cntr),'all','omitnan')/sum(apo_seeds_final_mask,'all');
    fiber_similarity_tpca_39_reshape = reshape(fiber_similarity_tpca_39(:,:,trial_cntr),[numel(fiber_similarity_tpca_39(:,:,trial_cntr)) 1]);
    fiber_similarity_tpca_39_reshape(fiber_similarity_tpca_39_reshape == 0) = [];
    fiber_similarity_tpca_39_reshape(isnan(fiber_similarity_tpca_39_reshape)) = 0;
    fiber_similarity_tpca_39_std(trial_cntr,1) = std(fiber_similarity_tpca_39_reshape);
    fiber_similarity_tpca_smooth_39_mean(trial_cntr,1) = sum(fiber_similarity_tpca_39_smooth(:,:,trial_cntr),'all','omitnan')/sum(apo_seeds_final_mask,'all');
    fiber_similarity_tpca_smooth_39_reshape = reshape(fiber_similarity_tpca_39_smooth(:,:,trial_cntr),[numel(fiber_similarity_tpca_39_smooth(:,:,trial_cntr)) 1]);
    fiber_similarity_tpca_smooth_39_reshape(fiber_similarity_tpca_smooth_39_reshape == 0) = [];
    fiber_similarity_tpca_smooth_39_reshape(isnan(fiber_similarity_tpca_smooth_39_reshape)) = 0;
    fiber_similarity_tpca_smooth_39_std(trial_cntr,1) = std(fiber_similarity_tpca_smooth_39_reshape);
    fiber_similarity_54_mean(trial_cntr,1) = sum(fiber_similarity_54(:,:,trial_cntr),'all','omitnan')/sum(apo_seeds_final_mask,'all');
    fiber_similarity_54_reshape = reshape(fiber_similarity_54(:,:,trial_cntr),[numel(fiber_similarity_54(:,:,trial_cntr)) 1]);
    fiber_similarity_54_reshape(fiber_similarity_54_reshape == 0) = [];
    fiber_similarity_54_reshape(isnan(fiber_similarity_54_reshape)) = 0;
    fiber_similarity_54_std(trial_cntr,1) = std(fiber_similarity_54_reshape);
    fiber_similarity_tpca_54_mean(trial_cntr,1) = sum(fiber_similarity_tpca_54(:,:,trial_cntr),'all','omitnan')/sum(apo_seeds_final_mask,'all');
    fiber_similarity_tpca_54_reshape = reshape(fiber_similarity_tpca_54(:,:,trial_cntr),[numel(fiber_similarity_tpca_54(:,:,trial_cntr)) 1]);
    fiber_similarity_tpca_54_reshape(fiber_similarity_tpca_54_reshape == 0) = [];
    fiber_similarity_tpca_54_reshape(isnan(fiber_similarity_tpca_54_reshape)) = 0;
    fiber_similarity_tpca_54_std(trial_cntr,1) = std(fiber_similarity_tpca_54_reshape);
    fiber_similarity_tpca_smooth_54_mean(trial_cntr,1) = sum(fiber_similarity_tpca_54_smooth(:,:,trial_cntr),'all','omitnan')/sum(apo_seeds_final_mask,'all');
    fiber_similarity_tpca_smooth_54_reshape = reshape(fiber_similarity_tpca_54_smooth(:,:,trial_cntr),[numel(fiber_similarity_tpca_54_smooth(:,:,trial_cntr)) 1]);
    fiber_similarity_tpca_smooth_54_reshape(fiber_similarity_tpca_smooth_54_reshape == 0) = [];
    fiber_similarity_tpca_smooth_54_reshape(isnan(fiber_similarity_tpca_smooth_54_reshape)) = 0;
    fiber_similarity_tpca_smooth_54_std(trial_cntr,1) = std(fiber_similarity_tpca_smooth_54_reshape);

    % Mean euclidean distance of all fiber tracts -  noisy tracts that didn't
    % start have fiber similarity = 0
    euclid_distance_39_mean(trial_cntr,1) = sum(euclid_distance_39(:,:,trial_cntr),'all','omitnan')/sum(apo_seeds_final_mask,'all');
    euclid_distance_39_reshape = reshape(euclid_distance_39(:,:,trial_cntr),[numel(euclid_distance_39(:,:,trial_cntr)) 1]);
    euclid_distance_39_reshape(euclid_distance_39_reshape == 0) = [];
    euclid_distance_39_reshape(isnan(euclid_distance_39_reshape)) = 0;
    euclid_distance_39_std(trial_cntr,1) = std(euclid_distance_39_reshape);
    euclid_distance_tpca_39_mean(trial_cntr,1) = sum(euclid_distance_tpca_39(:,:,trial_cntr),'all','omitnan')/sum(apo_seeds_final_mask,'all');
    euclid_distance_tpca_39_reshape = reshape(euclid_distance_tpca_39(:,:,trial_cntr),[numel(euclid_distance_tpca_39(:,:,trial_cntr)) 1]);
    euclid_distance_tpca_39_reshape(euclid_distance_tpca_39_reshape == 0) = [];
    euclid_distance_tpca_39_reshape(isnan(euclid_distance_tpca_39_reshape)) = 0;
    euclid_distance_tpca_39_std(trial_cntr,1) = std(euclid_distance_tpca_39_reshape);
    euclid_distance_tpca_smooth_39_mean(trial_cntr,1) = sum(euclid_distance_tpca_39_smooth(:,:,trial_cntr),'all','omitnan')/sum(apo_seeds_final_mask,'all');
    euclid_distance_tpca_smooth_39_reshape = reshape(euclid_distance_tpca_39_smooth(:,:,trial_cntr),[numel(euclid_distance_tpca_39_smooth(:,:,trial_cntr)) 1]);
    euclid_distance_tpca_smooth_39_reshape(euclid_distance_tpca_smooth_39_reshape == 0) = [];
    euclid_distance_tpca_smooth_39_reshape(isnan(euclid_distance_tpca_smooth_39_reshape)) = 0;
    euclid_distance_tpca_smooth_39_std(trial_cntr,1) = std(euclid_distance_tpca_smooth_39_reshape);
    euclid_distance_54_mean(trial_cntr,1) = sum(euclid_distance_54(:,:,trial_cntr),'all','omitnan')/sum(apo_seeds_final_mask,'all');
    euclid_distance_54_reshape = reshape(euclid_distance_54(:,:,trial_cntr),[numel(euclid_distance_54(:,:,trial_cntr)) 1]);
    euclid_distance_54_reshape(euclid_distance_54_reshape == 0) = [];
    euclid_distance_54_reshape(isnan(euclid_distance_54_reshape)) = 0;
    euclid_distance_54_std(trial_cntr,1) = std(euclid_distance_54_reshape);
    euclid_distance_tpca_54_mean(trial_cntr,1) = sum(euclid_distance_tpca_54(:,:,trial_cntr),'all','omitnan')/sum(apo_seeds_final_mask,'all');
    euclid_distance_tpca_54_reshape = reshape(euclid_distance_tpca_54(:,:,trial_cntr),[numel(euclid_distance_tpca_54(:,:,trial_cntr)) 1]);
    euclid_distance_tpca_54_reshape(euclid_distance_tpca_54_reshape == 0) = [];
    euclid_distance_tpca_54_reshape(isnan(euclid_distance_tpca_54_reshape)) = 0;
    euclid_distance_tpca_54_std(trial_cntr,1) = std(euclid_distance_tpca_54_reshape);
    euclid_distance_tpca_smooth_54_mean(trial_cntr,1) = sum(euclid_distance_tpca_54_smooth(:,:,trial_cntr),'all','omitnan')/sum(apo_seeds_final_mask,'all');
    euclid_distance_tpca_smooth_54_reshape = reshape(euclid_distance_tpca_54_smooth(:,:,trial_cntr),[numel(euclid_distance_tpca_54_smooth(:,:,trial_cntr)) 1]);
    euclid_distance_tpca_smooth_54_reshape(euclid_distance_tpca_smooth_54_reshape == 0) = [];
    euclid_distance_tpca_smooth_54_reshape(isnan(euclid_distance_tpca_smooth_54_reshape)) = 0;
    euclid_distance_tpca_smooth_54_std(trial_cntr,1) = std(euclid_distance_tpca_smooth_54_reshape);
    
end

elapsed_time = toc;

%% Calculate summary data for fiber similarity metrics

percent_empty_fiber_tracts_39_mean = mean(percent_empty_fiber_tracts_39);
percent_empty_fiber_tracts_39_std = std(percent_empty_fiber_tracts_39);
percent_empty_fiber_tracts_tpca_39_mean = mean(percent_empty_fiber_tracts_tpca_39);
percent_empty_fiber_tracts_tpca_39_std = std(percent_empty_fiber_tracts_tpca_39);
percent_empty_fiber_tracts_tpca_smooth_39_mean = mean(percent_empty_fiber_tracts_tpca_smooth_39);
percent_empty_fiber_tracts_tpca_smooth_39_std = std(percent_empty_fiber_tracts_tpca_smooth_39);
percent_empty_fiber_tracts_54_mean = mean(percent_empty_fiber_tracts_54);
percent_empty_fiber_tracts_54_std = std(percent_empty_fiber_tracts_54);
percent_empty_fiber_tracts_tpca_54_mean = mean(percent_empty_fiber_tracts_tpca_54);
percent_empty_fiber_tracts_tpca_54_std = std(percent_empty_fiber_tracts_tpca_54);
percent_empty_fiber_tracts_tpca_smooth_54_mean = mean(percent_empty_fiber_tracts_tpca_smooth_54);
percent_empty_fiber_tracts_tpca_smooth_54_std = std(percent_empty_fiber_tracts_tpca_smooth_54);

fiber_similarity_39_mean_all=mean(fiber_similarity_39_mean);
fiber_similarity_39_std_all=mean(fiber_similarity_39_std);
fiber_similarity_tpca_39_mean_all=mean(fiber_similarity_tpca_39_mean);
fiber_similarity_tpca_39_std_all=mean(fiber_similarity_tpca_39_std);
fiber_similarity_tpca_smooth_39_mean_all=mean(fiber_similarity_tpca_smooth_39_mean);
fiber_similarity_tpca_smooth_39_std_all=mean(fiber_similarity_tpca_smooth_39_std);
fiber_similarity_54_mean_all=mean(fiber_similarity_54_mean);
fiber_similarity_54_std_all=mean(fiber_similarity_54_std);
fiber_similarity_tpca_54_mean_all=mean(fiber_similarity_tpca_54_mean);
fiber_similarity_tpca_54_std_all=mean(fiber_similarity_tpca_54_std);
fiber_similarity_tpca_smooth_54_mean_all=mean(fiber_similarity_tpca_smooth_54_mean);
fiber_similarity_tpca_smooth_54_std_all=mean(fiber_similarity_tpca_smooth_54_std);

euclid_distance_39_mean_all=mean(euclid_distance_39_mean);
euclid_distance_39_std_all=mean(euclid_distance_39_std);
euclid_distance_tpca_39_mean_all=mean(euclid_distance_tpca_39_mean);
euclid_distance_tpca_39_std_all=mean(euclid_distance_tpca_39_std);
euclid_distance_tpca_smooth_39_mean_all=mean(euclid_distance_tpca_smooth_39_mean);
euclid_distance_tpca_smooth_39_std_all=mean(euclid_distance_tpca_smooth_39_std);
euclid_distance_54_mean_all=mean(euclid_distance_54_mean);
euclid_distance_54_std_all=mean(euclid_distance_54_std);
euclid_distance_tpca_54_mean_all=mean(euclid_distance_tpca_54_mean);
euclid_distance_tpca_54_std_all=mean(euclid_distance_tpca_54_std);
euclid_distance_tpca_smooth_54_mean_all=mean(euclid_distance_tpca_smooth_54_mean);
euclid_distance_tpca_smooth_54_std_all=mean(euclid_distance_tpca_smooth_54_std);

%% Create table to plot data on R

SNR_Cat = repelem(["SNR_39";"SNR_54"],[3000,3000]);
Data_Cat = repelem(["Raw";"tpca";"tpca_Smoothn"],[1000,1000,1000]);
Data_Cat = [Data_Cat;Data_Cat];
Curvature_m1 = [all_mean_curvature_m1_noise_39;all_mean_curvature_m1_noise_tpca_39;...
    all_mean_curvature_m1_noise_tpca_39_smooth;...
    all_mean_curvature_m1_noise_54;all_mean_curvature_m1_noise_tpca_54;...
    all_mean_curvature_m1_noise_tpca_54_smooth];
Fiber_length_mm = [all_total_distance_mm_noise_39;all_total_distance_mm_noise_tpca_39;...
    all_total_distance_mm_noise_tpca_39_smooth;...
    all_total_distance_mm_noise_54;all_total_distance_mm_noise_tpca_54;...
    all_total_distance_mm_noise_tpca_54_smooth];
Pennation_angle_deg = [all_mean_alpha_deg_noise_39;all_mean_alpha_deg_noise_tpca_39;...
    all_mean_alpha_deg_noise_tpca_39_smooth;...
    all_mean_alpha_deg_noise_54;all_mean_alpha_deg_noise_tpca_54;...
    all_mean_alpha_deg_noise_tpca_54_smooth];
Fiber_similarity = [fiber_similarity_39_mean;fiber_similarity_tpca_39_mean;...
    fiber_similarity_tpca_smooth_39_mean;
    fiber_similarity_54_mean;fiber_similarity_tpca_54_mean;...
    fiber_similarity_tpca_smooth_54_mean];
Simulation_Table = table(SNR_Cat,Data_Cat,Fiber_length_mm,Pennation_angle_deg,...
    Curvature_m1,Fiber_similarity); 
writetable(Simulation_Table,'Sim_Results_TPCA_APO3_SNR39_54.csv')

%% 
% Results for the manuscript

fiber_length_39_mean=mean(all_total_distance_mm_noise_39)
fiber_length_39_ci=quantile(all_total_distance_mm_noise_39,[0.025 0.975])

fiber_length_tpca_39_mean=mean(all_total_distance_mm_noise_tpca_39)
fiber_length_tpca_39_ci=quantile(all_total_distance_mm_noise_tpca_39,[0.025 0.975])

fiber_length_tpca_smooth_39_mean=mean(all_total_distance_mm_noise_tpca_39_smooth)
fiber_length_tpca_smooth_39_ci=quantile(all_total_distance_mm_noise_tpca_39_smooth,[0.025 0.975])

fiber_length_54_mean=mean(all_total_distance_mm_noise_54)
fiber_length_54_ci=quantile(all_total_distance_mm_noise_54,[0.025 0.975])

fiber_length_tpca_54_mean=mean(all_total_distance_mm_noise_tpca_54)
fiber_length_tpca_54_ci=quantile(all_total_distance_mm_noise_tpca_54,[0.025 0.975])

fiber_length_tpca_smooth_54_mean=mean(all_total_distance_mm_noise_tpca_54_smooth)
fiber_length_tpca_smooth_54_ci=quantile(all_total_distance_mm_noise_tpca_54_smooth,[0.025 0.975])

pen_angle_39_mean=mean(all_mean_alpha_deg_noise_39)
pen_angle_39_ci=quantile(all_mean_alpha_deg_noise_39,[0.025 0.975])

pen_angle_tpca_39_mean=mean(all_mean_alpha_deg_noise_tpca_39)
pen_angle_tpca_39_ci=quantile(all_mean_alpha_deg_noise_tpca_39,[0.025 0.975])

pen_angle_tpca_smooth_39_mean=mean(all_mean_alpha_deg_noise_tpca_39_smooth)
pen_angle_tpca_smooth_39_ci=quantile(all_mean_alpha_deg_noise_tpca_39_smooth,[0.025 0.975])

pen_angle_54_mean=mean(all_mean_alpha_deg_noise_54)
pen_angle_54_ci=quantile(all_mean_alpha_deg_noise_54,[0.025 0.975])

pen_angle_tpca_54_mean=mean(all_mean_alpha_deg_noise_tpca_54)
pen_angle_tpca_54_ci=quantile(all_mean_alpha_deg_noise_tpca_54,[0.025 0.975])

pen_angle_tpca_smooth_54_mean=mean(all_mean_alpha_deg_noise_tpca_54_smooth)
pen_angle_tpca_smooth_54_ci=quantile(all_mean_alpha_deg_noise_tpca_54_smooth,[0.025 0.975])

curvature_39_mean=mean(all_mean_curvature_m1_noise_39)
curvature_39_ci=quantile(all_mean_curvature_m1_noise_39,[0.025 0.975])

curvature_tpca_39_mean=mean(all_mean_curvature_m1_noise_tpca_39)
curvature_tpca_39_ci=quantile(all_mean_curvature_m1_noise_tpca_39,[0.025 0.975])

curvature_tpca_smooth_39_mean=mean(all_mean_curvature_m1_noise_tpca_39_smooth)
curvature_tpca_smooth_39_ci=quantile(all_mean_curvature_m1_noise_tpca_39_smooth,[0.025 0.975])

curvature_54_mean=mean(all_mean_curvature_m1_noise_54)
curvature_54_ci=quantile(all_mean_curvature_m1_noise_54,[0.025 0.975])

curvature_tpca_54_mean=mean(all_mean_curvature_m1_noise_tpca_54)
curvature_tpca_54_ci=quantile(all_mean_curvature_m1_noise_tpca_54,[0.025 0.975])

curvature_tpca_smooth_54_mean=mean(all_mean_curvature_m1_noise_tpca_54_smooth)
curvature_tpca_smooth_54_ci=quantile(all_mean_curvature_m1_noise_tpca_54_smooth,[0.025 0.975])

%% Plot 1st eigenvector maps

clim=[-0.5 0.5];

figure
t = tiledlayout(2,5);

nexttile(1,[2 2]);
imagesc(squeeze(squeeze(e1_map(:,:,21,1))).*muscle_mask(:,:,21),clim); pbaspect([1 1 1])
title('Noise-free map')

nexttile(3)
imagesc(squeeze(squeeze(e1map_39(:,:,21,1))).*muscle_mask(:,:,21),clim); pbaspect([1 1 1])
title('SNR = 39')

nexttile(4)
imagesc(squeeze(squeeze(e1map_tpca39(:,:,21,1))).*muscle_mask(:,:,21),clim); pbaspect([1 1 1])
title('SNR = 39, TPCA')

nexttile(5)
imagesc(squeeze(squeeze(e1map_39_tpca_smooth(:,:,21,1))).*muscle_mask(:,:,21),clim); pbaspect([1 1 1])
title('SNR = 39, TPCA + smoothn')

nexttile(8)
imagesc(squeeze(squeeze(e1map_54(:,:,21,1))).*muscle_mask(:,:,21),clim); pbaspect([1 1 1])
title('SNR = 54')

nexttile(9)
imagesc(squeeze(squeeze(e1map_tpca54(:,:,21,1))).*muscle_mask(:,:,21),clim); pbaspect([1 1 1])
title('SNR = 54, TPCA')

nexttile(10)
imagesc(squeeze(squeeze(e1map_54_tpca_smooth(:,:,21,1))).*muscle_mask(:,:,21),clim); pbaspect([1 1 1])
title('SNR = 54, TPCA + smoothn')


cb = colorbar;
cb.Layout.Tile = 'east';
title(t,'1st eigenvector - x component')

clim=[-0.5 0.5];

figure
t = tiledlayout(2,5);

nexttile(1,[2 2]);
imagesc(squeeze(squeeze(e1_map(:,:,21,2))).*muscle_mask(:,:,21),clim); pbaspect([1 1 1])
title('Noise-free map')

nexttile(3)
imagesc(squeeze(squeeze(e1map_39(:,:,21,2))).*muscle_mask(:,:,21),clim); pbaspect([1 1 1])
title('SNR = 39')

nexttile(4)
imagesc(squeeze(squeeze(e1map_tpca39(:,:,21,2))).*muscle_mask(:,:,21),clim); pbaspect([1 1 1])
title('SNR = 39, TPCA')

nexttile(5)
imagesc(squeeze(squeeze(e1map_39_tpca_smooth(:,:,21,2))).*muscle_mask(:,:,21),clim); pbaspect([1 1 1])
title('SNR = 39, TPCA + smoothn')

nexttile(8)
imagesc(squeeze(squeeze(e1map_54(:,:,21,2))).*muscle_mask(:,:,21),clim); pbaspect([1 1 1])
title('SNR = 54')

nexttile(9)
imagesc(squeeze(squeeze(e1map_tpca54(:,:,21,2))).*muscle_mask(:,:,21),clim); pbaspect([1 1 1])
title('SNR = 54, TPCA')

nexttile(10)
imagesc(squeeze(squeeze(e1map_54_tpca_smooth(:,:,21,2))).*muscle_mask(:,:,21),clim); pbaspect([1 1 1])
title('SNR = 54, TPCA + smoothn')


cb = colorbar;
cb.Layout.Tile = 'east';
title(t,'1st eigenvector - y component')

figure
t = tiledlayout(2,5);
clim=[0.7 1];

nexttile(1,[2 2]);
imagesc(squeeze(squeeze(e1_map(:,:,21,3))).*muscle_mask(:,:,21),clim); pbaspect([1 1 1])
title('Noise-free map')

nexttile(3)
imagesc(squeeze(squeeze(e1map_39(:,:,21,3))).*muscle_mask(:,:,21),clim); pbaspect([1 1 1])
title('SNR = 39')

nexttile(4)
imagesc(squeeze(squeeze(e1map_tpca39(:,:,21,3))).*muscle_mask(:,:,21),clim); pbaspect([1 1 1])
title('SNR = 39, TPCA')

nexttile(5)
imagesc(squeeze(squeeze(e1map_39_tpca_smooth(:,:,21,3))).*muscle_mask(:,:,21),clim); pbaspect([1 1 1])
title('SNR = 39, TPCA + smoothn')

nexttile(8)
imagesc(squeeze(squeeze(e1map_54(:,:,21,3))).*muscle_mask(:,:,21),clim); pbaspect([1 1 1])
title('SNR = 54')

nexttile(9)
imagesc(squeeze(squeeze(e1map_tpca54(:,:,21,3))).*muscle_mask(:,:,21),clim); pbaspect([1 1 1])
title('SNR = 54, TPCA')

nexttile(10)
imagesc(squeeze(squeeze(e1map_54_tpca_smooth(:,:,21,3))).*muscle_mask(:,:,21),clim); pbaspect([1 1 1])
title('SNR = 54, TPCA + smoothn')


cb = colorbar;
cb.Layout.Tile = 'east';
title(t,'1st eigenvector - z component')

%% Plot angle of deviation maps

clim=[0 15];

angle_noise_39(angle_noise_39 ==0) = NaN;
angle_noise_tpca_39(angle_noise_tpca_39 ==0) = NaN;
angle_noise_tpca_39_smooth(angle_noise_tpca_39_smooth ==0) = NaN;
angle_noise_54(angle_noise_54 ==0) = NaN;
angle_noise_tpca_54(angle_noise_tpca_54 ==0) = NaN;
angle_noise_tpca_54_smooth(angle_noise_tpca_54_smooth ==0) = NaN;

figure
t = tiledlayout(2,3);
nexttile
imagesc(squeeze(squeeze(angle_noise_39(:,:,21,end))),'AlphaData',~isnan(squeeze(squeeze(angle_noise_39(:,:,21,end)))),clim);set(gca,'color',[0.7 0.7 0.7]); pbaspect([1 1 1])
title('SNR = 39')

nexttile
imagesc(squeeze(squeeze(angle_noise_tpca_39(:,:,21,end))),'AlphaData',~isnan(squeeze(squeeze(angle_noise_39(:,:,21,end)))),clim);set(gca,'color',[0.7 0.7 0.7]); pbaspect([1 1 1])
title('SNR = 39, TPCA')

nexttile
imagesc(squeeze(squeeze(angle_noise_tpca_39_smooth(:,:,21,end))),'AlphaData',~isnan(squeeze(squeeze(angle_noise_39(:,:,21,end)))),clim);set(gca,'color',[0.7 0.7 0.7]); pbaspect([1 1 1])
title('SNR = 39, TPCA + smoothn')

nexttile
imagesc(squeeze(squeeze(angle_noise_54(:,:,21,end))),'AlphaData',~isnan(squeeze(squeeze(angle_noise_39(:,:,21,end)))),clim);set(gca,'color',[0.7 0.7 0.7]); pbaspect([1 1 1])
title('SNR = 54')

nexttile
imagesc(squeeze(squeeze(angle_noise_tpca_54(:,:,21,end))),'AlphaData',~isnan(squeeze(squeeze(angle_noise_39(:,:,21,end)))),clim);set(gca,'color',[0.7 0.7 0.7]); pbaspect([1 1 1])
title('SNR = 54, TPCA')

nexttile
imagesc(squeeze(squeeze(angle_noise_tpca_54_smooth(:,:,21,end))),'AlphaData',~isnan(squeeze(squeeze(angle_noise_39(:,:,21,end)))),clim);set(gca,'color',[0.7 0.7 0.7]); pbaspect([1 1 1])
title('SNR = 54, TPCA + smoothn')

cb = colorbar;
cb.Layout.Tile = 'east';
title(t,'Angular deviation \theta wrt 1st eigenvector')

cd('S:\Muscle_DTI\Roberto\DTI_Muscle_Analysis\Processed_Data\Aim1E')
print(gcf,'Angular_Deviation_Images_SNR39_54_TPCA.png','-dpng','-r1000'); 

figure
imagesc(squeeze(squeeze(angle_noise_39(:,:,21,end))),'AlphaData',~isnan(squeeze(squeeze(angle_noise_39(:,:,21,end)))),clim);set(gca,'color',[0.7 0.7 0.7]); pbaspect([1 1 1])
title('SNR = 39')
%print(gcf,'Angular_Deviation_Images_SNR39_Raw.tif','-dtiffn','-r300'); 

figure
imagesc(squeeze(squeeze(angle_noise_tpca_39(:,:,21,end))),'AlphaData',~isnan(squeeze(squeeze(angle_noise_39(:,:,21,end)))),clim);set(gca,'color',[0.7 0.7 0.7]); pbaspect([1 1 1])
title('SNR = 39, TPCA')
%print(gcf,'Angular_Deviation_Images_SNR39_tpca.png','-dpng','-r2000'); 

figure
imagesc(squeeze(squeeze(angle_noise_tpca_39_smooth(:,:,21,end))),'AlphaData',~isnan(squeeze(squeeze(angle_noise_39(:,:,21,end)))),clim);set(gca,'color',[0.7 0.7 0.7]); pbaspect([1 1 1])
title('SNR = 39, TPCA + smoothn')
%print(gcf,'Angular_Deviation_Images_SNR39_tpca_Smoothn.png','-dpng','-r2000'); 

figure
imagesc(squeeze(squeeze(angle_noise_54(:,:,21,end))),'AlphaData',~isnan(squeeze(squeeze(angle_noise_39(:,:,21,end)))),clim);set(gca,'color',[0.7 0.7 0.7]); pbaspect([1 1 1])
title('SNR = 54')
%print(gcf,'Angular_Deviation_Images_SNR54_Raw.png','-dpng','-r2000'); 

figure
imagesc(squeeze(squeeze(angle_noise_tpca_54(:,:,21,end))),'AlphaData',~isnan(squeeze(squeeze(angle_noise_39(:,:,21,end)))),clim);set(gca,'color',[0.7 0.7 0.7]); pbaspect([1 1 1])
title('SNR = 54, TPCA')
%print(gcf,'Angular_Deviation_Images_SNR54_tpca.png','-dpng','-r2000'); 

figure
imagesc(squeeze(squeeze(angle_noise_tpca_54_smooth(:,:,21,end))),'AlphaData',~isnan(squeeze(squeeze(angle_noise_39(:,:,21,end)))),clim);set(gca,'color',[0.7 0.7 0.7]); pbaspect([1 1 1])
title('SNR = 54, TPCA + smoothn')
%print(gcf,'Angular_Deviation_Images_SNR54_tpca_Smoothn.png','-dpng','-r2000'); 

%% 
% Results for the manuscript

angle_noise_39_iterations=squeeze(mean(angle_noise_39,[1 2 3],"omitnan"));
angle_noise_39_iterations_mean = mean(angle_noise_39_iterations)
angle_noise_39_ci=quantile(angle_noise_39_iterations,[0.025 0.975])

angle_noise_tpca_39_iterations = squeeze(mean(angle_noise_tpca_39,[1 2 3],"omitnan"));
angle_noise_tpca_39_iterations_mean = mean(angle_noise_tpca_39_iterations)
angle_noise_tpca_39_ci=quantile(angle_noise_tpca_39_iterations,[0.025 0.975])

angle_noise_tpca_39_smooth_iterations = squeeze(mean(angle_noise_tpca_39_smooth,[1 2 3],"omitnan"));
angle_noise_tpca_39_smooth_iterations_mean = mean(angle_noise_tpca_39_smooth_iterations)
angle_noise_tpca_39_smooth_ci=quantile(angle_noise_tpca_39_smooth_iterations,[0.025 0.975])

angle_noise_54_iterations=squeeze(mean(angle_noise_54,[1 2 3],"omitnan"));
angle_noise_54_iterations_mean = mean(angle_noise_54_iterations)
angle_noise_54_ci=quantile(angle_noise_54_iterations,[0.025 0.975])

angle_noise_tpca_54_iterations = squeeze(mean(angle_noise_tpca_54,[1 2 3],"omitnan"));
angle_noise_tpca_54_iterations_mean = mean(angle_noise_tpca_54_iterations)
angle_noise_tpca_54_ci=quantile(angle_noise_tpca_54_iterations,[0.025 0.975])

angle_noise_tpca_54_smooth_iterations = squeeze(mean(angle_noise_tpca_54_smooth,[1 2 3],"omitnan"));
angle_noise_tpca_54_smooth_iterations_mean = mean(angle_noise_tpca_54_smooth_iterations)
angle_noise_tpca_54_smooth_ci=quantile(angle_noise_tpca_54_smooth_iterations,[0.025 0.975])