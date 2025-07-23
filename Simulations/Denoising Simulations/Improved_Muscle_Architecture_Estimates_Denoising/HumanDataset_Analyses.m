%% This code contains the analyses made on the human in vivo dataset presented
% in the manuscript 'Improved DTI-based skeletal muscle architecture
% estimation via diffusion weighted image denoising and first eigenvector
% map smoothing'

% Written by: Roberto Pineda Guzman, Carle Foundation Hospital

% Supporting functions:
% smoothn: https://www.mathworks.com/matlabcentral/fileexchange/25634-smoothn/

close all
clear 
clc

%% Open pre-processed data from the human dataset

% Dataset available upon request
load("Dataset.mat")

% Consistent variable names with simulations

muscle_mask = TA_mask;
roi_mesh_rstrct = roi_mesh;
apo_mask = roi_mask;

composite_mask(:,:,:,1) = muscle_mask - apo_mask;
composite_mask(:,:,:,2) = apo_mask;
composite_mask(composite_mask<0)=0;

% Obtain DTI data using only some of the excitations to decrease SNR

dti_img_4d_reg_nonavg = dti_img_raw_reg(:,:,:,1:11); % DTI data using only 1 excitation
dti_img_b0_reg_2avg = (dti_img_raw_reg(:,:,:,1)+dti_img_raw_reg(:,:,:,12))/2; % b = 0 image using 2 excitations
dti_img_4d_reg_b0_2avg = cat(4,dti_img_b0_reg_2avg,dti_img_raw_reg(:,:,:,2:11)); % DTI data using 2 b = 0 excitations and 1 excitation per diffusion encoding direction
dti_img_4d_reg_b0_4avg = cat(4,dti_img_4d_reg(:,:,:,1),dti_img_raw_reg(:,:,:,2:11)); % DTI data using 4 b = 0 excitations
dti_img_4d_reg_2avg = cat(4,dti_img_b0_reg_2avg,dti_img_4d_reg(:,:,:,2:11)); % DTI data using 2 excitations

% Obtain 4D array with only b = 0 values
dti_img_4d_reg_b0 = cat(4,dti_img_raw_reg(:,:,:,1),dti_img_raw_reg(:,:,:,12),dti_img_raw_reg(:,:,:,23),dti_img_raw_reg(:,:,:,24));

% Compute variance of b = 0 images
img_b0_std = std(dti_img_4d_reg_b0,0,4);
img_b0_var = img_b0_std.^2;

%% Anisotropic smoothing on the DTI data 

%set parameters based on Buck et al
noise = 5;
sigma = noise/100;
rho = 2*sigma;
delta_t = noise*3/44;
schemetype = 'implicit_multisplitting';
isfasteig  = true;
isnormg = false;
dti_res = [dim_param_dti_1.voxel_dim(1) dim_param_dti_1.voxel_dim(2) dim_param_dti_1.voxel_dim(3)];

tic %measure computation time
dti_img_4d_smooth = aniso4D_smoothing(dti_img_4d_reg, sigma, rho, delta_t, dti_res, schemetype, isnormg, isfasteig);
dti_img_4d_smooth_b0_4avg = aniso4D_smoothing(dti_img_4d_reg_b0_4avg, sigma, rho, delta_t, dti_res, schemetype, isnormg, isfasteig);
dti_img_4d_smooth_2avg = aniso4D_smoothing(dti_img_4d_reg_2avg, sigma, rho, delta_t, dti_res, schemetype, isnormg, isfasteig);
dti_img_4d_smooth_b0_2avg = aniso4D_smoothing(dti_img_4d_reg_b0_2avg, sigma, rho, delta_t, dti_res, schemetype, isnormg, isfasteig);
dti_img_4d_smooth_nonavg = aniso4D_smoothing(dti_img_4d_reg_nonavg, sigma, rho, delta_t, dti_res, schemetype, isnormg, isfasteig);
toc

%% Apply TPCA smoothing on the data

tic %measure computation time
dti_img_4d_tpca = TPCA_denoising(dti_img_4d_reg,muscle_mask,[9 9 1],img_b0_var);
dti_img_4d_tpca_b0_4avg = TPCA_denoising(dti_img_4d_reg_b0_4avg,muscle_mask,[9 9 1],img_b0_var);
dti_img_4d_tpca_2avg = TPCA_denoising(dti_img_4d_reg_2avg,muscle_mask,[9 9 1],img_b0_var); 
dti_img_4d_tpca_b0_2avg = TPCA_denoising(dti_img_4d_reg_b0_2avg,muscle_mask,[9 9 1],img_b0_var);
dti_img_4d_tpca_nonavg = TPCA_denoising(dti_img_4d_reg_nonavg,muscle_mask,[9 9 1],img_b0_var);
toc

%% Analyze b = 0 (4 avg), b = 450 (2 avg) data

%% 
% Compute FA and first eigenvector map of DTI data

%Correct diffusion directions fror FFS scans
diff_dir=diff_param_1.bvec; %diff_param1 and diff_param2 are the same

b_val=mean(diff_param_1.bval(2:end));

%Tensor computation

FA_raw = zeros(size(dti_img_4d_reg,1),size(dti_img_4d_reg,2),size(dti_img_4d_reg,3));
E1_450_raw = zeros(size(dti_img_4d_reg,1),size(dti_img_4d_reg,2),size(dti_img_4d_reg,3),3);
FA_aniso = zeros(size(dti_img_4d_reg,1),size(dti_img_4d_reg,2),size(dti_img_4d_reg,3));
E1_450_aniso = zeros(size(dti_img_4d_reg,1),size(dti_img_4d_reg,2),size(dti_img_4d_reg,3),3);
FA_tpca = zeros(size(dti_img_4d_reg,1),size(dti_img_4d_reg,2),size(dti_img_4d_reg,3));
E1_450_tpca = zeros(size(dti_img_4d_reg,1),size(dti_img_4d_reg,2),size(dti_img_4d_reg,3),3);

for s = 1:size(dti_img_4d_reg,3)
    for r = 1:size(dti_img_4d_reg,1)
        for c = 1:size(dti_img_4d_reg,2)
                dir_m = diff_dir;  %diffusion encoding matrix
                dir_m = dir_m(2:end,:);

                %estimate tensor
                %use TA mask to make process faster
                if muscle_mask(r,c,s)>0
                    % Raw data
                    sig_v = squeeze(dti_img_4d_reg(r,c,s,:));    %signal vector
                    D = signal2tensor2(sig_v, dir_m, b_val);
                    [E,L] = svd(D);
    
                    L_v = diag(L);
                    md = mean(L_v);
                    fa = sqrt(3/2) * sqrt(((L_v(1) - md)^2 + (L_v(2) - md)^2 + (L_v(3) - md)^2) / ...
                        (L_v(1)^2 + L_v(2)^2 + L_v(3)^2));

                    loop_E1 = E(:,L_v==max(L_v));
                    if loop_E1(3)<0
                        loop_E1 = -loop_E1;
                    end
                    
                    E1_450_raw(r,c,s,:) = loop_E1;
                    FA_raw(r,c,s,:) = fa;

                    % AIS
                    sig_v = squeeze(dti_img_4d_smooth(r,c,s,:));    %signal vector
                    D = signal2tensor2(sig_v, dir_m, b_val);
                    [E,L] = svd(D);
    
                    L_v = diag(L);
                    md = mean(L_v);
                    fa = sqrt(3/2) * sqrt(((L_v(1) - md)^2 + (L_v(2) - md)^2 + (L_v(3) - md)^2) / ...
                        (L_v(1)^2 + L_v(2)^2 + L_v(3)^2));

                    loop_E1 = E(:,L_v==max(L_v));
                    if loop_E1(3)<0
                        loop_E1 = -loop_E1;
                    end
                    
                    E1_450_aniso(r,c,s,:) = loop_E1;
                    FA_aniso(r,c,s,:) = fa;

                    % TPCA
                    sig_v = squeeze(dti_img_4d_tpca(r,c,s,:));    %signal vector
                    D = signal2tensor2(sig_v, dir_m, b_val);
                    [E,L] = svd(D);
    
                    L_v = diag(L);
                    md = mean(L_v);
                    fa = sqrt(3/2) * sqrt(((L_v(1) - md)^2 + (L_v(2) - md)^2 + (L_v(3) - md)^2) / ...
                        (L_v(1)^2 + L_v(2)^2 + L_v(3)^2));

                    loop_E1 = E(:,L_v==max(L_v));
                    if loop_E1(3)<0
                        loop_E1 = -loop_E1;
                    end
                    
                    E1_450_tpca(r,c,s,:) = loop_E1;
                    FA_tpca(r,c,s,:) = fa;

                end
        end
    end
end

%% 
% Smooth first eigenvector field of TA muscle

[E1_450_raw_smooth,smooth_parameter_raw]=smoothn({E1_450_raw(:,:,:,1),E1_450_raw(:,:,:,2),E1_450_raw(:,:,:,3)});
[E1_450_aniso_smooth,smooth_parameter_aniso]=smoothn({E1_450_aniso(:,:,:,1),E1_450_aniso(:,:,:,2),E1_450_aniso(:,:,:,3)});
[E1_450_tpca_smooth,smooth_parameter_tpca]=smoothn({E1_450_tpca(:,:,:,1),E1_450_tpca(:,:,:,2),E1_450_tpca(:,:,:,3)});


E1_450_raw_smooth=cat(4,E1_450_raw_smooth{1},E1_450_raw_smooth{2},E1_450_raw_smooth{3});
E1_450_aniso_smooth=cat(4,E1_450_aniso_smooth{1},E1_450_aniso_smooth{2},E1_450_aniso_smooth{3});
E1_450_tpca_smooth=cat(4,E1_450_tpca_smooth{1},E1_450_tpca_smooth{2},E1_450_tpca_smooth{3});


for r_img=1:size(E1_450_aniso,1)
    for c_img=1:size(E1_450_aniso,2)
        for s_img=1:size(E1_450_aniso,3)
            if muscle_mask(r_img,c_img,s_img)>0

                E1_raw_smooth = squeeze(E1_450_raw_smooth(r_img,c_img,s_img,:));
                E1_aniso_smooth = squeeze(E1_450_aniso_smooth(r_img,c_img,s_img,:));
                E1_tpca_smooth = squeeze(E1_450_tpca_smooth(r_img,c_img,s_img,:));

                E1_raw_smooth = E1_raw_smooth/norm(E1_raw_smooth);
                E1_aniso_smooth = E1_aniso_smooth/norm(E1_aniso_smooth);
                E1_tpca_smooth = E1_tpca_smooth/norm(E1_tpca_smooth);

                E1_450_raw_smooth(r_img,c_img,s_img,:)=E1_raw_smooth;
                E1_450_aniso_smooth(r_img,c_img,s_img,:)=E1_aniso_smooth;
                E1_450_tpca_smooth(r_img,c_img,s_img,:)=E1_tpca_smooth;
            end
        end
    end
end

for i = 1:3
    E1_450_raw_smooth(:,:,:,i)=squeeze(E1_450_raw_smooth(:,:,:,i)).*muscle_mask;
    E1_450_aniso_smooth(:,:,:,i)=squeeze(E1_450_aniso_smooth(:,:,:,i)).*muscle_mask;
    E1_450_tpca_smooth(:,:,:,i)=squeeze(E1_450_tpca_smooth(:,:,:,i)).*muscle_mask;
end

%% 
% Create E1-FA maps for fiber tractography

e1fa_raw = cat(4,E1_450_raw,FA_raw);
e1fa_aniso = cat(4,E1_450_aniso,FA_aniso);
e1fa_tpca = cat(4,E1_450_tpca,FA_tpca);
e1fa_raw_smooth = cat(4,E1_450_raw_smooth,FA_raw);
e1fa_aniso_smooth = cat(4,E1_450_aniso_smooth,FA_aniso);
e1fa_tpca_smooth = cat(4,E1_450_tpca_smooth,FA_tpca);

%%
% Visualize first eigenvector maps

rect_crop = [39, 22, 41, 41];

E1_x_raw=squeeze(E1_450_raw(:,:,:,1));
E1_y_raw=squeeze(E1_450_raw(:,:,:,2));
E1_z_raw=squeeze(E1_450_raw(:,:,:,3));

E1_x_raw_Crop=imcrop3(E1_x_raw,[rect_crop(1),rect_crop(2),1,rect_crop(3),rect_crop(4),size(E1_450_raw,3)-1]);
E1_y_raw_Crop=imcrop3(E1_y_raw,[rect_crop(1),rect_crop(2),1,rect_crop(3),rect_crop(4),size(E1_450_raw,3)-1]);
E1_z_raw_Crop=imcrop3(E1_z_raw,[rect_crop(1),rect_crop(2),1,rect_crop(3),rect_crop(4),size(E1_450_raw,3)-1]);

E1_x_raw_Crop(E1_x_raw_Crop==0)=NaN;
E1_y_raw_Crop(E1_y_raw_Crop==0)=NaN;
E1_z_raw_Crop(E1_z_raw_Crop==0)=NaN;

E1_x_aniso=squeeze(E1_450_aniso(:,:,:,1));
E1_y_aniso=squeeze(E1_450_aniso(:,:,:,2));
E1_z_aniso=squeeze(E1_450_aniso(:,:,:,3));

E1_x_aniso_Crop=imcrop3(E1_x_aniso,[rect_crop(1),rect_crop(2),1,rect_crop(3),rect_crop(4),size(E1_450_raw,3)-1]);
E1_y_aniso_Crop=imcrop3(E1_y_aniso,[rect_crop(1),rect_crop(2),1,rect_crop(3),rect_crop(4),size(E1_450_raw,3)-1]);
E1_z_aniso_Crop=imcrop3(E1_z_aniso,[rect_crop(1),rect_crop(2),1,rect_crop(3),rect_crop(4),size(E1_450_raw,3)-1]);

E1_x_aniso_Crop(E1_x_aniso_Crop==0)=NaN;
E1_y_aniso_Crop(E1_y_aniso_Crop==0)=NaN;
E1_z_aniso_Crop(E1_z_aniso_Crop==0)=NaN;

E1_x_raw_smooth=squeeze(E1_450_raw_smooth(:,:,:,1));
E1_y_raw_smooth=squeeze(E1_450_raw_smooth(:,:,:,2));
E1_z_raw_smooth=squeeze(E1_450_raw_smooth(:,:,:,3));

E1_x_raw_smooth_Crop=imcrop3(E1_x_raw_smooth,[rect_crop(1),rect_crop(2),1,rect_crop(3),rect_crop(4),size(E1_450_raw,3)-1]);
E1_y_raw_smooth_Crop=imcrop3(E1_y_raw_smooth,[rect_crop(1),rect_crop(2),1,rect_crop(3),rect_crop(4),size(E1_450_raw,3)-1]);
E1_z_raw_smooth_Crop=imcrop3(E1_z_raw_smooth,[rect_crop(1),rect_crop(2),1,rect_crop(3),rect_crop(4),size(E1_450_raw,3)-1]);

E1_x_raw_smooth_Crop(E1_x_raw_smooth_Crop==0)=NaN;
E1_y_raw_smooth_Crop(E1_y_raw_smooth_Crop==0)=NaN;
E1_z_raw_smooth_Crop(E1_z_raw_smooth_Crop==0)=NaN;

E1_x_aniso_smooth=squeeze(E1_450_aniso_smooth(:,:,:,1));
E1_y_aniso_smooth=squeeze(E1_450_aniso_smooth(:,:,:,2));
E1_z_aniso_smooth=squeeze(E1_450_aniso_smooth(:,:,:,3));

E1_x_aniso_smooth_Crop=imcrop3(E1_x_aniso_smooth,[rect_crop(1),rect_crop(2),1,rect_crop(3),rect_crop(4),size(E1_450_raw,3)-1]);
E1_y_aniso_smooth_Crop=imcrop3(E1_y_aniso_smooth,[rect_crop(1),rect_crop(2),1,rect_crop(3),rect_crop(4),size(E1_450_raw,3)-1]);
E1_z_aniso_smooth_Crop=imcrop3(E1_z_aniso_smooth,[rect_crop(1),rect_crop(2),1,rect_crop(3),rect_crop(4),size(E1_450_raw,3)-1]);

E1_x_aniso_smooth_Crop(E1_x_aniso_smooth_Crop==0)=NaN;
E1_y_aniso_smooth_Crop(E1_y_aniso_smooth_Crop==0)=NaN;
E1_z_aniso_smooth_Crop(E1_z_aniso_smooth_Crop==0)=NaN;

E1_x_tpca=squeeze(E1_450_tpca(:,:,:,1));
E1_y_tpca=squeeze(E1_450_tpca(:,:,:,2));
E1_z_tpca=squeeze(E1_450_tpca(:,:,:,3));

E1_x_tpca_Crop=imcrop3(E1_x_tpca,[rect_crop(1),rect_crop(2),1,rect_crop(3),rect_crop(4),size(E1_450_raw,3)-1]);
E1_y_tpca_Crop=imcrop3(E1_y_tpca,[rect_crop(1),rect_crop(2),1,rect_crop(3),rect_crop(4),size(E1_450_raw,3)-1]);
E1_z_tpca_Crop=imcrop3(E1_z_tpca,[rect_crop(1),rect_crop(2),1,rect_crop(3),rect_crop(4),size(E1_450_raw,3)-1]);

E1_x_tpca_Crop(E1_x_tpca_Crop==0)=NaN;
E1_y_tpca_Crop(E1_y_tpca_Crop==0)=NaN;
E1_z_tpca_Crop(E1_z_tpca_Crop==0)=NaN;

E1_x_tpca_smooth=squeeze(E1_450_tpca_smooth(:,:,:,1));
E1_y_tpca_smooth=squeeze(E1_450_tpca_smooth(:,:,:,2));
E1_z_tpca_smooth=squeeze(E1_450_tpca_smooth(:,:,:,3));

E1_x_tpca_smooth_Crop=imcrop3(E1_x_tpca_smooth,[rect_crop(1),rect_crop(2),1,rect_crop(3),rect_crop(4),size(E1_450_raw,3)-1]);
E1_y_tpca_smooth_Crop=imcrop3(E1_y_tpca_smooth,[rect_crop(1),rect_crop(2),1,rect_crop(3),rect_crop(4),size(E1_450_raw,3)-1]);
E1_z_tpca_smooth_Crop=imcrop3(E1_z_tpca_smooth,[rect_crop(1),rect_crop(2),1,rect_crop(3),rect_crop(4),size(E1_450_raw,3)-1]);

E1_x_tpca_smooth_Crop(E1_x_tpca_smooth_Crop==0)=NaN;
E1_y_tpca_smooth_Crop(E1_y_tpca_smooth_Crop==0)=NaN;
E1_z_tpca_smooth_Crop(E1_z_tpca_smooth_Crop==0)=NaN;

%% 
% Obtain image of first eigenvector map of muscle cross-section

figure
t = tiledlayout(3,6);

nexttile
imagesc(E1_x_raw_Crop(:,:,30),'AlphaData',~isnan(E1_x_raw_Crop(:,:,30)),[-0.3 0.3]); pbaspect([1 1 1]); colorbar
title('E1_x - TA - Raw')

nexttile
imagesc(E1_x_aniso_Crop(:,:,30),'AlphaData',~isnan(E1_x_aniso_Crop(:,:,30)),[-0.3 0.3]); pbaspect([1 1 1])
title('E1_x - TA - AnisoSmooth')

nexttile
imagesc(E1_x_tpca_Crop(:,:,30),'AlphaData',~isnan(E1_x_tpca_Crop(:,:,30)),[-0.3 0.3]); pbaspect([1 1 1])
title('E1_x - TA - TPCA')

nexttile
imagesc(E1_x_raw_smooth_Crop(:,:,30),'AlphaData',~isnan(E1_x_raw_smooth_Crop(:,:,30)),[-0.3 0.3]); pbaspect([1 1 1])
title('E1_x - TA - Smoothn')

nexttile
imagesc(E1_x_aniso_smooth_Crop(:,:,30),'AlphaData',~isnan(E1_x_aniso_smooth_Crop(:,:,30)),[-0.3 0.3]); pbaspect([1 1 1])
title('E1_x - TA - AnisoSmooth + Smoothn')

nexttile
imagesc(E1_x_tpca_smooth_Crop(:,:,30),'AlphaData',~isnan(E1_x_tpca_smooth_Crop(:,:,30)),[-0.3 0.3]); pbaspect([1 1 1])
title('E1_x - TA - TPCA + Smoothn')

nexttile
imagesc(E1_y_raw_Crop(:,:,30),'AlphaData',~isnan(E1_y_raw_Crop(:,:,30)),[-0.5 0.5]); pbaspect([1 1 1]); colorbar
title('E1_y - TA - Raw')

nexttile
imagesc(E1_y_aniso_Crop(:,:,30),'AlphaData',~isnan(E1_y_aniso_Crop(:,:,30)),[-0.5 0.5]); pbaspect([1 1 1])
title('E1_y - TA - AnisoSmooth')

nexttile
imagesc(E1_y_tpca_Crop(:,:,30),'AlphaData',~isnan(E1_y_tpca_Crop(:,:,30)),[-0.5 0.5]); pbaspect([1 1 1])
title('E1_y - TA - TPCA')

nexttile
imagesc(E1_y_raw_smooth_Crop(:,:,30),'AlphaData',~isnan(E1_y_raw_smooth_Crop(:,:,30)),[-0.5 0.5]); pbaspect([1 1 1])
title('E1_y - TA - Smoothn')

nexttile
imagesc(E1_y_aniso_smooth_Crop(:,:,30),'AlphaData',~isnan(E1_y_aniso_smooth_Crop(:,:,30)),[-0.5 0.5]); pbaspect([1 1 1])
title('E1_y - TA - AnisoSmooth + Smoothn')

nexttile
imagesc(E1_y_tpca_smooth_Crop(:,:,30),'AlphaData',~isnan(E1_y_tpca_smooth_Crop(:,:,30)),[-0.5 0.5]); pbaspect([1 1 1])
title('E1_y - TA - TPCA + Smoothn')

nexttile
imagesc(E1_z_raw_Crop(:,:,30),'AlphaData',~isnan(E1_z_raw_Crop(:,:,30)),[0.9 1]); pbaspect([1 1 1]); colorbar
title('E1_z - TA - Raw')

nexttile
imagesc(E1_z_aniso_Crop(:,:,30),'AlphaData',~isnan(E1_z_aniso_Crop(:,:,30)),[0.9 1]); pbaspect([1 1 1])
title('E1_z - TA - AnisoSmooth')

nexttile
imagesc(E1_z_tpca_Crop(:,:,30),'AlphaData',~isnan(E1_z_tpca_Crop(:,:,30)),[0.9 1]); pbaspect([1 1 1])
title('E1_z - TA - TPCA')

nexttile
imagesc(E1_z_raw_smooth_Crop(:,:,30),'AlphaData',~isnan(E1_z_raw_smooth_Crop(:,:,30)),[0.9 1]); pbaspect([1 1 1])
title('E1_z - TA - Smoothn')

nexttile
imagesc(E1_z_aniso_smooth_Crop(:,:,30),'AlphaData',~isnan(E1_z_aniso_smooth_Crop(:,:,30)),[0.9 1]); pbaspect([1 1 1])
title('E1_z - TA - AnisoSmooth + Smoothn')

nexttile
imagesc(E1_z_tpca_smooth_Crop(:,:,30),'AlphaData',~isnan(E1_z_tpca_smooth_Crop(:,:,30)),[0.9 1]); pbaspect([1 1 1])
title('E1_z - TA - TPCA + Smoothn')

cd('S:\Muscle_DTI\Roberto\DTI_Muscle_Analysis\Processed_Data\Aim1E')
print(gcf,'TA_eigenvector_Images_b0_4_b450_2.png','-dpng','-r1500'); 
savefig('TA_eigenvector_Images_b0_4_b450_2.fig')

%% 
% Set fiber-tractography parameters 

% set baseline visualization options 
fv_options.anat_dims = [200 dim_param_dti_1.voxel_dim(3)];
fv_options.anat_slices = 15:20:55;
fv_options.anat_alpha = 0;
fv_options.plot_mesh = 1;
fv_options.plot_mask = 0;
fv_options.plot_fibers = 1;
fv_options.mask_size = [144 144];
fv_options.mask_dims = [dim_param_dti_1.voxel_dim(1)*144 dim_param_dti_1.voxel_dim(3)];
fv_options.mask_color = [.5 0 0];
fv_options.mesh_size = [144 144];
fv_options.mesh_dims = [dim_param_dti_1.voxel_dim(1)*144 dim_param_dti_1.voxel_dim(3)];
fv_options.mesh_color = [.7 .7 .7];
fv_options.dti_size = [144 144];
fv_options.dti_dims = [dim_param_dti_1.voxel_dim(1)*144 dim_param_dti_1.voxel_dim(3)];
fv_options.mesh_dist = 0;

% set baseline fiber-tracking options
ft_options.ref_frame = 'LPS';                                               %LPS frame of reference
ft_options.image_orient = 'AL';                                             %anterior/left image orientation
ft_options.mesh_dist = 0;                                                   %no shift in mesh position
ft_options.prop_algo = 'eul';                                               %integrate E1 using Euler
ft_options.step_size = 1;                                                   %step size of 1 voxel width
ft_options.term_mthd = 'bin2';                                              %BIN2 requires two points to meet stop criteria
ft_options.angle_thrsh = [30 2];                                            %30 degree angle between current step and the step 2 points earlier
ft_options.fa_thrsh = [.1 .4];                                              %FA limits
ft_options.depth_ratio =  dim_param_dti_1.voxel_dim(3)/dim_param_dti_1.voxel_dim(1);                  

% set baseline fiber smoothing options
fs_options.dwi_res = [dim_param_dti_1.voxel_dim(1)*144 144 dim_param_dti_1.voxel_dim(3)];
fs_options.interpolation_step = 1;
fs_options.p_order = [3 3 3];
fs_options.tract_units = 'vx';
fs_options.seed_point_weight = 20;

% set baseline fiber quantifier options
fq_options.dwi_res = [dim_param_dti_1.voxel_dim(1)*144 144 dim_param_dti_1.voxel_dim(3)];
fq_options.filt_kernel = 3;
fq_options.mesh_units = 'vx';
fq_options.tract_units = 'vx';

fv_options.fiber_color = [1 0 1];

allapo_rstrct_ft_options = ft_options;
allapo_rstrct_ft_options.seed_method = 'apo';

allapo_rstrct_fq_options = fq_options;
allapo_rstrct_fq_options.seed_method = 'apo';

%% Fiber-tracking - APO3 - Aniso4D + Smoothn

% track, smooth, and quantify fibers
[allapo_rstrct_fiber_all, ~, allapo_rstrct_stop_list_aniso_smoothn, ~] = ...
    fiber_track_v20(allapo_rstrct_ft_options, e1fa_aniso_smooth, muscle_mask, roi_mesh_rstrct);
[allapo_rstrct_smoothed_fiber_all, allapo_rstrct_fiber_all_mm, allapo_rstrct_smoothed_fiber_all_mm, ~, ~, ~, ~, ~, ~] = ...
    fiber_smoother_v14a(allapo_rstrct_fiber_all, fs_options);

% Filter fiber-tracts starting in the same voxel
allapo_rstrct_seeds_mm_rounded = roi_mesh_rstrct;
allapo_rstrct_seeds_mm_rounded(:,:,1) = allapo_rstrct_seeds_mm_rounded(:,:,1)*(fv_options.dti_dims(1)/fv_options.dti_size(1));
allapo_rstrct_seeds_mm_rounded(:,:,2) = allapo_rstrct_seeds_mm_rounded(:,:,2)*(fv_options.dti_dims(1)/fv_options.dti_size(1));
allapo_rstrct_seeds_mm_rounded(:,:,3) = allapo_rstrct_seeds_mm_rounded(:,:,3)*fv_options.dti_dims(2);
allapo_rstrct_seeds_mm_rounded = round(allapo_rstrct_seeds_mm_rounded);

% Loop to find indices of tracts with identical voxel coordinates for seed points
allapo_rstrct_smoothed_fiber_all_mm_final = zeros(size(allapo_rstrct_smoothed_fiber_all_mm));
allapo_rstrct_smoothed_fiber_all_final = zeros(size(allapo_rstrct_smoothed_fiber_all));
allapo_rstrct_seeds_final = zeros(10,10,5);

n=1;
for r=min(min(allapo_rstrct_seeds_mm_rounded(:,:,1))):max(max(allapo_rstrct_seeds_mm_rounded(:,:,1)))
    for c=min(min(allapo_rstrct_seeds_mm_rounded(:,:,2))):max(max(allapo_rstrct_seeds_mm_rounded(:,:,2)))
        for s=min(min(allapo_rstrct_seeds_mm_rounded(:,:,3))):max(max(allapo_rstrct_seeds_mm_rounded(:,:,3)))
            idx = find(allapo_rstrct_seeds_mm_rounded(:,:,1)==r...
                & allapo_rstrct_seeds_mm_rounded(:,:,2)==c...
                & allapo_rstrct_seeds_mm_rounded(:,:,3)==s);
            if isempty(idx)
                continue
            else
                [loop_r, loop_c] = ind2sub(size(allapo_rstrct_seeds_mm_rounded(:,:,1)), idx);
                randm_v = randperm(length(loop_r));
                loop_r = loop_r(randm_v);
                loop_c = loop_c(randm_v);

                for k=1:length(loop_r)
                    if allapo_rstrct_stop_list_aniso_smoothn(loop_r(k), loop_c(k))==4 && ...
                                length(nonzeros(allapo_rstrct_fiber_all_mm(loop_r(k),loop_c(k),:,1)))>7 && ...
                                length(nonzeros(allapo_rstrct_fiber_all_mm(loop_r(k),loop_c(k),:,1)))<80
                        allapo_rstrct_seeds_final(loop_r(k),loop_c(k),:) = [r c s loop_r(k) loop_c(k)];
                        allapo_rstrct_smoothed_fiber_all_mm_final(loop_r(k),loop_c(k),:,1:3) = ...
                                allapo_rstrct_smoothed_fiber_all_mm(loop_r(k),loop_c(k),:,:);
                        allapo_rstrct_smoothed_fiber_all_final(loop_r(k),loop_c(k),:,1:3) = ...
                                allapo_rstrct_smoothed_fiber_all(loop_r(k),loop_c(k),:,:);
                        break;
                    end
                end

                n=n+1;
            end
        end
    end
end

[allapo_rstrct_angle_list_final, allapo_rstrct_distance_list_final, allapo_rstrct_curvature_list_final, ~, allapo_rstrct_n_points_final, allapo_rstrct_area_final] = ...
    fiber_quantifier_v20a(allapo_rstrct_fq_options, allapo_rstrct_smoothed_fiber_all_final, roi_mesh_rstrct, muscle_mask);

% calculate summary statistics
allapo_rstrct_total_distance_final = squeeze(max(allapo_rstrct_distance_list_final, [], 3));
allapo_rstrct_total_distance_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final = sum(squeeze(allapo_rstrct_angle_list_final(:,:,:,1)), 3)./allapo_rstrct_n_points_final(:,:,2);
allapo_rstrct_mean_alpha_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_curvature_final = sum(allapo_rstrct_curvature_list_final, 3)./allapo_rstrct_n_points_final(:,:,3);
allapo_rstrct_mean_curvature_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final(isnan(allapo_rstrct_mean_alpha_final)) = 0;
allapo_rstrct_mean_alpha_final = nonzeros(allapo_rstrct_mean_alpha_final);
allapo_rstrct_mean_curvature_final(isnan(allapo_rstrct_mean_curvature_final))=0;
allapo_rstrct_mean_curvature_final = nonzeros(allapo_rstrct_mean_curvature_final);
allapo_rstrct_total_distance_final(isnan(allapo_rstrct_total_distance_final))=0;
allapo_rstrct_total_distance_final = nonzeros(allapo_rstrct_total_distance_final);

mean_fiber_length_aniso_smoothn = mean(allapo_rstrct_total_distance_final);
mean_pennation_angle_aniso_smoothn = mean(allapo_rstrct_mean_alpha_final);
mean_curvature_aniso_smoothn = mean(allapo_rstrct_mean_curvature_final);

fiber_visualizer_v11(muscle_mask, fv_options, roi_mesh_rstrct, anat_img, allapo_rstrct_smoothed_fiber_all_final);
set(gcf, 'position', [1200 50 500 800])
set(gca, 'fontsize', 14)
xlabel('Column', 'fontsize', 14)
ylabel('Row', 'fontsize', 14)
zlabel('Slice', 'fontsize', 14)
title('Aponeurosis Seeding-Restricted Mesh, Aniso + Smoothn')

allapo_rstrct_seeds_final_mask = allapo_rstrct_seeds_final(:,:,4); % Create mask of apo_seeds_final to conserve same seed points on other datasets
allapo_rstrct_seeds_final_mask (allapo_rstrct_seeds_final_mask > 0) = 1;

stop_pct_mask_aniso_smoothn = numel(find(allapo_rstrct_stop_list_aniso_smoothn == 4))/numel(allapo_rstrct_stop_list_aniso_smoothn);
stop_pct_fa_aniso_smoothn = numel(find(allapo_rstrct_stop_list_aniso_smoothn == 2))/numel(allapo_rstrct_stop_list_aniso_smoothn);
stop_pct_ang_aniso_smoothn = numel(find(allapo_rstrct_stop_list_aniso_smoothn == 3))/numel(allapo_rstrct_stop_list_aniso_smoothn);

%% Fiber-tracking - APO3 - Raw

% track, smooth, and quantify fibers
[allapo_rstrct_fiber_all, allapo_rstrct_roi_flag, allapo_rstrct_stop_list_raw, allapo_rstrct_fiber_len] = ...
    fiber_track_v20(allapo_rstrct_ft_options, e1fa_raw, muscle_mask, roi_mesh_rstrct);
[allapo_rstrct_smoothed_fiber_all, ~, allapo_rstrct_smoothed_fiber_all_mm, allapo_rstrct_pcoeff_r, allapo_rstrct_pcoeff_c, allapo_rstrct_pcoeff_s, allapo_rstrct_n_points_smoothed, allapo_rstrct_residuals, allapo_rstrct_residuals_mm] = ...
    fiber_smoother_v14a(allapo_rstrct_fiber_all, fs_options);

% Use same fiber-tract seed points in all the conditions evaluated
allapo_rstrct_smoothed_fiber_all_mm_final = zeros(size(allapo_rstrct_smoothed_fiber_all_mm));   %will hold tracts after seed-location filtering
allapo_rstrct_smoothed_fiber_all_final = zeros(size(allapo_rstrct_smoothed_fiber_all));

for step_cntr = 1:size(allapo_rstrct_smoothed_fiber_all,3)
    for dim_cntr = 1:size(allapo_rstrct_smoothed_fiber_all,4)
            allapo_rstrct_smoothed_fiber_all_final(:,:,step_cntr,dim_cntr) = ...
                squeeze(allapo_rstrct_smoothed_fiber_all(:,:,step_cntr,dim_cntr)).*allapo_rstrct_seeds_final_mask;
            allapo_rstrct_smoothed_fiber_all_mm_final (:,:,step_cntr,dim_cntr) = ...
                squeeze(allapo_rstrct_smoothed_fiber_all_mm(:,:,step_cntr,dim_cntr)).*allapo_rstrct_seeds_final_mask;
    end
end

[allapo_rstrct_angle_list_final, allapo_rstrct_distance_list_final, allapo_rstrct_curvature_list_final, ~, allapo_rstrct_n_points_final, allapo_rstrct_area_final] = ...
    fiber_quantifier_v20a(allapo_rstrct_fq_options, allapo_rstrct_smoothed_fiber_all_final, roi_mesh_rstrct, muscle_mask);

% calculate summary statistics
allapo_rstrct_total_distance_final = squeeze(max(allapo_rstrct_distance_list_final, [], 3));
allapo_rstrct_total_distance_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final = sum(squeeze(allapo_rstrct_angle_list_final(:,:,:,1)), 3)./allapo_rstrct_n_points_final(:,:,2);
allapo_rstrct_mean_alpha_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_curvature_final = sum(allapo_rstrct_curvature_list_final, 3)./allapo_rstrct_n_points_final(:,:,3);
allapo_rstrct_mean_curvature_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final(isnan(allapo_rstrct_mean_alpha_final)) = 0;
allapo_rstrct_mean_alpha_final = nonzeros(allapo_rstrct_mean_alpha_final);
allapo_rstrct_mean_curvature_final(isnan(allapo_rstrct_mean_curvature_final))=0;
allapo_rstrct_mean_curvature_final = nonzeros(allapo_rstrct_mean_curvature_final);
allapo_rstrct_total_distance_final(isnan(allapo_rstrct_total_distance_final))=0;
allapo_rstrct_total_distance_final = nonzeros(allapo_rstrct_total_distance_final);

mean_fiber_length_raw = mean(allapo_rstrct_total_distance_final);
mean_pennation_angle_raw = mean(allapo_rstrct_mean_alpha_final);
mean_curvature_raw = mean(allapo_rstrct_mean_curvature_final);

fiber_visualizer_v11(muscle_mask, fv_options, roi_mesh_rstrct, anat_img, allapo_rstrct_smoothed_fiber_all_final);
set(gcf, 'position', [1200 50 500 800])
set(gca, 'fontsize', 14)
xlabel('Column', 'fontsize', 14)
ylabel('Row', 'fontsize', 14)
zlabel('Slice', 'fontsize', 14)
title('Aponeurosis Seeding-Restricted Mesh, Raw')

stop_pct_mask_raw = numel(find(allapo_rstrct_stop_list_raw == 4))/numel(allapo_rstrct_stop_list_raw);
stop_pct_fa_raw = numel(find(allapo_rstrct_stop_list_raw == 2))/numel(allapo_rstrct_stop_list_raw);
stop_pct_ang_raw = numel(find(allapo_rstrct_stop_list_raw == 3))/numel(allapo_rstrct_stop_list_raw);

%% Fiber-tracking - APO3 - Smoothn

% track, smooth, and quantify fibers
[allapo_rstrct_fiber_all, ~, allapo_rstrct_stop_list_smoothn, ~] = ...
    fiber_track_v20(allapo_rstrct_ft_options, e1fa_raw_smooth, muscle_mask, roi_mesh_rstrct);
[allapo_rstrct_smoothed_fiber_all, ~, allapo_rstrct_smoothed_fiber_all_mm, ~, ~, ~, ~, ~, ~] = ...
    fiber_smoother_v14a(allapo_rstrct_fiber_all, fs_options);

% Use same fiber-tract seed points in all the conditions evaluated
allapo_rstrct_smoothed_fiber_all_mm_final = zeros(size(allapo_rstrct_smoothed_fiber_all_mm));   %will hold tracts after seed-location filtering
allapo_rstrct_smoothed_fiber_all_final = zeros(size(allapo_rstrct_smoothed_fiber_all));

for step_cntr = 1:size(allapo_rstrct_smoothed_fiber_all,3)
    for dim_cntr = 1:size(allapo_rstrct_smoothed_fiber_all,4)
            allapo_rstrct_smoothed_fiber_all_final(:,:,step_cntr,dim_cntr) = ...
                squeeze(allapo_rstrct_smoothed_fiber_all(:,:,step_cntr,dim_cntr)).*allapo_rstrct_seeds_final_mask;
            allapo_rstrct_smoothed_fiber_all_mm_final (:,:,step_cntr,dim_cntr) = ...
                squeeze(allapo_rstrct_smoothed_fiber_all_mm(:,:,step_cntr,dim_cntr)).*allapo_rstrct_seeds_final_mask;
    end
end

[allapo_rstrct_angle_list_final, allapo_rstrct_distance_list_final, allapo_rstrct_curvature_list_final, ~, allapo_rstrct_n_points_final, allapo_rstrct_area_final] = ...
    fiber_quantifier_v20a(allapo_rstrct_fq_options, allapo_rstrct_smoothed_fiber_all_final, roi_mesh_rstrct, muscle_mask);

% calculate summary statistics
allapo_rstrct_total_distance_final = squeeze(max(allapo_rstrct_distance_list_final, [], 3));
allapo_rstrct_total_distance_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final = sum(squeeze(allapo_rstrct_angle_list_final(:,:,:,1)), 3)./allapo_rstrct_n_points_final(:,:,2);
allapo_rstrct_mean_alpha_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_curvature_final = sum(allapo_rstrct_curvature_list_final, 3)./allapo_rstrct_n_points_final(:,:,3);
allapo_rstrct_mean_curvature_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final(isnan(allapo_rstrct_mean_alpha_final)) = 0;
allapo_rstrct_mean_alpha_final = nonzeros(allapo_rstrct_mean_alpha_final);
allapo_rstrct_mean_curvature_final(isnan(allapo_rstrct_mean_curvature_final))=0;
allapo_rstrct_mean_curvature_final = nonzeros(allapo_rstrct_mean_curvature_final);
allapo_rstrct_total_distance_final(isnan(allapo_rstrct_total_distance_final))=0;
allapo_rstrct_total_distance_final = nonzeros(allapo_rstrct_total_distance_final);

mean_fiber_length_smoothn = mean(allapo_rstrct_total_distance_final);
mean_pennation_angle_smoothn = mean(allapo_rstrct_mean_alpha_final);
mean_curvature_smoothn = mean(allapo_rstrct_mean_curvature_final);

fiber_visualizer_v11(muscle_mask, fv_options, roi_mesh_rstrct, anat_img, allapo_rstrct_smoothed_fiber_all_final);
set(gcf, 'position', [1200 50 500 800])
set(gca, 'fontsize', 14)
xlabel('Column', 'fontsize', 14)
ylabel('Row', 'fontsize', 14)
zlabel('Slice', 'fontsize', 14)
title('Aponeurosis Seeding-Restricted Mesh, Smoothn')

stop_pct_mask_smoothn = numel(find(allapo_rstrct_stop_list_smoothn == 4))/numel(allapo_rstrct_stop_list_smoothn);
stop_pct_fa_smoothn = numel(find(allapo_rstrct_stop_list_smoothn == 2))/numel(allapo_rstrct_stop_list_smoothn);
stop_pct_ang_smoothn = numel(find(allapo_rstrct_stop_list_smoothn == 3))/numel(allapo_rstrct_stop_list_smoothn);

%% Fiber-tracking - APO3 - Aniso

% track, smooth, and quantify fibers
[allapo_rstrct_fiber_all, ~, allapo_rstrct_stop_list_aniso, ~] = ...
    fiber_track_v20(allapo_rstrct_ft_options, e1fa_aniso, muscle_mask, roi_mesh_rstrct);
[allapo_rstrct_smoothed_fiber_all, ~, allapo_rstrct_smoothed_fiber_all_mm, ~, ~, ~, ~, ~, ~] = ...
    fiber_smoother_v14a(allapo_rstrct_fiber_all, fs_options);

% Use same fiber-tract seed points in all the conditions evaluated
allapo_rstrct_smoothed_fiber_all_mm_final = zeros(size(allapo_rstrct_smoothed_fiber_all_mm));   %will hold tracts after seed-location filtering
allapo_rstrct_smoothed_fiber_all_final = zeros(size(allapo_rstrct_smoothed_fiber_all));

for step_cntr = 1:size(allapo_rstrct_smoothed_fiber_all,3)
    for dim_cntr = 1:size(allapo_rstrct_smoothed_fiber_all,4)
            allapo_rstrct_smoothed_fiber_all_final(:,:,step_cntr,dim_cntr) = ...
                squeeze(allapo_rstrct_smoothed_fiber_all(:,:,step_cntr,dim_cntr)).*allapo_rstrct_seeds_final_mask;
            allapo_rstrct_smoothed_fiber_all_mm_final (:,:,step_cntr,dim_cntr) = ...
                squeeze(allapo_rstrct_smoothed_fiber_all_mm(:,:,step_cntr,dim_cntr)).*allapo_rstrct_seeds_final_mask;
    end
end

[allapo_rstrct_angle_list_final, allapo_rstrct_distance_list_final, allapo_rstrct_curvature_list_final, ~, allapo_rstrct_n_points_final, allapo_rstrct_area_final] = ...
    fiber_quantifier_v20a(allapo_rstrct_fq_options, allapo_rstrct_smoothed_fiber_all_final, roi_mesh_rstrct, muscle_mask);

% calculate summary statistics
allapo_rstrct_total_distance_final = squeeze(max(allapo_rstrct_distance_list_final, [], 3));
allapo_rstrct_total_distance_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final = sum(squeeze(allapo_rstrct_angle_list_final(:,:,:,1)), 3)./allapo_rstrct_n_points_final(:,:,2);
allapo_rstrct_mean_alpha_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_curvature_final = sum(allapo_rstrct_curvature_list_final, 3)./allapo_rstrct_n_points_final(:,:,3);
allapo_rstrct_mean_curvature_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final(isnan(allapo_rstrct_mean_alpha_final)) = 0;
allapo_rstrct_mean_alpha_final = nonzeros(allapo_rstrct_mean_alpha_final);
allapo_rstrct_mean_curvature_final(isnan(allapo_rstrct_mean_curvature_final))=0;
allapo_rstrct_mean_curvature_final = nonzeros(allapo_rstrct_mean_curvature_final);
allapo_rstrct_total_distance_final(isnan(allapo_rstrct_total_distance_final))=0;
allapo_rstrct_total_distance_final = nonzeros(allapo_rstrct_total_distance_final);

mean_fiber_length_aniso = mean(allapo_rstrct_total_distance_final)
mean_pennation_angle_aniso = mean(allapo_rstrct_mean_alpha_final)
mean_curvature_aniso = mean(allapo_rstrct_mean_curvature_final)

fiber_visualizer_v11(muscle_mask, fv_options, roi_mesh_rstrct, anat_img, allapo_rstrct_smoothed_fiber_all_final);
set(gcf, 'position', [1200 50 500 800])
set(gca, 'fontsize', 14)
xlabel('Column', 'fontsize', 14)
ylabel('Row', 'fontsize', 14)
zlabel('Slice', 'fontsize', 14)
title('Aponeurosis Seeding-Restricted Mesh, Aniso')

stop_pct_mask_aniso = numel(find(allapo_rstrct_stop_list_aniso == 4))/numel(allapo_rstrct_stop_list_aniso)
stop_pct_fa_aniso = numel(find(allapo_rstrct_stop_list_aniso == 2))/numel(allapo_rstrct_stop_list_aniso)
stop_pct_ang_aniso = numel(find(allapo_rstrct_stop_list_aniso == 3))/numel(allapo_rstrct_stop_list_aniso)

%% Fiber-tracking - APO3 - TPCA

% track, smooth, and quantify fibers
[allapo_rstrct_fiber_all, ~, allapo_rstrct_stop_list_tpca, ~] = ...
    fiber_track_v20(allapo_rstrct_ft_options, e1fa_tpca, muscle_mask, roi_mesh_rstrct);
[allapo_rstrct_smoothed_fiber_all, ~, allapo_rstrct_smoothed_fiber_all_mm, ~, ~, ~, ~, ~, ~] = ...
    fiber_smoother_v14a(allapo_rstrct_fiber_all, fs_options);

% Use same fiber-tract seed points in all the conditions evaluated
allapo_rstrct_smoothed_fiber_all_mm_final = zeros(size(allapo_rstrct_smoothed_fiber_all_mm));   %will hold tracts after seed-location filtering
allapo_rstrct_smoothed_fiber_all_final = zeros(size(allapo_rstrct_smoothed_fiber_all));

for step_cntr = 1:size(allapo_rstrct_smoothed_fiber_all,3)
    for dim_cntr = 1:size(allapo_rstrct_smoothed_fiber_all,4)
            allapo_rstrct_smoothed_fiber_all_final(:,:,step_cntr,dim_cntr) = ...
                squeeze(allapo_rstrct_smoothed_fiber_all(:,:,step_cntr,dim_cntr)).*allapo_rstrct_seeds_final_mask;
            allapo_rstrct_smoothed_fiber_all_mm_final (:,:,step_cntr,dim_cntr) = ...
                squeeze(allapo_rstrct_smoothed_fiber_all_mm(:,:,step_cntr,dim_cntr)).*allapo_rstrct_seeds_final_mask;
    end
end

[allapo_rstrct_angle_list_final, allapo_rstrct_distance_list_final, allapo_rstrct_curvature_list_final, ~, allapo_rstrct_n_points_final, allapo_rstrct_area_final] = ...
    fiber_quantifier_v20a(allapo_rstrct_fq_options, allapo_rstrct_smoothed_fiber_all_final, roi_mesh_rstrct, muscle_mask);

% calculate summary statistics
allapo_rstrct_total_distance_final = squeeze(max(allapo_rstrct_distance_list_final, [], 3));
allapo_rstrct_total_distance_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final = sum(squeeze(allapo_rstrct_angle_list_final(:,:,:,1)), 3)./allapo_rstrct_n_points_final(:,:,2);
allapo_rstrct_mean_alpha_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_curvature_final = sum(allapo_rstrct_curvature_list_final, 3)./allapo_rstrct_n_points_final(:,:,3);
allapo_rstrct_mean_curvature_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final(isnan(allapo_rstrct_mean_alpha_final)) = 0;
allapo_rstrct_mean_alpha_final = nonzeros(allapo_rstrct_mean_alpha_final);
allapo_rstrct_mean_curvature_final(isnan(allapo_rstrct_mean_curvature_final))=0;
allapo_rstrct_mean_curvature_final = nonzeros(allapo_rstrct_mean_curvature_final);
allapo_rstrct_total_distance_final(isnan(allapo_rstrct_total_distance_final))=0;
allapo_rstrct_total_distance_final = nonzeros(allapo_rstrct_total_distance_final);

mean_fiber_length_tpca = mean(allapo_rstrct_total_distance_final)
mean_pennation_angle_tpca = mean(allapo_rstrct_mean_alpha_final)
mean_curvature_tpca = mean(allapo_rstrct_mean_curvature_final)

fiber_visualizer_v11(muscle_mask, fv_options, roi_mesh_rstrct, anat_img, allapo_rstrct_smoothed_fiber_all_final);
set(gcf, 'position', [1200 50 500 800])
set(gca, 'fontsize', 14)
xlabel('Column', 'fontsize', 14)
ylabel('Row', 'fontsize', 14)
zlabel('Slice', 'fontsize', 14)
title('Aponeurosis Seeding-Restricted Mesh, TPCA')

stop_pct_mask_tpca = numel(find(allapo_rstrct_stop_list_tpca == 4))/numel(allapo_rstrct_stop_list_tpca)
stop_pct_fa_tpca = numel(find(allapo_rstrct_stop_list_tpca == 2))/numel(allapo_rstrct_stop_list_tpca)
stop_pct_ang_tpca = numel(find(allapo_rstrct_stop_list_tpca == 3))/numel(allapo_rstrct_stop_list_tpca)

%% Fiber-tracking - APO3 - TPCA + smoothn

% track, smooth, and quantify fibers
[allapo_rstrct_fiber_all, ~, allapo_rstrct_stop_list_tpca_smooth, ~] = ...
    fiber_track_v20(allapo_rstrct_ft_options, e1fa_tpca_smooth, muscle_mask, roi_mesh_rstrct);
[allapo_rstrct_smoothed_fiber_all, ~, allapo_rstrct_smoothed_fiber_all_mm, ~, ~, ~, ~, ~, ~] = ...
    fiber_smoother_v14a(allapo_rstrct_fiber_all, fs_options);

% Use same fiber-tract seed points in all the conditions evaluated
allapo_rstrct_smoothed_fiber_all_mm_final = zeros(size(allapo_rstrct_smoothed_fiber_all_mm));   %will hold tracts after seed-location filtering
allapo_rstrct_smoothed_fiber_all_final = zeros(size(allapo_rstrct_smoothed_fiber_all));

for step_cntr = 1:size(allapo_rstrct_smoothed_fiber_all,3)
    for dim_cntr = 1:size(allapo_rstrct_smoothed_fiber_all,4)
            allapo_rstrct_smoothed_fiber_all_final(:,:,step_cntr,dim_cntr) = ...
                squeeze(allapo_rstrct_smoothed_fiber_all(:,:,step_cntr,dim_cntr)).*allapo_rstrct_seeds_final_mask;
            allapo_rstrct_smoothed_fiber_all_mm_final (:,:,step_cntr,dim_cntr) = ...
                squeeze(allapo_rstrct_smoothed_fiber_all_mm(:,:,step_cntr,dim_cntr)).*allapo_rstrct_seeds_final_mask;
    end
end

[allapo_rstrct_angle_list_final, allapo_rstrct_distance_list_final, allapo_rstrct_curvature_list_final, ~, allapo_rstrct_n_points_final, allapo_rstrct_area_final] = ...
    fiber_quantifier_v20a(allapo_rstrct_fq_options, allapo_rstrct_smoothed_fiber_all_final, roi_mesh_rstrct, muscle_mask);

% calculate summary statistics
allapo_rstrct_total_distance_final = squeeze(max(allapo_rstrct_distance_list_final, [], 3));
allapo_rstrct_total_distance_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final = sum(squeeze(allapo_rstrct_angle_list_final(:,:,:,1)), 3)./allapo_rstrct_n_points_final(:,:,2);
allapo_rstrct_mean_alpha_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_curvature_final = sum(allapo_rstrct_curvature_list_final, 3)./allapo_rstrct_n_points_final(:,:,3);
allapo_rstrct_mean_curvature_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final(isnan(allapo_rstrct_mean_alpha_final)) = 0;
allapo_rstrct_mean_alpha_final = nonzeros(allapo_rstrct_mean_alpha_final);
allapo_rstrct_mean_curvature_final(isnan(allapo_rstrct_mean_curvature_final))=0;
allapo_rstrct_mean_curvature_final = nonzeros(allapo_rstrct_mean_curvature_final);
allapo_rstrct_total_distance_final(isnan(allapo_rstrct_total_distance_final))=0;
allapo_rstrct_total_distance_final = nonzeros(allapo_rstrct_total_distance_final);

mean_fiber_length_tpca_smooth = mean(allapo_rstrct_total_distance_final)
mean_pennation_angle_tpca_smooth = mean(allapo_rstrct_mean_alpha_final)
mean_curvature_tpca_smooth = mean(allapo_rstrct_mean_curvature_final)

fiber_visualizer_v11(muscle_mask, fv_options, roi_mesh_rstrct, anat_img, allapo_rstrct_smoothed_fiber_all_final);
set(gcf, 'position', [1200 50 500 800])
set(gca, 'fontsize', 14)
xlabel('Column', 'fontsize', 14)
ylabel('Row', 'fontsize', 14)
zlabel('Slice', 'fontsize', 14)
title('Aponeurosis Seeding-Restricted Mesh, TPCA + smoothn')

stop_pct_mask_tpca_smooth = numel(find(allapo_rstrct_stop_list_tpca_smooth == 4))/numel(allapo_rstrct_stop_list_tpca_smooth)
stop_pct_fa_tpca_smooth = numel(find(allapo_rstrct_stop_list_tpca_smooth == 2))/numel(allapo_rstrct_stop_list_tpca_smooth)
stop_pct_ang_tpca_smooth = numel(find(allapo_rstrct_stop_list_tpca_smooth == 3))/numel(allapo_rstrct_stop_list_tpca_smooth)

%% Analyze b = 0 (4 avg), b = 450 (1 avg) data

%% 
% Compute FA and first eigenvector map of DTI data

%Correct diffusion directions fror FFS scans
diff_dir=diff_param_1.bvec; %diff_param1 and diff_param2 are the same

b_val=mean(diff_param_1.bval(2:end));

%Tensor computation

FA_raw = zeros(size(dti_img_4d_reg_b0_4avg,1),size(dti_img_4d_reg_b0_4avg,2),size(dti_img_4d_reg_b0_4avg,3));
E1_450_raw = zeros(size(dti_img_4d_reg_b0_4avg,1),size(dti_img_4d_reg_b0_4avg,2),size(dti_img_4d_reg_b0_4avg,3),3);
FA_aniso = zeros(size(dti_img_4d_reg_b0_4avg,1),size(dti_img_4d_reg_b0_4avg,2),size(dti_img_4d_reg_b0_4avg,3));
E1_450_aniso = zeros(size(dti_img_4d_reg_b0_4avg,1),size(dti_img_4d_reg_b0_4avg,2),size(dti_img_4d_reg_b0_4avg,3),3);
FA_tpca = zeros(size(dti_img_4d_reg_b0_4avg,1),size(dti_img_4d_reg_b0_4avg,2),size(dti_img_4d_reg_b0_4avg,3));
E1_450_tpca = zeros(size(dti_img_4d_reg_b0_4avg,1),size(dti_img_4d_reg_b0_4avg,2),size(dti_img_4d_reg_b0_4avg,3),3);

for s = 1:size(dti_img_4d_reg_b0_4avg,3)
    for r = 1:size(dti_img_4d_reg_b0_4avg,1)
        for c = 1:size(dti_img_4d_reg_b0_4avg,2)
                dir_m = diff_dir;  %diffusion encoding matrix
                dir_m = dir_m(2:end,:);

                %estimate tensor
                %use TA mask to make process faster
                if muscle_mask(r,c,s)>0
                    % Raw data
                    sig_v = squeeze(dti_img_4d_reg_b0_4avg(r,c,s,:));    %signal vector
                    D = signal2tensor2(sig_v, dir_m, b_val);
                    [E,L] = svd(D);
    
                    L_v = diag(L);
                    md = mean(L_v);
                    fa = sqrt(3/2) * sqrt(((L_v(1) - md)^2 + (L_v(2) - md)^2 + (L_v(3) - md)^2) / ...
                        (L_v(1)^2 + L_v(2)^2 + L_v(3)^2));

                    loop_E1 = E(:,L_v==max(L_v));
                    if loop_E1(3)<0
                        loop_E1 = -loop_E1;
                    end
                    
                    E1_450_raw(r,c,s,:) = loop_E1;
                    FA_raw(r,c,s,:) = fa;

                    % AIS
                    sig_v = squeeze(dti_img_4d_smooth_b0_4avg(r,c,s,:));    %signal vector
                    D = signal2tensor2(sig_v, dir_m, b_val);
                    [E,L] = svd(D);
    
                    L_v = diag(L);
                    md = mean(L_v);
                    fa = sqrt(3/2) * sqrt(((L_v(1) - md)^2 + (L_v(2) - md)^2 + (L_v(3) - md)^2) / ...
                        (L_v(1)^2 + L_v(2)^2 + L_v(3)^2));

                    loop_E1 = E(:,L_v==max(L_v));
                    if loop_E1(3)<0
                        loop_E1 = -loop_E1;
                    end
                    
                    E1_450_aniso(r,c,s,:) = loop_E1;
                    FA_aniso(r,c,s,:) = fa;

                    % TPCA
                    sig_v = squeeze(dti_img_4d_tpca_b0_4avg(r,c,s,:));    %signal vector
                    D = signal2tensor2(sig_v, dir_m, b_val);
                    [E,L] = svd(D);
    
                    L_v = diag(L);
                    md = mean(L_v);
                    fa = sqrt(3/2) * sqrt(((L_v(1) - md)^2 + (L_v(2) - md)^2 + (L_v(3) - md)^2) / ...
                        (L_v(1)^2 + L_v(2)^2 + L_v(3)^2));

                    loop_E1 = E(:,L_v==max(L_v));
                    if loop_E1(3)<0
                        loop_E1 = -loop_E1;
                    end
                    
                    E1_450_tpca(r,c,s,:) = loop_E1;
                    FA_tpca(r,c,s,:) = fa;

                end
        end
    end
end

%% 
% Smooth first eigenvector field of TA muscle

[E1_450_raw_smooth,smooth_parameter_raw]=smoothn({E1_450_raw(:,:,:,1),E1_450_raw(:,:,:,2),E1_450_raw(:,:,:,3)});
[E1_450_aniso_smooth,smooth_parameter_aniso]=smoothn({E1_450_aniso(:,:,:,1),E1_450_aniso(:,:,:,2),E1_450_aniso(:,:,:,3)});
[E1_450_tpca_smooth,smooth_parameter_tpca]=smoothn({E1_450_tpca(:,:,:,1),E1_450_tpca(:,:,:,2),E1_450_tpca(:,:,:,3)});

E1_450_raw_smooth=cat(4,E1_450_raw_smooth{1},E1_450_raw_smooth{2},E1_450_raw_smooth{3});
E1_450_aniso_smooth=cat(4,E1_450_aniso_smooth{1},E1_450_aniso_smooth{2},E1_450_aniso_smooth{3});
E1_450_tpca_smooth=cat(4,E1_450_tpca_smooth{1},E1_450_tpca_smooth{2},E1_450_tpca_smooth{3});

for r_img=1:size(E1_450_aniso,1)
    for c_img=1:size(E1_450_aniso,2)
        for s_img=1:size(E1_450_aniso,3)
            if muscle_mask(r_img,c_img,s_img)>0

                E1_raw_smooth = squeeze(E1_450_raw_smooth(r_img,c_img,s_img,:));
                E1_aniso_smooth = squeeze(E1_450_aniso_smooth(r_img,c_img,s_img,:));
                E1_tpca_smooth = squeeze(E1_450_tpca_smooth(r_img,c_img,s_img,:));

                E1_raw_smooth = E1_raw_smooth/norm(E1_raw_smooth);
                E1_aniso_smooth = E1_aniso_smooth/norm(E1_aniso_smooth);
                E1_tpca_smooth = E1_tpca_smooth/norm(E1_tpca_smooth);

                E1_450_raw_smooth(r_img,c_img,s_img,:)=E1_raw_smooth;
                E1_450_aniso_smooth(r_img,c_img,s_img,:)=E1_aniso_smooth;
                E1_450_tpca_smooth(r_img,c_img,s_img,:)=E1_tpca_smooth;
            end
        end
    end
end

for i = 1:3
    E1_450_raw_smooth(:,:,:,i)=squeeze(E1_450_raw_smooth(:,:,:,i)).*muscle_mask;
    E1_450_aniso_smooth(:,:,:,i)=squeeze(E1_450_aniso_smooth(:,:,:,i)).*muscle_mask;
    E1_450_tpca_smooth(:,:,:,i)=squeeze(E1_450_tpca_smooth(:,:,:,i)).*muscle_mask;
end

%% 
% Create E1-FA maps for fiber tractography

e1fa_raw = cat(4,E1_450_raw,FA_raw);
e1fa_aniso = cat(4,E1_450_aniso,FA_aniso);
e1fa_tpca = cat(4,E1_450_tpca,FA_tpca);
e1fa_raw_smooth = cat(4,E1_450_raw_smooth,FA_raw);
e1fa_aniso_smooth = cat(4,E1_450_aniso_smooth,FA_aniso);
e1fa_tpca_smooth = cat(4,E1_450_tpca_smooth,FA_tpca);

%% Fiber-tracking - APO3 - Aniso4D + Smoothn

% track, smooth, and quantify fibers
[allapo_rstrct_fiber_all, ~, allapo_rstrct_stop_list_aniso_smoothn, ~] = ...
    fiber_track_v20(allapo_rstrct_ft_options, e1fa_aniso_smooth, muscle_mask, roi_mesh_rstrct);
[allapo_rstrct_smoothed_fiber_all, allapo_rstrct_fiber_all_mm, allapo_rstrct_smoothed_fiber_all_mm, ~, ~, ~, ~, ~, ~] = ...
    fiber_smoother_v14a(allapo_rstrct_fiber_all, fs_options);

% Use same fiber-tract seed points in all the conditions evaluated
allapo_rstrct_smoothed_fiber_all_mm_final = zeros(size(allapo_rstrct_smoothed_fiber_all_mm));   %will hold tracts after seed-location filtering
allapo_rstrct_smoothed_fiber_all_final = zeros(size(allapo_rstrct_smoothed_fiber_all));

for step_cntr = 1:size(allapo_rstrct_smoothed_fiber_all,3)
    for dim_cntr = 1:size(allapo_rstrct_smoothed_fiber_all,4)
            allapo_rstrct_smoothed_fiber_all_final(:,:,step_cntr,dim_cntr) = ...
                squeeze(allapo_rstrct_smoothed_fiber_all(:,:,step_cntr,dim_cntr)).*allapo_rstrct_seeds_final_mask;
            allapo_rstrct_smoothed_fiber_all_mm_final (:,:,step_cntr,dim_cntr) = ...
                squeeze(allapo_rstrct_smoothed_fiber_all_mm(:,:,step_cntr,dim_cntr)).*allapo_rstrct_seeds_final_mask;
    end
end

[allapo_rstrct_angle_list_final, allapo_rstrct_distance_list_final, allapo_rstrct_curvature_list_final, ~, allapo_rstrct_n_points_final, allapo_rstrct_area_final] = ...
    fiber_quantifier_v20a(allapo_rstrct_fq_options, allapo_rstrct_smoothed_fiber_all_final, roi_mesh_rstrct, muscle_mask);

% calculate summary statistics
allapo_rstrct_total_distance_final = squeeze(max(allapo_rstrct_distance_list_final, [], 3));
allapo_rstrct_total_distance_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final = sum(squeeze(allapo_rstrct_angle_list_final(:,:,:,1)), 3)./allapo_rstrct_n_points_final(:,:,2);
allapo_rstrct_mean_alpha_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_curvature_final = sum(allapo_rstrct_curvature_list_final, 3)./allapo_rstrct_n_points_final(:,:,3);
allapo_rstrct_mean_curvature_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final(isnan(allapo_rstrct_mean_alpha_final)) = 0;
allapo_rstrct_mean_alpha_final = nonzeros(allapo_rstrct_mean_alpha_final);
allapo_rstrct_mean_curvature_final(isnan(allapo_rstrct_mean_curvature_final))=0;
allapo_rstrct_mean_curvature_final = nonzeros(allapo_rstrct_mean_curvature_final);
allapo_rstrct_total_distance_final(isnan(allapo_rstrct_total_distance_final))=0;
allapo_rstrct_total_distance_final = nonzeros(allapo_rstrct_total_distance_final);

mean_fiber_length_aniso_smoothn = mean(allapo_rstrct_total_distance_final)
mean_pennation_angle_aniso_smoothn = mean(allapo_rstrct_mean_alpha_final)
mean_curvature_aniso_smoothn = mean(allapo_rstrct_mean_curvature_final)

fiber_visualizer_v11(muscle_mask, fv_options, roi_mesh_rstrct, anat_img, allapo_rstrct_smoothed_fiber_all_final);
set(gcf, 'position', [1200 50 500 800])
set(gca, 'fontsize', 14)
xlabel('Column', 'fontsize', 14)
ylabel('Row', 'fontsize', 14)
zlabel('Slice', 'fontsize', 14)
title('Aponeurosis Seeding-Restricted Mesh, Aniso + Smoothn')

stop_pct_mask_aniso_smoothn = numel(find(allapo_rstrct_stop_list_aniso_smoothn == 4))/numel(allapo_rstrct_stop_list_aniso_smoothn)
stop_pct_fa_aniso_smoothn = numel(find(allapo_rstrct_stop_list_aniso_smoothn == 2))/numel(allapo_rstrct_stop_list_aniso_smoothn)
stop_pct_ang_aniso_smoothn = numel(find(allapo_rstrct_stop_list_aniso_smoothn == 3))/numel(allapo_rstrct_stop_list_aniso_smoothn)

%% Fiber-tracking - APO3 - Raw

% track, smooth, and quantify fibers
[allapo_rstrct_fiber_all, allapo_rstrct_roi_flag, allapo_rstrct_stop_list_raw, allapo_rstrct_fiber_len] = ...
    fiber_track_v20(allapo_rstrct_ft_options, e1fa_raw, muscle_mask, roi_mesh_rstrct);
[allapo_rstrct_smoothed_fiber_all, ~, allapo_rstrct_smoothed_fiber_all_mm, allapo_rstrct_pcoeff_r, allapo_rstrct_pcoeff_c, allapo_rstrct_pcoeff_s, allapo_rstrct_n_points_smoothed, allapo_rstrct_residuals, allapo_rstrct_residuals_mm] = ...
    fiber_smoother_v14a(allapo_rstrct_fiber_all, fs_options);

% Use same fiber-tract seed points in all the conditions evaluated
allapo_rstrct_smoothed_fiber_all_mm_final = zeros(size(allapo_rstrct_smoothed_fiber_all_mm));   %will hold tracts after seed-location filtering
allapo_rstrct_smoothed_fiber_all_final = zeros(size(allapo_rstrct_smoothed_fiber_all));

for step_cntr = 1:size(allapo_rstrct_smoothed_fiber_all,3)
    for dim_cntr = 1:size(allapo_rstrct_smoothed_fiber_all,4)
            allapo_rstrct_smoothed_fiber_all_final(:,:,step_cntr,dim_cntr) = ...
                squeeze(allapo_rstrct_smoothed_fiber_all(:,:,step_cntr,dim_cntr)).*allapo_rstrct_seeds_final_mask;
            allapo_rstrct_smoothed_fiber_all_mm_final (:,:,step_cntr,dim_cntr) = ...
                squeeze(allapo_rstrct_smoothed_fiber_all_mm(:,:,step_cntr,dim_cntr)).*allapo_rstrct_seeds_final_mask;
    end
end

[allapo_rstrct_angle_list_final, allapo_rstrct_distance_list_final, allapo_rstrct_curvature_list_final, ~, allapo_rstrct_n_points_final, allapo_rstrct_area_final] = ...
    fiber_quantifier_v20a(allapo_rstrct_fq_options, allapo_rstrct_smoothed_fiber_all_final, roi_mesh_rstrct, muscle_mask);

% calculate summary statistics
allapo_rstrct_total_distance_final = squeeze(max(allapo_rstrct_distance_list_final, [], 3));
allapo_rstrct_total_distance_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final = sum(squeeze(allapo_rstrct_angle_list_final(:,:,:,1)), 3)./allapo_rstrct_n_points_final(:,:,2);
allapo_rstrct_mean_alpha_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_curvature_final = sum(allapo_rstrct_curvature_list_final, 3)./allapo_rstrct_n_points_final(:,:,3);
allapo_rstrct_mean_curvature_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final(isnan(allapo_rstrct_mean_alpha_final)) = 0;
allapo_rstrct_mean_alpha_final = nonzeros(allapo_rstrct_mean_alpha_final);
allapo_rstrct_mean_curvature_final(isnan(allapo_rstrct_mean_curvature_final))=0;
allapo_rstrct_mean_curvature_final = nonzeros(allapo_rstrct_mean_curvature_final);
allapo_rstrct_total_distance_final(isnan(allapo_rstrct_total_distance_final))=0;
allapo_rstrct_total_distance_final = nonzeros(allapo_rstrct_total_distance_final);

mean_fiber_length_raw = mean(allapo_rstrct_total_distance_final)
mean_pennation_angle_raw = mean(allapo_rstrct_mean_alpha_final)
mean_curvature_raw = mean(allapo_rstrct_mean_curvature_final)

fiber_visualizer_v11(muscle_mask, fv_options, roi_mesh_rstrct, anat_img, allapo_rstrct_smoothed_fiber_all_final);
set(gcf, 'position', [1200 50 500 800])
set(gca, 'fontsize', 14)
xlabel('Column', 'fontsize', 14)
ylabel('Row', 'fontsize', 14)
zlabel('Slice', 'fontsize', 14)
title('Aponeurosis Seeding-Restricted Mesh, Raw')

stop_pct_mask_raw = numel(find(allapo_rstrct_stop_list_raw == 4))/numel(allapo_rstrct_stop_list_raw)
stop_pct_fa_raw = numel(find(allapo_rstrct_stop_list_raw == 2))/numel(allapo_rstrct_stop_list_raw)
stop_pct_ang_raw = numel(find(allapo_rstrct_stop_list_raw == 3))/numel(allapo_rstrct_stop_list_raw)

%% Fiber-tracking - APO3 - Smoothn

% track, smooth, and quantify fibers
[allapo_rstrct_fiber_all, ~, allapo_rstrct_stop_list_smoothn, ~] = ...
    fiber_track_v20(allapo_rstrct_ft_options, e1fa_raw_smooth, muscle_mask, roi_mesh_rstrct);
[allapo_rstrct_smoothed_fiber_all, ~, allapo_rstrct_smoothed_fiber_all_mm, ~, ~, ~, ~, ~, ~] = ...
    fiber_smoother_v14a(allapo_rstrct_fiber_all, fs_options);

% Use same fiber-tract seed points in all the conditions evaluated
allapo_rstrct_smoothed_fiber_all_mm_final = zeros(size(allapo_rstrct_smoothed_fiber_all_mm));   %will hold tracts after seed-location filtering
allapo_rstrct_smoothed_fiber_all_final = zeros(size(allapo_rstrct_smoothed_fiber_all));

for step_cntr = 1:size(allapo_rstrct_smoothed_fiber_all,3)
    for dim_cntr = 1:size(allapo_rstrct_smoothed_fiber_all,4)
            allapo_rstrct_smoothed_fiber_all_final(:,:,step_cntr,dim_cntr) = ...
                squeeze(allapo_rstrct_smoothed_fiber_all(:,:,step_cntr,dim_cntr)).*allapo_rstrct_seeds_final_mask;
            allapo_rstrct_smoothed_fiber_all_mm_final (:,:,step_cntr,dim_cntr) = ...
                squeeze(allapo_rstrct_smoothed_fiber_all_mm(:,:,step_cntr,dim_cntr)).*allapo_rstrct_seeds_final_mask;
    end
end

[allapo_rstrct_angle_list_final, allapo_rstrct_distance_list_final, allapo_rstrct_curvature_list_final, ~, allapo_rstrct_n_points_final, allapo_rstrct_area_final] = ...
    fiber_quantifier_v20a(allapo_rstrct_fq_options, allapo_rstrct_smoothed_fiber_all_final, roi_mesh_rstrct, muscle_mask);

% calculate summary statistics
allapo_rstrct_total_distance_final = squeeze(max(allapo_rstrct_distance_list_final, [], 3));
allapo_rstrct_total_distance_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final = sum(squeeze(allapo_rstrct_angle_list_final(:,:,:,1)), 3)./allapo_rstrct_n_points_final(:,:,2);
allapo_rstrct_mean_alpha_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_curvature_final = sum(allapo_rstrct_curvature_list_final, 3)./allapo_rstrct_n_points_final(:,:,3);
allapo_rstrct_mean_curvature_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final(isnan(allapo_rstrct_mean_alpha_final)) = 0;
allapo_rstrct_mean_alpha_final = nonzeros(allapo_rstrct_mean_alpha_final);
allapo_rstrct_mean_curvature_final(isnan(allapo_rstrct_mean_curvature_final))=0;
allapo_rstrct_mean_curvature_final = nonzeros(allapo_rstrct_mean_curvature_final);
allapo_rstrct_total_distance_final(isnan(allapo_rstrct_total_distance_final))=0;
allapo_rstrct_total_distance_final = nonzeros(allapo_rstrct_total_distance_final);

mean_fiber_length_smoothn = mean(allapo_rstrct_total_distance_final)
mean_pennation_angle_smoothn = mean(allapo_rstrct_mean_alpha_final)
mean_curvature_smoothn = mean(allapo_rstrct_mean_curvature_final)

fiber_visualizer_v11(muscle_mask, fv_options, roi_mesh_rstrct, anat_img, allapo_rstrct_smoothed_fiber_all_final);
set(gcf, 'position', [1200 50 500 800])
set(gca, 'fontsize', 14)
xlabel('Column', 'fontsize', 14)
ylabel('Row', 'fontsize', 14)
zlabel('Slice', 'fontsize', 14)
title('Aponeurosis Seeding-Restricted Mesh, Smoothn')

stop_pct_mask_smoothn = numel(find(allapo_rstrct_stop_list_smoothn == 4))/numel(allapo_rstrct_stop_list_smoothn)
stop_pct_fa_smoothn = numel(find(allapo_rstrct_stop_list_smoothn == 2))/numel(allapo_rstrct_stop_list_smoothn)
stop_pct_ang_smoothn = numel(find(allapo_rstrct_stop_list_smoothn == 3))/numel(allapo_rstrct_stop_list_smoothn)

%% Fiber-tracking - APO3 - Aniso

% track, smooth, and quantify fibers
[allapo_rstrct_fiber_all, ~, allapo_rstrct_stop_list_aniso, ~] = ...
    fiber_track_v20(allapo_rstrct_ft_options, e1fa_aniso, muscle_mask, roi_mesh_rstrct);
[allapo_rstrct_smoothed_fiber_all, ~, allapo_rstrct_smoothed_fiber_all_mm, ~, ~, ~, ~, ~, ~] = ...
    fiber_smoother_v14a(allapo_rstrct_fiber_all, fs_options);

% Use same fiber-tract seed points in all the conditions evaluated
allapo_rstrct_smoothed_fiber_all_mm_final = zeros(size(allapo_rstrct_smoothed_fiber_all_mm));   %will hold tracts after seed-location filtering
allapo_rstrct_smoothed_fiber_all_final = zeros(size(allapo_rstrct_smoothed_fiber_all));

for step_cntr = 1:size(allapo_rstrct_smoothed_fiber_all,3)
    for dim_cntr = 1:size(allapo_rstrct_smoothed_fiber_all,4)
            allapo_rstrct_smoothed_fiber_all_final(:,:,step_cntr,dim_cntr) = ...
                squeeze(allapo_rstrct_smoothed_fiber_all(:,:,step_cntr,dim_cntr)).*allapo_rstrct_seeds_final_mask;
            allapo_rstrct_smoothed_fiber_all_mm_final (:,:,step_cntr,dim_cntr) = ...
                squeeze(allapo_rstrct_smoothed_fiber_all_mm(:,:,step_cntr,dim_cntr)).*allapo_rstrct_seeds_final_mask;
    end
end

[allapo_rstrct_angle_list_final, allapo_rstrct_distance_list_final, allapo_rstrct_curvature_list_final, ~, allapo_rstrct_n_points_final, allapo_rstrct_area_final] = ...
    fiber_quantifier_v20a(allapo_rstrct_fq_options, allapo_rstrct_smoothed_fiber_all_final, roi_mesh_rstrct, muscle_mask);

% calculate summary statistics
allapo_rstrct_total_distance_final = squeeze(max(allapo_rstrct_distance_list_final, [], 3));
allapo_rstrct_total_distance_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final = sum(squeeze(allapo_rstrct_angle_list_final(:,:,:,1)), 3)./allapo_rstrct_n_points_final(:,:,2);
allapo_rstrct_mean_alpha_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_curvature_final = sum(allapo_rstrct_curvature_list_final, 3)./allapo_rstrct_n_points_final(:,:,3);
allapo_rstrct_mean_curvature_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final(isnan(allapo_rstrct_mean_alpha_final)) = 0;
allapo_rstrct_mean_alpha_final = nonzeros(allapo_rstrct_mean_alpha_final);
allapo_rstrct_mean_curvature_final(isnan(allapo_rstrct_mean_curvature_final))=0;
allapo_rstrct_mean_curvature_final = nonzeros(allapo_rstrct_mean_curvature_final);
allapo_rstrct_total_distance_final(isnan(allapo_rstrct_total_distance_final))=0;
allapo_rstrct_total_distance_final = nonzeros(allapo_rstrct_total_distance_final);

mean_fiber_length_aniso = mean(allapo_rstrct_total_distance_final)
mean_pennation_angle_aniso = mean(allapo_rstrct_mean_alpha_final)
mean_curvature_aniso = mean(allapo_rstrct_mean_curvature_final)

fiber_visualizer_v11(muscle_mask, fv_options, roi_mesh_rstrct, anat_img, allapo_rstrct_smoothed_fiber_all_final);
set(gcf, 'position', [1200 50 500 800])
set(gca, 'fontsize', 14)
xlabel('Column', 'fontsize', 14)
ylabel('Row', 'fontsize', 14)
zlabel('Slice', 'fontsize', 14)
title('Aponeurosis Seeding-Restricted Mesh, Aniso')

stop_pct_mask_aniso = numel(find(allapo_rstrct_stop_list_aniso == 4))/numel(allapo_rstrct_stop_list_aniso)
stop_pct_fa_aniso = numel(find(allapo_rstrct_stop_list_aniso == 2))/numel(allapo_rstrct_stop_list_aniso)
stop_pct_ang_aniso = numel(find(allapo_rstrct_stop_list_aniso == 3))/numel(allapo_rstrct_stop_list_aniso)

%% Fiber-tracking - APO3 - TPCA

% track, smooth, and quantify fibers
[allapo_rstrct_fiber_all, ~, allapo_rstrct_stop_list_tpca, ~] = ...
    fiber_track_v20(allapo_rstrct_ft_options, e1fa_tpca, muscle_mask, roi_mesh_rstrct);
[allapo_rstrct_smoothed_fiber_all, ~, allapo_rstrct_smoothed_fiber_all_mm, ~, ~, ~, ~, ~, ~] = ...
    fiber_smoother_v14a(allapo_rstrct_fiber_all, fs_options);

% Use same fiber-tract seed points in all the conditions evaluated
allapo_rstrct_smoothed_fiber_all_mm_final = zeros(size(allapo_rstrct_smoothed_fiber_all_mm));   %will hold tracts after seed-location filtering
allapo_rstrct_smoothed_fiber_all_final = zeros(size(allapo_rstrct_smoothed_fiber_all));

for step_cntr = 1:size(allapo_rstrct_smoothed_fiber_all,3)
    for dim_cntr = 1:size(allapo_rstrct_smoothed_fiber_all,4)
            allapo_rstrct_smoothed_fiber_all_final(:,:,step_cntr,dim_cntr) = ...
                squeeze(allapo_rstrct_smoothed_fiber_all(:,:,step_cntr,dim_cntr)).*allapo_rstrct_seeds_final_mask;
            allapo_rstrct_smoothed_fiber_all_mm_final (:,:,step_cntr,dim_cntr) = ...
                squeeze(allapo_rstrct_smoothed_fiber_all_mm(:,:,step_cntr,dim_cntr)).*allapo_rstrct_seeds_final_mask;
    end
end

[allapo_rstrct_angle_list_final, allapo_rstrct_distance_list_final, allapo_rstrct_curvature_list_final, ~, allapo_rstrct_n_points_final, allapo_rstrct_area_final] = ...
    fiber_quantifier_v20a(allapo_rstrct_fq_options, allapo_rstrct_smoothed_fiber_all_final, roi_mesh_rstrct, muscle_mask);

% calculate summary statistics
allapo_rstrct_total_distance_final = squeeze(max(allapo_rstrct_distance_list_final, [], 3));
allapo_rstrct_total_distance_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final = sum(squeeze(allapo_rstrct_angle_list_final(:,:,:,1)), 3)./allapo_rstrct_n_points_final(:,:,2);
allapo_rstrct_mean_alpha_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_curvature_final = sum(allapo_rstrct_curvature_list_final, 3)./allapo_rstrct_n_points_final(:,:,3);
allapo_rstrct_mean_curvature_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final(isnan(allapo_rstrct_mean_alpha_final)) = 0;
allapo_rstrct_mean_alpha_final = nonzeros(allapo_rstrct_mean_alpha_final);
allapo_rstrct_mean_curvature_final(isnan(allapo_rstrct_mean_curvature_final))=0;
allapo_rstrct_mean_curvature_final = nonzeros(allapo_rstrct_mean_curvature_final);
allapo_rstrct_total_distance_final(isnan(allapo_rstrct_total_distance_final))=0;
allapo_rstrct_total_distance_final = nonzeros(allapo_rstrct_total_distance_final);

mean_fiber_length_tpca = mean(allapo_rstrct_total_distance_final)
mean_pennation_angle_tpca = mean(allapo_rstrct_mean_alpha_final)
mean_curvature_tpca = mean(allapo_rstrct_mean_curvature_final)

fiber_visualizer_v11(muscle_mask, fv_options, roi_mesh_rstrct, anat_img, allapo_rstrct_smoothed_fiber_all_final);
set(gcf, 'position', [1200 50 500 800])
set(gca, 'fontsize', 14)
xlabel('Column', 'fontsize', 14)
ylabel('Row', 'fontsize', 14)
zlabel('Slice', 'fontsize', 14)
title('Aponeurosis Seeding-Restricted Mesh, TPCA')

stop_pct_mask_tpca = numel(find(allapo_rstrct_stop_list_tpca == 4))/numel(allapo_rstrct_stop_list_tpca)
stop_pct_fa_tpca = numel(find(allapo_rstrct_stop_list_tpca == 2))/numel(allapo_rstrct_stop_list_tpca)
stop_pct_ang_tpca = numel(find(allapo_rstrct_stop_list_tpca == 3))/numel(allapo_rstrct_stop_list_tpca)

%% Fiber-tracking - APO3 - TPCA + smoothn

% track, smooth, and quantify fibers
[allapo_rstrct_fiber_all, ~, allapo_rstrct_stop_list_tpca_smooth, ~] = ...
    fiber_track_v20(allapo_rstrct_ft_options, e1fa_tpca_smooth, muscle_mask, roi_mesh_rstrct);
[allapo_rstrct_smoothed_fiber_all, ~, allapo_rstrct_smoothed_fiber_all_mm, ~, ~, ~, ~, ~, ~] = ...
    fiber_smoother_v14a(allapo_rstrct_fiber_all, fs_options);

% Use same fiber-tract seed points in all the conditions evaluated
allapo_rstrct_smoothed_fiber_all_mm_final = zeros(size(allapo_rstrct_smoothed_fiber_all_mm));   %will hold tracts after seed-location filtering
allapo_rstrct_smoothed_fiber_all_final = zeros(size(allapo_rstrct_smoothed_fiber_all));

for step_cntr = 1:size(allapo_rstrct_smoothed_fiber_all,3)
    for dim_cntr = 1:size(allapo_rstrct_smoothed_fiber_all,4)
            allapo_rstrct_smoothed_fiber_all_final(:,:,step_cntr,dim_cntr) = ...
                squeeze(allapo_rstrct_smoothed_fiber_all(:,:,step_cntr,dim_cntr)).*allapo_rstrct_seeds_final_mask;
            allapo_rstrct_smoothed_fiber_all_mm_final (:,:,step_cntr,dim_cntr) = ...
                squeeze(allapo_rstrct_smoothed_fiber_all_mm(:,:,step_cntr,dim_cntr)).*allapo_rstrct_seeds_final_mask;
    end
end

[allapo_rstrct_angle_list_final, allapo_rstrct_distance_list_final, allapo_rstrct_curvature_list_final, ~, allapo_rstrct_n_points_final, allapo_rstrct_area_final] = ...
    fiber_quantifier_v20a(allapo_rstrct_fq_options, allapo_rstrct_smoothed_fiber_all_final, roi_mesh_rstrct, muscle_mask);

% calculate summary statistics
allapo_rstrct_total_distance_final = squeeze(max(allapo_rstrct_distance_list_final, [], 3));
allapo_rstrct_total_distance_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final = sum(squeeze(allapo_rstrct_angle_list_final(:,:,:,1)), 3)./allapo_rstrct_n_points_final(:,:,2);
allapo_rstrct_mean_alpha_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_curvature_final = sum(allapo_rstrct_curvature_list_final, 3)./allapo_rstrct_n_points_final(:,:,3);
allapo_rstrct_mean_curvature_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final(isnan(allapo_rstrct_mean_alpha_final)) = 0;
allapo_rstrct_mean_alpha_final = nonzeros(allapo_rstrct_mean_alpha_final);
allapo_rstrct_mean_curvature_final(isnan(allapo_rstrct_mean_curvature_final))=0;
allapo_rstrct_mean_curvature_final = nonzeros(allapo_rstrct_mean_curvature_final);
allapo_rstrct_total_distance_final(isnan(allapo_rstrct_total_distance_final))=0;
allapo_rstrct_total_distance_final = nonzeros(allapo_rstrct_total_distance_final);

mean_fiber_length_tpca_smooth = mean(allapo_rstrct_total_distance_final)
mean_pennation_angle_tpca_smooth = mean(allapo_rstrct_mean_alpha_final)
mean_curvature_tpca_smooth = mean(allapo_rstrct_mean_curvature_final)

fiber_visualizer_v11(muscle_mask, fv_options, roi_mesh_rstrct, anat_img, allapo_rstrct_smoothed_fiber_all_final);
set(gcf, 'position', [1200 50 500 800])
set(gca, 'fontsize', 14)
xlabel('Column', 'fontsize', 14)
ylabel('Row', 'fontsize', 14)
zlabel('Slice', 'fontsize', 14)
title('Aponeurosis Seeding-Restricted Mesh, TPCA + smoothn')

stop_pct_mask_tpca_smooth = numel(find(allapo_rstrct_stop_list_tpca_smooth == 4))/numel(allapo_rstrct_stop_list_tpca_smooth)
stop_pct_fa_tpca_smooth = numel(find(allapo_rstrct_stop_list_tpca_smooth == 2))/numel(allapo_rstrct_stop_list_tpca_smooth)
stop_pct_ang_tpca_smooth = numel(find(allapo_rstrct_stop_list_tpca_smooth == 3))/numel(allapo_rstrct_stop_list_tpca_smooth)

%% Analyze b = 0 (2 avg), b = 450 (2 avg) data

%% 
% Compute FA and first eigenvector map of DTI data

%Correct diffusion directions fror FFS scans
diff_dir=diff_param_1.bvec; %diff_param1 and diff_param2 are the same

b_val=mean(diff_param_1.bval(2:end));

%Tensor computation

FA_raw = zeros(size(dti_img_4d_reg_2avg,1),size(dti_img_4d_reg_2avg,2),size(dti_img_4d_reg_2avg,3));
E1_450_raw = zeros(size(dti_img_4d_reg_2avg,1),size(dti_img_4d_reg_2avg,2),size(dti_img_4d_reg_2avg,3),3);
FA_aniso = zeros(size(dti_img_4d_reg_2avg,1),size(dti_img_4d_reg_2avg,2),size(dti_img_4d_reg_2avg,3));
E1_450_aniso = zeros(size(dti_img_4d_reg_2avg,1),size(dti_img_4d_reg_2avg,2),size(dti_img_4d_reg_2avg,3),3);
FA_tpca = zeros(size(dti_img_4d_reg_2avg,1),size(dti_img_4d_reg_2avg,2),size(dti_img_4d_reg_2avg,3));
E1_450_tpca = zeros(size(dti_img_4d_reg_2avg,1),size(dti_img_4d_reg_2avg,2),size(dti_img_4d_reg_2avg,3),3);

for s = 1:size(dti_img_4d_reg_2avg,3)
    for r = 1:size(dti_img_4d_reg_2avg,1)
        for c = 1:size(dti_img_4d_reg_2avg,2)
                dir_m = diff_dir;  %diffusion encoding matrix
                dir_m = dir_m(2:end,:);

                %estimate tensor
                %use TA mask to make process faster
                if muscle_mask(r,c,s)>0
                    % Raw data
                    sig_v = squeeze(dti_img_4d_reg_2avg(r,c,s,:));    %signal vector
                    D = signal2tensor2(sig_v, dir_m, b_val);
                    [E,L] = svd(D);
    
                    L_v = diag(L);
                    md = mean(L_v);
                    fa = sqrt(3/2) * sqrt(((L_v(1) - md)^2 + (L_v(2) - md)^2 + (L_v(3) - md)^2) / ...
                        (L_v(1)^2 + L_v(2)^2 + L_v(3)^2));

                    loop_E1 = E(:,L_v==max(L_v));
                    if loop_E1(3)<0
                        loop_E1 = -loop_E1;
                    end
                    
                    E1_450_raw(r,c,s,:) = loop_E1;
                    FA_raw(r,c,s,:) = fa;

                    % AIS
                    sig_v = squeeze(dti_img_4d_smooth_2avg(r,c,s,:));    %signal vector
                    D = signal2tensor2(sig_v, dir_m, b_val);
                    [E,L] = svd(D);
    
                    L_v = diag(L);
                    md = mean(L_v);
                    fa = sqrt(3/2) * sqrt(((L_v(1) - md)^2 + (L_v(2) - md)^2 + (L_v(3) - md)^2) / ...
                        (L_v(1)^2 + L_v(2)^2 + L_v(3)^2));

                    loop_E1 = E(:,L_v==max(L_v));
                    if loop_E1(3)<0
                        loop_E1 = -loop_E1;
                    end
                    
                    E1_450_aniso(r,c,s,:) = loop_E1;
                    FA_aniso(r,c,s,:) = fa;

                    % TPCA
                    sig_v = squeeze(dti_img_4d_tpca_2avg(r,c,s,:));    %signal vector
                    D = signal2tensor2(sig_v, dir_m, b_val);
                    [E,L] = svd(D);
    
                    L_v = diag(L);
                    md = mean(L_v);
                    fa = sqrt(3/2) * sqrt(((L_v(1) - md)^2 + (L_v(2) - md)^2 + (L_v(3) - md)^2) / ...
                        (L_v(1)^2 + L_v(2)^2 + L_v(3)^2));

                    loop_E1 = E(:,L_v==max(L_v));
                    if loop_E1(3)<0
                        loop_E1 = -loop_E1;
                    end
                    
                    E1_450_tpca(r,c,s,:) = loop_E1;
                    FA_tpca(r,c,s,:) = fa;

                end
        end
    end
end

%% 
% Smooth first eigenvector field of TA muscle

[E1_450_raw_smooth,smooth_parameter_raw]=smoothn({E1_450_raw(:,:,:,1),E1_450_raw(:,:,:,2),E1_450_raw(:,:,:,3)});
[E1_450_aniso_smooth,smooth_parameter_aniso]=smoothn({E1_450_aniso(:,:,:,1),E1_450_aniso(:,:,:,2),E1_450_aniso(:,:,:,3)});
[E1_450_tpca_smooth,smooth_parameter_tpca]=smoothn({E1_450_tpca(:,:,:,1),E1_450_tpca(:,:,:,2),E1_450_tpca(:,:,:,3)});

E1_450_raw_smooth=cat(4,E1_450_raw_smooth{1},E1_450_raw_smooth{2},E1_450_raw_smooth{3});
E1_450_aniso_smooth=cat(4,E1_450_aniso_smooth{1},E1_450_aniso_smooth{2},E1_450_aniso_smooth{3});
E1_450_tpca_smooth=cat(4,E1_450_tpca_smooth{1},E1_450_tpca_smooth{2},E1_450_tpca_smooth{3});

for r_img=1:size(E1_450_aniso,1)
    for c_img=1:size(E1_450_aniso,2)
        for s_img=1:size(E1_450_aniso,3)
            if muscle_mask(r_img,c_img,s_img)>0

                E1_raw_smooth = squeeze(E1_450_raw_smooth(r_img,c_img,s_img,:));
                E1_aniso_smooth = squeeze(E1_450_aniso_smooth(r_img,c_img,s_img,:));
                E1_tpca_smooth = squeeze(E1_450_tpca_smooth(r_img,c_img,s_img,:));

                E1_raw_smooth = E1_raw_smooth/norm(E1_raw_smooth);
                E1_aniso_smooth = E1_aniso_smooth/norm(E1_aniso_smooth);
                E1_tpca_smooth = E1_tpca_smooth/norm(E1_tpca_smooth);

                E1_450_raw_smooth(r_img,c_img,s_img,:)=E1_raw_smooth;
                E1_450_aniso_smooth(r_img,c_img,s_img,:)=E1_aniso_smooth;
                E1_450_tpca_smooth(r_img,c_img,s_img,:)=E1_tpca_smooth;
            end
        end
    end
end

for i = 1:3
    E1_450_raw_smooth(:,:,:,i)=squeeze(E1_450_raw_smooth(:,:,:,i)).*muscle_mask;
    E1_450_aniso_smooth(:,:,:,i)=squeeze(E1_450_aniso_smooth(:,:,:,i)).*muscle_mask;
    E1_450_tpca_smooth(:,:,:,i)=squeeze(E1_450_tpca_smooth(:,:,:,i)).*muscle_mask;
end

%% 
% Create E1-FA maps for fiber tractography

e1fa_raw = cat(4,E1_450_raw,FA_raw);
e1fa_aniso = cat(4,E1_450_aniso,FA_aniso);
e1fa_tpca = cat(4,E1_450_tpca,FA_tpca);
e1fa_raw_smooth = cat(4,E1_450_raw_smooth,FA_raw);
e1fa_aniso_smooth = cat(4,E1_450_aniso_smooth,FA_aniso);
e1fa_tpca_smooth = cat(4,E1_450_tpca_smooth,FA_tpca);

%% Fiber-tracking - APO3 - Aniso4D + Smoothn

% track, smooth, and quantify fibers
[allapo_rstrct_fiber_all, ~, allapo_rstrct_stop_list_aniso_smoothn, ~] = ...
    fiber_track_v20(allapo_rstrct_ft_options, e1fa_aniso_smooth, muscle_mask, roi_mesh_rstrct);
[allapo_rstrct_smoothed_fiber_all, allapo_rstrct_fiber_all_mm, allapo_rstrct_smoothed_fiber_all_mm, ~, ~, ~, ~, ~, ~] = ...
    fiber_smoother_v14a(allapo_rstrct_fiber_all, fs_options);

% Use same fiber-tract seed points in all the conditions evaluated
allapo_rstrct_smoothed_fiber_all_mm_final = zeros(size(allapo_rstrct_smoothed_fiber_all_mm));   %will hold tracts after seed-location filtering
allapo_rstrct_smoothed_fiber_all_final = zeros(size(allapo_rstrct_smoothed_fiber_all));

for step_cntr = 1:size(allapo_rstrct_smoothed_fiber_all,3)
    for dim_cntr = 1:size(allapo_rstrct_smoothed_fiber_all,4)
            allapo_rstrct_smoothed_fiber_all_final(:,:,step_cntr,dim_cntr) = ...
                squeeze(allapo_rstrct_smoothed_fiber_all(:,:,step_cntr,dim_cntr)).*allapo_rstrct_seeds_final_mask;
            allapo_rstrct_smoothed_fiber_all_mm_final (:,:,step_cntr,dim_cntr) = ...
                squeeze(allapo_rstrct_smoothed_fiber_all_mm(:,:,step_cntr,dim_cntr)).*allapo_rstrct_seeds_final_mask;
    end
end

[allapo_rstrct_angle_list_final, allapo_rstrct_distance_list_final, allapo_rstrct_curvature_list_final, ~, allapo_rstrct_n_points_final, allapo_rstrct_area_final] = ...
    fiber_quantifier_v20a(allapo_rstrct_fq_options, allapo_rstrct_smoothed_fiber_all_final, roi_mesh_rstrct, muscle_mask);

% calculate summary statistics
allapo_rstrct_total_distance_final = squeeze(max(allapo_rstrct_distance_list_final, [], 3));
allapo_rstrct_total_distance_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final = sum(squeeze(allapo_rstrct_angle_list_final(:,:,:,1)), 3)./allapo_rstrct_n_points_final(:,:,2);
allapo_rstrct_mean_alpha_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_curvature_final = sum(allapo_rstrct_curvature_list_final, 3)./allapo_rstrct_n_points_final(:,:,3);
allapo_rstrct_mean_curvature_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final(isnan(allapo_rstrct_mean_alpha_final)) = 0;
allapo_rstrct_mean_alpha_final = nonzeros(allapo_rstrct_mean_alpha_final);
allapo_rstrct_mean_curvature_final(isnan(allapo_rstrct_mean_curvature_final))=0;
allapo_rstrct_mean_curvature_final = nonzeros(allapo_rstrct_mean_curvature_final);
allapo_rstrct_total_distance_final(isnan(allapo_rstrct_total_distance_final))=0;
allapo_rstrct_total_distance_final = nonzeros(allapo_rstrct_total_distance_final);

mean_fiber_length_aniso_smoothn = mean(allapo_rstrct_total_distance_final)
mean_pennation_angle_aniso_smoothn = mean(allapo_rstrct_mean_alpha_final)
mean_curvature_aniso_smoothn = mean(allapo_rstrct_mean_curvature_final)

fiber_visualizer_v11(muscle_mask, fv_options, roi_mesh_rstrct, anat_img, allapo_rstrct_smoothed_fiber_all_final);
set(gcf, 'position', [1200 50 500 800])
set(gca, 'fontsize', 14)
xlabel('Column', 'fontsize', 14)
ylabel('Row', 'fontsize', 14)
zlabel('Slice', 'fontsize', 14)
title('Aponeurosis Seeding-Restricted Mesh, Aniso + Smoothn')

stop_pct_mask_aniso_smoothn = numel(find(allapo_rstrct_stop_list_aniso_smoothn == 4))/numel(allapo_rstrct_stop_list_aniso_smoothn)
stop_pct_fa_aniso_smoothn = numel(find(allapo_rstrct_stop_list_aniso_smoothn == 2))/numel(allapo_rstrct_stop_list_aniso_smoothn)
stop_pct_ang_aniso_smoothn = numel(find(allapo_rstrct_stop_list_aniso_smoothn == 3))/numel(allapo_rstrct_stop_list_aniso_smoothn)

%% Fiber-tracking - APO3 - Raw

% track, smooth, and quantify fibers
[allapo_rstrct_fiber_all, allapo_rstrct_roi_flag, allapo_rstrct_stop_list_raw, allapo_rstrct_fiber_len] = ...
    fiber_track_v20(allapo_rstrct_ft_options, e1fa_raw, muscle_mask, roi_mesh_rstrct);
[allapo_rstrct_smoothed_fiber_all, ~, allapo_rstrct_smoothed_fiber_all_mm, allapo_rstrct_pcoeff_r, allapo_rstrct_pcoeff_c, allapo_rstrct_pcoeff_s, allapo_rstrct_n_points_smoothed, allapo_rstrct_residuals, allapo_rstrct_residuals_mm] = ...
    fiber_smoother_v14a(allapo_rstrct_fiber_all, fs_options);

% Use same fiber-tract seed points in all the conditions evaluated
allapo_rstrct_smoothed_fiber_all_mm_final = zeros(size(allapo_rstrct_smoothed_fiber_all_mm));   %will hold tracts after seed-location filtering
allapo_rstrct_smoothed_fiber_all_final = zeros(size(allapo_rstrct_smoothed_fiber_all));

for step_cntr = 1:size(allapo_rstrct_smoothed_fiber_all,3)
    for dim_cntr = 1:size(allapo_rstrct_smoothed_fiber_all,4)
            allapo_rstrct_smoothed_fiber_all_final(:,:,step_cntr,dim_cntr) = ...
                squeeze(allapo_rstrct_smoothed_fiber_all(:,:,step_cntr,dim_cntr)).*allapo_rstrct_seeds_final_mask;
            allapo_rstrct_smoothed_fiber_all_mm_final (:,:,step_cntr,dim_cntr) = ...
                squeeze(allapo_rstrct_smoothed_fiber_all_mm(:,:,step_cntr,dim_cntr)).*allapo_rstrct_seeds_final_mask;
    end
end

[allapo_rstrct_angle_list_final, allapo_rstrct_distance_list_final, allapo_rstrct_curvature_list_final, ~, allapo_rstrct_n_points_final, allapo_rstrct_area_final] = ...
    fiber_quantifier_v20a(allapo_rstrct_fq_options, allapo_rstrct_smoothed_fiber_all_final, roi_mesh_rstrct, muscle_mask);

% calculate summary statistics
allapo_rstrct_total_distance_final = squeeze(max(allapo_rstrct_distance_list_final, [], 3));
allapo_rstrct_total_distance_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final = sum(squeeze(allapo_rstrct_angle_list_final(:,:,:,1)), 3)./allapo_rstrct_n_points_final(:,:,2);
allapo_rstrct_mean_alpha_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_curvature_final = sum(allapo_rstrct_curvature_list_final, 3)./allapo_rstrct_n_points_final(:,:,3);
allapo_rstrct_mean_curvature_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final(isnan(allapo_rstrct_mean_alpha_final)) = 0;
allapo_rstrct_mean_alpha_final = nonzeros(allapo_rstrct_mean_alpha_final);
allapo_rstrct_mean_curvature_final(isnan(allapo_rstrct_mean_curvature_final))=0;
allapo_rstrct_mean_curvature_final = nonzeros(allapo_rstrct_mean_curvature_final);
allapo_rstrct_total_distance_final(isnan(allapo_rstrct_total_distance_final))=0;
allapo_rstrct_total_distance_final = nonzeros(allapo_rstrct_total_distance_final);

mean_fiber_length_raw = mean(allapo_rstrct_total_distance_final)
mean_pennation_angle_raw = mean(allapo_rstrct_mean_alpha_final)
mean_curvature_raw = mean(allapo_rstrct_mean_curvature_final)

fiber_visualizer_v11(muscle_mask, fv_options, roi_mesh_rstrct, anat_img, allapo_rstrct_smoothed_fiber_all_final);
set(gcf, 'position', [1200 50 500 800])
set(gca, 'fontsize', 14)
xlabel('Column', 'fontsize', 14)
ylabel('Row', 'fontsize', 14)
zlabel('Slice', 'fontsize', 14)
title('Aponeurosis Seeding-Restricted Mesh, Raw')

stop_pct_mask_raw = numel(find(allapo_rstrct_stop_list_raw == 4))/numel(allapo_rstrct_stop_list_raw)
stop_pct_fa_raw = numel(find(allapo_rstrct_stop_list_raw == 2))/numel(allapo_rstrct_stop_list_raw)
stop_pct_ang_raw = numel(find(allapo_rstrct_stop_list_raw == 3))/numel(allapo_rstrct_stop_list_raw)

%% Fiber-tracking - APO3 - Smoothn

% track, smooth, and quantify fibers
[allapo_rstrct_fiber_all, ~, allapo_rstrct_stop_list_smoothn, ~] = ...
    fiber_track_v20(allapo_rstrct_ft_options, e1fa_raw_smooth, muscle_mask, roi_mesh_rstrct);
[allapo_rstrct_smoothed_fiber_all, ~, allapo_rstrct_smoothed_fiber_all_mm, ~, ~, ~, ~, ~, ~] = ...
    fiber_smoother_v14a(allapo_rstrct_fiber_all, fs_options);

% Use same fiber-tract seed points in all the conditions evaluated
allapo_rstrct_smoothed_fiber_all_mm_final = zeros(size(allapo_rstrct_smoothed_fiber_all_mm));   %will hold tracts after seed-location filtering
allapo_rstrct_smoothed_fiber_all_final = zeros(size(allapo_rstrct_smoothed_fiber_all));

for step_cntr = 1:size(allapo_rstrct_smoothed_fiber_all,3)
    for dim_cntr = 1:size(allapo_rstrct_smoothed_fiber_all,4)
            allapo_rstrct_smoothed_fiber_all_final(:,:,step_cntr,dim_cntr) = ...
                squeeze(allapo_rstrct_smoothed_fiber_all(:,:,step_cntr,dim_cntr)).*allapo_rstrct_seeds_final_mask;
            allapo_rstrct_smoothed_fiber_all_mm_final (:,:,step_cntr,dim_cntr) = ...
                squeeze(allapo_rstrct_smoothed_fiber_all_mm(:,:,step_cntr,dim_cntr)).*allapo_rstrct_seeds_final_mask;
    end
end

[allapo_rstrct_angle_list_final, allapo_rstrct_distance_list_final, allapo_rstrct_curvature_list_final, ~, allapo_rstrct_n_points_final, allapo_rstrct_area_final] = ...
    fiber_quantifier_v20a(allapo_rstrct_fq_options, allapo_rstrct_smoothed_fiber_all_final, roi_mesh_rstrct, muscle_mask);

% calculate summary statistics
allapo_rstrct_total_distance_final = squeeze(max(allapo_rstrct_distance_list_final, [], 3));
allapo_rstrct_total_distance_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final = sum(squeeze(allapo_rstrct_angle_list_final(:,:,:,1)), 3)./allapo_rstrct_n_points_final(:,:,2);
allapo_rstrct_mean_alpha_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_curvature_final = sum(allapo_rstrct_curvature_list_final, 3)./allapo_rstrct_n_points_final(:,:,3);
allapo_rstrct_mean_curvature_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final(isnan(allapo_rstrct_mean_alpha_final)) = 0;
allapo_rstrct_mean_alpha_final = nonzeros(allapo_rstrct_mean_alpha_final);
allapo_rstrct_mean_curvature_final(isnan(allapo_rstrct_mean_curvature_final))=0;
allapo_rstrct_mean_curvature_final = nonzeros(allapo_rstrct_mean_curvature_final);
allapo_rstrct_total_distance_final(isnan(allapo_rstrct_total_distance_final))=0;
allapo_rstrct_total_distance_final = nonzeros(allapo_rstrct_total_distance_final);

mean_fiber_length_smoothn = mean(allapo_rstrct_total_distance_final)
mean_pennation_angle_smoothn = mean(allapo_rstrct_mean_alpha_final)
mean_curvature_smoothn = mean(allapo_rstrct_mean_curvature_final)

fiber_visualizer_v11(muscle_mask, fv_options, roi_mesh_rstrct, anat_img, allapo_rstrct_smoothed_fiber_all_final);
set(gcf, 'position', [1200 50 500 800])
set(gca, 'fontsize', 14)
xlabel('Column', 'fontsize', 14)
ylabel('Row', 'fontsize', 14)
zlabel('Slice', 'fontsize', 14)
title('Aponeurosis Seeding-Restricted Mesh, Smoothn')

stop_pct_mask_smoothn = numel(find(allapo_rstrct_stop_list_smoothn == 4))/numel(allapo_rstrct_stop_list_smoothn)
stop_pct_fa_smoothn = numel(find(allapo_rstrct_stop_list_smoothn == 2))/numel(allapo_rstrct_stop_list_smoothn)
stop_pct_ang_smoothn = numel(find(allapo_rstrct_stop_list_smoothn == 3))/numel(allapo_rstrct_stop_list_smoothn)

%% Fiber-tracking - APO3 - Aniso

% track, smooth, and quantify fibers
[allapo_rstrct_fiber_all, ~, allapo_rstrct_stop_list_aniso, ~] = ...
    fiber_track_v20(allapo_rstrct_ft_options, e1fa_aniso, muscle_mask, roi_mesh_rstrct);
[allapo_rstrct_smoothed_fiber_all, ~, allapo_rstrct_smoothed_fiber_all_mm, ~, ~, ~, ~, ~, ~] = ...
    fiber_smoother_v14a(allapo_rstrct_fiber_all, fs_options);

% Use same fiber-tract seed points in all the conditions evaluated
allapo_rstrct_smoothed_fiber_all_mm_final = zeros(size(allapo_rstrct_smoothed_fiber_all_mm));   %will hold tracts after seed-location filtering
allapo_rstrct_smoothed_fiber_all_final = zeros(size(allapo_rstrct_smoothed_fiber_all));

for step_cntr = 1:size(allapo_rstrct_smoothed_fiber_all,3)
    for dim_cntr = 1:size(allapo_rstrct_smoothed_fiber_all,4)
            allapo_rstrct_smoothed_fiber_all_final(:,:,step_cntr,dim_cntr) = ...
                squeeze(allapo_rstrct_smoothed_fiber_all(:,:,step_cntr,dim_cntr)).*allapo_rstrct_seeds_final_mask;
            allapo_rstrct_smoothed_fiber_all_mm_final (:,:,step_cntr,dim_cntr) = ...
                squeeze(allapo_rstrct_smoothed_fiber_all_mm(:,:,step_cntr,dim_cntr)).*allapo_rstrct_seeds_final_mask;
    end
end

[allapo_rstrct_angle_list_final, allapo_rstrct_distance_list_final, allapo_rstrct_curvature_list_final, ~, allapo_rstrct_n_points_final, allapo_rstrct_area_final] = ...
    fiber_quantifier_v20a(allapo_rstrct_fq_options, allapo_rstrct_smoothed_fiber_all_final, roi_mesh_rstrct, muscle_mask);

% calculate summary statistics
allapo_rstrct_total_distance_final = squeeze(max(allapo_rstrct_distance_list_final, [], 3));
allapo_rstrct_total_distance_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final = sum(squeeze(allapo_rstrct_angle_list_final(:,:,:,1)), 3)./allapo_rstrct_n_points_final(:,:,2);
allapo_rstrct_mean_alpha_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_curvature_final = sum(allapo_rstrct_curvature_list_final, 3)./allapo_rstrct_n_points_final(:,:,3);
allapo_rstrct_mean_curvature_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final(isnan(allapo_rstrct_mean_alpha_final)) = 0;
allapo_rstrct_mean_alpha_final = nonzeros(allapo_rstrct_mean_alpha_final);
allapo_rstrct_mean_curvature_final(isnan(allapo_rstrct_mean_curvature_final))=0;
allapo_rstrct_mean_curvature_final = nonzeros(allapo_rstrct_mean_curvature_final);
allapo_rstrct_total_distance_final(isnan(allapo_rstrct_total_distance_final))=0;
allapo_rstrct_total_distance_final = nonzeros(allapo_rstrct_total_distance_final);

mean_fiber_length_aniso = mean(allapo_rstrct_total_distance_final)
mean_pennation_angle_aniso = mean(allapo_rstrct_mean_alpha_final)
mean_curvature_aniso = mean(allapo_rstrct_mean_curvature_final)

fiber_visualizer_v11(muscle_mask, fv_options, roi_mesh_rstrct, anat_img, allapo_rstrct_smoothed_fiber_all_final);
set(gcf, 'position', [1200 50 500 800])
set(gca, 'fontsize', 14)
xlabel('Column', 'fontsize', 14)
ylabel('Row', 'fontsize', 14)
zlabel('Slice', 'fontsize', 14)
title('Aponeurosis Seeding-Restricted Mesh, Aniso')

stop_pct_mask_aniso = numel(find(allapo_rstrct_stop_list_aniso == 4))/numel(allapo_rstrct_stop_list_aniso)
stop_pct_fa_aniso = numel(find(allapo_rstrct_stop_list_aniso == 2))/numel(allapo_rstrct_stop_list_aniso)
stop_pct_ang_aniso = numel(find(allapo_rstrct_stop_list_aniso == 3))/numel(allapo_rstrct_stop_list_aniso)

%% Fiber-tracking - APO3 - TPCA

% track, smooth, and quantify fibers
[allapo_rstrct_fiber_all, ~, allapo_rstrct_stop_list_tpca, ~] = ...
    fiber_track_v20(allapo_rstrct_ft_options, e1fa_tpca, muscle_mask, roi_mesh_rstrct);
[allapo_rstrct_smoothed_fiber_all, ~, allapo_rstrct_smoothed_fiber_all_mm, ~, ~, ~, ~, ~, ~] = ...
    fiber_smoother_v14a(allapo_rstrct_fiber_all, fs_options);

% Use same fiber-tract seed points in all the conditions evaluated
allapo_rstrct_smoothed_fiber_all_mm_final = zeros(size(allapo_rstrct_smoothed_fiber_all_mm));   %will hold tracts after seed-location filtering
allapo_rstrct_smoothed_fiber_all_final = zeros(size(allapo_rstrct_smoothed_fiber_all));

for step_cntr = 1:size(allapo_rstrct_smoothed_fiber_all,3)
    for dim_cntr = 1:size(allapo_rstrct_smoothed_fiber_all,4)
            allapo_rstrct_smoothed_fiber_all_final(:,:,step_cntr,dim_cntr) = ...
                squeeze(allapo_rstrct_smoothed_fiber_all(:,:,step_cntr,dim_cntr)).*allapo_rstrct_seeds_final_mask;
            allapo_rstrct_smoothed_fiber_all_mm_final (:,:,step_cntr,dim_cntr) = ...
                squeeze(allapo_rstrct_smoothed_fiber_all_mm(:,:,step_cntr,dim_cntr)).*allapo_rstrct_seeds_final_mask;
    end
end

[allapo_rstrct_angle_list_final, allapo_rstrct_distance_list_final, allapo_rstrct_curvature_list_final, ~, allapo_rstrct_n_points_final, allapo_rstrct_area_final] = ...
    fiber_quantifier_v20a(allapo_rstrct_fq_options, allapo_rstrct_smoothed_fiber_all_final, roi_mesh_rstrct, muscle_mask);

% calculate summary statistics
allapo_rstrct_total_distance_final = squeeze(max(allapo_rstrct_distance_list_final, [], 3));
allapo_rstrct_total_distance_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final = sum(squeeze(allapo_rstrct_angle_list_final(:,:,:,1)), 3)./allapo_rstrct_n_points_final(:,:,2);
allapo_rstrct_mean_alpha_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_curvature_final = sum(allapo_rstrct_curvature_list_final, 3)./allapo_rstrct_n_points_final(:,:,3);
allapo_rstrct_mean_curvature_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final(isnan(allapo_rstrct_mean_alpha_final)) = 0;
allapo_rstrct_mean_alpha_final = nonzeros(allapo_rstrct_mean_alpha_final);
allapo_rstrct_mean_curvature_final(isnan(allapo_rstrct_mean_curvature_final))=0;
allapo_rstrct_mean_curvature_final = nonzeros(allapo_rstrct_mean_curvature_final);
allapo_rstrct_total_distance_final(isnan(allapo_rstrct_total_distance_final))=0;
allapo_rstrct_total_distance_final = nonzeros(allapo_rstrct_total_distance_final);

mean_fiber_length_tpca = mean(allapo_rstrct_total_distance_final)
mean_pennation_angle_tpca = mean(allapo_rstrct_mean_alpha_final)
mean_curvature_tpca = mean(allapo_rstrct_mean_curvature_final)

fiber_visualizer_v11(muscle_mask, fv_options, roi_mesh_rstrct, anat_img, allapo_rstrct_smoothed_fiber_all_final);
set(gcf, 'position', [1200 50 500 800])
set(gca, 'fontsize', 14)
xlabel('Column', 'fontsize', 14)
ylabel('Row', 'fontsize', 14)
zlabel('Slice', 'fontsize', 14)
title('Aponeurosis Seeding-Restricted Mesh, TPCA')

stop_pct_mask_tpca = numel(find(allapo_rstrct_stop_list_tpca == 4))/numel(allapo_rstrct_stop_list_tpca)
stop_pct_fa_tpca = numel(find(allapo_rstrct_stop_list_tpca == 2))/numel(allapo_rstrct_stop_list_tpca)
stop_pct_ang_tpca = numel(find(allapo_rstrct_stop_list_tpca == 3))/numel(allapo_rstrct_stop_list_tpca)

%% Fiber-tracking - APO3 - TPCA + smoothn

% track, smooth, and quantify fibers
[allapo_rstrct_fiber_all, ~, allapo_rstrct_stop_list_tpca_smooth, ~] = ...
    fiber_track_v20(allapo_rstrct_ft_options, e1fa_tpca_smooth, muscle_mask, roi_mesh_rstrct);
[allapo_rstrct_smoothed_fiber_all, ~, allapo_rstrct_smoothed_fiber_all_mm, ~, ~, ~, ~, ~, ~] = ...
    fiber_smoother_v14a(allapo_rstrct_fiber_all, fs_options);

% Use same fiber-tract seed points in all the conditions evaluated
allapo_rstrct_smoothed_fiber_all_mm_final = zeros(size(allapo_rstrct_smoothed_fiber_all_mm));   %will hold tracts after seed-location filtering
allapo_rstrct_smoothed_fiber_all_final = zeros(size(allapo_rstrct_smoothed_fiber_all));

for step_cntr = 1:size(allapo_rstrct_smoothed_fiber_all,3)
    for dim_cntr = 1:size(allapo_rstrct_smoothed_fiber_all,4)
            allapo_rstrct_smoothed_fiber_all_final(:,:,step_cntr,dim_cntr) = ...
                squeeze(allapo_rstrct_smoothed_fiber_all(:,:,step_cntr,dim_cntr)).*allapo_rstrct_seeds_final_mask;
            allapo_rstrct_smoothed_fiber_all_mm_final (:,:,step_cntr,dim_cntr) = ...
                squeeze(allapo_rstrct_smoothed_fiber_all_mm(:,:,step_cntr,dim_cntr)).*allapo_rstrct_seeds_final_mask;
    end
end

[allapo_rstrct_angle_list_final, allapo_rstrct_distance_list_final, allapo_rstrct_curvature_list_final, ~, allapo_rstrct_n_points_final, allapo_rstrct_area_final] = ...
    fiber_quantifier_v20a(allapo_rstrct_fq_options, allapo_rstrct_smoothed_fiber_all_final, roi_mesh_rstrct, muscle_mask);

% calculate summary statistics
allapo_rstrct_total_distance_final = squeeze(max(allapo_rstrct_distance_list_final, [], 3));
allapo_rstrct_total_distance_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final = sum(squeeze(allapo_rstrct_angle_list_final(:,:,:,1)), 3)./allapo_rstrct_n_points_final(:,:,2);
allapo_rstrct_mean_alpha_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_curvature_final = sum(allapo_rstrct_curvature_list_final, 3)./allapo_rstrct_n_points_final(:,:,3);
allapo_rstrct_mean_curvature_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final(isnan(allapo_rstrct_mean_alpha_final)) = 0;
allapo_rstrct_mean_alpha_final = nonzeros(allapo_rstrct_mean_alpha_final);
allapo_rstrct_mean_curvature_final(isnan(allapo_rstrct_mean_curvature_final))=0;
allapo_rstrct_mean_curvature_final = nonzeros(allapo_rstrct_mean_curvature_final);
allapo_rstrct_total_distance_final(isnan(allapo_rstrct_total_distance_final))=0;
allapo_rstrct_total_distance_final = nonzeros(allapo_rstrct_total_distance_final);

mean_fiber_length_tpca_smooth = mean(allapo_rstrct_total_distance_final)
mean_pennation_angle_tpca_smooth = mean(allapo_rstrct_mean_alpha_final)
mean_curvature_tpca_smooth = mean(allapo_rstrct_mean_curvature_final)

fiber_visualizer_v11(muscle_mask, fv_options, roi_mesh_rstrct, anat_img, allapo_rstrct_smoothed_fiber_all_final);
set(gcf, 'position', [1200 50 500 800])
set(gca, 'fontsize', 14)
xlabel('Column', 'fontsize', 14)
ylabel('Row', 'fontsize', 14)
zlabel('Slice', 'fontsize', 14)
title('Aponeurosis Seeding-Restricted Mesh, TPCA + smoothn')

stop_pct_mask_tpca_smooth = numel(find(allapo_rstrct_stop_list_tpca_smooth == 4))/numel(allapo_rstrct_stop_list_tpca_smooth)
stop_pct_fa_tpca_smooth = numel(find(allapo_rstrct_stop_list_tpca_smooth == 2))/numel(allapo_rstrct_stop_list_tpca_smooth)
stop_pct_ang_tpca_smooth = numel(find(allapo_rstrct_stop_list_tpca_smooth == 3))/numel(allapo_rstrct_stop_list_tpca_smooth)

%% Analyze b = 0 (2 avg), b = 450 (1 avg) data

%% 
% Compute FA and first eigenvector map of DTI data

%Correct diffusion directions fror FFS scans
diff_dir=diff_param_1.bvec; %diff_param1 and diff_param2 are the same

b_val=mean(diff_param_1.bval(2:end));

%Tensor computation

FA_raw = zeros(size(dti_img_4d_reg_b0_2avg,1),size(dti_img_4d_reg_b0_2avg,2),size(dti_img_4d_reg_b0_2avg,3));
E1_450_raw = zeros(size(dti_img_4d_reg_b0_2avg,1),size(dti_img_4d_reg_b0_2avg,2),size(dti_img_4d_reg_b0_2avg,3),3);
FA_aniso = zeros(size(dti_img_4d_reg_b0_2avg,1),size(dti_img_4d_reg_b0_2avg,2),size(dti_img_4d_reg_b0_2avg,3));
E1_450_aniso = zeros(size(dti_img_4d_reg_b0_2avg,1),size(dti_img_4d_reg_b0_2avg,2),size(dti_img_4d_reg_b0_2avg,3),3);
FA_tpca = zeros(size(dti_img_4d_reg_b0_2avg,1),size(dti_img_4d_reg_b0_2avg,2),size(dti_img_4d_reg_b0_2avg,3));
E1_450_tpca = zeros(size(dti_img_4d_reg_b0_2avg,1),size(dti_img_4d_reg_b0_2avg,2),size(dti_img_4d_reg_b0_2avg,3),3);

for s = 1:size(dti_img_4d_reg_b0_2avg,3)
    for r = 1:size(dti_img_4d_reg_b0_2avg,1)
        for c = 1:size(dti_img_4d_reg_b0_2avg,2)
                dir_m = diff_dir;  %diffusion encoding matrix
                dir_m = dir_m(2:end,:);

                %estimate tensor
                %use TA mask to make process faster
                if muscle_mask(r,c,s)>0
                    % Raw data
                    sig_v = squeeze(dti_img_4d_reg_b0_2avg(r,c,s,:));    %signal vector
                    D = signal2tensor2(sig_v, dir_m, b_val);
                    [E,L] = svd(D);
    
                    L_v = diag(L);
                    md = mean(L_v);
                    fa = sqrt(3/2) * sqrt(((L_v(1) - md)^2 + (L_v(2) - md)^2 + (L_v(3) - md)^2) / ...
                        (L_v(1)^2 + L_v(2)^2 + L_v(3)^2));

                    loop_E1 = E(:,L_v==max(L_v));
                    if loop_E1(3)<0
                        loop_E1 = -loop_E1;
                    end
                    
                    E1_450_raw(r,c,s,:) = loop_E1;
                    FA_raw(r,c,s,:) = fa;

                    % AIS
                    sig_v = squeeze(dti_img_4d_smooth_b0_2avg(r,c,s,:));    %signal vector
                    D = signal2tensor2(sig_v, dir_m, b_val);
                    [E,L] = svd(D);
    
                    L_v = diag(L);
                    md = mean(L_v);
                    fa = sqrt(3/2) * sqrt(((L_v(1) - md)^2 + (L_v(2) - md)^2 + (L_v(3) - md)^2) / ...
                        (L_v(1)^2 + L_v(2)^2 + L_v(3)^2));

                    loop_E1 = E(:,L_v==max(L_v));
                    if loop_E1(3)<0
                        loop_E1 = -loop_E1;
                    end
                    
                    E1_450_aniso(r,c,s,:) = loop_E1;
                    FA_aniso(r,c,s,:) = fa;

                    % TPCA
                    sig_v = squeeze(dti_img_4d_tpca_b0_2avg(r,c,s,:));    %signal vector
                    D = signal2tensor2(sig_v, dir_m, b_val);
                    [E,L] = svd(D);
    
                    L_v = diag(L);
                    md = mean(L_v);
                    fa = sqrt(3/2) * sqrt(((L_v(1) - md)^2 + (L_v(2) - md)^2 + (L_v(3) - md)^2) / ...
                        (L_v(1)^2 + L_v(2)^2 + L_v(3)^2));

                    loop_E1 = E(:,L_v==max(L_v));
                    if loop_E1(3)<0
                        loop_E1 = -loop_E1;
                    end
                    
                    E1_450_tpca(r,c,s,:) = loop_E1;
                    FA_tpca(r,c,s,:) = fa;

                end
        end
    end
end

%% 
% Smooth first eigenvector field of TA muscle

[E1_450_raw_smooth,smooth_parameter_raw]=smoothn({E1_450_raw(:,:,:,1),E1_450_raw(:,:,:,2),E1_450_raw(:,:,:,3)});
[E1_450_aniso_smooth,smooth_parameter_aniso]=smoothn({E1_450_aniso(:,:,:,1),E1_450_aniso(:,:,:,2),E1_450_aniso(:,:,:,3)});
[E1_450_tpca_smooth,smooth_parameter_tpca]=smoothn({E1_450_tpca(:,:,:,1),E1_450_tpca(:,:,:,2),E1_450_tpca(:,:,:,3)});

E1_450_raw_smooth=cat(4,E1_450_raw_smooth{1},E1_450_raw_smooth{2},E1_450_raw_smooth{3});
E1_450_aniso_smooth=cat(4,E1_450_aniso_smooth{1},E1_450_aniso_smooth{2},E1_450_aniso_smooth{3});
E1_450_tpca_smooth=cat(4,E1_450_tpca_smooth{1},E1_450_tpca_smooth{2},E1_450_tpca_smooth{3});

for r_img=1:size(E1_450_aniso,1)
    for c_img=1:size(E1_450_aniso,2)
        for s_img=1:size(E1_450_aniso,3)
            if muscle_mask(r_img,c_img,s_img)>0

                E1_raw_smooth = squeeze(E1_450_raw_smooth(r_img,c_img,s_img,:));
                E1_aniso_smooth = squeeze(E1_450_aniso_smooth(r_img,c_img,s_img,:));
                E1_tpca_smooth = squeeze(E1_450_tpca_smooth(r_img,c_img,s_img,:));

                E1_raw_smooth = E1_raw_smooth/norm(E1_raw_smooth);
                E1_aniso_smooth = E1_aniso_smooth/norm(E1_aniso_smooth);
                E1_tpca_smooth = E1_tpca_smooth/norm(E1_tpca_smooth);

                E1_450_raw_smooth(r_img,c_img,s_img,:)=E1_raw_smooth;
                E1_450_aniso_smooth(r_img,c_img,s_img,:)=E1_aniso_smooth;
                E1_450_tpca_smooth(r_img,c_img,s_img,:)=E1_tpca_smooth;
            end
        end
    end
end

for i = 1:3
    E1_450_raw_smooth(:,:,:,i)=squeeze(E1_450_raw_smooth(:,:,:,i)).*muscle_mask;
    E1_450_aniso_smooth(:,:,:,i)=squeeze(E1_450_aniso_smooth(:,:,:,i)).*muscle_mask;
    E1_450_tpca_smooth(:,:,:,i)=squeeze(E1_450_tpca_smooth(:,:,:,i)).*muscle_mask;
end

%% 
% Create E1-FA maps for fiber tractography

e1fa_raw = cat(4,E1_450_raw,FA_raw);
e1fa_aniso = cat(4,E1_450_aniso,FA_aniso);
e1fa_tpca = cat(4,E1_450_tpca,FA_tpca);
e1fa_raw_smooth = cat(4,E1_450_raw_smooth,FA_raw);
e1fa_aniso_smooth = cat(4,E1_450_aniso_smooth,FA_aniso);
e1fa_tpca_smooth = cat(4,E1_450_tpca_smooth,FA_tpca);

%% Fiber-tracking - APO3 - Aniso4D + Smoothn

% track, smooth, and quantify fibers
[allapo_rstrct_fiber_all, ~, allapo_rstrct_stop_list_aniso_smoothn, ~] = ...
    fiber_track_v20(allapo_rstrct_ft_options, e1fa_aniso_smooth, muscle_mask, roi_mesh_rstrct);
[allapo_rstrct_smoothed_fiber_all, allapo_rstrct_fiber_all_mm, allapo_rstrct_smoothed_fiber_all_mm, ~, ~, ~, ~, ~, ~] = ...
    fiber_smoother_v14a(allapo_rstrct_fiber_all, fs_options);

% Use same fiber-tract seed points in all the conditions evaluated
allapo_rstrct_smoothed_fiber_all_mm_final = zeros(size(allapo_rstrct_smoothed_fiber_all_mm));   %will hold tracts after seed-location filtering
allapo_rstrct_smoothed_fiber_all_final = zeros(size(allapo_rstrct_smoothed_fiber_all));

for step_cntr = 1:size(allapo_rstrct_smoothed_fiber_all,3)
    for dim_cntr = 1:size(allapo_rstrct_smoothed_fiber_all,4)
            allapo_rstrct_smoothed_fiber_all_final(:,:,step_cntr,dim_cntr) = ...
                squeeze(allapo_rstrct_smoothed_fiber_all(:,:,step_cntr,dim_cntr)).*allapo_rstrct_seeds_final_mask;
            allapo_rstrct_smoothed_fiber_all_mm_final (:,:,step_cntr,dim_cntr) = ...
                squeeze(allapo_rstrct_smoothed_fiber_all_mm(:,:,step_cntr,dim_cntr)).*allapo_rstrct_seeds_final_mask;
    end
end

[allapo_rstrct_angle_list_final, allapo_rstrct_distance_list_final, allapo_rstrct_curvature_list_final, ~, allapo_rstrct_n_points_final, allapo_rstrct_area_final] = ...
    fiber_quantifier_v20a(allapo_rstrct_fq_options, allapo_rstrct_smoothed_fiber_all_final, roi_mesh_rstrct, muscle_mask);

% calculate summary statistics
allapo_rstrct_total_distance_final = squeeze(max(allapo_rstrct_distance_list_final, [], 3));
allapo_rstrct_total_distance_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final = sum(squeeze(allapo_rstrct_angle_list_final(:,:,:,1)), 3)./allapo_rstrct_n_points_final(:,:,2);
allapo_rstrct_mean_alpha_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_curvature_final = sum(allapo_rstrct_curvature_list_final, 3)./allapo_rstrct_n_points_final(:,:,3);
allapo_rstrct_mean_curvature_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final(isnan(allapo_rstrct_mean_alpha_final)) = 0;
allapo_rstrct_mean_alpha_final = nonzeros(allapo_rstrct_mean_alpha_final);
allapo_rstrct_mean_curvature_final(isnan(allapo_rstrct_mean_curvature_final))=0;
allapo_rstrct_mean_curvature_final = nonzeros(allapo_rstrct_mean_curvature_final);
allapo_rstrct_total_distance_final(isnan(allapo_rstrct_total_distance_final))=0;
allapo_rstrct_total_distance_final = nonzeros(allapo_rstrct_total_distance_final);

mean_fiber_length_aniso_smoothn = mean(allapo_rstrct_total_distance_final)
mean_pennation_angle_aniso_smoothn = mean(allapo_rstrct_mean_alpha_final)
mean_curvature_aniso_smoothn = mean(allapo_rstrct_mean_curvature_final)

fiber_visualizer_v11(muscle_mask, fv_options, roi_mesh_rstrct, anat_img, allapo_rstrct_smoothed_fiber_all_final);
set(gcf, 'position', [1200 50 500 800])
set(gca, 'fontsize', 14)
xlabel('Column', 'fontsize', 14)
ylabel('Row', 'fontsize', 14)
zlabel('Slice', 'fontsize', 14)
title('Aponeurosis Seeding-Restricted Mesh, Aniso + Smoothn')

stop_pct_mask_aniso_smoothn = numel(find(allapo_rstrct_stop_list_aniso_smoothn == 4))/numel(allapo_rstrct_stop_list_aniso_smoothn)
stop_pct_fa_aniso_smoothn = numel(find(allapo_rstrct_stop_list_aniso_smoothn == 2))/numel(allapo_rstrct_stop_list_aniso_smoothn)
stop_pct_ang_aniso_smoothn = numel(find(allapo_rstrct_stop_list_aniso_smoothn == 3))/numel(allapo_rstrct_stop_list_aniso_smoothn)

%% Fiber-tracking - APO3 - Raw

% track, smooth, and quantify fibers
[allapo_rstrct_fiber_all, allapo_rstrct_roi_flag, allapo_rstrct_stop_list_raw, allapo_rstrct_fiber_len] = ...
    fiber_track_v20(allapo_rstrct_ft_options, e1fa_raw, muscle_mask, roi_mesh_rstrct);
[allapo_rstrct_smoothed_fiber_all, ~, allapo_rstrct_smoothed_fiber_all_mm, allapo_rstrct_pcoeff_r, allapo_rstrct_pcoeff_c, allapo_rstrct_pcoeff_s, allapo_rstrct_n_points_smoothed, allapo_rstrct_residuals, allapo_rstrct_residuals_mm] = ...
    fiber_smoother_v14a(allapo_rstrct_fiber_all, fs_options);

% Use same fiber-tract seed points in all the conditions evaluated
allapo_rstrct_smoothed_fiber_all_mm_final = zeros(size(allapo_rstrct_smoothed_fiber_all_mm));   %will hold tracts after seed-location filtering
allapo_rstrct_smoothed_fiber_all_final = zeros(size(allapo_rstrct_smoothed_fiber_all));

for step_cntr = 1:size(allapo_rstrct_smoothed_fiber_all,3)
    for dim_cntr = 1:size(allapo_rstrct_smoothed_fiber_all,4)
            allapo_rstrct_smoothed_fiber_all_final(:,:,step_cntr,dim_cntr) = ...
                squeeze(allapo_rstrct_smoothed_fiber_all(:,:,step_cntr,dim_cntr)).*allapo_rstrct_seeds_final_mask;
            allapo_rstrct_smoothed_fiber_all_mm_final (:,:,step_cntr,dim_cntr) = ...
                squeeze(allapo_rstrct_smoothed_fiber_all_mm(:,:,step_cntr,dim_cntr)).*allapo_rstrct_seeds_final_mask;
    end
end

[allapo_rstrct_angle_list_final, allapo_rstrct_distance_list_final, allapo_rstrct_curvature_list_final, ~, allapo_rstrct_n_points_final, allapo_rstrct_area_final] = ...
    fiber_quantifier_v20a(allapo_rstrct_fq_options, allapo_rstrct_smoothed_fiber_all_final, roi_mesh_rstrct, muscle_mask);

% calculate summary statistics
allapo_rstrct_total_distance_final = squeeze(max(allapo_rstrct_distance_list_final, [], 3));
allapo_rstrct_total_distance_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final = sum(squeeze(allapo_rstrct_angle_list_final(:,:,:,1)), 3)./allapo_rstrct_n_points_final(:,:,2);
allapo_rstrct_mean_alpha_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_curvature_final = sum(allapo_rstrct_curvature_list_final, 3)./allapo_rstrct_n_points_final(:,:,3);
allapo_rstrct_mean_curvature_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final(isnan(allapo_rstrct_mean_alpha_final)) = 0;
allapo_rstrct_mean_alpha_final = nonzeros(allapo_rstrct_mean_alpha_final);
allapo_rstrct_mean_curvature_final(isnan(allapo_rstrct_mean_curvature_final))=0;
allapo_rstrct_mean_curvature_final = nonzeros(allapo_rstrct_mean_curvature_final);
allapo_rstrct_total_distance_final(isnan(allapo_rstrct_total_distance_final))=0;
allapo_rstrct_total_distance_final = nonzeros(allapo_rstrct_total_distance_final);

mean_fiber_length_raw = mean(allapo_rstrct_total_distance_final)
mean_pennation_angle_raw = mean(allapo_rstrct_mean_alpha_final)
mean_curvature_raw = mean(allapo_rstrct_mean_curvature_final)

fiber_visualizer_v11(muscle_mask, fv_options, roi_mesh_rstrct, anat_img, allapo_rstrct_smoothed_fiber_all_final);
set(gcf, 'position', [1200 50 500 800])
set(gca, 'fontsize', 14)
xlabel('Column', 'fontsize', 14)
ylabel('Row', 'fontsize', 14)
zlabel('Slice', 'fontsize', 14)
title('Aponeurosis Seeding-Restricted Mesh, Raw')

stop_pct_mask_raw = numel(find(allapo_rstrct_stop_list_raw == 4))/numel(allapo_rstrct_stop_list_raw)
stop_pct_fa_raw = numel(find(allapo_rstrct_stop_list_raw == 2))/numel(allapo_rstrct_stop_list_raw)
stop_pct_ang_raw = numel(find(allapo_rstrct_stop_list_raw == 3))/numel(allapo_rstrct_stop_list_raw)

%% Fiber-tracking - APO3 - Smoothn

% track, smooth, and quantify fibers
[allapo_rstrct_fiber_all, ~, allapo_rstrct_stop_list_smoothn, ~] = ...
    fiber_track_v20(allapo_rstrct_ft_options, e1fa_raw_smooth, muscle_mask, roi_mesh_rstrct);
[allapo_rstrct_smoothed_fiber_all, ~, allapo_rstrct_smoothed_fiber_all_mm, ~, ~, ~, ~, ~, ~] = ...
    fiber_smoother_v14a(allapo_rstrct_fiber_all, fs_options);

% Use same fiber-tract seed points in all the conditions evaluated
allapo_rstrct_smoothed_fiber_all_mm_final = zeros(size(allapo_rstrct_smoothed_fiber_all_mm));   %will hold tracts after seed-location filtering
allapo_rstrct_smoothed_fiber_all_final = zeros(size(allapo_rstrct_smoothed_fiber_all));

for step_cntr = 1:size(allapo_rstrct_smoothed_fiber_all,3)
    for dim_cntr = 1:size(allapo_rstrct_smoothed_fiber_all,4)
            allapo_rstrct_smoothed_fiber_all_final(:,:,step_cntr,dim_cntr) = ...
                squeeze(allapo_rstrct_smoothed_fiber_all(:,:,step_cntr,dim_cntr)).*allapo_rstrct_seeds_final_mask;
            allapo_rstrct_smoothed_fiber_all_mm_final (:,:,step_cntr,dim_cntr) = ...
                squeeze(allapo_rstrct_smoothed_fiber_all_mm(:,:,step_cntr,dim_cntr)).*allapo_rstrct_seeds_final_mask;
    end
end

[allapo_rstrct_angle_list_final, allapo_rstrct_distance_list_final, allapo_rstrct_curvature_list_final, ~, allapo_rstrct_n_points_final, allapo_rstrct_area_final] = ...
    fiber_quantifier_v20a(allapo_rstrct_fq_options, allapo_rstrct_smoothed_fiber_all_final, roi_mesh_rstrct, muscle_mask);

% calculate summary statistics
allapo_rstrct_total_distance_final = squeeze(max(allapo_rstrct_distance_list_final, [], 3));
allapo_rstrct_total_distance_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final = sum(squeeze(allapo_rstrct_angle_list_final(:,:,:,1)), 3)./allapo_rstrct_n_points_final(:,:,2);
allapo_rstrct_mean_alpha_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_curvature_final = sum(allapo_rstrct_curvature_list_final, 3)./allapo_rstrct_n_points_final(:,:,3);
allapo_rstrct_mean_curvature_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final(isnan(allapo_rstrct_mean_alpha_final)) = 0;
allapo_rstrct_mean_alpha_final = nonzeros(allapo_rstrct_mean_alpha_final);
allapo_rstrct_mean_curvature_final(isnan(allapo_rstrct_mean_curvature_final))=0;
allapo_rstrct_mean_curvature_final = nonzeros(allapo_rstrct_mean_curvature_final);
allapo_rstrct_total_distance_final(isnan(allapo_rstrct_total_distance_final))=0;
allapo_rstrct_total_distance_final = nonzeros(allapo_rstrct_total_distance_final);

mean_fiber_length_smoothn = mean(allapo_rstrct_total_distance_final)
mean_pennation_angle_smoothn = mean(allapo_rstrct_mean_alpha_final)
mean_curvature_smoothn = mean(allapo_rstrct_mean_curvature_final)

fiber_visualizer_v11(muscle_mask, fv_options, roi_mesh_rstrct, anat_img, allapo_rstrct_smoothed_fiber_all_final);
set(gcf, 'position', [1200 50 500 800])
set(gca, 'fontsize', 14)
xlabel('Column', 'fontsize', 14)
ylabel('Row', 'fontsize', 14)
zlabel('Slice', 'fontsize', 14)
title('Aponeurosis Seeding-Restricted Mesh, Smoothn')

stop_pct_mask_smoothn = numel(find(allapo_rstrct_stop_list_smoothn == 4))/numel(allapo_rstrct_stop_list_smoothn)
stop_pct_fa_smoothn = numel(find(allapo_rstrct_stop_list_smoothn == 2))/numel(allapo_rstrct_stop_list_smoothn)
stop_pct_ang_smoothn = numel(find(allapo_rstrct_stop_list_smoothn == 3))/numel(allapo_rstrct_stop_list_smoothn)

%% Fiber-tracking - APO3 - Aniso

% track, smooth, and quantify fibers
[allapo_rstrct_fiber_all, ~, allapo_rstrct_stop_list_aniso, ~] = ...
    fiber_track_v20(allapo_rstrct_ft_options, e1fa_aniso, muscle_mask, roi_mesh_rstrct);
[allapo_rstrct_smoothed_fiber_all, ~, allapo_rstrct_smoothed_fiber_all_mm, ~, ~, ~, ~, ~, ~] = ...
    fiber_smoother_v14a(allapo_rstrct_fiber_all, fs_options);

% Use same fiber-tract seed points in all the conditions evaluated
allapo_rstrct_smoothed_fiber_all_mm_final = zeros(size(allapo_rstrct_smoothed_fiber_all_mm));   %will hold tracts after seed-location filtering
allapo_rstrct_smoothed_fiber_all_final = zeros(size(allapo_rstrct_smoothed_fiber_all));

for step_cntr = 1:size(allapo_rstrct_smoothed_fiber_all,3)
    for dim_cntr = 1:size(allapo_rstrct_smoothed_fiber_all,4)
            allapo_rstrct_smoothed_fiber_all_final(:,:,step_cntr,dim_cntr) = ...
                squeeze(allapo_rstrct_smoothed_fiber_all(:,:,step_cntr,dim_cntr)).*allapo_rstrct_seeds_final_mask;
            allapo_rstrct_smoothed_fiber_all_mm_final (:,:,step_cntr,dim_cntr) = ...
                squeeze(allapo_rstrct_smoothed_fiber_all_mm(:,:,step_cntr,dim_cntr)).*allapo_rstrct_seeds_final_mask;
    end
end

[allapo_rstrct_angle_list_final, allapo_rstrct_distance_list_final, allapo_rstrct_curvature_list_final, ~, allapo_rstrct_n_points_final, allapo_rstrct_area_final] = ...
    fiber_quantifier_v20a(allapo_rstrct_fq_options, allapo_rstrct_smoothed_fiber_all_final, roi_mesh_rstrct, muscle_mask);

% calculate summary statistics
allapo_rstrct_total_distance_final = squeeze(max(allapo_rstrct_distance_list_final, [], 3));
allapo_rstrct_total_distance_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final = sum(squeeze(allapo_rstrct_angle_list_final(:,:,:,1)), 3)./allapo_rstrct_n_points_final(:,:,2);
allapo_rstrct_mean_alpha_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_curvature_final = sum(allapo_rstrct_curvature_list_final, 3)./allapo_rstrct_n_points_final(:,:,3);
allapo_rstrct_mean_curvature_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final(isnan(allapo_rstrct_mean_alpha_final)) = 0;
allapo_rstrct_mean_alpha_final = nonzeros(allapo_rstrct_mean_alpha_final);
allapo_rstrct_mean_curvature_final(isnan(allapo_rstrct_mean_curvature_final))=0;
allapo_rstrct_mean_curvature_final = nonzeros(allapo_rstrct_mean_curvature_final);
allapo_rstrct_total_distance_final(isnan(allapo_rstrct_total_distance_final))=0;
allapo_rstrct_total_distance_final = nonzeros(allapo_rstrct_total_distance_final);

mean_fiber_length_aniso = mean(allapo_rstrct_total_distance_final)
mean_pennation_angle_aniso = mean(allapo_rstrct_mean_alpha_final)
mean_curvature_aniso = mean(allapo_rstrct_mean_curvature_final)

fiber_visualizer_v11(muscle_mask, fv_options, roi_mesh_rstrct, anat_img, allapo_rstrct_smoothed_fiber_all_final);
set(gcf, 'position', [1200 50 500 800])
set(gca, 'fontsize', 14)
xlabel('Column', 'fontsize', 14)
ylabel('Row', 'fontsize', 14)
zlabel('Slice', 'fontsize', 14)
title('Aponeurosis Seeding-Restricted Mesh, Aniso')

stop_pct_mask_aniso = numel(find(allapo_rstrct_stop_list_aniso == 4))/numel(allapo_rstrct_stop_list_aniso)
stop_pct_fa_aniso = numel(find(allapo_rstrct_stop_list_aniso == 2))/numel(allapo_rstrct_stop_list_aniso)
stop_pct_ang_aniso = numel(find(allapo_rstrct_stop_list_aniso == 3))/numel(allapo_rstrct_stop_list_aniso)

%% Fiber-tracking - APO3 - TPCA

% track, smooth, and quantify fibers
[allapo_rstrct_fiber_all, ~, allapo_rstrct_stop_list_tpca, ~] = ...
    fiber_track_v20(allapo_rstrct_ft_options, e1fa_tpca, muscle_mask, roi_mesh_rstrct);
[allapo_rstrct_smoothed_fiber_all, ~, allapo_rstrct_smoothed_fiber_all_mm, ~, ~, ~, ~, ~, ~] = ...
    fiber_smoother_v14a(allapo_rstrct_fiber_all, fs_options);

% Use same fiber-tract seed points in all the conditions evaluated
allapo_rstrct_smoothed_fiber_all_mm_final = zeros(size(allapo_rstrct_smoothed_fiber_all_mm));   %will hold tracts after seed-location filtering
allapo_rstrct_smoothed_fiber_all_final = zeros(size(allapo_rstrct_smoothed_fiber_all));

for step_cntr = 1:size(allapo_rstrct_smoothed_fiber_all,3)
    for dim_cntr = 1:size(allapo_rstrct_smoothed_fiber_all,4)
            allapo_rstrct_smoothed_fiber_all_final(:,:,step_cntr,dim_cntr) = ...
                squeeze(allapo_rstrct_smoothed_fiber_all(:,:,step_cntr,dim_cntr)).*allapo_rstrct_seeds_final_mask;
            allapo_rstrct_smoothed_fiber_all_mm_final (:,:,step_cntr,dim_cntr) = ...
                squeeze(allapo_rstrct_smoothed_fiber_all_mm(:,:,step_cntr,dim_cntr)).*allapo_rstrct_seeds_final_mask;
    end
end

[allapo_rstrct_angle_list_final, allapo_rstrct_distance_list_final, allapo_rstrct_curvature_list_final, ~, allapo_rstrct_n_points_final, allapo_rstrct_area_final] = ...
    fiber_quantifier_v20a(allapo_rstrct_fq_options, allapo_rstrct_smoothed_fiber_all_final, roi_mesh_rstrct, muscle_mask);

% calculate summary statistics
allapo_rstrct_total_distance_final = squeeze(max(allapo_rstrct_distance_list_final, [], 3));
allapo_rstrct_total_distance_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final = sum(squeeze(allapo_rstrct_angle_list_final(:,:,:,1)), 3)./allapo_rstrct_n_points_final(:,:,2);
allapo_rstrct_mean_alpha_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_curvature_final = sum(allapo_rstrct_curvature_list_final, 3)./allapo_rstrct_n_points_final(:,:,3);
allapo_rstrct_mean_curvature_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final(isnan(allapo_rstrct_mean_alpha_final)) = 0;
allapo_rstrct_mean_alpha_final = nonzeros(allapo_rstrct_mean_alpha_final);
allapo_rstrct_mean_curvature_final(isnan(allapo_rstrct_mean_curvature_final))=0;
allapo_rstrct_mean_curvature_final = nonzeros(allapo_rstrct_mean_curvature_final);
allapo_rstrct_total_distance_final(isnan(allapo_rstrct_total_distance_final))=0;
allapo_rstrct_total_distance_final = nonzeros(allapo_rstrct_total_distance_final);

mean_fiber_length_tpca = mean(allapo_rstrct_total_distance_final)
mean_pennation_angle_tpca = mean(allapo_rstrct_mean_alpha_final)
mean_curvature_tpca = mean(allapo_rstrct_mean_curvature_final)

fiber_visualizer_v11(muscle_mask, fv_options, roi_mesh_rstrct, anat_img, allapo_rstrct_smoothed_fiber_all_final);
set(gcf, 'position', [1200 50 500 800])
set(gca, 'fontsize', 14)
xlabel('Column', 'fontsize', 14)
ylabel('Row', 'fontsize', 14)
zlabel('Slice', 'fontsize', 14)
title('Aponeurosis Seeding-Restricted Mesh, TPCA')

stop_pct_mask_tpca = numel(find(allapo_rstrct_stop_list_tpca == 4))/numel(allapo_rstrct_stop_list_tpca)
stop_pct_fa_tpca = numel(find(allapo_rstrct_stop_list_tpca == 2))/numel(allapo_rstrct_stop_list_tpca)
stop_pct_ang_tpca = numel(find(allapo_rstrct_stop_list_tpca == 3))/numel(allapo_rstrct_stop_list_tpca)

%% Fiber-tracking - APO3 - TPCA + smoothn

% track, smooth, and quantify fibers
[allapo_rstrct_fiber_all, ~, allapo_rstrct_stop_list_tpca_smooth, ~] = ...
    fiber_track_v20(allapo_rstrct_ft_options, e1fa_tpca_smooth, muscle_mask, roi_mesh_rstrct);
[allapo_rstrct_smoothed_fiber_all, ~, allapo_rstrct_smoothed_fiber_all_mm, ~, ~, ~, ~, ~, ~] = ...
    fiber_smoother_v14a(allapo_rstrct_fiber_all, fs_options);

% Use same fiber-tract seed points in all the conditions evaluated
allapo_rstrct_smoothed_fiber_all_mm_final = zeros(size(allapo_rstrct_smoothed_fiber_all_mm));   %will hold tracts after seed-location filtering
allapo_rstrct_smoothed_fiber_all_final = zeros(size(allapo_rstrct_smoothed_fiber_all));

for step_cntr = 1:size(allapo_rstrct_smoothed_fiber_all,3)
    for dim_cntr = 1:size(allapo_rstrct_smoothed_fiber_all,4)
            allapo_rstrct_smoothed_fiber_all_final(:,:,step_cntr,dim_cntr) = ...
                squeeze(allapo_rstrct_smoothed_fiber_all(:,:,step_cntr,dim_cntr)).*allapo_rstrct_seeds_final_mask;
            allapo_rstrct_smoothed_fiber_all_mm_final (:,:,step_cntr,dim_cntr) = ...
                squeeze(allapo_rstrct_smoothed_fiber_all_mm(:,:,step_cntr,dim_cntr)).*allapo_rstrct_seeds_final_mask;
    end
end

[allapo_rstrct_angle_list_final, allapo_rstrct_distance_list_final, allapo_rstrct_curvature_list_final, ~, allapo_rstrct_n_points_final, allapo_rstrct_area_final] = ...
    fiber_quantifier_v20a(allapo_rstrct_fq_options, allapo_rstrct_smoothed_fiber_all_final, roi_mesh_rstrct, muscle_mask);

% calculate summary statistics
allapo_rstrct_total_distance_final = squeeze(max(allapo_rstrct_distance_list_final, [], 3));
allapo_rstrct_total_distance_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final = sum(squeeze(allapo_rstrct_angle_list_final(:,:,:,1)), 3)./allapo_rstrct_n_points_final(:,:,2);
allapo_rstrct_mean_alpha_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_curvature_final = sum(allapo_rstrct_curvature_list_final, 3)./allapo_rstrct_n_points_final(:,:,3);
allapo_rstrct_mean_curvature_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final(isnan(allapo_rstrct_mean_alpha_final)) = 0;
allapo_rstrct_mean_alpha_final = nonzeros(allapo_rstrct_mean_alpha_final);
allapo_rstrct_mean_curvature_final(isnan(allapo_rstrct_mean_curvature_final))=0;
allapo_rstrct_mean_curvature_final = nonzeros(allapo_rstrct_mean_curvature_final);
allapo_rstrct_total_distance_final(isnan(allapo_rstrct_total_distance_final))=0;
allapo_rstrct_total_distance_final = nonzeros(allapo_rstrct_total_distance_final);

mean_fiber_length_tpca_smooth = mean(allapo_rstrct_total_distance_final)
mean_pennation_angle_tpca_smooth = mean(allapo_rstrct_mean_alpha_final)
mean_curvature_tpca_smooth = mean(allapo_rstrct_mean_curvature_final)

fiber_visualizer_v11(muscle_mask, fv_options, roi_mesh_rstrct, anat_img, allapo_rstrct_smoothed_fiber_all_final);
set(gcf, 'position', [1200 50 500 800])
set(gca, 'fontsize', 14)
xlabel('Column', 'fontsize', 14)
ylabel('Row', 'fontsize', 14)
zlabel('Slice', 'fontsize', 14)
title('Aponeurosis Seeding-Restricted Mesh, TPCA + smoothn')

stop_pct_mask_tpca_smooth = numel(find(allapo_rstrct_stop_list_tpca_smooth == 4))/numel(allapo_rstrct_stop_list_tpca_smooth)
stop_pct_fa_tpca_smooth = numel(find(allapo_rstrct_stop_list_tpca_smooth == 2))/numel(allapo_rstrct_stop_list_tpca_smooth)
stop_pct_ang_tpca_smooth = numel(find(allapo_rstrct_stop_list_tpca_smooth == 3))/numel(allapo_rstrct_stop_list_tpca_smooth)

%% Analyze b = 0 (1 avg), b = 450 (1 avg) data

%% 
% Compute FA and first eigenvector map of DTI data

%Correct diffusion directions fror FFS scans
diff_dir=diff_param_1.bvec; %diff_param1 and diff_param2 are the same

b_val=mean(diff_param_1.bval(2:end));

%Tensor computation

FA_raw = zeros(size(dti_img_4d_reg_nonavg,1),size(dti_img_4d_reg_nonavg,2),size(dti_img_4d_reg_nonavg,3));
E1_450_raw = zeros(size(dti_img_4d_reg_nonavg,1),size(dti_img_4d_reg_nonavg,2),size(dti_img_4d_reg_nonavg,3),3);
FA_aniso = zeros(size(dti_img_4d_reg_nonavg,1),size(dti_img_4d_reg_nonavg,2),size(dti_img_4d_reg_nonavg,3));
E1_450_aniso = zeros(size(dti_img_4d_reg_nonavg,1),size(dti_img_4d_reg_nonavg,2),size(dti_img_4d_reg_nonavg,3),3);
FA_tpca = zeros(size(dti_img_4d_reg_nonavg,1),size(dti_img_4d_reg_nonavg,2),size(dti_img_4d_reg_nonavg,3));
E1_450_tpca = zeros(size(dti_img_4d_reg_nonavg,1),size(dti_img_4d_reg_nonavg,2),size(dti_img_4d_reg_nonavg,3),3);

for s = 1:size(dti_img_4d_reg_nonavg,3)
    for r = 1:size(dti_img_4d_reg_nonavg,1)
        for c = 1:size(dti_img_4d_reg_nonavg,2)
                dir_m = diff_dir;  %diffusion encoding matrix
                dir_m = dir_m(2:end,:);

                %estimate tensor
                %use TA mask to make process faster
                if muscle_mask(r,c,s)>0
                    % Raw data
                    sig_v = squeeze(dti_img_4d_reg_nonavg(r,c,s,:));    %signal vector
                    D = signal2tensor2(sig_v, dir_m, b_val);
                    [E,L] = svd(D);
    
                    L_v = diag(L);
                    md = mean(L_v);
                    fa = sqrt(3/2) * sqrt(((L_v(1) - md)^2 + (L_v(2) - md)^2 + (L_v(3) - md)^2) / ...
                        (L_v(1)^2 + L_v(2)^2 + L_v(3)^2));

                    loop_E1 = E(:,L_v==max(L_v));
                    if loop_E1(3)<0
                        loop_E1 = -loop_E1;
                    end
                    
                    E1_450_raw(r,c,s,:) = loop_E1;
                    FA_raw(r,c,s,:) = fa;

                    % AIS
                    sig_v = squeeze(dti_img_4d_smooth_nonavg(r,c,s,:));    %signal vector
                    D = signal2tensor2(sig_v, dir_m, b_val);
                    [E,L] = svd(D);
    
                    L_v = diag(L);
                    md = mean(L_v);
                    fa = sqrt(3/2) * sqrt(((L_v(1) - md)^2 + (L_v(2) - md)^2 + (L_v(3) - md)^2) / ...
                        (L_v(1)^2 + L_v(2)^2 + L_v(3)^2));

                    loop_E1 = E(:,L_v==max(L_v));
                    if loop_E1(3)<0
                        loop_E1 = -loop_E1;
                    end
                    
                    E1_450_aniso(r,c,s,:) = loop_E1;
                    FA_aniso(r,c,s,:) = fa;

                    % TPCA
                    sig_v = squeeze(dti_img_4d_tpca_nonavg(r,c,s,:));    %signal vector
                    D = signal2tensor2(sig_v, dir_m, b_val);
                    [E,L] = svd(D);
    
                    L_v = diag(L);
                    md = mean(L_v);
                    fa = sqrt(3/2) * sqrt(((L_v(1) - md)^2 + (L_v(2) - md)^2 + (L_v(3) - md)^2) / ...
                        (L_v(1)^2 + L_v(2)^2 + L_v(3)^2));

                    loop_E1 = E(:,L_v==max(L_v));
                    if loop_E1(3)<0
                        loop_E1 = -loop_E1;
                    end
                    
                    E1_450_tpca(r,c,s,:) = loop_E1;
                    FA_tpca(r,c,s,:) = fa;

                end
        end
    end
end

%% 
% Smooth first eigenvector field of TA muscle

[E1_450_raw_smooth,smooth_parameter_raw]=smoothn({E1_450_raw(:,:,:,1),E1_450_raw(:,:,:,2),E1_450_raw(:,:,:,3)});
[E1_450_aniso_smooth,smooth_parameter_aniso]=smoothn({E1_450_aniso(:,:,:,1),E1_450_aniso(:,:,:,2),E1_450_aniso(:,:,:,3)});
[E1_450_tpca_smooth,smooth_parameter_tpca]=smoothn({E1_450_tpca(:,:,:,1),E1_450_tpca(:,:,:,2),E1_450_tpca(:,:,:,3)});

E1_450_raw_smooth=cat(4,E1_450_raw_smooth{1},E1_450_raw_smooth{2},E1_450_raw_smooth{3});
E1_450_aniso_smooth=cat(4,E1_450_aniso_smooth{1},E1_450_aniso_smooth{2},E1_450_aniso_smooth{3});
E1_450_tpca_smooth=cat(4,E1_450_tpca_smooth{1},E1_450_tpca_smooth{2},E1_450_tpca_smooth{3});

for r_img=1:size(E1_450_aniso,1)
    for c_img=1:size(E1_450_aniso,2)
        for s_img=1:size(E1_450_aniso,3)
            if muscle_mask(r_img,c_img,s_img)>0

                E1_raw_smooth = squeeze(E1_450_raw_smooth(r_img,c_img,s_img,:));
                E1_aniso_smooth = squeeze(E1_450_aniso_smooth(r_img,c_img,s_img,:));
                E1_tpca_smooth = squeeze(E1_450_tpca_smooth(r_img,c_img,s_img,:));

                E1_raw_smooth = E1_raw_smooth/norm(E1_raw_smooth);
                E1_aniso_smooth = E1_aniso_smooth/norm(E1_aniso_smooth);
                E1_tpca_smooth = E1_tpca_smooth/norm(E1_tpca_smooth);

                E1_450_raw_smooth(r_img,c_img,s_img,:)=E1_raw_smooth;
                E1_450_aniso_smooth(r_img,c_img,s_img,:)=E1_aniso_smooth;
                E1_450_tpca_smooth(r_img,c_img,s_img,:)=E1_tpca_smooth;
            end
        end
    end
end

for i = 1:3
    E1_450_raw_smooth(:,:,:,i)=squeeze(E1_450_raw_smooth(:,:,:,i)).*muscle_mask;
    E1_450_aniso_smooth(:,:,:,i)=squeeze(E1_450_aniso_smooth(:,:,:,i)).*muscle_mask;
    E1_450_tpca_smooth(:,:,:,i)=squeeze(E1_450_tpca_smooth(:,:,:,i)).*muscle_mask;
end

%% 
% Create E1-FA maps for fiber tractography

e1fa_raw = cat(4,E1_450_raw,FA_raw);
e1fa_aniso = cat(4,E1_450_aniso,FA_aniso);
e1fa_tpca = cat(4,E1_450_tpca,FA_tpca);
e1fa_raw_smooth = cat(4,E1_450_raw_smooth,FA_raw);
e1fa_aniso_smooth = cat(4,E1_450_aniso_smooth,FA_aniso);
e1fa_tpca_smooth = cat(4,E1_450_tpca_smooth,FA_tpca);

%%
% Visualize first eigenvector maps

rect_crop = [39, 22, 41, 41];

E1_x_raw=squeeze(E1_450_raw(:,:,:,1));
E1_y_raw=squeeze(E1_450_raw(:,:,:,2));
E1_z_raw=squeeze(E1_450_raw(:,:,:,3));

E1_x_raw_Crop=imcrop3(E1_x_raw,[rect_crop(1),rect_crop(2),1,rect_crop(3),rect_crop(4),size(E1_450_raw,3)-1]);
E1_y_raw_Crop=imcrop3(E1_y_raw,[rect_crop(1),rect_crop(2),1,rect_crop(3),rect_crop(4),size(E1_450_raw,3)-1]);
E1_z_raw_Crop=imcrop3(E1_z_raw,[rect_crop(1),rect_crop(2),1,rect_crop(3),rect_crop(4),size(E1_450_raw,3)-1]);

E1_x_raw_Crop(E1_x_raw_Crop==0)=NaN;
E1_y_raw_Crop(E1_y_raw_Crop==0)=NaN;
E1_z_raw_Crop(E1_z_raw_Crop==0)=NaN;

E1_x_aniso=squeeze(E1_450_aniso(:,:,:,1));
E1_y_aniso=squeeze(E1_450_aniso(:,:,:,2));
E1_z_aniso=squeeze(E1_450_aniso(:,:,:,3));

E1_x_aniso_Crop=imcrop3(E1_x_aniso,[rect_crop(1),rect_crop(2),1,rect_crop(3),rect_crop(4),size(E1_450_raw,3)-1]);
E1_y_aniso_Crop=imcrop3(E1_y_aniso,[rect_crop(1),rect_crop(2),1,rect_crop(3),rect_crop(4),size(E1_450_raw,3)-1]);
E1_z_aniso_Crop=imcrop3(E1_z_aniso,[rect_crop(1),rect_crop(2),1,rect_crop(3),rect_crop(4),size(E1_450_raw,3)-1]);

E1_x_aniso_Crop(E1_x_aniso_Crop==0)=NaN;
E1_y_aniso_Crop(E1_y_aniso_Crop==0)=NaN;
E1_z_aniso_Crop(E1_z_aniso_Crop==0)=NaN;

E1_x_raw_smooth=squeeze(E1_450_raw_smooth(:,:,:,1));
E1_y_raw_smooth=squeeze(E1_450_raw_smooth(:,:,:,2));
E1_z_raw_smooth=squeeze(E1_450_raw_smooth(:,:,:,3));

E1_x_raw_smooth_Crop=imcrop3(E1_x_raw_smooth,[rect_crop(1),rect_crop(2),1,rect_crop(3),rect_crop(4),size(E1_450_raw,3)-1]);
E1_y_raw_smooth_Crop=imcrop3(E1_y_raw_smooth,[rect_crop(1),rect_crop(2),1,rect_crop(3),rect_crop(4),size(E1_450_raw,3)-1]);
E1_z_raw_smooth_Crop=imcrop3(E1_z_raw_smooth,[rect_crop(1),rect_crop(2),1,rect_crop(3),rect_crop(4),size(E1_450_raw,3)-1]);

E1_x_raw_smooth_Crop(E1_x_raw_smooth_Crop==0)=NaN;
E1_y_raw_smooth_Crop(E1_y_raw_smooth_Crop==0)=NaN;
E1_z_raw_smooth_Crop(E1_z_raw_smooth_Crop==0)=NaN;

E1_x_aniso_smooth=squeeze(E1_450_aniso_smooth(:,:,:,1));
E1_y_aniso_smooth=squeeze(E1_450_aniso_smooth(:,:,:,2));
E1_z_aniso_smooth=squeeze(E1_450_aniso_smooth(:,:,:,3));

E1_x_aniso_smooth_Crop=imcrop3(E1_x_aniso_smooth,[rect_crop(1),rect_crop(2),1,rect_crop(3),rect_crop(4),size(E1_450_raw,3)-1]);
E1_y_aniso_smooth_Crop=imcrop3(E1_y_aniso_smooth,[rect_crop(1),rect_crop(2),1,rect_crop(3),rect_crop(4),size(E1_450_raw,3)-1]);
E1_z_aniso_smooth_Crop=imcrop3(E1_z_aniso_smooth,[rect_crop(1),rect_crop(2),1,rect_crop(3),rect_crop(4),size(E1_450_raw,3)-1]);

E1_x_aniso_smooth_Crop(E1_x_aniso_smooth_Crop==0)=NaN;
E1_y_aniso_smooth_Crop(E1_y_aniso_smooth_Crop==0)=NaN;
E1_z_aniso_smooth_Crop(E1_z_aniso_smooth_Crop==0)=NaN;

E1_x_tpca=squeeze(E1_450_tpca(:,:,:,1));
E1_y_tpca=squeeze(E1_450_tpca(:,:,:,2));
E1_z_tpca=squeeze(E1_450_tpca(:,:,:,3));

E1_x_tpca_Crop=imcrop3(E1_x_tpca,[rect_crop(1),rect_crop(2),1,rect_crop(3),rect_crop(4),size(E1_450_raw,3)-1]);
E1_y_tpca_Crop=imcrop3(E1_y_tpca,[rect_crop(1),rect_crop(2),1,rect_crop(3),rect_crop(4),size(E1_450_raw,3)-1]);
E1_z_tpca_Crop=imcrop3(E1_z_tpca,[rect_crop(1),rect_crop(2),1,rect_crop(3),rect_crop(4),size(E1_450_raw,3)-1]);

E1_x_tpca_Crop(E1_x_tpca_Crop==0)=NaN;
E1_y_tpca_Crop(E1_y_tpca_Crop==0)=NaN;
E1_z_tpca_Crop(E1_z_tpca_Crop==0)=NaN;

E1_x_tpca_smooth=squeeze(E1_450_tpca_smooth(:,:,:,1));
E1_y_tpca_smooth=squeeze(E1_450_tpca_smooth(:,:,:,2));
E1_z_tpca_smooth=squeeze(E1_450_tpca_smooth(:,:,:,3));

E1_x_tpca_smooth_Crop=imcrop3(E1_x_tpca_smooth,[rect_crop(1),rect_crop(2),1,rect_crop(3),rect_crop(4),size(E1_450_raw,3)-1]);
E1_y_tpca_smooth_Crop=imcrop3(E1_y_tpca_smooth,[rect_crop(1),rect_crop(2),1,rect_crop(3),rect_crop(4),size(E1_450_raw,3)-1]);
E1_z_tpca_smooth_Crop=imcrop3(E1_z_tpca_smooth,[rect_crop(1),rect_crop(2),1,rect_crop(3),rect_crop(4),size(E1_450_raw,3)-1]);

E1_x_tpca_smooth_Crop(E1_x_tpca_smooth_Crop==0)=NaN;
E1_y_tpca_smooth_Crop(E1_y_tpca_smooth_Crop==0)=NaN;
E1_z_tpca_smooth_Crop(E1_z_tpca_smooth_Crop==0)=NaN;

%% 
% Obtain image of first eigenvector map of muscle cross-section

figure
t = tiledlayout(3,6);

nexttile
imagesc(E1_x_raw_Crop(:,:,30),'AlphaData',~isnan(E1_x_raw_Crop(:,:,30)),[-0.3 0.3]); pbaspect([1 1 1]); colorbar
title('E1_x - TA - Raw')

nexttile
imagesc(E1_x_aniso_Crop(:,:,30),'AlphaData',~isnan(E1_x_aniso_Crop(:,:,30)),[-0.3 0.3]); pbaspect([1 1 1])
title('E1_x - TA - AnisoSmooth')

nexttile
imagesc(E1_x_tpca_Crop(:,:,30),'AlphaData',~isnan(E1_x_tpca_Crop(:,:,30)),[-0.3 0.3]); pbaspect([1 1 1])
title('E1_x - TA - TPCA')

nexttile
imagesc(E1_x_raw_smooth_Crop(:,:,30),'AlphaData',~isnan(E1_x_raw_smooth_Crop(:,:,30)),[-0.3 0.3]); pbaspect([1 1 1])
title('E1_x - TA - Smoothn')

nexttile
imagesc(E1_x_aniso_smooth_Crop(:,:,30),'AlphaData',~isnan(E1_x_aniso_smooth_Crop(:,:,30)),[-0.3 0.3]); pbaspect([1 1 1])
title('E1_x - TA - AnisoSmooth + Smoothn')

nexttile
imagesc(E1_x_tpca_smooth_Crop(:,:,30),'AlphaData',~isnan(E1_x_tpca_smooth_Crop(:,:,30)),[-0.3 0.3]); pbaspect([1 1 1])
title('E1_x - TA - TPCA + Smoothn')

nexttile
imagesc(E1_y_raw_Crop(:,:,30),'AlphaData',~isnan(E1_y_raw_Crop(:,:,30)),[-0.5 0.5]); pbaspect([1 1 1]); colorbar
title('E1_y - TA - Raw')

nexttile
imagesc(E1_y_aniso_Crop(:,:,30),'AlphaData',~isnan(E1_y_aniso_Crop(:,:,30)),[-0.5 0.5]); pbaspect([1 1 1])
title('E1_y - TA - AnisoSmooth')

nexttile
imagesc(E1_y_tpca_Crop(:,:,30),'AlphaData',~isnan(E1_y_tpca_Crop(:,:,30)),[-0.5 0.5]); pbaspect([1 1 1])
title('E1_y - TA - TPCA')

nexttile
imagesc(E1_y_raw_smooth_Crop(:,:,30),'AlphaData',~isnan(E1_y_raw_smooth_Crop(:,:,30)),[-0.5 0.5]); pbaspect([1 1 1])
title('E1_y - TA - Smoothn')

nexttile
imagesc(E1_y_aniso_smooth_Crop(:,:,30),'AlphaData',~isnan(E1_y_aniso_smooth_Crop(:,:,30)),[-0.5 0.5]); pbaspect([1 1 1])
title('E1_y - TA - AnisoSmooth + Smoothn')

nexttile
imagesc(E1_y_tpca_smooth_Crop(:,:,30),'AlphaData',~isnan(E1_y_tpca_smooth_Crop(:,:,30)),[-0.5 0.5]); pbaspect([1 1 1])
title('E1_y - TA - TPCA + Smoothn')

nexttile
imagesc(E1_z_raw_Crop(:,:,30),'AlphaData',~isnan(E1_z_raw_Crop(:,:,30)),[0.9 1]); pbaspect([1 1 1]); colorbar
title('E1_z - TA - Raw')

nexttile
imagesc(E1_z_aniso_Crop(:,:,30),'AlphaData',~isnan(E1_z_aniso_Crop(:,:,30)),[0.9 1]); pbaspect([1 1 1])
title('E1_z - TA - AnisoSmooth')

nexttile
imagesc(E1_z_tpca_Crop(:,:,30),'AlphaData',~isnan(E1_z_tpca_Crop(:,:,30)),[0.9 1]); pbaspect([1 1 1])
title('E1_z - TA - TPCA')

nexttile
imagesc(E1_z_raw_smooth_Crop(:,:,30),'AlphaData',~isnan(E1_z_raw_smooth_Crop(:,:,30)),[0.9 1]); pbaspect([1 1 1])
title('E1_z - TA - Smoothn')

nexttile
imagesc(E1_z_aniso_smooth_Crop(:,:,30),'AlphaData',~isnan(E1_z_aniso_smooth_Crop(:,:,30)),[0.9 1]); pbaspect([1 1 1])
title('E1_z - TA - AnisoSmooth + Smoothn')

nexttile
imagesc(E1_z_tpca_smooth_Crop(:,:,30),'AlphaData',~isnan(E1_z_tpca_smooth_Crop(:,:,30)),[0.9 1]); pbaspect([1 1 1])
title('E1_z - TA - TPCA + Smoothn')

cd('S:\Muscle_DTI\Roberto\DTI_Muscle_Analysis\Processed_Data\Aim1E')
print(gcf,'TA_eigenvector_Images_b0_1_b450_1.png','-dpng','-r1500'); 
savefig('TA_eigenvector_Images_b0_1_b450_1.fig')

%% Fiber-tracking - APO3 - Aniso4D + Smoothn

% track, smooth, and quantify fibers
[allapo_rstrct_fiber_all, ~, allapo_rstrct_stop_list_aniso_smoothn, ~] = ...
    fiber_track_v20(allapo_rstrct_ft_options, e1fa_aniso_smooth, muscle_mask, roi_mesh_rstrct);
[allapo_rstrct_smoothed_fiber_all, allapo_rstrct_fiber_all_mm, allapo_rstrct_smoothed_fiber_all_mm, ~, ~, ~, ~, ~, ~] = ...
    fiber_smoother_v14a(allapo_rstrct_fiber_all, fs_options);

% Use same fiber-tract seed points in all the conditions evaluated
allapo_rstrct_smoothed_fiber_all_mm_final = zeros(size(allapo_rstrct_smoothed_fiber_all_mm));   %will hold tracts after seed-location filtering
allapo_rstrct_smoothed_fiber_all_final = zeros(size(allapo_rstrct_smoothed_fiber_all));

for step_cntr = 1:size(allapo_rstrct_smoothed_fiber_all,3)
    for dim_cntr = 1:size(allapo_rstrct_smoothed_fiber_all,4)
            allapo_rstrct_smoothed_fiber_all_final(:,:,step_cntr,dim_cntr) = ...
                squeeze(allapo_rstrct_smoothed_fiber_all(:,:,step_cntr,dim_cntr)).*allapo_rstrct_seeds_final_mask;
            allapo_rstrct_smoothed_fiber_all_mm_final (:,:,step_cntr,dim_cntr) = ...
                squeeze(allapo_rstrct_smoothed_fiber_all_mm(:,:,step_cntr,dim_cntr)).*allapo_rstrct_seeds_final_mask;
    end
end

[allapo_rstrct_angle_list_final, allapo_rstrct_distance_list_final, allapo_rstrct_curvature_list_final, ~, allapo_rstrct_n_points_final, allapo_rstrct_area_final] = ...
    fiber_quantifier_v20a(allapo_rstrct_fq_options, allapo_rstrct_smoothed_fiber_all_final, roi_mesh_rstrct, muscle_mask);

% calculate summary statistics
allapo_rstrct_total_distance_final = squeeze(max(allapo_rstrct_distance_list_final, [], 3));
allapo_rstrct_total_distance_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final = sum(squeeze(allapo_rstrct_angle_list_final(:,:,:,1)), 3)./allapo_rstrct_n_points_final(:,:,2);
allapo_rstrct_mean_alpha_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_curvature_final = sum(allapo_rstrct_curvature_list_final, 3)./allapo_rstrct_n_points_final(:,:,3);
allapo_rstrct_mean_curvature_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final(isnan(allapo_rstrct_mean_alpha_final)) = 0;
allapo_rstrct_mean_alpha_final = nonzeros(allapo_rstrct_mean_alpha_final);
allapo_rstrct_mean_curvature_final(isnan(allapo_rstrct_mean_curvature_final))=0;
allapo_rstrct_mean_curvature_final = nonzeros(allapo_rstrct_mean_curvature_final);
allapo_rstrct_total_distance_final(isnan(allapo_rstrct_total_distance_final))=0;
allapo_rstrct_total_distance_final = nonzeros(allapo_rstrct_total_distance_final);

mean_fiber_length_aniso_smoothn = mean(allapo_rstrct_total_distance_final)
mean_pennation_angle_aniso_smoothn = mean(allapo_rstrct_mean_alpha_final)
mean_curvature_aniso_smoothn = mean(allapo_rstrct_mean_curvature_final)

fiber_visualizer_v11(muscle_mask, fv_options, roi_mesh_rstrct, anat_img, allapo_rstrct_smoothed_fiber_all_final);
set(gcf, 'position', [1200 50 500 800])
set(gca, 'fontsize', 14)
xlabel('Column', 'fontsize', 14)
ylabel('Row', 'fontsize', 14)
zlabel('Slice', 'fontsize', 14)
title('Aponeurosis Seeding-Restricted Mesh, Aniso + Smoothn')

stop_pct_mask_aniso_smoothn = numel(find(allapo_rstrct_stop_list_aniso_smoothn == 4))/numel(allapo_rstrct_stop_list_aniso_smoothn)
stop_pct_fa_aniso_smoothn = numel(find(allapo_rstrct_stop_list_aniso_smoothn == 2))/numel(allapo_rstrct_stop_list_aniso_smoothn)
stop_pct_ang_aniso_smoothn = numel(find(allapo_rstrct_stop_list_aniso_smoothn == 3))/numel(allapo_rstrct_stop_list_aniso_smoothn)

%% Fiber-tracking - APO3 - Raw

% track, smooth, and quantify fibers
[allapo_rstrct_fiber_all, allapo_rstrct_roi_flag, allapo_rstrct_stop_list_raw, allapo_rstrct_fiber_len] = ...
    fiber_track_v20(allapo_rstrct_ft_options, e1fa_raw, muscle_mask, roi_mesh_rstrct);
[allapo_rstrct_smoothed_fiber_all, ~, allapo_rstrct_smoothed_fiber_all_mm, allapo_rstrct_pcoeff_r, allapo_rstrct_pcoeff_c, allapo_rstrct_pcoeff_s, allapo_rstrct_n_points_smoothed, allapo_rstrct_residuals, allapo_rstrct_residuals_mm] = ...
    fiber_smoother_v14a(allapo_rstrct_fiber_all, fs_options);

% Use same fiber-tract seed points in all the conditions evaluated
allapo_rstrct_smoothed_fiber_all_mm_final = zeros(size(allapo_rstrct_smoothed_fiber_all_mm));   %will hold tracts after seed-location filtering
allapo_rstrct_smoothed_fiber_all_final = zeros(size(allapo_rstrct_smoothed_fiber_all));

for step_cntr = 1:size(allapo_rstrct_smoothed_fiber_all,3)
    for dim_cntr = 1:size(allapo_rstrct_smoothed_fiber_all,4)
            allapo_rstrct_smoothed_fiber_all_final(:,:,step_cntr,dim_cntr) = ...
                squeeze(allapo_rstrct_smoothed_fiber_all(:,:,step_cntr,dim_cntr)).*allapo_rstrct_seeds_final_mask;
            allapo_rstrct_smoothed_fiber_all_mm_final (:,:,step_cntr,dim_cntr) = ...
                squeeze(allapo_rstrct_smoothed_fiber_all_mm(:,:,step_cntr,dim_cntr)).*allapo_rstrct_seeds_final_mask;
    end
end

[allapo_rstrct_angle_list_final, allapo_rstrct_distance_list_final, allapo_rstrct_curvature_list_final, ~, allapo_rstrct_n_points_final, allapo_rstrct_area_final] = ...
    fiber_quantifier_v20a(allapo_rstrct_fq_options, allapo_rstrct_smoothed_fiber_all_final, roi_mesh_rstrct, muscle_mask);

% calculate summary statistics
allapo_rstrct_total_distance_final = squeeze(max(allapo_rstrct_distance_list_final, [], 3));
allapo_rstrct_total_distance_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final = sum(squeeze(allapo_rstrct_angle_list_final(:,:,:,1)), 3)./allapo_rstrct_n_points_final(:,:,2);
allapo_rstrct_mean_alpha_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_curvature_final = sum(allapo_rstrct_curvature_list_final, 3)./allapo_rstrct_n_points_final(:,:,3);
allapo_rstrct_mean_curvature_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final(isnan(allapo_rstrct_mean_alpha_final)) = 0;
allapo_rstrct_mean_alpha_final = nonzeros(allapo_rstrct_mean_alpha_final);
allapo_rstrct_mean_curvature_final(isnan(allapo_rstrct_mean_curvature_final))=0;
allapo_rstrct_mean_curvature_final = nonzeros(allapo_rstrct_mean_curvature_final);
allapo_rstrct_total_distance_final(isnan(allapo_rstrct_total_distance_final))=0;
allapo_rstrct_total_distance_final = nonzeros(allapo_rstrct_total_distance_final);

mean_fiber_length_raw = mean(allapo_rstrct_total_distance_final)
mean_pennation_angle_raw = mean(allapo_rstrct_mean_alpha_final)
mean_curvature_raw = mean(allapo_rstrct_mean_curvature_final)

fiber_visualizer_v11(muscle_mask, fv_options, roi_mesh_rstrct, anat_img, allapo_rstrct_smoothed_fiber_all_final);
set(gcf, 'position', [1200 50 500 800])
set(gca, 'fontsize', 14)
xlabel('Column', 'fontsize', 14)
ylabel('Row', 'fontsize', 14)
zlabel('Slice', 'fontsize', 14)
title('Aponeurosis Seeding-Restricted Mesh, Raw')

stop_pct_mask_raw = numel(find(allapo_rstrct_stop_list_raw == 4))/numel(allapo_rstrct_stop_list_raw)
stop_pct_fa_raw = numel(find(allapo_rstrct_stop_list_raw == 2))/numel(allapo_rstrct_stop_list_raw)
stop_pct_ang_raw = numel(find(allapo_rstrct_stop_list_raw == 3))/numel(allapo_rstrct_stop_list_raw)

%% Fiber-tracking - APO3 - Smoothn

% track, smooth, and quantify fibers
[allapo_rstrct_fiber_all, ~, allapo_rstrct_stop_list_smoothn, ~] = ...
    fiber_track_v20(allapo_rstrct_ft_options, e1fa_raw_smooth, muscle_mask, roi_mesh_rstrct);
[allapo_rstrct_smoothed_fiber_all, ~, allapo_rstrct_smoothed_fiber_all_mm, ~, ~, ~, ~, ~, ~] = ...
    fiber_smoother_v14a(allapo_rstrct_fiber_all, fs_options);

% Use same fiber-tract seed points in all the conditions evaluated
allapo_rstrct_smoothed_fiber_all_mm_final = zeros(size(allapo_rstrct_smoothed_fiber_all_mm));   %will hold tracts after seed-location filtering
allapo_rstrct_smoothed_fiber_all_final = zeros(size(allapo_rstrct_smoothed_fiber_all));

for step_cntr = 1:size(allapo_rstrct_smoothed_fiber_all,3)
    for dim_cntr = 1:size(allapo_rstrct_smoothed_fiber_all,4)
            allapo_rstrct_smoothed_fiber_all_final(:,:,step_cntr,dim_cntr) = ...
                squeeze(allapo_rstrct_smoothed_fiber_all(:,:,step_cntr,dim_cntr)).*allapo_rstrct_seeds_final_mask;
            allapo_rstrct_smoothed_fiber_all_mm_final (:,:,step_cntr,dim_cntr) = ...
                squeeze(allapo_rstrct_smoothed_fiber_all_mm(:,:,step_cntr,dim_cntr)).*allapo_rstrct_seeds_final_mask;
    end
end

[allapo_rstrct_angle_list_final, allapo_rstrct_distance_list_final, allapo_rstrct_curvature_list_final, ~, allapo_rstrct_n_points_final, allapo_rstrct_area_final] = ...
    fiber_quantifier_v20a(allapo_rstrct_fq_options, allapo_rstrct_smoothed_fiber_all_final, roi_mesh_rstrct, muscle_mask);

% calculate summary statistics
allapo_rstrct_total_distance_final = squeeze(max(allapo_rstrct_distance_list_final, [], 3));
allapo_rstrct_total_distance_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final = sum(squeeze(allapo_rstrct_angle_list_final(:,:,:,1)), 3)./allapo_rstrct_n_points_final(:,:,2);
allapo_rstrct_mean_alpha_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_curvature_final = sum(allapo_rstrct_curvature_list_final, 3)./allapo_rstrct_n_points_final(:,:,3);
allapo_rstrct_mean_curvature_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final(isnan(allapo_rstrct_mean_alpha_final)) = 0;
allapo_rstrct_mean_alpha_final = nonzeros(allapo_rstrct_mean_alpha_final);
allapo_rstrct_mean_curvature_final(isnan(allapo_rstrct_mean_curvature_final))=0;
allapo_rstrct_mean_curvature_final = nonzeros(allapo_rstrct_mean_curvature_final);
allapo_rstrct_total_distance_final(isnan(allapo_rstrct_total_distance_final))=0;
allapo_rstrct_total_distance_final = nonzeros(allapo_rstrct_total_distance_final);

mean_fiber_length_smoothn = mean(allapo_rstrct_total_distance_final)
mean_pennation_angle_smoothn = mean(allapo_rstrct_mean_alpha_final)
mean_curvature_smoothn = mean(allapo_rstrct_mean_curvature_final)

fiber_visualizer_v11(muscle_mask, fv_options, roi_mesh_rstrct, anat_img, allapo_rstrct_smoothed_fiber_all_final);
set(gcf, 'position', [1200 50 500 800])
set(gca, 'fontsize', 14)
xlabel('Column', 'fontsize', 14)
ylabel('Row', 'fontsize', 14)
zlabel('Slice', 'fontsize', 14)
title('Aponeurosis Seeding-Restricted Mesh, Smoothn')

stop_pct_mask_smoothn = numel(find(allapo_rstrct_stop_list_smoothn == 4))/numel(allapo_rstrct_stop_list_smoothn)
stop_pct_fa_smoothn = numel(find(allapo_rstrct_stop_list_smoothn == 2))/numel(allapo_rstrct_stop_list_smoothn)
stop_pct_ang_smoothn = numel(find(allapo_rstrct_stop_list_smoothn == 3))/numel(allapo_rstrct_stop_list_smoothn)

%% Fiber-tracking - APO3 - Aniso

% track, smooth, and quantify fibers
[allapo_rstrct_fiber_all, ~, allapo_rstrct_stop_list_aniso, ~] = ...
    fiber_track_v20(allapo_rstrct_ft_options, e1fa_aniso, muscle_mask, roi_mesh_rstrct);
[allapo_rstrct_smoothed_fiber_all, ~, allapo_rstrct_smoothed_fiber_all_mm, ~, ~, ~, ~, ~, ~] = ...
    fiber_smoother_v14a(allapo_rstrct_fiber_all, fs_options);

% Use same fiber-tract seed points in all the conditions evaluated
allapo_rstrct_smoothed_fiber_all_mm_final = zeros(size(allapo_rstrct_smoothed_fiber_all_mm));   %will hold tracts after seed-location filtering
allapo_rstrct_smoothed_fiber_all_final = zeros(size(allapo_rstrct_smoothed_fiber_all));

for step_cntr = 1:size(allapo_rstrct_smoothed_fiber_all,3)
    for dim_cntr = 1:size(allapo_rstrct_smoothed_fiber_all,4)
            allapo_rstrct_smoothed_fiber_all_final(:,:,step_cntr,dim_cntr) = ...
                squeeze(allapo_rstrct_smoothed_fiber_all(:,:,step_cntr,dim_cntr)).*allapo_rstrct_seeds_final_mask;
            allapo_rstrct_smoothed_fiber_all_mm_final (:,:,step_cntr,dim_cntr) = ...
                squeeze(allapo_rstrct_smoothed_fiber_all_mm(:,:,step_cntr,dim_cntr)).*allapo_rstrct_seeds_final_mask;
    end
end

[allapo_rstrct_angle_list_final, allapo_rstrct_distance_list_final, allapo_rstrct_curvature_list_final, ~, allapo_rstrct_n_points_final, allapo_rstrct_area_final] = ...
    fiber_quantifier_v20a(allapo_rstrct_fq_options, allapo_rstrct_smoothed_fiber_all_final, roi_mesh_rstrct, muscle_mask);

% calculate summary statistics
allapo_rstrct_total_distance_final = squeeze(max(allapo_rstrct_distance_list_final, [], 3));
allapo_rstrct_total_distance_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final = sum(squeeze(allapo_rstrct_angle_list_final(:,:,:,1)), 3)./allapo_rstrct_n_points_final(:,:,2);
allapo_rstrct_mean_alpha_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_curvature_final = sum(allapo_rstrct_curvature_list_final, 3)./allapo_rstrct_n_points_final(:,:,3);
allapo_rstrct_mean_curvature_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final(isnan(allapo_rstrct_mean_alpha_final)) = 0;
allapo_rstrct_mean_alpha_final = nonzeros(allapo_rstrct_mean_alpha_final);
allapo_rstrct_mean_curvature_final(isnan(allapo_rstrct_mean_curvature_final))=0;
allapo_rstrct_mean_curvature_final = nonzeros(allapo_rstrct_mean_curvature_final);
allapo_rstrct_total_distance_final(isnan(allapo_rstrct_total_distance_final))=0;
allapo_rstrct_total_distance_final = nonzeros(allapo_rstrct_total_distance_final);

mean_fiber_length_aniso = mean(allapo_rstrct_total_distance_final)
mean_pennation_angle_aniso = mean(allapo_rstrct_mean_alpha_final)
mean_curvature_aniso = mean(allapo_rstrct_mean_curvature_final)

fiber_visualizer_v11(muscle_mask, fv_options, roi_mesh_rstrct, anat_img, allapo_rstrct_smoothed_fiber_all_final);
set(gcf, 'position', [1200 50 500 800])
set(gca, 'fontsize', 14)
xlabel('Column', 'fontsize', 14)
ylabel('Row', 'fontsize', 14)
zlabel('Slice', 'fontsize', 14)
title('Aponeurosis Seeding-Restricted Mesh, Aniso')

stop_pct_mask_aniso = numel(find(allapo_rstrct_stop_list_aniso == 4))/numel(allapo_rstrct_stop_list_aniso)
stop_pct_fa_aniso = numel(find(allapo_rstrct_stop_list_aniso == 2))/numel(allapo_rstrct_stop_list_aniso)
stop_pct_ang_aniso = numel(find(allapo_rstrct_stop_list_aniso == 3))/numel(allapo_rstrct_stop_list_aniso)

%% Fiber-tracking - APO3 - TPCA

% track, smooth, and quantify fibers
[allapo_rstrct_fiber_all, ~, allapo_rstrct_stop_list_tpca, ~] = ...
    fiber_track_v20(allapo_rstrct_ft_options, e1fa_tpca, muscle_mask, roi_mesh_rstrct);
[allapo_rstrct_smoothed_fiber_all, ~, allapo_rstrct_smoothed_fiber_all_mm, ~, ~, ~, ~, ~, ~] = ...
    fiber_smoother_v14a(allapo_rstrct_fiber_all, fs_options);

% Use same fiber-tract seed points in all the conditions evaluated
allapo_rstrct_smoothed_fiber_all_mm_final = zeros(size(allapo_rstrct_smoothed_fiber_all_mm));   %will hold tracts after seed-location filtering
allapo_rstrct_smoothed_fiber_all_final = zeros(size(allapo_rstrct_smoothed_fiber_all));

for step_cntr = 1:size(allapo_rstrct_smoothed_fiber_all,3)
    for dim_cntr = 1:size(allapo_rstrct_smoothed_fiber_all,4)
            allapo_rstrct_smoothed_fiber_all_final(:,:,step_cntr,dim_cntr) = ...
                squeeze(allapo_rstrct_smoothed_fiber_all(:,:,step_cntr,dim_cntr)).*allapo_rstrct_seeds_final_mask;
            allapo_rstrct_smoothed_fiber_all_mm_final (:,:,step_cntr,dim_cntr) = ...
                squeeze(allapo_rstrct_smoothed_fiber_all_mm(:,:,step_cntr,dim_cntr)).*allapo_rstrct_seeds_final_mask;
    end
end

[allapo_rstrct_angle_list_final, allapo_rstrct_distance_list_final, allapo_rstrct_curvature_list_final, ~, allapo_rstrct_n_points_final, allapo_rstrct_area_final] = ...
    fiber_quantifier_v20a(allapo_rstrct_fq_options, allapo_rstrct_smoothed_fiber_all_final, roi_mesh_rstrct, muscle_mask);

% calculate summary statistics
allapo_rstrct_total_distance_final = squeeze(max(allapo_rstrct_distance_list_final, [], 3));
allapo_rstrct_total_distance_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final = sum(squeeze(allapo_rstrct_angle_list_final(:,:,:,1)), 3)./allapo_rstrct_n_points_final(:,:,2);
allapo_rstrct_mean_alpha_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_curvature_final = sum(allapo_rstrct_curvature_list_final, 3)./allapo_rstrct_n_points_final(:,:,3);
allapo_rstrct_mean_curvature_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final(isnan(allapo_rstrct_mean_alpha_final)) = 0;
allapo_rstrct_mean_alpha_final = nonzeros(allapo_rstrct_mean_alpha_final);
allapo_rstrct_mean_curvature_final(isnan(allapo_rstrct_mean_curvature_final))=0;
allapo_rstrct_mean_curvature_final = nonzeros(allapo_rstrct_mean_curvature_final);
allapo_rstrct_total_distance_final(isnan(allapo_rstrct_total_distance_final))=0;
allapo_rstrct_total_distance_final = nonzeros(allapo_rstrct_total_distance_final);

mean_fiber_length_tpca = mean(allapo_rstrct_total_distance_final)
mean_pennation_angle_tpca = mean(allapo_rstrct_mean_alpha_final)
mean_curvature_tpca = mean(allapo_rstrct_mean_curvature_final)

fiber_visualizer_v11(muscle_mask, fv_options, roi_mesh_rstrct, anat_img, allapo_rstrct_smoothed_fiber_all_final);
set(gcf, 'position', [1200 50 500 800])
set(gca, 'fontsize', 14)
xlabel('Column', 'fontsize', 14)
ylabel('Row', 'fontsize', 14)
zlabel('Slice', 'fontsize', 14)
title('Aponeurosis Seeding-Restricted Mesh, TPCA')

stop_pct_mask_tpca = numel(find(allapo_rstrct_stop_list_tpca == 4))/numel(allapo_rstrct_stop_list_tpca)
stop_pct_fa_tpca = numel(find(allapo_rstrct_stop_list_tpca == 2))/numel(allapo_rstrct_stop_list_tpca)
stop_pct_ang_tpca = numel(find(allapo_rstrct_stop_list_tpca == 3))/numel(allapo_rstrct_stop_list_tpca)

%% Fiber-tracking - APO3 - TPCA + smoothn

% track, smooth, and quantify fibers
[allapo_rstrct_fiber_all, ~, allapo_rstrct_stop_list_tpca_smooth, ~] = ...
    fiber_track_v20(allapo_rstrct_ft_options, e1fa_tpca_smooth, muscle_mask, roi_mesh_rstrct);
[allapo_rstrct_smoothed_fiber_all, ~, allapo_rstrct_smoothed_fiber_all_mm, ~, ~, ~, ~, ~, ~] = ...
    fiber_smoother_v14a(allapo_rstrct_fiber_all, fs_options);

% Use same fiber-tract seed points in all the conditions evaluated
allapo_rstrct_smoothed_fiber_all_mm_final = zeros(size(allapo_rstrct_smoothed_fiber_all_mm));   %will hold tracts after seed-location filtering
allapo_rstrct_smoothed_fiber_all_final = zeros(size(allapo_rstrct_smoothed_fiber_all));

for step_cntr = 1:size(allapo_rstrct_smoothed_fiber_all,3)
    for dim_cntr = 1:size(allapo_rstrct_smoothed_fiber_all,4)
            allapo_rstrct_smoothed_fiber_all_final(:,:,step_cntr,dim_cntr) = ...
                squeeze(allapo_rstrct_smoothed_fiber_all(:,:,step_cntr,dim_cntr)).*allapo_rstrct_seeds_final_mask;
            allapo_rstrct_smoothed_fiber_all_mm_final (:,:,step_cntr,dim_cntr) = ...
                squeeze(allapo_rstrct_smoothed_fiber_all_mm(:,:,step_cntr,dim_cntr)).*allapo_rstrct_seeds_final_mask;
    end
end

[allapo_rstrct_angle_list_final, allapo_rstrct_distance_list_final, allapo_rstrct_curvature_list_final, ~, allapo_rstrct_n_points_final, allapo_rstrct_area_final] = ...
    fiber_quantifier_v20a(allapo_rstrct_fq_options, allapo_rstrct_smoothed_fiber_all_final, roi_mesh_rstrct, muscle_mask);

% calculate summary statistics
allapo_rstrct_total_distance_final = squeeze(max(allapo_rstrct_distance_list_final, [], 3));
allapo_rstrct_total_distance_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final = sum(squeeze(allapo_rstrct_angle_list_final(:,:,:,1)), 3)./allapo_rstrct_n_points_final(:,:,2);
allapo_rstrct_mean_alpha_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_curvature_final = sum(allapo_rstrct_curvature_list_final, 3)./allapo_rstrct_n_points_final(:,:,3);
allapo_rstrct_mean_curvature_final(allapo_rstrct_total_distance_final<10)=0;

allapo_rstrct_mean_alpha_final(isnan(allapo_rstrct_mean_alpha_final)) = 0;
allapo_rstrct_mean_alpha_final = nonzeros(allapo_rstrct_mean_alpha_final);
allapo_rstrct_mean_curvature_final(isnan(allapo_rstrct_mean_curvature_final))=0;
allapo_rstrct_mean_curvature_final = nonzeros(allapo_rstrct_mean_curvature_final);
allapo_rstrct_total_distance_final(isnan(allapo_rstrct_total_distance_final))=0;
allapo_rstrct_total_distance_final = nonzeros(allapo_rstrct_total_distance_final);

mean_fiber_length_tpca_smooth = mean(allapo_rstrct_total_distance_final)
mean_pennation_angle_tpca_smooth = mean(allapo_rstrct_mean_alpha_final)
mean_curvature_tpca_smooth = mean(allapo_rstrct_mean_curvature_final)

fiber_visualizer_v11(muscle_mask, fv_options, roi_mesh_rstrct, anat_img, allapo_rstrct_smoothed_fiber_all_final);
set(gcf, 'position', [1200 50 500 800])
set(gca, 'fontsize', 14)
xlabel('Column', 'fontsize', 14)
ylabel('Row', 'fontsize', 14)
zlabel('Slice', 'fontsize', 14)
title('Aponeurosis Seeding-Restricted Mesh, TPCA + smoothn')

stop_pct_mask_tpca_smooth = numel(find(allapo_rstrct_stop_list_tpca_smooth == 4))/numel(allapo_rstrct_stop_list_tpca_smooth)
stop_pct_fa_tpca_smooth = numel(find(allapo_rstrct_stop_list_tpca_smooth == 2))/numel(allapo_rstrct_stop_list_tpca_smooth)
stop_pct_ang_tpca_smooth = numel(find(allapo_rstrct_stop_list_tpca_smooth == 3))/numel(allapo_rstrct_stop_list_tpca_smooth)
