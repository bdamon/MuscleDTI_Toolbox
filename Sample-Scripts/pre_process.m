%% pre_process.m
% This file was used to process the data published in the manuscript.
% Completed on 17 Jan 2021, Bruce Damon


%% File input

% start with a clean slate:
clear
close all
clc

% select images
[dixon_file, dixon_directory] = uigetfile('*.DCM', 'Select Dixon Image');   %Dixon image
[dti_distal_file, dti_distal_directory] = uigetfile('*.DCM', 'Select Distal Stack of DTMRI Image'); %Distal DTI data
[dti_proximal_file, dti_proximal_directory] = uigetfile('*.DCM', 'Select Proximal Stack of DTMRI Image'); %Proximal DTI data

% open Dixon images
cd(dixon_directory);
dixon_img = dicomread(dixon_file);                                          %open dixon image data
dixon_img = double(dixon_img);                                              %convert to double precision
dixon_info = dicominfo(dixon_file);                                         %open dixon header

% open Dixon file header/get image geometry information
dixon_numslcs = double(dixon_info.Private_2001_1018);                       %number of slices
dixon_numrows = double(dixon_info.Rows);                                    %dixon matrix size, row dimension
dixon_numcols = double(dixon_info.Columns);                                 %dixon matrix size, column dimension
dixon_pixelspace = (double(dixon_info.PerFrameFunctionalGroupsSequence.Item_1.PixelMeasuresSequence.Item_1.PixelSpacing))'; %in-plane resolution
dixon_fov = dixon_pixelspace*dixon_numrows;                                 %in-plane resolution x number of rows
dixon_slcthick =  double(dixon_info.PerFrameFunctionalGroupsSequence.Item_1.Private_2005_140f.Item_1.SpacingBetweenSlices); %original Dixon slice thickness

% open distal DTI images
cd(dti_distal_directory);                                                   %start with distal stack
dti_distal_data = dicomread(dti_distal_file);                               %rows x columns x image number
dti_distal_info = dicominfo(dti_distal_file);                               %get image information from DTI header/

% open distal DTI file header/get image geometry, diffusion encoding information
dti_numrows = double(dti_distal_info.Rows);                                 %DTI matrix size, row dimension
dti_numcols = double(dti_distal_info.Columns);                              %DTI matrix size, column dimension
dti_slcthick = double(dti_distal_info.PerFrameFunctionalGroupsSequence.Item_1.Private_2005_140f.Item_1.SpacingBetweenSlices);
dti_numslcs = double(dti_distal_info.Private_2001_1018);                    %DTI matrix size, slice direction
dti_numdir = length(dti_distal_data(1,1,:))/dti_numslcs - 1;                %calculate number of directions
dti_pixelspace = (double(dti_distal_info.PerFrameFunctionalGroupsSequence.Item_1.PixelMeasuresSequence.Item_1.PixelSpacing))'; %in-plane resolution
dti_fov = dti_pixelspace*dti_numrows;                                       %in-plane resolution x number of rows
bval = dti_distal_info.PerFrameFunctionalGroupsSequence.Item_2.MRDiffusionSequence.Item_1.DiffusionBValue; %get the single b-value
for d=1:dti_numdir                                                          %loop to get directions
    eval(['bvect(d,:) = dti_distal_info.PerFrameFunctionalGroupsSequence.Item_' num2str(d+1) '.MRDiffusionSequence.Item_1.DiffusionGradientDirectionSequence.Item_1.DiffusionGradientOrientation;']);
end

% reshape DTI data
dti_distal_data = reshape(dti_distal_data, dti_numrows, dti_numcols, (dti_numdir+1), dti_numslcs); %data are now rows x columns x direction x slice number
dti_distal_data = permute(dti_distal_data, [1 2 4 3]);                      %data are now rows x columns x slice number x direction
dti_distal_data = double(dti_distal_data);                                  %convert to double precision

% open proximal images
cd(dti_proximal_directory);                                                 %repeat for proximal stack
dti_proximal_data = dicomread(dti_proximal_file);
dti_proximal_data = reshape(dti_proximal_data, dti_numrows, dti_numcols, (dti_numdir+1), dti_numslcs);
dti_proximal_data = permute(dti_proximal_data, [1 2 4 3]);
dti_proximal_data = double(dti_proximal_data);


%% Adjust for image geometry differences
% warning - some hard-coding in this section, specific to the acquisition scheme

% extract water image; use as anatomical image
anat_image_orig = squeeze(dixon_img(:,:,((2*dixon_numslcs + 1):(3*dixon_numslcs))));
anat_image = imresize3(anat_image_orig, [dixon_numrows dixon_numcols dixon_numslcs/4]); %DTI and water image sizes now match (DTI slice thickness=4*water image slice thickness)
anat_numslcs = dixon_numslcs/4;                                             %will resize slice dimension down by a factor of four
anat_numrows = dixon_numrows;                                               %anatomical image matrix size, row dimension
anat_numcols = dixon_numcols;                                               %anatomical image, column dimension
anat_pixelspace = dixon_pixelspace;                                         %anatomical image, in-plane resolution
anat_fov = dixon_fov;                                                       %anatomical image, in-plane resolution x number of rows
anat_slcthick = dixon_slcthick*4;                                           %will resize slice dimension down by a factor of four

% concatenate the DTI data files
dti_all_unreg(:,:,1:24,:) = dti_distal_data;                                %put distal slices into the dti_all matrix
dti_all_unreg(:,:,21:44,:)= dti_proximal_data;                              %put proximal slices into the dti_all matrix; account for overlap
dti_all_numslcs = 44;


%% Form whole-leg mask for registration and diffusion tensor calculations

% need to remove skin pixels from water image before registration
fat_image = squeeze(dixon_img(:,:,(1:dixon_numslcs)));                      %get the fat image
pd_image = fat_image + anat_image_orig;                                     %make a fat+water (pseudo proton-density) image
pd_image = imresize3(pd_image, [anat_numrows anat_numcols anat_numslcs]);   %resize slice dimension to match DTI images

% pre-allocate memory for mask
pd_mask = zeros(size(pd_image));

% mask formation loop
for dw_slc=1:length(pd_image(1,1,:))
    
    % get the slice/normalize
    loop_image = pd_image(:,:,dw_slc);
    loop_image = loop_image/max(max(loop_image));                           %normalize signals
    
    %form mask/do some processing
    loop_anat_mask = zeros(anat_numrows, anat_numcols);
    loop_threshold = graythresh(loop_image);                                %signal threshold using Otsu's method
    loop_anat_mask(loop_image>loop_threshold) = 1;                          %form the mask
    loop_anat_mask = imfill(loop_anat_mask, 'holes');                       %fill holes
    loop_anat_mask=bwmorph(loop_anat_mask, 'erode', 2);                     %erode twice to eliminate skin pixels
    
    pd_mask(:,:,dw_slc) = loop_anat_mask;                                   %put loop mask into the big matrix
    
end


%% Register images

% pre-allocate memory
dti_all_reg = zeros(size(dti_all_unreg));

% indices for defining the four quadrants
idx_r = [1 anat_numrows/2; 1 anat_numrows/2; (1+anat_numrows/2) anat_numrows; (1+anat_numrows/2) anat_numrows];
idx_c = [1 anat_numcols/2; (1+anat_numcols/2) anat_numcols; 1 anat_numcols/2; (1+anat_numcols/2) anat_numcols];

for dw_slc=1:dti_all_numslcs
    
    loop_anat_mask = squeeze(pd_mask(:,:,dw_slc));
    
    loop_anat_img = anat_image(:,:,dw_slc);
    loop_anat_img = loop_anat_img.*loop_anat_mask;
    loop_anat_img_norm = loop_anat_img/max(max(loop_anat_img));
    
    for d = 1:(dti_numdir+1)
        
        %get dti image for the loop
        loop_dti_unreg = squeeze(dti_all_unreg(:,:,dw_slc,d));
        loop_dti_unreg_norm = loop_dti_unreg/max(max(loop_dti_unreg));
        
        % form mask for each DW image, based on regional gray level threshold:
        loop_dwi_mask = zeros(anat_numrows,anat_numcols,4);
        for k=1:4
            
            % divide unregisteredimage into four quadrants
            quarter_img = zeros(192);
            quarter_img(idx_r(k,1):idx_r(k,2), idx_c(k,1):idx_c(k,2)) = loop_dti_unreg(idx_r(k,1):idx_r(k,2), idx_c(k,1):idx_c(k,2));
            quarter_img = quarter_img/max(max(quarter_img));
            
            %create mask for teh quadrant
            quarter_mask = zeros(anat_numcols);
            threshold = graythresh(quarter_img);                            %use Otsu's method
            quarter_mask(quarter_img>threshold) = 1;
            loop_dwi_mask(:,:,k) = quarter_mask;
            
        end
        loop_dwi_mask = sum(loop_dwi_mask,3);
        
        %some morphological processing on the mask
        loop_dwi_mask = imfill(loop_dwi_mask, 'holes');
        loop_dwi_mask = bwmorph(loop_dwi_mask, 'open');
        loop_dwi_mask = bwmorph(loop_dwi_mask, 'dilate');
        
        % extra processing of normalized b=0 image to match signal patterns to water image
        if d==1
            loop_dti_unreg_norm(loop_dti_unreg_norm>0.5) = 0.5;             %in b=0 image, find blood vessels and then flatten out their signal 
            loop_dti_unreg_norm = 2*loop_dti_unreg_norm;                    %match signal to anatomical image
        end
        loop_dti_unreg_norm = loop_dti_unreg_norm.*loop_dwi_mask;           %mask the normalized image

        %registration
        [disp_field, temp] = imregdemons(loop_dti_unreg_norm, loop_anat_img_norm); %calculate displacement field using normalized image
        loop_dti_reg = imwarp(loop_dti_unreg, disp_field);                  %apply displacement field to non-normalized image
        dti_all_reg(:,:,dw_slc,d) = loop_dti_reg;
        
    end                                                                     %of diffusion-weighted image loop
end                                                                         %of slice loop

% inspect results - comment out for unsupervised execution
% close all
% for dw_slc=1:4:44
%     for d=1:25
%         subplot(5,5,d)
%         imshowpair(anat_image(:,:,dw_slc), dti_all_reg(:,:,dw_slc,d))
%         axis([50 150 1 100])
%     end
%     pause
% end

%% De-noise images

%set parameters based on Buck et al
noise = 5;
sigma = noise/100;
rho = 2*sigma;
delta_t = noise*3/44;
schemetype = 'implicit_multisplitting';
isfasteig  = true;
isnormg = false;
dti_res = [dti_pixelspace dti_slcthick];

dti_all_reg_smooth = aniso4D_smoothing(dti_all_reg, sigma, rho, delta_t, dti_res, schemetype, isnormg, isfasteig);

%% Estimate the diffusion tensor

tensor_m = zeros(dti_numrows, dti_numcols, dti_all_numslcs, 3, 3);
for dw_row = 1:dti_numrows
    for dw_col = 1:dti_numcols
        for dw_slc = 1:dti_all_numslcs
            if pd_mask(dw_row,dw_col,dw_slc)==1
                
                signal_v = squeeze(dti_all_reg_smooth(dw_row,dw_col,dw_slc,:));
                tensor_m(dw_row,dw_col,dw_slc,:,:) = signal2tensor2(signal_v, bvect, bval);
                
            end
        end
    end
end

% inspect results - comment out for unsupervised execution
% close all
% for dw_slc=1:4:44
%     for k=1:3
%         subplot(1,3,k)
%         imagesc(squeeze(tensor_m(:,:,dw_slc,k,k)))
%         axis image
%     end
%     pause
% end

% save all data
save all_preprocessed_data

% save data needed for fiber tracking
save data_for_fibertracking tensor_m anat_image anat_numslcs anat_numrows anat_numcols anat_pixelspace anat_fov anat_slcthick ...
    dti_numrows dti_numcols dti_slcthick dti_all_numslcs dti_numdir dti_pixelspace dti_fov



