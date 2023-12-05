function output_structure = form_model_muscle(dti_options, simulation_vars, isdisp)
%FUNCTION form_model_muscle
%  output_structure = form_model_muscle(dti_options, simulation_vars, isdisp)
%
%USAGE
%    The function form_model_muscle is used to simulate a stack of noise-free
%  DTI images. Muscle fiber orientations are modeled, along with adipose
%  tissue outside of it. The images are formed at 1 mm isotropic resolution using
%  architecture and diffusion-encoding parameters defined by the user. The
%  results are output as a structure. The muscle eigenvalues are assumed to
%  be L1=2x10^-5 cm^2/s; L2=1.6x10^-5 cm^2/s; and L3=1.4x10^-5 cm^2/s.
%  Adipose eigenvalues are assumed to be isotropic with L1=L2=L3=
%  0.7x10^-5 cm^2/s; slight variability is added to facilitate solution of
%  the diffusion tensor in later steps. Plots and images are generated at
%  selected major steps.
%
%INPUT ARGUMENTS
%   dti_options: Structure to describe DTI scan setup
%       -dti_dir: Diffusion encoding direction matrix
%       -b_val: b-value
%
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
%OUTPUT ARGUMENTS
%   output_structure: A structure containing fields:
%
%   diff_img: noise-free diffusion images with FOV = 42x42 mm and 120 1 mm
%   thick slices
%
%   tensor_m: spatial maps of the diffusion tensor
%
%   slice_mask: spatial maps of the distribution of tissue types
%
%   mask_muscle: spatial maps of the distribution of the muscle
%
%   model_fiber: the fiber used as a model for forming the images.  Note
%   that because of curvature, aponeurosis orientation, and/or excessive
%   length, the model fiber may not match all of the fibers in the muscles.
%
%EXAMPLE USAGE
% Set a 24-direction encoding scheme, a b-value of 450 mm^2/s, and fiber
% variability in the Y and Z directions:
%
% y_params = [0.05 0 0];
% z_params = [0.05 2 0];
% dti_dir = [0    0   0
%     0.1890    0.6940    0.6940
%     0.6940    0.6940    0.1890
%     0.1890    0.6940   -0.6940
%    -0.6940    0.6940    0.1890
%    -0.6940    0.1890   -0.6940
%     0.6940    0.1890   -0.6940
%    -0.6340    0.4420    0.6340
%     0.6340    0.4420    0.6340
%     0.6340    0.6340   -0.4420
%    -0.4420    0.6340   -0.6340
%     0.4420    0.6340    0.6340
%     0.4420    0.6340   -0.6340
%    -0.6940    0.6940   -0.1890
%     0.6340    0.6340    0.4420
%    -0.6940    0.1890    0.6940
%    -0.6340    0.4420   -0.6340
%     0.6940    0.1890    0.6940
%     0.6340    0.4420   -0.6340
%    -0.1890    0.6940   -0.6940
%    -0.4420    0.6340    0.6340
%     0.6940    0.6940   -0.1890
%    -0.1890    0.6940    0.6940
%    -0.6340    0.6340   -0.4420
%    -0.6340    0.6340    0.4420];
% b_val=450;
%
% The function is then called as:
%   output_structure = form_model_muscle(x_params, y_params, z_params, dti_dir, b_val);
%
%ACKNOWLEDGEMENTS
%  Grant support: NIH/NIAMS R01 AR073831

%% Initialization
poly_params = simulation_vars.poly_params;
voxel_dim  = simulation_vars.voxel_dim;
SNR         = simulation_vars.SNR;
dti_dir     = dti_options.dti_dir;
b_val       = dti_options.b_val;

%% form and visualize a model fiber
Nsample = 100;

dim_xy = 30;
dim_z  = 120;
FOV_core_xy_vx = ceil(dim_xy / voxel_dim(1));
FOV_core_z_vx  = ceil(dim_z / voxel_dim(3));
% roi_width   = round(full_dim_xy / voxel_dim(1));
% roi_thick   = round(full_dim_z / voxel_dim(3));

FOV_core_xy_mm = FOV_core_xy_vx * voxel_dim(1);
FOV_core_z_mm  = FOV_core_z_vx * voxel_dim(3);

h_v = linspace(0, 1.5, Nsample)';
model_fiber = zeros(FOV_core_xy_vx + 1,3);

fiber_x_v = polyval(poly_params(1, : ), h_v);
fiber_y_v = polyval(poly_params(2, : ), h_v);
fiber_z_v = polyval(poly_params(3, : ), h_v);

model_fiber( : , 2) = (0 : voxel_dim(1) : FOV_core_xy_mm + 1e-3)';
model_fiber( : , 1) = interp1(fiber_y_v, fiber_x_v, model_fiber( : , 2));
model_fiber( : , 3) = interp1(fiber_y_v, fiber_z_v, model_fiber( : , 2));

if isdisp == 1
    %view model fiber in 3D and its X, Y, and Z components
    figure('Name', 'Fiber Geometry');
    subplot(2, 3, [1, 4]);
    plot3(model_fiber( : , 1), model_fiber( : , 2), model_fiber( : , 3));
    xlabel('X Position (mm)')
    ylabel('Y Position (mm)')
    zlabel('Z Position (mm)')

    subplot(2, 3, 2);
    plot(model_fiber( : , 1), model_fiber( : , 2));
    xlabel('X Position (mm)')
    ylabel('Y Position (mm)')

    subplot(2, 3, [3, 6]);
    plot(model_fiber( : , 2), model_fiber( : , 3));
    xlabel('Y Position (mm)')
    ylabel('Z Position (mm)')

    subplot(2, 3, 5);
    plot(model_fiber( : , 3), model_fiber( : , 1));
    xlabel('Z Position (mm)')
    ylabel('X Position (mm)')
    set(gcf, 'position', [100 50 560 420])
end

%% synthesize the model muscle
% fiber_cum_dist = cumsum(sqrt(sum(diff(model_fiber) .^ 2, 2)));
% fiber_len = fiber_cum_dist(end);
% fiber_orientation_matrix = zeros(19,ceil(fiber_len),3);                  %matrix to hold muscle fiber directions

fiber_orient_m = diff(model_fiber);
fiber_orient_m = fiber_orient_m ./ sqrt(sum(fiber_orient_m .^ 2, 2));
fiber_orient_m = repmat(reshape(fiber_orient_m, [FOV_core_xy_vx, 1, 3]), ...
    [1, FOV_core_xy_vx, 1]);

% visualize
if isdisp == 1
    figure('Name', 'Fiber Orientation');
    imagesc(squeeze(fiber_orient_m( : , : , : )));
    title('(-)<--- X --->(+)')
    ylabel('(-)<--- Y --->(+)', 'fontweight', 'bold')
    axis image
    set(gcf, 'position', [100 50 560 420])
end

%% create an image containing the model muscle.
img_n_row = FOV_core_xy_vx + 4;
img_n_col = FOV_core_xy_vx + 4;
img_n_slc = FOV_core_z_vx + 4;

% form the diffusion tensor for the model muscles
tensor_muscle = zeros(FOV_core_xy_vx, FOV_core_xy_vx, FOV_core_z_vx, 3, 3);

for row_cntr = 1 : FOV_core_xy_vx
    for col_cntr = 1 : FOV_core_xy_vx
        %E1 is the fiber orientation
        E1_m1 = squeeze(fiber_orient_m(row_cntr, col_cntr, : ));

        %find a temporary vector that's orthogonal to E1 and [0 1 0];
        temp_vector_m1 = cross(E1_m1, [0, 1, 0]');
        temp_vector_m1 = temp_vector_m1 / norm(temp_vector_m1);

        %define E3 as being orthogonal to E1 and temp_v
        E3_m1 = cross(E1_m1, temp_vector_m1);

        %define E2 as being orthogonal to E1 and E3
        E2_m1 = cross(E1_m1, E3_m1);
        E_m1 = [E1_m1, E2_m1, E3_m1];

        %form the tensor
        L_m1 = [2, 0, 0; 0, 1.6, 0; 0, 0, 1.4] * 1E-5;                     %assume diffusivities
        D_m1 = E_m1 * L_m1 * E_m1';

        tensor_muscle(row_cntr, col_cntr, : , : , : ) = ...
            repmat(reshape(D_m1, [1, 1, 1, 3, 3]), [1, 1, FOV_core_z_vx]);
    end
end

slice_mask = zeros(img_n_row, img_n_col, img_n_slc);

%%%%%%%%%%%%%% CAUTION: Redifinition of mask value %%%%%%%%%%%%%%
% 1 - muscle
% 0 - surrounding adipose tissue
[~, Ym, Zm] = meshgrid(1 : FOV_core_xy_vx, ...
                       1 : FOV_core_xy_vx, ...
                       1 : voxel_dim(3) : FOV_core_z_mm);
model_fiber_z = model_fiber( : , 3);
slice_mask_core = (Zm > model_fiber_z(Ym));
slice_mask(3 : img_n_row - 2, ...
           3 : img_n_col - 2, ...
           3 : img_n_slc - 2) = slice_mask_core;

%around the muscle is inter-muscular adipose tissue - encode an isotropic tensor
tensor_whole_slice = (repmat(reshape(eye(3), [1, 1, 1, 3, 3]), ...
    [img_n_row, img_n_col, img_n_slc]) ...
    .* randn(img_n_row, img_n_col, img_n_slc, 3, 3) ...
    .* 0.01 ...
    + reshape(eye(3) .* 0.7, [1, 1, 1, 3, 3])) .* 1e-5;

tensor_core_vol = tensor_whole_slice(3 : img_n_row - 2, ...
                                     3 : img_n_col - 2, ...
                                     3 : img_n_slc - 2, ...
                                     : , : ) .* (~slice_mask_core) ...
                  + tensor_muscle .* slice_mask_core;

tensor_whole_slice(3 : img_n_row - 2, ...
                   3 : img_n_col - 2, ...
                   3 : img_n_slc - 2, ...
                   : , : ) = tensor_core_vol;

if isdisp == 1
    %visualize masks of tissue type
    figure('Name', 'Tissue Masks');
    n = 1;
    for k = round(linspace(2, FOV_core_z_vx, 18))
        subplot(3, 6, n);
        imagesc(slice_mask( : , : , k))
        axis image;
        axis off;
        caxis([0, 1]);
        title(['Slice ' num2str(k)])
        n=n+1;
    end
    set(gcf, 'position', [100 50 560 420])

    % plot tensors for two representative slices
    plot_titles = 'XYZ';
    figure('Name', 'Diffusion Tensors');
    slc_demo_1 = ceil(5 / voxel_dim(3)) + 1;
    slc_demo_2 = ceil(70 / voxel_dim(3)) + 1;
    for p = 1 : 3
        subplot(2, 3, p);
        imagesc(1e5 * squeeze(tensor_whole_slice( : , : , slc_demo_1, p, p)));
        caxis([0, 2.01]);
        cb = colorbar;
        set(cb, 'ticks', 0 : 0.4 : 2);
        title(['D_' plot_titles(p) '_' plot_titles(p) ...
            ' (x10^{-5} cm^2/s), Slice ' int2str(slc_demo_1)]);
        axis image;
        axis off;

        subplot(2, 3, p + 3);
        imagesc(1e5 * squeeze(tensor_whole_slice( : , : , slc_demo_2, p, p)));
        caxis([0, 2.01]);
        cb = colorbar;
        set(cb, 'ticks', 0 : 0.4 : 2);
        title(['D_' plot_titles(p) '_' plot_titles(p) ...
            ' (x10^{-5} cm^2/s), Slice ' int2str(slc_demo_2)]);
        axis image;
        axis off;
    end
    set(gcf, 'position', [100 50 560 420])
end

%% simulate images
Ndir = size(dti_dir, 1);
diff_img = zeros(img_n_row, img_n_col, img_n_slc, Ndir);
tensor_m = zeros(size(tensor_whole_slice));

for slc_cntr = 1 : img_n_slc
    for row_cntr = 1 : img_n_row
        for col_cntr = 1 : img_n_col
            if slice_mask(row_cntr, col_cntr, slc_cntr) == 0
                %scale fat signals by 0.05 to account for fat suppression
                scale_coef = 0.05;
            elseif slice_mask(row_cntr, col_cntr, slc_cntr) == 1
                scale_coef = 1.0;
            end

            %convert from cm^2/s to mm^2/s
            d_tensor = squeeze(tensor_whole_slice(row_cntr, col_cntr, slc_cntr, : , : )) .* 100;

            for diff_dir = 1 : Ndir
                diff_img(row_cntr, col_cntr, slc_cntr, diff_dir) = ...
                    scale_coef .* ...
                    exp(-b_val * squeeze(dti_dir(diff_dir, : )) * d_tensor * squeeze(dti_dir(diff_dir, : ))');
            end
        end
    end
end

% Add noise
if SNR ~= 999
    diff_img = abs(diff_img + randn(size(diff_img)) ./ SNR);
end

% Compute tensor
for slc_cntr = 1 : img_n_slc
    for row_cntr = 1 : img_n_row
        for col_cntr = 1 : img_n_col
            signal_v = squeeze(diff_img(row_cntr, col_cntr, slc_cntr, : ));
            tensor_m(row_cntr, col_cntr, slc_cntr, : , : ) = ...
                signal2tensor2(signal_v, dti_dir(2 : end, : ), b_val);
        end %of column loop
    end %or row loop
end %of slice loop

if isdisp == 1
    eigvalue1_max = 0;
    eigvector1_m = zeros(img_n_row, img_n_col, img_n_slc, 3);
    for slc_cntr = 1 : img_n_slc
        for row_cntr = 1 : img_n_row
            for col_cntr = 1 : img_n_col
                [eigvector_m, eigvalue_m] = ...
                    eigs(squeeze(tensor_m(row_cntr,col_cntr,slc_cntr,:,:)));
                
                if max(max(eigvalue_m))>eigvalue1_max
                    eigvalue1_max = max(max(eigvalue_m));
                end
                
                eigvalue_m = diag(eigvalue_m);
                eigvector1_m(row_cntr, col_cntr, slc_cntr, : ) = ...
                    eigvector_m( : , find(eigvalue_m == max(max(eigvalue_m))));
                if eigvector1_m(row_cntr, col_cntr, slc_cntr, 3) < 0
                    eigvector1_m(row_cntr, col_cntr, slc_cntr, : ) = ...
                        -eigvector1_m(row_cntr, col_cntr, slc_cntr, : );
                end
            end
        end
    end
    % plot diffusion images
    figure('Name', 'Diffusion Images');
    im2show1=cat(2, ...
        diff_img( : , : , slc_demo_1, 1), ...
        diff_img( : , : , slc_demo_1, 2), ...
        diff_img( : , : , slc_demo_1, 3), ...
        diff_img( : , : , slc_demo_1, 4), ...
        diff_img( : , : , slc_demo_1, 5));
    im2show2=cat(2, ...
        diff_img( : , : , slc_demo_1, 6), ...
        diff_img( : , : , slc_demo_1, 7), ...
        diff_img( : , : , slc_demo_1, 8), ...
        diff_img( : , : , slc_demo_1, 9), ...
        diff_img( : , : , slc_demo_1, 10));
    im2show3=cat(2, ...
        diff_img( : , : , slc_demo_1, 11), ...
        diff_img( : , : , slc_demo_1, 12), ...
        diff_img( : , : , slc_demo_1, 13), ...
        diff_img( : , : , slc_demo_1, 14), ...
        diff_img( : , : , slc_demo_1, 15));
    im2show4=cat(2, ...
        diff_img( : , : , slc_demo_1, 16), ...
        diff_img( : , : , slc_demo_1, 17), ...
        diff_img( : , : , slc_demo_1, 18), ...
        diff_img( : , : , slc_demo_1, 19), ...
        diff_img( : , : , slc_demo_1, 20));
    im2show5=cat(2, ...
        diff_img( : , : , slc_demo_1, 21), ...
        diff_img( : , : , slc_demo_1, 22), ...
        diff_img( : , : , slc_demo_1, 23), ...
        diff_img( : , : , slc_demo_1, 24), ...
        diff_img( : , : , slc_demo_1, 25));
    im2show=cat(1, im2show1, im2show2, im2show3, im2show4, im2show5);
    imagesc(im2show)
    colormap gray
    axis image, axis off
    title(['Slice ' int2str(slc_demo_1)]);
    set(gcf, 'position', [100 50 560 420])

    figure('Name', 'Diffusion Images');
    im2show1=cat(2, ...
        diff_img( : , : , slc_demo_2, 1), ...
        diff_img( : , : , slc_demo_2, 2), ...
        diff_img( : , : , slc_demo_2, 3), ...
        diff_img( : , : , slc_demo_2, 4), ...
        diff_img( : , : , slc_demo_2, 5));
    im2show2=cat(2, ...
        diff_img( : , : , slc_demo_2, 6), ...
        diff_img( : , : , slc_demo_2, 7), ...
        diff_img( : , : , slc_demo_2, 8), ...
        diff_img( : , : , slc_demo_2, 9), ...
        diff_img( : , : , slc_demo_2, 10));
    im2show3=cat(2, ...
        diff_img( : , : , slc_demo_2, 11), ...
        diff_img( : , : , slc_demo_2, 12), ...
        diff_img( : , : , slc_demo_2, 13), ...
        diff_img( : , : , slc_demo_2, 14), ...
        diff_img( : , : , slc_demo_2, 15));
    im2show4=cat(2, ...
        diff_img( : , : , slc_demo_2, 16), ...
        diff_img( : , : , slc_demo_2, 17), ...
        diff_img( : , : , slc_demo_2, 18), ...
        diff_img( : , : , slc_demo_2, 19), ...
        diff_img( : , : , slc_demo_2, 20));
    im2show5=cat(2, ...
        diff_img( : , : , slc_demo_2, 21), ...
        diff_img( : , : , slc_demo_2, 22), ...
        diff_img( : , : , slc_demo_2, 23), ...
        diff_img( : , : , slc_demo_2, 24), ...
        diff_img( : , : , slc_demo_2, 25));
    im2show=cat(1, im2show1, im2show2, im2show3, im2show4, im2show5);
    imagesc(im2show)
    colormap gray
    axis image, axis off
    title(['Slice ' int2str(slc_demo_2)]);
    set(gcf, 'position', [100 50 560 420])

    %plot two representative eigenvector 1 maps
    figure('Name', 'E1 Maps');
    subplot(1,2,1)
    imagesc(squeeze(eigvector1_m(:,:,slc_demo_1,:)))
    title(['\epsilon_1, Slice ' int2str(slc_demo_1)]);
    axis image, axis off
    set(gcf, 'position', [100 50 560 420])

    subplot(1,2,2)
    imagesc(squeeze(eigvector1_m(:,:,slc_demo_2,:)))
    title(['\epsilon_1, Slice ' int2str(slc_demo_2)]);
    axis image, axis off
    set(gcf, 'position', [100 50 560 420])
end
%% prepare output structure
% mask_muscle = slice_mask;
% slice_mask(slice_mask == 0) = 2;
% mask_muscle(slice_mask == 2) = 0;

output_structure.diff_img    = single(diff_img);
output_structure.tensor_m    = single(tensor_m);
% output_structure.tissue_mask = slice_mask;
% output_structure.mask_muscle = mask_muscle;
output_structure.mask_muscle = slice_mask;
output_structure.x_params    = single(poly_params(1, : ));
output_structure.y_params    = single(poly_params(2, : ));
output_structure.z_params    = single(poly_params(3, : ));
output_structure.model_fiber = single(model_fiber);

%% end function
return
