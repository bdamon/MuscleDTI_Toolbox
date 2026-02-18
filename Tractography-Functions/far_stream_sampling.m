function [sampled_fiber_all_mm,sampled_fiber_all,sampled_fiber_all_idx] = ...
    far_stream_sampling(fiber_all_mm,fiber_all,sample_size,resampled_points)
%
%FUNCTION far_stream_sampling
%  sampled_fiber_all_mm = ...
%     far_stream_sampling(fiber_all_mm,sample_size,resampled_points)
%
%USAGE
%  The function far_stream_sampling is used to uniformly sample
%  fiber-tracts using the farthest streamline sampling method proposed by
%  Li et al., IEEE Trans Biomed Eng, 2024. Briefly, this method resamples
%  each fiber tract to N points and then finds the set of N tracts that
%  result in the largest distance between all pairs of tracts.
%
%  The user supplies the fiber_all_mm and fiber_all matrices, the desired
%  final number of tracts K, and the desired number of points N to use when
%  resampling the tracts. The resampled fiber_all_mm and fiber_all matrices
%  and the indices into the fiber_all_mm and fiber_all matrices are returned.
%
%
%INPUT ARGUMENTS
%  fiber_all_mm: The fiber tracts. The rows and columns correspond to locations
%    on the roi_mesh or the different seeding locations if dimension 2 is equal 
%    to 1. Dimension 3 gives point numbers on the tract, and the fourth dimension
%    has row, column, and slice coordinates, with units of mm.
%
%  fiber_all: The fiber tracts. The rows and columns correspond to locations
%    on the roi_mesh or the different seeding locations if dimension 2 is equal 
%    to 1. Dimension 3 gives point numbers on the tract, and the fourth dimension
%    has row, column, and slice coordinates, with units of voxels.
%
%  sample_size: The number of fiber-tracts contained in the output fiber
%    tract matrix
%
%  resampled_points: Number of points to which the fiber-tracts will be
%    resampled prior to the computation of the resampling distance; the default
%    is 12
%
%OUTPUT ARGUMENTS
%  sampled_fiber_all_mm: The matrix including the sampled fiber tracts
%
%  sampled_fiber_all: The matrix including the sampled fiber tracts, with
%  units of voxels.
%
%  sampled_fiber_all_idx: The indices of the sampled fiber-tracts, 1D array for
%    voxel-seeded tracts, 2D array for roi mesh-seeded tracts
%
%OTHER FUNCTIONS IN THE MUSCLE DTI FIBER-TRACKING TOOLBOX
%  For help with anisotropic smoothing, see <a href="matlab: help aniso4D_smoothing">aniso4D_smoothing</a>.
%  For help with threshold PCA denoising, see <a href="matlab: help TPCA_denoising">TPCA_denoising</a>.
%  For help with eddy current correction, see <a href="matlab: help eddy_correct">eddy_correct</a>.
%  For help calculating the diffusion tensor, see <a href="matlab: help signal2tensor2">signal2tensor2</a>.
%  For help defining the muscle mask, see <a href="matlab: help define_muscle">define_muscle</a>.
%  For help defining the aponeurosis ROI, see <a href="matlab: help define_roi_v11">define_roi_v11</a>.
%  For help with fiber tracking, see <a href="matlab: help fiber_track_v20">fiber_track_v20</a>.
%  For help smoothing fiber tracts, see <a href="matlab: help fiber_smoother_v14">fiber_smoother_v12</a>.
%  For help uniformly sampling the fiber tracts prior their quantification, see <a href="matlab: help far_stream_sampling">far_stream_sampling</a>.
%  For help quantifying fiber tracts, see <a href="matlab: help fiber_quantifier_v20">fiber_quantifier_v20</a>.
%  For help selecting fiber tracts following their quantification, see <a href="matlab: help fiber_goodness">fiber_goodness</a>.
%  For help visualizing fiber tracts and other structures, see <a href="matlab: help fiber_visualizer_v12">fiber_visualizer_v12</a>.
%
%VERSION INFORMATION
%  v. 1.0.0 (initial release), June 13, 2025, Roberto Pineda Guzman
%  v. 1.1.0 October 19, 2025, added fiber_all input and sampled_fiber_all
%  output, Roberto Pineda Guzman
%
%ACKNOWLEDGEMENTS
%  Grant support: NIH/NIAMS R01 AR073831

%% 
if nargin<4
    resampled_points = 12;
end

% Resize to a 4D array with dimensions #fiber-tracts x 1 x fiber-tract
% points x 3
n_rows_fiber_all_mm = size(fiber_all_mm,1);
n_cols_fiber_all_mm = size(fiber_all_mm,2);
if n_cols_fiber_all_mm > 1
    fiber_all_mm_orig = fiber_all_mm; %Preserve original fiber_all array
    fiber_all_orig = fiber_all; %Preserve original fiber_all array
    fiber_all_mm = reshape(fiber_all_mm,[size(fiber_all_mm,1)*size(fiber_all_mm,2) 1 ...
        size(fiber_all_mm,3) size(fiber_all_mm,4)]);
    fiber_all = reshape(fiber_all,[size(fiber_all,1)*size(fiber_all,2) 1 ...
        size(fiber_all,3) size(fiber_all,4)]);
end

% Resample all fiber-tracts to the specified number of resampled_points
resampled_fiber_all_mm = zeros([size(fiber_all_mm,1) size(fiber_all_mm,2) resampled_points 3]);

for n_row = 1:size(resampled_fiber_all_mm ,1)
    for n_col = 1:size(resampled_fiber_all_mm,2)
        fiber_loop = squeeze(fiber_all_mm(n_row,n_col,:,:));
        if sum(fiber_loop,'all') > 0
            zero_idx = find(fiber_loop(:,1)>0,1,'last');
            fiber_loop_r = fiber_loop(1:zero_idx,1);
            fiber_loop_c = fiber_loop(1:zero_idx,2);
            fiber_loop_s = fiber_loop(1:zero_idx,3);
            resampled_fiber_loop_r = interp1(1:length(fiber_loop_r),fiber_loop_r,...
                linspace(1,length(fiber_loop_r),resampled_points));
            resampled_fiber_loop_c = interp1(1:length(fiber_loop_c),fiber_loop_c,...
                linspace(1,length(fiber_loop_c),resampled_points));
            resampled_fiber_loop_s = interp1(1:length(fiber_loop_s),fiber_loop_s,...
                linspace(1,length(fiber_loop_s),resampled_points));
        else
            resampled_fiber_loop_r = zeros(1,resampled_points);
            resampled_fiber_loop_c = zeros(1,resampled_points);
            resampled_fiber_loop_s = zeros(1,resampled_points);
        end

        resampled_fiber_all_mm(n_row,n_col,:,1) = resampled_fiber_loop_r;
        resampled_fiber_all_mm(n_row,n_col,:,2) = resampled_fiber_loop_c;
        resampled_fiber_all_mm(n_row,n_col,:,3) = resampled_fiber_loop_s;
    end
end

s = zeros([sample_size,1]); % store indices of preserved fiber-tracts

% Select random streamline
s(1) = randi(size(resampled_fiber_all_mm,1));

%Dobuble-check random streamline is not empty
while sum(squeeze(resampled_fiber_all_mm(s(1),1,:,:)),'all')==0
    s(1) = randi(size(resampled_fiber_all_mm,1));
end

dist = zeros([size(resampled_fiber_all_mm,1),1]); % store MDF between fiber-tracts

%Compute distance between each streamline and first sampled streamline
for i = 1:size(resampled_fiber_all_mm,1)
    dist(i) = min_direct_flip(squeeze(resampled_fiber_all_mm(i,1,:,:)), ...
        squeeze(resampled_fiber_all_mm(s(1),1,:,:)), resampled_points);
end

for t = 2:sample_size
    [~,k_idx] = max(dist);
    s(t) = k_idx;

    for i = 1:size(resampled_fiber_all_mm,1)
        if dist(i) > min_direct_flip(squeeze(resampled_fiber_all_mm(i,1,:,:)), ...
        squeeze(resampled_fiber_all_mm(s(t),1,:,:)), resampled_points)
            dist(i) = min_direct_flip(squeeze(resampled_fiber_all_mm(i,1,:,:)), ...
                squeeze(resampled_fiber_all_mm(s(t),1,:,:)), resampled_points);
        end
    end
end

sampled_fiber_all_mm = zeros([sample_size size(fiber_all_mm,2) size(fiber_all_mm,3) 3]);
sampled_fiber_all = zeros([sample_size size(fiber_all,2) size(fiber_all,3) 3]);

for t = 1:sample_size
    sampled_fiber_all_mm(t,1,:,:) = squeeze(fiber_all_mm(s(t),1,:,:));
    sampled_fiber_all(t,1,:,:) = squeeze(fiber_all(s(t),1,:,:));
end

sampled_fiber_all_idx = s; %indices of sampled fiber-tracts

% Convert sampled_fiber_all_mm to 4D array with the same size of the input
% fiber-tracts and convert sampled_fiber_all_idx  to 2D roi_mesh indices
if n_cols_fiber_all_mm > 1
    s_roi_mesh = zeros([sample_size,2]);
    for t = 1:sample_size
        [row_roi_mesh, col_roi_mesh] = ind2sub([n_rows_fiber_all_mm n_cols_fiber_all_mm],s(t));
        s_roi_mesh(t,1) = row_roi_mesh;
        s_roi_mesh(t,2) = col_roi_mesh;
    end
    sampled_fiber_all_idx = s_roi_mesh;

    sampled_fiber_all_mm = zeros(size(fiber_all_mm_orig));
    sampled_fiber_all = zeros(size(fiber_all_orig));
    for k = 1:sample_size
        sampled_fiber_all_mm(sampled_fiber_all_idx(k,1),sampled_fiber_all_idx(k,2),:,:) = fiber_all_mm_orig(sampled_fiber_all_idx(k,1),sampled_fiber_all_idx(k,2),:,:);
        sampled_fiber_all(sampled_fiber_all_idx(k,1),sampled_fiber_all_idx(k,2),:,:) = fiber_all_orig(sampled_fiber_all_idx(k,1),sampled_fiber_all_idx(k,2),:,:);
    end

end

end

% Function that computes minimum average direct flip distance
function distance = min_direct_flip(u, v ,n_points) 
    d_direct = 1/n_points*sum(vecnorm((u-v)'));
    v_flipped = flip(v,1);
    d_flipped = 1/n_points*sum(vecnorm((u-v_flipped)'));
    distance = min(d_direct,d_flipped);
end
