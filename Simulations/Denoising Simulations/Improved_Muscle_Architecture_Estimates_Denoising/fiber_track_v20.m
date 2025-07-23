function [fiber_all, roi_flag, stop_list, fiber_len, seed_points] = ...
    fiber_track_v20(ft_options, e1fa, mask, roi_mesh, fv_options, anat_image)
%FUNCTION fiber_track
%  [fiber_all, roi_flag, stop_list, fiber_len, seed_points] = ...
%    fiber_track_v20(ft_options, e1fa, mask, roi_mesh, fv_options, anat_image);
%
%USAGE
%  The function fiber_track_v2, beta is used to fiber-track a muscle DTI dataset.
%
%  The required inputs include a structure defining the fiber-tracking
%  options; a 4D matrix with [row column slice] dimensions matching those
%  of the DTMRI data and holding tracking-relevant diffusion data at each
%  voxel; and the muscle mask, output from define_muscle or other program.
%  The structure allows the user to set options such as the tracking algorithm,
%  step size, laboratory frame of reference,  image orientation, and tract
%  termination method. There are also several optional arguments, including
%  the aponeurosis mesh (output from define_roi; necessary for aponeurosis-
%  based  tracking methods) and variables related to visualization.
%
%  Fibers are tracked from the seed points according to the selected
%  propagation algorithm. The seed points may be defined according to a
%  mesh reconstruction of the aponeurosis of muscle fiber insertion, or
%  they may be defined in one or more image planes or in each voxel. Each
%  fiber tract is propagated until it reaches the edge of the mask or meets
%  another stop criterion, such as an excessive inter-segment angle or an
%  out-of-bounds value for fractional anisotropy (FA).  See the description
%  of the input arguments for additional information on these criteria.
%
%  The outputs include the fiber tracts, variables describing the outcomes
%  of the tracking, and selected data about the tracts.
%
%  The fiber tracts may be viewed using fiber_visualizer_v11, either as part of
%  the function call to fiber_track or directly from the command line.
%
%INPUT ARGUMENTS
%  The following input arguments are required:
%
%  ft_options: A structure containing the following fields:
%    ref_frame: The frame of reference in which the diffusion directions
%      are specified. For example, set ft_options.ref_frame='LPS'; if the
%      left, posterior, and superior anatomical positions are (+).
%
%    image_orient: The orientation of the images. Specify the anatomical
%      directions at the top and right edges of the image as A (anterior) or
%      P (posterior) and right (R) or left (L).  Input the result as a 2-element
%      string variable (ft_options.image_orient='RA', 'AL', etc.).
%
%   seed_method: The method used to seed the points. The method must be
%     specified as a string. The options are:
%       -apo: A mesh reconstruction of the aponeurosis of muscle fiber
%        insertion. If this method is used, the input argument roi_mesh
%        must be included.
%       -column: Tracts are seeded in a user-specified column or columns
%        of the image passing through the muscle mask. The selected
%        column(s) should be included in the field ft_options.planes.
%       -row: Tracts are seeded in a user-specified row or rows of the
%        image passing through the muscle mask. The selected rows should
%        be included in the field ft_options.planes.
%       -slice: Tracts are seeded in a user-specified slice or slices of
%        the image passing through the muscle mask. The selected slices
%        should be included in the field ft_options.planes.
%       -voxel: Tracts are seeded in individual voxels. By default, tracts
%        are seeded in every voxel in the image mask. To seed every Nth voxel,
%        set ft_options.skip_vxl to N.
%
%    mesh_dist: The number of pixels to shift the mesh into the muscle,
%      prior to fiber tracking. This can be a (+) or (-) number, depending
%      on the desired direction of the shift. This field is not used for
%      planar or voxel seeding.
%
%    prop_algo: A string variable that specifies the method for determining
%      the direction of fiber tract propagation. The available options
%      include:
%        -euler: Diagonalization of the observed diffusion tensor D at the current
%          fiber tracking point, followed by Euler integration of the first
%          eigenvector.
%        -rk4: 4th order Runge-Kutta integration of the first eigenvector.
%          Note that if the 2nd, 3rd, or 4th order points fall outside of the
%          mask, the propogation algorithm automatically changes to Euler
%          integration.
%
%    step_size: The fiber-tracking step size, in pixels. A step size of 1
%      indicates a step size of one voxel width
%
%    term_mthd: A string variable that specifies the method for determining
%      whether to terminate a fiber tract or not. Any fiber tracking point
%      that falls outside of the image mask will terminate the tract.
%      Other criteria using the inter-segment angle and the FA, processed
%      according to either of two algorithms:
%        -bin1: Inter- segment angle and FA are used as binary criteria to
%         decide whether to terminate the tract. The angle used is the angle
%         formed by two fiber tracking steps. The user can decide whether to
%         calculate this angle between the current step and its immediate
%         predecessor (1-back) or between the current step and a step that
%         looks back by M points. Using the lookback option allows a tract to
%         correct its direction following an initially aberrant result. When
%         bin1 is used the tract terminates if either the angle or FA value
%         is disallowed for a single point;
%        -bin2, At each point, the inter-segment angle and FA data are
%         treated as for bin1, but two consecutive points must fail in order
%         to terminate the tract.
%
%    angle_thrsh: A two-element vector containing the angle threshold in
%      degrees and the number of look-back steps.
%
%    fa_thrsh: a two-element vector containing the lower and upper bounds
%      of allowable FA values.
%
%    depth_ratio: ratio of slice thickness/in-plane resolution. Note that
%      the function assumes equal in-plane voxel dimensions.
%
% -e1fa: A 4D matrix, with the first-third dimensions matching the
%  [row column slice] size of the DTI images and the fourth dimension
%  holding the X, Y, and Z components of the first eigenvector and the FA.
%
%  mask: The mask delimiting the muscle to be fiber-tracked.
%
%
%  The following input argument is only required for aponeurosis seeding:
%
%  roi_mesh: The mesh reconstruction of the aponeurosis of muscle fiber
%    insertion, output from define_roi.
%
%
%  The following input arguments are only required if the user wishes to
%  plot the plot tracts from within the fiber_track function.
%
%  fv_options: If specified, this calls the fiber_visualizer_v11 function to
%    plot the fiber, mask, and roi mesh.
%
%  anat_image: The structural images.
%
%OUTPUT ARGUMENTS
%  fiber_all: The fiber tracts. The rows and columns correspond to locations
%    on the roi_mesh. Dimension 3 gives point numbers on the tract, and the
%    fourth dimension has row, column, and slice coordinates, with units of
%    voxels.
%
%  roi_flag: A matrix indicating the presence of fiber tracts that propagated
%    at least 1 point
%
%  stop_list: matrix containing the reason for fiber tract termination
%    (4=mask, 3=curvature, 2=FA, or R (the last for FACT only))
%
%  fiber_len: the length, in points, of each fiber tract.
%
% seed_points: for planar and voxel seeding methods, a variable called
%   seed_points is created to hold the seed points. Including this in the list
%   of output arguments allows the user to identify the seed point locations.
%
%OTHER FUNCTIONS IN THE MUSCLE DTI FIBER-TRACKING TOOLBOX
%  For help with anisotropic smoothing, see <a href="matlab: help aniso4D_smoothing">aniso4D_smoothing</a>.
%  For help calculating the diffusion tensor, see <a href="matlab: help signal2tensor2">signal2tensor2</a>.
%  For help defining the muscle mask, see <a href="matlab: help define_muscle">define_muscle</a>.
%  For help defining the aponeurosis ROI, see <a href="matlab: help define_roi">define_roi</a>.
%  For help smoothing fiber tracts, see <a href="matlab: help fiber_smoother">fiber_smoother</a>.
%  For help quantifying fiber tracts, see <a href="matlab: help fiber_quantifier">fiber_quantifier</a>.
%  For help selecting fiber tracts following their quantification, see <a href="matlab: help fiber_goodness">fiber_goodness</a>.
%  For help visualizing fiber tracts and other structures, see <a href="matlab: help fiber_visualizer_v11">fiber_visualizer_v11</a>.
%
%VERSION INFORMATION
%  v. 1.0.0 (initial release), 15 Jan 2021, Bruce Damon
%  v. 2.0.0 (updated to allow planar or voxel seeding; remove FACT; and speed
%   tract propagation time by inputing the e1fa matrix instead of the
%   diffusion tensor matrix), 7 June 2024, Bruce Damon. Not backwards
%   compatible to v. 1
%
%ACKNOWLEDGEMENTS
%  People: Zhaohua Ding, Adam Anderson, Amanda Buck, Anneriet Heemskerk,
%    and Justin Montenegro
%  Grant support: NIH/NIAMS R01 AR050101, NIH/NIAMS R01 AR073831

%% get user options for fiber tracking
%find seeding method/get method-specific options
seed_method = ft_options.seed_method(1);

% additions for v. 2.0
% mask definitions
if length(size(mask))==4                                                    % if the input mask contains both an internal mask and a boundary mask
    internal_mask = mask(:,:,:,2);
    diffusion_mask = mask(:,:,:,1);
elseif length(size(mask))==3                                                % if the input mask contains only a boundary mask
    internal_mask = zeros(size(mask));                                      % dummy variable for later subtraction step
    diffusion_mask = mask;
else
    error('Exiting due to inappropriate mask dimensions.')
end

% seed point definitions - additions for v. 2.0
switch seed_method
    case{'a'}
        if ~exist('roi_mesh', 'var')
            error('The input argument roi_mesh was not found.  Exiting fiber tracking.')
        else
            mesh_dist = ft_options.mesh_dist;
        end

    case{'r'}
        if ~isfield(ft_options, 'planes')                                    %if the user specificies the mesh plane location(s)
            error('The field ft_options.planes was not found.  Exiting fiber tracking.');
        end

        % get the planes of interest
        planes = ft_options.planes;
        seed_points = zeros(2, length(planes), 6);                          %initialize roi_mesh for the planar case, pad with extra zeros in dim. 3
        for i=1:length(planes)                                          %define the mesh as the row, column, and slice coordinates
            loop_plane = planes(i);
            loop_mask = squeeze(diffusion_mask(loop_plane,:,:));
            [seed_points_c, seed_points_s] = find(loop_mask);
            seed_points_r = loop_plane*ones(size(seed_points_c));
            seed_points(1:length(seed_points_r),i,1) = seed_points_r;
            seed_points(1:length(seed_points_r),i,2) = seed_points_c;
            seed_points(1:length(seed_points_r),i,3) = seed_points_s;
        end

    case{'c'}
        if ~isfield(ft_options, 'planes')                                    %if the user specificies the mesh plane location(s)
            error('The field ft_options.planes was not found.  Exiting fiber tracking.');
        end

        % get the planes of interest
        planes = ft_options.planes;
        seed_points = zeros(2, length(planes), 6);                          %initialize roi_mesh for the planar case, pad with extra zeros in dim. 3
        for i=1:length(planes)                                          %define the mesh as the row, column, and slice coordinates
            loop_plane = planes(i);
            loop_mask = squeeze(diffusion_mask(:,loop_plane,:));
            [seed_points_r, seed_points_s] = find(loop_mask);
            seed_points_c = loop_plane*ones(size(seed_points_r));
            seed_points(1:length(seed_points_r),i,1) = seed_points_r;
            seed_points(1:length(seed_points_r),i,2) = seed_points_c;
            seed_points(1:length(seed_points_r),i,3) = seed_points_s;
        end

    case{'s'}
        if ~isfield(ft_options, 'planes')                                    %if the user specificies the mesh plane location(s)
            error('The field ft_options.planes was not found.  Exiting fiber tracking.');
        end

        % get the planes of interest
        planes = ft_options.planes;
        seed_points = zeros(2, length(planes), 6);                          %initialize roi_mesh for the planar case, pad with extra zeros in dim. 3
        for i=1:length(planes)                                          %define the mesh as the row, column, and slice coordinates
            loop_plane = planes(i);
            loop_mask = squeeze(diffusion_mask(:,:,loop_plane));
            [seed_points_r, seed_points_c] = find(loop_mask);
            seed_points_s = loop_plane*ones(size(seed_points_r));
            seed_points(1:length(seed_points_r),i,1) = seed_points_r;
            seed_points(1:length(seed_points_r),i,2) = seed_points_c;
            seed_points(1:length(seed_points_r),i,3) = seed_points_s;
        end

    case{'v'}

        mask_seed=zeros(size(diffusion_mask));

        % get skip setting from options structure; if missing, default is 1
        if isfield(ft_options, 'skip_vxl')
            skip_vxl = ft_options.skip_vxl;
        else
            skip_vxl=1;
        end

        % define seed points at skip_vxl spacing; alternate starting point by slice
        if skip_vxl==0

            for s=1:length(diffusion_mask(1,1,:))

                loop_mask_orig = diffusion_mask(:,:,s);
                %                 retained_fraction = min(1, 1.15 - 0.0015*length(find(loop_mask_orig)));
                %                 retained_fraction = max(retained_fraction, 0.2);
                %                 loop_indices = find(loop_mask_orig);
                %                 random_indices = randperm(length(loop_indices), round(retained_fraction*length(find(loop_mask_orig))));
                %                 retained_indices = loop_indices(random_indices);
                %
                %                 loop_mask = zeros(size(loop_mask_orig));
                %                 loop_mask(retained_indices) = 1;
                %                 loop_mask(loop_mask>1) = 1;

                loop_mask_size = length(find(loop_mask_orig));
                skip_vxl = min(ceil(loop_mask_size/100), 6);
                loop_mask = zeros(size(loop_mask_orig));
                loop_offset = mod(s, 2);
                loop_mask((1+loop_offset):end, (1+loop_offset):skip_vxl:end) = 1;
                mask_seed(:,:,s) = loop_mask.*bwmorph(diffusion_mask(:,:,s), 'erode') - internal_mask(:,:,s);
            end
            mask_seed(mask_seed<0) = 0;

        elseif skip_vxl==1

            mask_seed = diffusion_mask;

        elseif skip_vxl>1

            for s=1:length(diffusion_mask(1,1,:))
                loop_offset = mod(s, skip_vxl);
                mask_seed((1+loop_offset):skip_vxl:end, (1+loop_offset):skip_vxl:end, s) = 1;
                mask_seed(:,:,s) = mask_seed(:,:,s).*bwmorph(diffusion_mask(:,:,s), 'erode');
            end

        end

        mask_seed = mask_seed.*diffusion_mask;
        mask_idx = find(mask_seed);
        [seed_points_r, seed_points_c, seed_points_s] = ind2sub(size(mask_seed), mask_idx);
        seed_points = zeros(length(seed_points_r), 1, 6);
        seed_points(:,1,1) = seed_points_r;
        seed_points(:,1,2) = seed_points_c;
        seed_points(:,1,3) = seed_points_s;

    case{'e'}

        mask_seed=zeros(size(diffusion_mask));

        for s=1:length(diffusion_mask(1,1,:))

            loop_mask_orig = diffusion_mask(:,:,s);
            loop_edge = bwmorph(loop_mask_orig, 'erode') - bwmorph(loop_mask_orig, 'erode', 2);
            loop_mask = zeros(size(loop_mask_orig));
            loop_mask = loop_mask + loop_edge;
            loop_mask(loop_mask>1) = 1;

            mask_seed(:,:,s) = loop_mask;
        end
        mask_seed = mask_seed.*diffusion_mask;

        mask_idx = find(mask_seed);
        [seed_points_r, seed_points_c, seed_points_s] = ind2sub(size(mask_seed), mask_idx);
        seed_points = zeros(length(seed_points_r), 1, 6);
        seed_points(:,1,1) = seed_points_r;
        seed_points(:,1,2) = seed_points_c;
        seed_points(:,1,3) = seed_points_s;
end

%basic tracking options
depth_ratio = ft_options.depth_ratio;

%get the frame of reference and image orientation to set conventions for
%propagating points
ref_frame = ft_options.ref_frame;
image_orient = ft_options.image_orient;

switch image_orient
    case{'AL'}                                                              %image top is anatomical anterior, image right is anatomical left
        switch ref_frame
            case{'LPS'}                                                     %image left, posterior, and superior are (+) directions ()
                e1_order = [2 1 3];
                e1r_sign = 1;                                               %sign of E1 component in row direction
                e1c_sign = 1;                                               %sign of E1 component in column direction
                e1s_sign = 1;                                               %sign of E1 component in slice direction

            case{'LAS'}                                                     %image left, anterior, and superior are (+) directions
                e1_order = [2 1 3];
                e1r_sign = -1;
                e1c_sign = 1;
                e1s_sign = 1;

            case{'RAS'}                                                     %image right, anterior, and superior are (+) directions
                e1_order = [2 1 3];
                e1r_sign = -1;
                e1c_sign = -1;
                e1s_sign = 1;
        end
    case {'RA'}                                                              %image top is anatomical right, image right is anatomical anterior
        switch ref_frame
            case{'LPS'}
                e1_order = [1 2 3];
                e1r_sign = 1;
                e1c_sign = -1;
                e1s_sign = 1;

            case{'LAS'}
                e1_order = [1 2 3];
                e1r_sign = 1;
                e1c_sign = 1;
                e1s_sign = 1;

            case{'RAS'}                                                     %image left, anterior, and superior are (+) directions (NIFTI)
                e1_order = [1 2 3];
                e1r_sign = -1;
                e1c_sign = 1;
                e1s_sign = 1;

        end
end

% get the tract propagation options
prop_algo = ft_options.prop_algo(1:2);

switch prop_algo
    case {'eu'}
        % tracking options
        step_incr = ft_options.step_size;

        %tract termination options
        term_mthd = ft_options.term_mthd;
        angle_thrsh=ft_options.angle_thrsh(1);
        angle_lookback = ft_options.angle_thrsh(2);
        fa_low=ft_options.fa_thrsh(1);
        fa_high=ft_options.fa_thrsh(2);

    case {'rk'}
        % tracking options
        step_incr = ft_options.step_size;

        %tract termination options
        term_mthd = ft_options.term_mthd;
        angle_thrsh=ft_options.angle_thrsh(1);
        angle_lookback = ft_options.angle_thrsh(2);
        fa_low=ft_options.fa_thrsh(1);
        fa_high=ft_options.fa_thrsh(2);
end


%% fiber track - aponeurosis method

% seed_method = aponeurosis is new to v. 2.0
if seed_method(1)=='a'

    % initialize matrix to hold the x, y, and z coordinates of the fiber tracts, with units of pixels
    fiber_all = zeros(length(roi_mesh(:, 1, 1)), length(roi_mesh(1, :, 1)), 50, 3);    %initialize as roi mesh size by up to 50 points long; can grow, as needed
    fa_all = zeros(length(roi_mesh(:, 1, 1)), length(roi_mesh(1, :, 1)), 50);    %initialize as roi mesh size by up to 50 points long; can grow, as needed

    % initialize quality-related variables
    roi_flag = zeros(size(squeeze(roi_mesh(:, :, 1))));                         %will hold information about tract initiation
    stop_list=zeros(length(roi_mesh(:,1,1)),length(roi_mesh(1,:,1)));         %will hold information about the stop criteria
    fiber_len=roi_flag;

    % fiber track by passing through rows and columns of the mesh
    for row_cntr=1:length(roi_mesh(:, 1, 1))                                  %start of the row loop

        for col_cntr=1:length(roi_mesh(1, :, 1))                              %start of the column loop

            % initialize some variables
            fiber_cntr = 1;                                                     %counter for the fiber tract points; initialize at 1 for each tract
            angle_all = zeros(50,1);                                           %keep track of all inter-step angles
            eigvector1_all = zeros(50,3);                                      %matrix to hold all of the first eigenvectors (needed for checking inter-point angle stop criterion)

            % 1) locate the seed point on the roi_mesh
            seed_point(1:3) = [roi_mesh(row_cntr, col_cntr, 1)+mesh_dist*roi_mesh(row_cntr, col_cntr, 4) ...
                roi_mesh(row_cntr, col_cntr, 2)+mesh_dist*roi_mesh(row_cntr, col_cntr, 5) ...
                roi_mesh(row_cntr, col_cntr, 3)+mesh_dist*roi_mesh(row_cntr, col_cntr, 6)];

            % addition for v. 2.0: change to rounding procedures for roi_mesh points
            e1fa_rf = floor(seed_point(1));                                   %round the slice coordiante; but take floor and ceiling for the row and column values
            e1fa_rc = ceil(seed_point(1));
            e1fa_cf = floor(seed_point(2));
            e1fa_cc = ceil(seed_point(2));
            e1fa_s = max([round(seed_point(3)) 1]);                           %just in case the z coordinate rounds to zero, make it a one

            % 2) add to the fiber_all matrix

            % calculate closest point that has mask=1
            mask_distance = zeros(4, 4);
            seed_point_m = [e1fa_rf e1fa_cf; e1fa_rf e1fa_cc; e1fa_rc e1fa_cf; e1fa_rc e1fa_cc];
            for seed_cntr = 1:4
                mask_distance(seed_cntr, 1) = diffusion_mask(seed_point_m(seed_cntr, 1), seed_point_m(seed_cntr, 2), e1fa_s);
                mask_distance(seed_cntr, 2) = (sum([seed_point_m(seed_cntr, 1) seed_point_m(seed_cntr, 2)] - ...
                    [seed_point(1) seed_point(2)]).^2).^0.5;
                mask_distance(seed_cntr, 3:4) = [seed_point_m(seed_cntr, 1) seed_point_m(seed_cntr, 2)];
            end
            if sum(mask_distance(:,1))==0
                continue
            else
                mask_distance(:,2) = mask_distance(:,1).*mask_distance(:,2);
                mask_distance = mask_distance(mask_distance(:,2)>0, :);
                mask_distance = sortrows(mask_distance, 2);
            end

            % use this points to look up initial values for mask and diff. tensor
            e1fa_r = mask_distance(1,3);
            e1fa_c = mask_distance(1,4);

            % replaced with above in v 2.0
            %             e1fa_r = round(seed_point(1));                                    %round the coordiantes to get a row, column, and slice location in the tensor matrix
            %             e1fa_c = round(seed_point(2));                                    %round the coordiantes to get a row, column, and slice location in the tensor matrix
            %             e1fa_s = max([round(seed_point(3)) 1]);                           %just in case the z coordinate rounds to zero, make it a one

            fiber_all(row_cntr,col_cntr, fiber_cntr,:) = seed_point;

            % 3) get the initial direction
            switch prop_algo

                case {'eu'}

                    % get the eigenvector
                    eigvector1 = squeeze(e1fa(e1fa_r, e1fa_c, e1fa_s, 1:3));
                    eigvector1 = check_E1_sign(seed_method, eigvector1, 1);         %correct, as needed, for antipodal symmetry
                    eigvector1_all(fiber_cntr,:) = eigvector1;

                    % store the FA data for termination checking
                    fa_all(row_cntr, col_cntr,fiber_cntr) = e1fa(e1fa_r, e1fa_c, e1fa_s, 4);

                    step_dir = eigvector1;                                              %initially set the step direction as E1
                    step_dir = step_dir(e1_order);                                      %adjust for frame of reference and row/column indexing
                    step_dir = [e1r_sign*step_dir(1) e1c_sign*step_dir(2) e1s_sign*step_dir(3)]';
                    step_dir(3) = step_dir(3)/depth_ratio;                              %account for the slice aspect ratio

                case {'rk'} %Runge-Kutta integration (Richard Hamming. 1989 Introduction to Applied Numerical Analysis. Dover Publications, p. 212)

                    % set default order for RK4
                    rk_order=1;

                    % for point p1, get the direction
                    point1 = squeeze(fiber_all(row_cntr,col_cntr, fiber_cntr,:));

                    % use indices from above
                    eigvector1_point1 = squeeze(e1fa(e1fa_r, e1fa_c, e1fa_s, 1:3));
                    eigvector1_point1 = check_E1_sign(seed_method, eigvector1_point1, 1);         %correct, as needed, for antipodal symmetry

                    % store the FA data for termination checking
                    fa_all(row_cntr, col_cntr,fiber_cntr) = e1fa(e1fa_r, e1fa_c, e1fa_s, 4);

                    % get path direction as the first eigenvector
                    k1 = eigvector1_point1;                                         %initially set the step direction as E1
                    k1 = k1(e1_order);                                              %adjust for frame of reference and row/column indexing
                    k1 = [e1r_sign*k1(1) e1c_sign*k1(2) e1s_sign*k1(3)]';
                    k1(3) = k1(3)/depth_ratio;                                      %account for slice aspect ratio in the z direction.

                    %calculate next point, p2
                    point2 = point1 + (step_incr/2)*k1;

                    % if mask at point 2 is valid, get the diffusion tensor there and calculate point three
                    if diffusion_mask(round(point2(1)), round(point2(2)), max([round(point2(3)) 1]))==1

                        e1fa_r = round(point2(1));
                        e1fa_c = round(point2(2));
                        e1fa_s = max([round(point2(3)) 1]);
                        eigvector1_point2 = squeeze(e1fa(e1fa_r, e1fa_c, e1fa_s, 1:3));
                        eigvector1_point2 = check_E1_sign(seed_method, eigvector1_point2, 1);         %correct, as needed, for antipodal symmetry

                        % get path direction as the first eigenvector
                        k2 = eigvector1_point2;                                         %initially set the step direction as E1
                        k2 = k2(e1_order);                                              %adjust for frame of reference and row/column indexing
                        k2 = [e1r_sign*k2(1) e1c_sign*k2(2) e1s_sign*k2(3)]';
                        k2(3) = k2(3)/depth_ratio;                                      %account for slice aspect ratio in the z direction.

                        %calculate next point, p3
                        point3 = point2 + (step_incr/2)*k2;

                        % if mask at point 3 is valid, get the diffusion tensor there and calculate point 4
                        if diffusion_mask(round(point3(1)), round(point3(2)), max([round(point3(3)) 1]))==1

                            e1fa_r = round(point3(1));
                            e1fa_c = round(point3(2));
                            e1fa_s = max([round(point3(3)) 1]);
                            eigvector1_point3 = squeeze(e1fa(e1fa_r, e1fa_c, e1fa_s, 1:3));
                            eigvector1_point3 = check_E1_sign(seed_method, eigvector1_point3, 1);         %correct, as needed, for antipodal symmetry

                            % get path direction as the first eigenvector
                            k3 = eigvector1_point3;                                         %initially set the step direction as E1
                            k3 = k3(e1_order);                                              %adjust for frame of reference and row/column indexing
                            k3 = [e1r_sign*k3(1) e1c_sign*k3(2) e1s_sign*k3(3)]';
                            k3(3) = k3(3)/depth_ratio;                                      %account for slice aspect ratio in the z direction.

                            %calculate next point, p4
                            point4 = point3 + step_incr*k3;

                            % if mask at point 4 is valid, get the diffusion tensor there
                            if diffusion_mask(round(point4(1)), round(point4(2)), max([round(point4(3)) 1]))==1
                                rk_order=4;

                                e1fa_r = round(point4(1));
                                e1fa_c = round(point4(2));
                                e1fa_s = max([round(point4(3)) 1]);
                                eigvector1_point4 = squeeze(e1fa(e1fa_r, e1fa_c, e1fa_s, 1:3));
                                eigvector1_point4 = check_E1_sign(seed_method, eigvector1_point4, 1);         %correct, as needed, for antipodal symmetry

                            end                                                             %of point 4 is valid if statement

                        end                                                                 %of point 3 is valid if statement

                    end                                                                     %of point 2 is valid if statement

                    if rk_order==4

                        %get the step direction from the average of the four eigenvectors.  First, find the mean dyadic tensor
                        dyadic_tensor_1 = eigvector1_point1*eigvector1_point1';
                        dyadic_tensor_2 = eigvector1_point2*eigvector1_point2';
                        dyadic_tensor_3 = eigvector1_point3*eigvector1_point3';
                        dyadic_tensor_4 = eigvector1_point4*eigvector1_point4';
                        dyadic_tensor_mean = (dyadic_tensor_1 + 2*dyadic_tensor_2 + 2*dyadic_tensor_3 + dyadic_tensor_4)/6;
                        [dyadic_tensor_eigvector_m, dyadic_tensor_eigvalue_m] = eig(dyadic_tensor_mean);                     %diagonalize to get e-vectors and e-values

                        %get the derived indices
                        [~, sort_order] = sort(diag(dyadic_tensor_eigvalue_m), 'descend');
                        eigvector1_mean = dyadic_tensor_eigvector_m(:, sort_order(1));
                        eigvector1 = check_E1_sign(seed_method, eigvector1_mean, 1);           %correct, as needed, for antipodal symmetry

                        eigvector1_all(fiber_cntr,:) = eigvector1;
                        step_dir = eigvector1;                                              %initially set the step direction as E1
                        step_dir = step_dir(e1_order);                                      %adjust for frame of reference and row/column indexing
                        step_dir = [e1r_sign*step_dir(1) e1c_sign*step_dir(2) e1s_sign*step_dir(3)]';
                        step_dir(3) = step_dir(3)/depth_ratio;                              %account for the slice aspect ratio

                    elseif rk_order==1

                        eigvector1 = eigvector1_point1;
                        eigvector1_all(fiber_cntr,:) = eigvector1;

                        step_dir = eigvector1;                                              %initially set the step direction as E1
                        step_dir = step_dir(e1_order);                                      %adjust for frame of reference and row/column indexing
                        step_dir = [e1r_sign*step_dir(1) e1c_sign*step_dir(2) e1s_sign*step_dir(3)]';
                        step_dir(3) = step_dir(3)/depth_ratio;                              %account for the slice aspect ratio

                    else

                        break

                    end

            end

            % add the point
            next_point = squeeze(fiber_all(row_cntr,col_cntr, fiber_cntr,:)) + ...
                step_incr(1)*step_dir;                                              %multiply step_dir by the step size before adding the step to the preceding point
            fiber_cntr = fiber_cntr+1;                                              %increment the fiber counter

            % get indices into tensor matrix for the next point
            e1fa_r = round(next_point(1));
            e1fa_c = round(next_point(2));
            e1fa_s = max([round(next_point(3)) 1]);                               %take teh maximum with slice 1, in case of low seed point

            % v. 2.0: check teh dilated mask; terminate tract if out of bounds. record as stop criterion = 4 and continue to next point
%             if dilated_mask(e1fa_r,e1fa_c,e1fa_s)==0
            if diffusion_mask(e1fa_r,e1fa_c,e1fa_s)==0                        %removed in v 2.0
                stop_list(row_cntr,col_cntr) = 4;
                continue;
            else
                fiber_all(row_cntr,col_cntr, fiber_cntr,:) = next_point;
            end

            % 4) begin the fiber tracking loop
            while 1

                % 4a) get the step direction
                switch prop_algo

                    case {'eu'}                                                %if Euler integration is selected

                        % v. 2.0: if within the diffusion mask, look up new E1 and FA. Otherwise, use previous values
%                         if diffusion_mask(e1fa_r, e1fa_c, e1fa_s)

                            % get the eigenvector
                            eigvector1 = squeeze(e1fa(e1fa_r, e1fa_c, e1fa_s, 1:3));
                            eigvector1 = check_E1_sign(seed_method, eigvector1, 1);         %correct, as needed, for antipodal symmetry
                            eigvector1_all(fiber_cntr,:) = eigvector1;

                            % store the FA data for termination checking
                            fa_all(row_cntr, col_cntr, fiber_cntr) = e1fa(e1fa_r, e1fa_c, e1fa_s, 4);

%                         else
% 
%                             eigvector1 = eigvector1_old;
%                             eigvector1_all(fiber_cntr,:) = eigvector1;
%                             fa_all(row_cntr, col_cntr, fiber_cntr) = fa_all(row_cntr, col_cntr, fiber_cntr-1);
% 
%                         end

                        step_dir = eigvector1;                                              %initially set the step direction as E1
                        step_dir = step_dir(e1_order);                                      %adjust for frame of reference and row/column indexing
                        step_dir = [e1r_sign*step_dir(1) e1c_sign*step_dir(2) e1s_sign*step_dir(3)]';
                        step_dir(3) = step_dir(3)/depth_ratio;                              %account for the slice aspect ratio

                    case {'rk'}

                        %create a variable to keep track of the RK order; use order = 1 (Euler's method) unless all four points are inside the mask
                        rk_order=1;

                        % define point p1; if it lies within the mask, get the diffusion tensor by calling the retrieve_tensor function.
                        point1 = squeeze(fiber_all(row_cntr,col_cntr, fiber_cntr,:));

                        if round(point1(1))>length(e1fa(:,1,1,1,1)) || round(point1(2))>length(e1fa(1,:,1,1,1)) || ...
                                round(point1(3))>length(e1fa(1,1,:,1,1))
                            stop_list(row_cntr,col_cntr) = 4;
                            break;
                        end

                        if diffusion_mask(round(point2(1)), round(point2(2)), max([round(point2(3)) 1]))==1

                            e1fa_r = round(point1(1));
                            e1fa_c = round(point1(2));
                            e1fa_s = max([round(point1(3)) 1]);
                            eigvector1_point1 = squeeze(e1fa(e1fa_r, e1fa_c, e1fa_s, 1:3));
                            eigvector1_point1 = check_E1_sign(seed_method, eigvector1_point1, 1);         %correct, as needed, for antipodal symmetry

                            % store the FA data for termination checking
                            fa_all(row_cntr, col_cntr,fiber_cntr) = e1fa(e1fa_r, e1fa_c, e1fa_s, 4);

                            % get path direction as the first eigenvector
                            k1 = eigvector1_point1;                                         %initially set the step direction as E1
                            k1 = k1(e1_order);                                              %adjust for frame of reference and row/column indexing
                            k1 = [e1r_sign*k1(1) e1c_sign*k1(2) e1s_sign*k1(3)]';
                            k1(3) = k1(3)/depth_ratio;                                      %account for slice aspect ratio in the z direction.

                            %calculate next point, p2
                            point2 = point1 + (step_incr/2)*k1;

                        else
                            stop_list(row_cntr, col_cntr)=4;
                            break

                        end

                        % if mask at point 2 is valid, get the diffusion tensor there and calculate point three
                        if diffusion_mask(round(point2(1)), round(point2(2)), max([round(point2(3)) 1]))==1

                            e1fa_r = round(point2(1));
                            e1fa_c = round(point2(2));
                            e1fa_s = max([round(point2(3)) 1]);
                            eigvector1_point2 = squeeze(e1fa(e1fa_r, e1fa_c, e1fa_s, 1:3));
                            eigvector1_point2 = check_E1_sign(seed_method, eigvector1_point2, 1);         %correct, as needed, for antipodal symmetry

                            % get path direction as the first eigenvector
                            k2 = eigvector1_point2;                                         %initially set the step direction as E1
                            k2 = k2(e1_order);                                              %adjust for frame of reference and row/column indexing
                            k2 = [e1r_sign*k2(1) e1c_sign*k2(2) e1s_sign*k2(3)]';
                            k2(3) = k2(3)/depth_ratio;                                      %account for slice aspect ratio in the z direction.

                            %calculate next point, p3
                            point3 = point2 + (step_incr/2)*k2;

                            % if mask at point 3 is valid, get the diffusion tensor there and calculate point 4
                            if diffusion_mask(round(point3(1)), round(point3(2)), max([round(point3(3)) 1]))==1

                                e1fa_r = round(point3(1));
                                e1fa_c = round(point3(2));
                                e1fa_s = max([round(point3(3)) 1]);
                                eigvector1_point3 = squeeze(e1fa(e1fa_r, e1fa_c, e1fa_s, 1:3));
                                eigvector1_point3 = check_E1_sign(seed_method, eigvector1_point3, 1);         %correct, as needed, for antipodal symmetry

                                % get path direction as the first eigenvector
                                k3 = eigvector1_point3;                                         %initially set the step direction as E1
                                k3 = k3(e1_order);                                              %adjust for frame of reference and row/column indexing
                                k3 = [e1r_sign*k3(1) e1c_sign*k3(2) e1s_sign*k3(3)]';
                                k3(3) = k3(3)/depth_ratio;                                      %account for slice aspect ratio in the z direction.

                                %calculate next point, p4
                                point4 = point3 + step_incr*k3;

                                % if mask at point 4 is valid, get the diffusion tensor there
                                if diffusion_mask(round(point4(1)), round(point4(2)), max([round(point4(3)) 1]))==1
                                    rk_order=4;

                                    e1fa_r = round(point4(1));
                                    e1fa_c = round(point4(2));
                                    e1fa_s = max([round(point4(3)) 1]);
                                    eigvector1_point4 = squeeze(e1fa(e1fa_r, e1fa_c, e1fa_s, 1:3));
                                    eigvector1_point4 = check_E1_sign(seed_method, eigvector1_point4, 1);         %correct, as needed, for antipodal symmetry

                                end                                                             %of point 4 is valid if statement

                            end                                                                 %of point 3 is valid if statement

                        end                                                                     %of point 2 is valid if statement

                        if rk_order==4

                            %get the step direction from the average of the four eigenvectors.  First, find the mean dyadic tensor
                            dyadic_tensor_1 = eigvector1_point1*eigvector1_point1';
                            dyadic_tensor_2 = eigvector1_point2*eigvector1_point2';
                            dyadic_tensor_3 = eigvector1_point3*eigvector1_point3';
                            dyadic_tensor_4 = eigvector1_point4*eigvector1_point4';
                            dyadic_tensor_mean = (dyadic_tensor_1 + 2*dyadic_tensor_2 + 2*dyadic_tensor_3 + dyadic_tensor_4)/6;
                            [dyadic_tensor_eigvector_m, dyadic_tensor_eigvalue_m] = eig(dyadic_tensor_mean);                     %diagonalize to get e-vectors and e-values

                            %get the derived indices
                            [~, sort_order] = sort(diag(dyadic_tensor_eigvalue_m), 'descend');
                            eigvector1_mean = dyadic_tensor_eigvector_m(:, sort_order(1));
                            eigvector1 = check_E1_sign(seed_method, eigvector1_mean, 1);           %correct, as needed, for antipodal symmetry

                            eigvector1_all(fiber_cntr,:) = eigvector1;
                            step_dir = eigvector1;                                              %initially set the step direction as E1
                            step_dir = step_dir(e1_order);                                      %adjust for frame of reference and row/column indexing
                            step_dir = [e1r_sign*step_dir(1) e1c_sign*step_dir(2) e1s_sign*step_dir(3)]';
                            step_dir(3) = step_dir(3)/depth_ratio;                              %account for the slice aspect ratio

                        elseif rk_order==1

                            eigvector1 = eigvector1_point1;
                            eigvector1_all(fiber_cntr,:) = eigvector1;

                            step_dir = eigvector1;                                              %initially set the step direction as E1
                            step_dir = step_dir(e1_order);                                      %adjust for frame of reference and row/column indexing
                            step_dir = [e1r_sign*step_dir(1) e1c_sign*step_dir(2) e1s_sign*step_dir(3)]';
                            step_dir(3) = step_dir(3)/depth_ratio;                              %account for the slice aspect ratio

                        else

                            break

                        end

                end                                                             %end of step direction switch statement

                % 4b) decide whether or not to add the point

                % get indices into e1fa matrix for the next point
                e1fa_r = round(next_point(1));                                    %round the coordiantes to get a row, column, and slice location in the tensor matrix
                e1fa_c = round(next_point(2));
                e1fa_s = max([round(next_point(3)) 1]);                           %just in case the z coordinate rounds to zero, make it a one

                % v. 2.0: use the dilated mask to ensure that tracts definitely propagate to muscle boundary. terminate tract if out of bounds and record as stop criterion = 4
%                 if dilated_mask(e1fa_r,e1fa_c,e1fa_s)==0
                if diffusion_mask(e1fa_r,e1fa_c,e1fa_s)==0  % removed in v. 2.0
                    stop_list(row_cntr,col_cntr) = 4;
                    break;
                end

                %other criteria depend on termination algorithm selected
                switch term_mthd                                                %select between tract termination methods

                    case {'bin1'}                                                %binary OR score, 1 point required for termination

                        %check for out of bounds FA values
                        if fa_all(row_cntr, col_cntr,fiber_cntr)<fa_low || fa_all(row_cntr, col_cntr,fiber_cntr)>fa_high
                            stop_list(row_cntr, col_cntr) = 2;
                            fiber_len(row_cntr, col_cntr) = fiber_cntr;
                            break;
                        end

                        %check for out of bounds curvature values
                        if fiber_cntr > angle_lookback
                            eigvector1_old = squeeze(eigvector1_all((fiber_cntr-angle_lookback),:));
                            step_angle = abs(acosd(dot(eigvector1, eigvector1_old)));
                            angle_all(fiber_cntr) = step_angle;
                            if step_angle>angle_thrsh
                                stop_list(row_cntr, col_cntr) = 3;
                                fiber_len(row_cntr, col_cntr) = fiber_cntr;
                                break;
                            end
                        end

                    case {'bin2'}                                                %binary OR score, 2 consecutive points required

                        %check for out of bounds FA values
                        if fa_all(row_cntr, col_cntr,fiber_cntr)<fa_low || fa_all(row_cntr, col_cntr,fiber_cntr)>fa_high
                            if fa_all(row_cntr, col_cntr,fiber_cntr-1)<fa_low || fa_all(row_cntr, col_cntr,fiber_cntr-1)>fa_high
                                stop_list(row_cntr, col_cntr) = 2;
                                fiber_len(row_cntr, col_cntr) = fiber_cntr;
                                break;
                            end
                        end

                        %check for out of bounds curvature values
                        if fiber_cntr > (angle_lookback+1)
                            eigvector1_old = squeeze(eigvector1_all((fiber_cntr-angle_lookback),:));
                            step_angle_old = abs(acosd(dot(eigvector1, eigvector1_old)));

                            eigvector1_last = squeeze(eigvector1_all((fiber_cntr-1),:));
                            step_angle_current = abs(acosd(dot(eigvector1, eigvector1_last)));
                            angle_all(fiber_cntr) = step_angle_current;

                            if step_angle_current>angle_thrsh && step_angle_old>angle_thrsh
                                stop_list(row_cntr, col_cntr) = 3;
                                fiber_len(row_cntr, col_cntr) = fiber_cntr;
                                break;
                            end
                        end

                end                                                             %end of termination method switch statement


                % 5) add the point and continue tracking
                fiber_cntr = fiber_cntr + 1;                                     %increment the fiber counter
                fiber_len(row_cntr, col_cntr) = fiber_cntr;
                next_point = squeeze(fiber_all(row_cntr,col_cntr, fiber_cntr-1,:)) + ...
                    step_incr*step_dir;                                         %multiply step_dir by the step increment before adding the step to the preceding point
                fiber_all(row_cntr,col_cntr, fiber_cntr,:) = next_point;        %add the point

                roi_flag(row_cntr, col_cntr)=1;                                 %record as a tracked fiber

            end                                                                 %end of fiber tracking while loop

        end                                                                     %end of column for loop

    end                                                                         %end of row for loop

end


%% fiber track - non aponeurosis methods

% planar, voxel, and edge seeding methods are new for v. 2.0
if seed_method(1)~='a'

    % initialize matrix to hold the r, c, and s coordinates of the fiber tracts
    % create separate matrices for fibers to be tracked up and fibers to be tracked down; merge later.
    fiber_all_1 = zeros(length(seed_points(:,1,1)), length(seed_points(1,:,1)), 50, 3); %initialize as roi mesh size by up to 50 points per tract
    fa_all_1 = zeros(length(seed_points(:,1,1)), length(seed_points(1,:,1)), 50);
    fiber_all_2 = zeros(length(seed_points(:,1,1)), length(seed_points(1,:,1)), 50, 3);
    fa_all_2 = zeros(length(seed_points(:,1,1)), length(seed_points(1,:,1)), 50);

    % initialize quality-related variables
    roi_flag_1 = zeros(length(seed_points(:,1,1)), length(seed_points(1,:,1)));  %will hold information about tract initiation
    stop_list_1 = zeros(length(seed_points(:,1,1)), length(seed_points(1,:,1)));
    fiber_len_1 = zeros(length(seed_points(:,1,1)), length(seed_points(1,:,1)));
    roi_flag_2 = zeros(length(seed_points(:,1,1)), length(seed_points(1,:,1)));
    stop_list_2 = zeros(length(seed_points(:,1,1)), length(seed_points(1,:,1)));
    fiber_len_2 = zeros(length(seed_points(:,1,1)), length(seed_points(1,:,1)));

    % fiber track by passing through rows and columns of the mesh
    for plane_cntr=1:length(seed_points(1,:,1))

        for seed_cntr=1:length(seed_points(:,plane_cntr,1))  %start of the seed point loop

            %%%%%%%%%%%%%%% Track in direction 1 %%%%%%%%%%%%%%

            % initialize some variables
            fiber_cntr = 1;                                                     %counter for the fiber tract points; initialize at 1 for each tract
            angle_1_all = zeros(50,1);                                           %keep track of all inter-step angles
            eigvector1_1_all = zeros(50,3);                                      %matrix to hold all of the first eigenvectors (needed for checking inter-point angle stop criterion)

            % 1) locate the seed point on the seed_point_1 matrix. uses standard rounding criteria, because points will be guaranteed to be int eh mask
            seed_point_1(1:3) = [seed_points(seed_cntr, plane_cntr, 1)...
                seed_points(seed_cntr, plane_cntr, 2)...
                seed_points(seed_cntr, plane_cntr, 3)];
            e1fa_r_1 = round(seed_point_1(1));                                    %round the coordinates to get a row, column, and slice location in the tensor matrix
            e1fa_c_1 = round(seed_point_1(2));
            e1fa_s_1 = max([round(seed_point_1(3)) 1]);                      %just in case the z coordinate rounds to zero, make it a one

            % 2) add to the fiber_all matrix. does not use the dilated mask
            if diffusion_mask(e1fa_r_1,e1fa_c_1,e1fa_s_1)==0
                stop_list_1(seed_cntr,plane_cntr)=4;
                continue;
            end

            fiber_all_1(seed_cntr, plane_cntr, fiber_cntr,:) = seed_point_1;

            % 3) get the initial direction (NB initial direction methods not updated for v. 2.0 because the seed points are, by definition, within the mask

            switch prop_algo
                case {'eu'}

                    % get the eigenvector
                    eigvector1_1 = squeeze(e1fa(e1fa_r_1, e1fa_c_1, e1fa_s_1, 1:3));
                    eigvector1_1 = check_E1_sign(seed_method, eigvector1_1, 1);                        %correct, as needed, for antipodal symmetry
                    eigvector1_1_all(fiber_cntr,:) = eigvector1_1;
                    eigvector1_1_old = eigvector1_1;

                    % store the FA data for termination checking
                    fa_all_1(seed_cntr, plane_cntr, fiber_cntr) = e1fa(e1fa_r_1, e1fa_c_1, e1fa_s_1, 4);

                    step_dir_1 = eigvector1_1;                                              %initially set the step direction as E1
                    step_dir_1 = step_dir_1(e1_order);                                      %adjust for frame of reference and row/column indexing
                    step_dir_1 = [e1r_sign*step_dir_1(1) e1c_sign*step_dir_1(2) e1s_sign*step_dir_1(3)]';
                    step_dir_1(3) = step_dir_1(3)/depth_ratio;                              %account for the slice aspect ratio

                case {'rk'}

                    % for point p1, get the direction
                    point1 = squeeze(fiber_all_1(seed_cntr, plane_cntr, fiber_cntr,:));

                    e1fa_r_1 = round(point1(1));
                    e1fa_c_1 = round(point1(2));
                    e1fa_s_1 = max([round(point1(3)) 1]);
                    eigvector1_point1 = squeeze(e1fa(e1fa_r_1, e1fa_c_1, e1fa_s_1, 1:3));
                    eigvector1_point1 = check_E1_sign(seed_method, eigvector1_point1, 1);         %correct, as needed, for antipodal symmetry

                    % store the FA data for termination checking
                    fa_all_1(seed_cntr, plane_cntr,fiber_cntr) = e1fa(e1fa_r_1, e1fa_c_1, e1fa_s_1, 4);

                    % get path direction as the first eigenvector
                    k1 = eigvector1_point1;                                         %initially set the step direction as E1
                    k1 = k1(e1_order);                                              %adjust for frame of reference and row/column indexing
                    k1 = [e1r_sign*k1(1) e1c_sign*k1(2) e1s_sign*k1(3)]';
                    k1(3) = k1(3)/depth_ratio;                                      %account for slice aspect ratio in the z direction.

                    %calculate next point, p2
                    point2 = point1 + (step_incr/2)*k1;

                    % get the direction at point p2
                    e1fa_r_1 = round(point2(1));
                    e1fa_c_1 = round(point2(2));
                    e1fa_s_1 = max([round(point2(3)) 1]);
                    eigvector1_point2 = squeeze(e1fa(e1fa_r_1, e1fa_c_1, e1fa_s_1, 1:3));
                    eigvector1_point2 = check_E1_sign(seed_method, eigvector1_point2, eigvector1_point1);         %correct, as needed, for antipodal symmetry

                    % get path direction as the first eigenvector
                    k2 = eigvector1_point2;                                         %initially set the step direction as E1
                    k2 = k2(e1_order);                                              %adjust for frame of reference and row/column indexing
                    k2 = [e1r_sign*k2(1) e1c_sign*k2(2) e1s_sign*k2(3)]';
                    k2(3) = k2(3)/depth_ratio;                                      %account for slice aspect ratio in the z direction.

                    %calculate next point, p3
                    point3 = point2 + (step_incr/2)*k2;

                    % get the direction at point p3
                    e1fa_r_1 = round(point3(1));
                    e1fa_c_1 = round(point3(2));
                    e1fa_s_1 = max([round(point3(3)) 1]);
                    eigvector1_point3 = squeeze(e1fa(e1fa_r_1, e1fa_c_1, e1fa_s_1, 1:3));
                    eigvector1_point3 = check_E1_sign(seed_method, eigvector1_point3, eigvector1_point2);         %correct, as needed, for antipodal symmetry

                    % get path direction as the first eigenvector
                    k3 = eigvector1_point3;                                         %initially set the step direction as E1
                    k3 = k3(e1_order);                                              %adjust for frame of reference and row/column indexing
                    k3 = [e1r_sign*k3(1) e1c_sign*k3(2) e1s_sign*k3(3)]';
                    k3(3) = k3(3)/depth_ratio;                                      %account for slice aspect ratio in the z direction.

                    %calculate next point, p4
                    point4 = point3 + step_incr*k3;

                    % get the direction at point p4
                    e1fa_r = round(point4(1));
                    e1fa_c = round(point4(2));
                    e1fa_s = max([round(point4(3)) 1]);
                    eigvector1_point4 = squeeze(e1fa(e1fa_r, e1fa_c, e1fa_s, 1:3));
                    eigvector1_point4 = check_E1_sign(seed_method, eigvector1_point4, eigvector1_point3);         %correct, as needed, for antipodal symmetry

                    %get the step direction from the average of the four eigenvectors.  First, find the mean dyadic tensor
                    dyadic_tensor_1 = eigvector1_point1*eigvector1_point1';
                    dyadic_tensor_2 = eigvector1_point2*eigvector1_point2';
                    dyadic_tensor_3 = eigvector1_point3*eigvector1_point3';
                    dyadic_tensor_4 = eigvector1_point4*eigvector1_point4';
                    dyadic_tensor_mean = (dyadic_tensor_1 + 2*dyadic_tensor_2 + dyadic_tensor_3*2 + dyadic_tensor_4)/6;
                    [dyadic_tensor_eigvector_m, dyadic_tensor_eigvalue_m] = eig(dyadic_tensor_mean);                     %diagonalize to get e-vectors and e-values

                    %get the derived indices
                    [~, sort_order] = sort(diag(dyadic_tensor_eigvalue_m), 'descend');
                    eigvector1_mean = dyadic_tensor_eigvector_m(:, sort_order(1));
                    eigvector1_1 = check_E1_sign(seed_method, eigvector1_mean, 1);                  %correct, as needed, for antipodal symmetry
                    eigvector1_1_old = eigvector1_1;

                    eigvector1_1_all(fiber_cntr, :) = eigvector1_1;
                    step_dir_1 = eigvector1_1;                                              %initially set the step direction as E1
                    step_dir_1 = step_dir_1(e1_order);                                      %adjust for frame of reference and row/column indexing
                    step_dir_1 = [e1r_sign*step_dir_1(1) e1c_sign*step_dir_1(2) e1s_sign*step_dir_1(3)]';
                    step_dir_1(3) = step_dir_1(3)/depth_ratio;                              %account for the slice aspect ratio

            end % of switch tracking mehtod for direction 1

            % add the point
            next_point_1 = squeeze(fiber_all_1(seed_cntr, plane_cntr, fiber_cntr,:)) + ...
                step_incr(1)*step_dir_1;                                              %multiply step_dir by the step size before adding the step to the preceding point
            fiber_cntr = fiber_cntr+1;                                              %increment the fiber counter

            % get indices into tensor matrix from the current point
            e1fa_r_1 = round(next_point_1(1));
            e1fa_c_1 = round(next_point_1(2));
            e1fa_s_1 = max([round(next_point_1(3)) 1]);                               %take teh maximum with slice 1, in case of low seed point

            % check mask; terminate tract if out of bounds. record as stop criterion = 4 and continue to next point. can use the
            % diffusion mask (i.e. the original, undilated muscle mask)
            if diffusion_mask(e1fa_r_1,e1fa_c_1,e1fa_s_1)==0
                stop_list_1(seed_cntr) = 4;
                continue;
            else
                fiber_all_1(seed_cntr, plane_cntr, fiber_cntr,:) = next_point_1;
            end

            % 4) begin the fiber tracking loop
            while 1

                % 4a) get the step direction

                switch prop_algo
                    case {'eu'}

                        % if within the diffusion mask, look up new E1 and FA. Otherwise, use previous values
                        %                         if diffusion_mask(e1fa_r_1, e1fa_c_1, e1fa_s_1)

                        eigvector1_1 = squeeze(e1fa(e1fa_r_1, e1fa_c_1, e1fa_s_1, 1:3));
                        eigvector1_1 = check_E1_sign(seed_method, eigvector1_1, eigvector1_1_old);                        %correct, as needed, for antipodal symmetry
                        eigvector1_1_all(fiber_cntr,:) = eigvector1_1;

                        % store the FA data for termination checking
                        fa_all_1(seed_cntr, plane_cntr,fiber_cntr) = e1fa(e1fa_r_1, e1fa_c_1, e1fa_s_1, 4);
                        %
                        %                         else
                        %
                        %                             eigvector1_1 = eigvector1_1_old;
                        %                             eigvector1_1_all(fiber_cntr,:) = eigvector1_1;
                        %                             fa_all_1(seed_cntr, plane_cntr, fiber_cntr) = fa_all_1(seed_cntr, plane_cntr, fiber_cntr-1);
                        % %
                        %                         end

                        % define the step
                        step_dir_1 = eigvector1_1;                                              %initially set the step direction as E1
                        step_dir_1 = step_dir_1(e1_order);                                      %adjust for frame of reference and row/column indexing
                        step_dir_1 = [e1r_sign*step_dir_1(1) e1c_sign*step_dir_1(2) e1s_sign*step_dir_1(3)]';
                        step_dir_1(3) = step_dir_1(3)/depth_ratio;                              %account for the slice aspect ratio

                    case {'rk'}

                        %create a variable to keep track of the RK order; use order = 1 (Euler's method) unless all four points are inside the mask
                        rk_order=1;

                        % define point p1; if it lies within the mask, get the diffusion tensor by calling the retrieve_tensor function.
                        point1 = squeeze(fiber_all_1(seed_cntr, plane_cntr,fiber_cntr,:));

                        if round(point1(1))>length(e1fa(:,1,1,1)) || round(point1(2))>length(e1fa(1,:,1,1)) || ...
                                round(point1(3))>length(e1fa(1,1,:,1))
                            stop_list_1(seed_cntr, plane_cntr) = 4;
                            break;
                        end

                        if diffusion_mask(round(point1(1)), round(point1(2)), max([round(point1(3)) 1]))==1

                            % get the eigenvector at point 1
                            e1fa_r_1 = round(point1(1));
                            e1fa_c_1 = round(point1(2));
                            e1fa_s_1 = max([round(point1(3)) 1]);
                            eigvector1_point1 = squeeze(e1fa(e1fa_r_1, e1fa_c_1, e1fa_s_1, 1:3));
                            eigvector1_point1 = check_E1_sign(seed_method, eigvector1_point1, eigvector1_1_old);         %correct, as needed, for antipodal symmetry

                            % store the FA data for termination checking
                            fa_all_1(seed_cntr, plane_cntr,fiber_cntr) = e1fa(e1fa_r, e1fa_c, e1fa_s, 4);

                            % get path direction as the first eigenvector
                            k1 = eigvector1_point1;                                         %initially set the step direction as E1
                            k1 = k1(e1_order);                                              %adjust for frame of reference and row/column indexing
                            k1 = [e1r_sign*k1(1) e1c_sign*k1(2) e1s_sign*k1(3)]';
                            k1(3) = k1(3)/depth_ratio;                                      %account for slice aspect ratio in the z direction.

                            %calculate next point, p2
                            point2 = point1 + (step_incr/2)*k1;

                        else
                            stop_list_1(seed_cntr, plane_cntr)=4;
                            break

                        end %of if point 1 is valid if statement

                        % if mask at point 2 is valid, get the diffusion tensor there and calculate point three
                        if diffusion_mask(round(point2(1)), round(point2(2)), max([round(point2(3)) 1]))==1

                            % get teh eigenvector at point 1
                            e1fa_r_1 = round(point2(1));
                            e1fa_c_1 = round(point2(2));
                            e1fa_s_1 = max([round(point2(3)) 1]);
                            eigvector1_point2 = squeeze(e1fa(e1fa_r_1, e1fa_c_1, e1fa_s_1, 1:3));
                            eigvector1_point2 = check_E1_sign(seed_method, eigvector1_point2, eigvector1_point1);         %correct, as needed, for antipodal symmetry

                            % get the step
                            k2 = eigvector1_point2;                                         %initially set the step direction as E1
                            k2 = k2(e1_order);                                              %adjust for frame of reference and row/column indexing
                            k2 = [e1r_sign*k2(1) e1c_sign*k2(2) e1s_sign*k2(3)]';
                            k2(3) = k2(3)/depth_ratio;                                      %account for slice aspect ratio in the z direction.

                            %calculate next point, p3
                            point3 = point2 + (step_incr/2)*k2;

                            % if mask at point 3 is valid, get the diffusion tensor there and calculate point 4
                            if diffusion_mask(round(point3(1)), round(point3(2)), max([round(point3(3)) 1]))==1

                                % get the eigenvector
                                e1fa_r_1 = round(point3(1));
                                e1fa_c_1 = round(point3(2));
                                e1fa_s_1 = max([round(point3(3)) 1]);
                                eigvector1_point3 = squeeze(e1fa(e1fa_r_1, e1fa_c_1, e1fa_s_1, 1:3));
                                eigvector1_point3 = check_E1_sign(seed_method, eigvector1_point3, eigvector1_point2);         %correct, as needed, for antipodal symmetry

                                % get path direction as the first eigenvector
                                k3 = eigvector1_point3;                                         %initially set the step direction as E1
                                k3 = k3(e1_order);                                              %adjust for frame of reference and row/column indexing
                                k3 = [e1r_sign*k3(1) e1c_sign*k3(2) e1s_sign*k3(3)]';
                                k3(3) = k3(3)/depth_ratio;                                      %account for slice aspect ratio in the z direction.

                                %calculate next point, p4
                                point4 = point3 + step_incr*k3;

                                % if mask at point 4 is valid, get the diffusion tensor there
                                if diffusion_mask(round(point4(1)), round(point4(2)), max([round(point4(3)) 1]))==1
                                    rk_order=4;

                                    % get the eigenvector
                                    e1fa_r_1 = round(point4(1));
                                    e1fa_c_1 = round(point4(2));
                                    e1fa_s_1 = max([round(point4(3)) 1]);
                                    eigvector1_point4 = squeeze(e1fa(e1fa_r_1, e1fa_c_1, e1fa_s_1, 1:3));
                                    eigvector1_point4 = check_E1_sign(seed_method, eigvector1_point4, eigvector1_point3);         %correct, as needed, for antipodal symmetry

                                end                                                                 %of point 4 is valid if statement

                            end                                                                 %of point 3 is valid if statement

                        end                                                                 %of point 2 is valid if statement

                        if rk_order==4

                            %get the step direction from the average of the four eigenvectors.  First, find the mean dyadic tensor
                            dyadic_tensor_1 = eigvector1_point1*eigvector1_point1';
                            dyadic_tensor_2 = eigvector1_point2*eigvector1_point2';
                            dyadic_tensor_3 = eigvector1_point3*eigvector1_point3';
                            dyadic_tensor_4 = eigvector1_point4*eigvector1_point4';
                            dyadic_tensor_mean = (dyadic_tensor_1 + 2*dyadic_tensor_2 + 2*dyadic_tensor_3 + dyadic_tensor_4)/6;
                            [dyadic_tensor_eigvector_m, dyadic_tensor_eigvalue_m] = eig(dyadic_tensor_mean);                     %diagonalize to get e-vectors and e-values

                            %get the derived indices
                            [~, sort_order] = sort(diag(dyadic_tensor_eigvalue_m), 'descend');
                            eigvector1_mean = dyadic_tensor_eigvector_m(:, sort_order(1));
                            eigvector1_1 = check_E1_sign(seed_method, eigvector1_mean, eigvector1_1_old);                  %correct, as needed, for antipodal symmetry
                            eigvector1_1_old = eigvector1_mean;

                            eigvector1_1_all(fiber_cntr,:) = eigvector1_1;
                            step_dir_1 = eigvector1_1;                                              %initially set the step direction as E1
                            step_dir_1 = step_dir_1(e1_order);                                      %adjust for frame of reference and row/column indexing
                            step_dir_1 = [e1r_sign*step_dir_1(1) e1c_sign*step_dir_1(2) e1s_sign*step_dir_1(3)]';
                            step_dir_1(3) = step_dir_1(3)/depth_ratio;                              %account for the slice aspect ratio

                        elseif rk_order==1

                            eigvector1_1 = eigvector1_point1;
                            eigvector1_1_all(fiber_cntr,:) = eigvector1_1;
                            step_dir_1 = eigvector1_1;                                              %initially set the step direction as E1
                            step_dir_1 = step_dir_1(e1_order);                                      %adjust for frame of reference and row/column indexing
                            step_dir_1 = [e1r_sign*step_dir_1(1) e1c_sign*step_dir_1(2) e1s_sign*step_dir_1(3)]';
                            step_dir_1(3) = step_dir_1(3)/depth_ratio;                              %account for the slice aspect ratio

                        end


                end %of propagation algorithm switch statement


                % 4b) decide whether or not to add the point
                % get indices into tensor matrix from the current point
                e1fa_r_1 = round(next_point_1(1));
                e1fa_c_1 = round(next_point_1(2));
                e1fa_s_1 = max([round(next_point_1(3)) 1]);
                e1fa_s_1 = min([length(e1fa(1,1,:,1,1)) e1fa_s_1]);

                % first criterion - check mask; terminate tract if out of bounds and record as stop criterion = 4. intentionally uses the diffusion mask
                if diffusion_mask(e1fa_r_1,e1fa_c_1,e1fa_s_1)==0
                    stop_list_1(seed_cntr, plane_cntr) = 4;
                    break;
                end

                %other criteria depend on termination algorithm selected
                switch term_mthd                                                %select between tract termination methods
                    case {'bin1'}                                                %binary OR score, 1 point required for termination
                        %check for out of bounds FA values
                        if fa_all_1(seed_cntr, plane_cntr, fiber_cntr)<fa_low || fa_all_1(seed_cntr, plane_cntr, fiber_cntr)>fa_high
                            stop_list_1(seed_cntr, plane_cntr) = 2;
                            fiber_len_1(seed_cntr, plane_cntr) = fiber_cntr;
                            break;
                        end

                        %check for out of bounds curvature values
                        if fiber_cntr > angle_lookback
                            eigvector1_1_prev = squeeze(eigvector1_1_all((fiber_cntr-angle_lookback),:));
                            step_angle_1 = abs(acosd(dot(eigvector1_1, eigvector1_1_prev)));
                            angle_1_all(fiber_cntr) = step_angle_1;
                            if step_angle_1>angle_thrsh
                                stop_list_1(seed_cntr, plane_cntr) = 3;
                                fiber_len_1(seed_cntr, plane_cntr) = fiber_cntr;
                                break;
                            end
                        end

                    case {'bin2'}                                                %binary OR score, 2 consecutive points required

                        %check for out of bounds FA values
                        if fa_all_1(seed_cntr, plane_cntr, fiber_cntr)<fa_low || fa_all_1(seed_cntr, plane_cntr, fiber_cntr)>fa_high
                            if fa_all_1(seed_cntr, plane_cntr, fiber_cntr-1)<fa_low || fa_all_1(seed_cntr, plane_cntr, fiber_cntr-1)>fa_high
                                stop_list_1(seed_cntr, plane_cntr) = 2;
                                fiber_len_1(seed_cntr, plane_cntr) = fiber_cntr;
                                break;
                            end
                        end

                        %check for out of bounds curvature values
                        if fiber_cntr > (angle_lookback+1)
                            eigvector1_1_prev = squeeze(eigvector1_1_all((fiber_cntr-angle_lookback),:));
                            step_angle_old = abs(acosd(dot(eigvector1_1, eigvector1_1_prev)));

                            eigvector1_last = squeeze(eigvector1_1_all((fiber_cntr-1),:));
                            step_angle_current = abs(acosd(dot(eigvector1_1, eigvector1_last)));
                            angle_1_all(fiber_cntr) = step_angle_current;

                            if step_angle_current>angle_thrsh && step_angle_old>angle_thrsh
                                stop_list_1(seed_cntr, plane_cntr) = 3;
                                fiber_len_1(seed_cntr, plane_cntr) = fiber_cntr;
                                break;
                            end
                        end

                end                                                             %end of termination method switch statement

                % 5) add the point and continue tracking
                fiber_cntr = fiber_cntr + 1;                                     %increment the fiber counter
                fiber_len_1(seed_cntr, plane_cntr) = fiber_cntr;
                next_point_1 = squeeze(fiber_all_1(seed_cntr, plane_cntr,fiber_cntr-1,:)) + ...
                    step_incr*step_dir_1;                                         %multiply step_dir by the step increment before adding the step to the preceding point
                fiber_all_1(seed_cntr, plane_cntr,fiber_cntr,:) = next_point_1;        %add the point
                fa_all_1(seed_cntr, plane_cntr, fiber_cntr) = e1fa(e1fa_r_1, e1fa_c_1, e1fa_s_1);

                roi_flag_1(seed_cntr, plane_cntr)=1;                                 %record as a tracked fiber

            end                                                                 %end of fiber tracking while loop

            %%%%%%%%%%%%%%%%%%%%%% Track in direction 2  %%%%%%%%%%%%%%%%%%%%

            % initialize some variables
            fiber_cntr = 1;                                                     %counter for the fiber tract points; initialize at 1 for each tract
            angle_2_all = zeros(50,1);                                           %keep track of all inter-step angles
            eigvector1_2_all = zeros(50,3);                                      %matrix to hold all of the first eigenvectors (needed for checking inter-point angle stop criterion)

            % 1) locate the seed point on the seed_point_2 matrix
            seed_point_2(1:3) = [seed_points(seed_cntr, plane_cntr, 1)...
                seed_points(seed_cntr, plane_cntr, 2)...
                seed_points(seed_cntr, plane_cntr, 3)];
            e1fa_r_2 = round(seed_point_2(1));                                    %round the coordinates to get a row, column, and slice location in the tensor matrix
            e1fa_c_2 = round(seed_point_2(2));
            e1fa_s_2 = max([round(seed_point_2(3)) 1]);                      %just in case the z coordinate rounds to zero, make it a one

            % 2) add to the fiber_all matrix
            if diffusion_mask(e1fa_r_2,e1fa_c_2,e1fa_s_2)==0
                stop_list_2(seed_cntr,plane_cntr)=4;
                continue;
            end

            fiber_all_2(seed_cntr, plane_cntr, fiber_cntr,:) = seed_point_2;

            % 3) get the initial direction (NB initial direction methods not updated for v. 2.0 because the seed points
            %    are, by definition, within the mask

            switch prop_algo
                case {'eu'}

                    % get the eigenvector
                    eigvector1_2 = squeeze(e1fa(e1fa_r_2, e1fa_c_2, e1fa_s_2, 1:3));
                    eigvector1_2 = check_E1_sign(seed_method, eigvector1_2, -1);                        %correct, as needed, for antipodal symmetry. note change of sign from direction 1
                    eigvector1_2_all(fiber_cntr,:) = eigvector1_2;
                    eigvector1_2_old = eigvector1_2;

                    % store the FA data for termination checking
                    fa_all_2(seed_cntr, plane_cntr, fiber_cntr) = e1fa(e1fa_r_2, e1fa_c_2, e1fa_s_2, 4);

                    step_dir_2 = eigvector1_2;                                              %initially set the step direction as E1
                    step_dir_2 = step_dir_2(e1_order);                                      %adjust for frame of reference and row/column indexing
                    step_dir_2 = [e1r_sign*step_dir_2(1) e1c_sign*step_dir_2(2) e1s_sign*step_dir_2(3)]';
                    step_dir_2(3) = step_dir_2(3)/depth_ratio;                              %account for the slice aspect ratio

                case {'rk'}

                    % for point p1, get the direction
                    point1 = squeeze(fiber_all_2(seed_cntr, plane_cntr, fiber_cntr,:));

                    e1fa_r_2 = round(point1(1));
                    e1fa_c_2 = round(point1(2));
                    e1fa_s_2 = max([round(point1(3)) 1]);
                    eigvector1_point1 = squeeze(e1fa(e1fa_r_2, e1fa_c_2, e1fa_s_2, 1:3));
                    eigvector1_point1 = check_E1_sign(seed_method, eigvector1_point1, -1);         %correct, as needed, for antipodal symmetry

                    % store the FA data for termination checking
                    fa_all_2(seed_cntr, plane_cntr,fiber_cntr) = e1fa(e1fa_r_2, e1fa_c_2, e1fa_s_2, 4);

                    % get path direction as the first eigenvector
                    k1 = eigvector1_point1;                                         %initially set the step direction as E1
                    k1 = k1(e1_order);                                              %adjust for frame of reference and row/column indexing
                    k1 = [e1r_sign*k1(1) e1c_sign*k1(2) e1s_sign*k1(3)]';
                    k1(3) = k1(3)/depth_ratio;                                      %account for slice aspect ratio in the z direction.

                    %calculate next point, p2
                    point2 = point1 + (step_incr/2)*k1;

                    % get the direction at point p2
                    e1fa_r_2 = round(point2(1));
                    e1fa_c_2 = round(point2(2));
                    e1fa_s_2 = max([round(point2(3)) 1]);
                    eigvector1_point2 = squeeze(e1fa(e1fa_r_2, e1fa_c_2, e1fa_s_2, 1:3));
                    eigvector1_point2 = check_E1_sign(seed_method, eigvector1_point2, eigvector1_point1);         %correct, as needed, for antipodal symmetry

                    % get path direction as the first eigenvector
                    k2 = eigvector1_point2;                                         %initially set the step direction as E1
                    k2 = k2(e1_order);                                              %adjust for frame of reference and row/column indexing
                    k2 = [e1r_sign*k2(1) e1c_sign*k2(2) e1s_sign*k2(3)]';
                    k2(3) = k2(3)/depth_ratio;                                      %account for slice aspect ratio in the z direction.

                    %calculate next point, p3
                    point3 = point2 + (step_incr/2)*k2;

                    % get the direction at point p3
                    e1fa_r_2 = round(point3(1));
                    e1fa_c_2 = round(point3(2));
                    e1fa_s_2 = max([round(point3(3)) 1]);
                    eigvector1_point3 = squeeze(e1fa(e1fa_r_2, e1fa_c_2, e1fa_s_2, 1:3));
                    eigvector1_point3 = check_E1_sign(seed_method, eigvector1_point3, eigvector1_point2);         %correct, as needed, for antipodal symmetry

                    % get path direction as the first eigenvector
                    k3 = eigvector1_point3;                                         %initially set the step direction as E1
                    k3 = k3(e1_order);                                              %adjust for frame of reference and row/column indexing
                    k3 = [e1r_sign*k3(1) e1c_sign*k3(2) e1s_sign*k3(3)]';
                    k3(3) = k3(3)/depth_ratio;                                      %account for slice aspect ratio in the z direction.

                    %calculate next point, p4
                    point4 = point3 + step_incr*k3;

                    % get the direction at point p4
                    e1fa_r = round(point4(1));
                    e1fa_c = round(point4(2));
                    e1fa_s = max([round(point4(3)) 1]);
                    eigvector1_point4 = squeeze(e1fa(e1fa_r, e1fa_c, e1fa_s, 1:3));
                    eigvector1_point4 = check_E1_sign(seed_method, eigvector1_point4, eigvector1_point3);         %correct, as needed, for antipodal symmetry

                    %get the step direction from the average of the four eigenvectors.  First, find the mean dyadic tensor
                    dyadic_tensor_1 = eigvector1_point1*eigvector1_point1';
                    dyadic_tensor_2 = eigvector1_point2*eigvector1_point2';
                    dyadic_tensor_3 = eigvector1_point3*eigvector1_point3';
                    dyadic_tensor_4 = eigvector1_point4*eigvector1_point4';
                    dyadic_tensor_mean = (dyadic_tensor_1 + 2*dyadic_tensor_2 + dyadic_tensor_3*2 + dyadic_tensor_4)/6;
                    [dyadic_tensor_eigvector_m, dyadic_tensor_eigvalue_m] = eig(dyadic_tensor_mean);                     %diagonalize to get e-vectors and e-values

                    %get the derived indices
                    [~, sort_order] = sort(diag(dyadic_tensor_eigvalue_m), 'descend');
                    eigvector1_mean = dyadic_tensor_eigvector_m(:, sort_order(1));
                    eigvector1_2 = check_E1_sign(seed_method, eigvector1_mean, -1);                  %correct, as needed, for antipodal symmetry
                    eigvector1_2_old = eigvector1_2;

                    eigvector1_2_all(fiber_cntr, :) = eigvector1_2;
                    step_dir_2 = eigvector1_2;                                              %initially set the step direction as E1
                    step_dir_2 = step_dir_2(e1_order);                                      %adjust for frame of reference and row/column indexing
                    step_dir_2 = [e1r_sign*step_dir_2(1) e1c_sign*step_dir_2(2) e1s_sign*step_dir_2(3)]';
                    step_dir_2(3) = step_dir_2(3)/depth_ratio;                              %account for the slice aspect ratio

            end % of switch tracking mehtod for direction 2, initial point

            % add the point
            next_point_2 = squeeze(fiber_all_2(seed_cntr, plane_cntr, fiber_cntr,:)) + ...
                step_incr(1)*step_dir_2;                                              %multiply step_dir by the step size before adding the step to the preceding point
            fiber_cntr = fiber_cntr+1;                                              %increment the fiber counter

            % get indices into tensor matrix from the current point
            e1fa_r_2 = round(next_point_2(1));
            e1fa_c_2 = round(next_point_2(2));
            e1fa_s_2 = max([round(next_point_2(3)) 1]);                               %take teh maximum with slice 1, in case of low seed point

            % check mask; terminate tract if out of bounds. record as stop criterion = 4 and continue to next point
            if diffusion_mask(e1fa_r_2,e1fa_c_2,e1fa_s_2)==0
                stop_list_2(seed_cntr) = 4;
                continue;
            else
                fiber_all_2(seed_cntr, plane_cntr, fiber_cntr,:) = next_point_2;
            end

            % 4) begin the fiber tracking loop
            while 1

                % 4a) get the step direction

                switch prop_algo
                    case {'eu'}

                        % if within the diffusion mask, look up new E1 and FA. Otherwise, use previous values
                        %                         if diffusion_mask(e1fa_r_2, e1fa_c_2, e1fa_s_2)

                        eigvector1_2 = squeeze(e1fa(e1fa_r_2, e1fa_c_2, e1fa_s_2, 1:3));
                        eigvector1_2 = check_E1_sign(seed_method, eigvector1_2, eigvector1_2_old);                        %correct, as needed, for antipodal symmetry
                        eigvector1_2_all(fiber_cntr,:) = eigvector1_2;

                        % store the FA data for termination checking
                        fa_all_2(seed_cntr, plane_cntr,fiber_cntr) = e1fa(e1fa_r_2, e1fa_c_2, e1fa_s_2, 4);

                        %                         else
                        %
                        %                             eigvector1_2 = eigvector1_2_old;
                        %                             eigvector1_2_all(fiber_cntr,:) = eigvector1_2;
                        %                             fa_all_2(seed_cntr, plane_cntr, fiber_cntr) = fa_all_2(seed_cntr, plane_cntr, fiber_cntr-1);
                        %
                        %                         end

                        % define the step
                        step_dir_2 = eigvector1_2;                                              %initially set the step direction as E1
                        step_dir_2 = step_dir_2(e1_order);                                      %adjust for frame of reference and row/column indexing
                        step_dir_2 = [e1r_sign*step_dir_2(1) e1c_sign*step_dir_2(2) e1s_sign*step_dir_2(3)]';
                        step_dir_2(3) = step_dir_2(3)/depth_ratio;                              %account for the slice aspect ratio

                    case {'rk'}

                        %create a variable to keep track of the RK order; use order = 1 (Euler's method) unless all four points are inside the mask
                        rk_order=1;

                        % define point p1; if it lies within the mask, get the diffusion tensor by calling the retrieve_tensor function.
                        point1 = squeeze(fiber_all_2(seed_cntr, plane_cntr,fiber_cntr,:));

                        if round(point1(1))>length(e1fa(:,1,1,1)) || round(point1(2))>length(e1fa(1,:,1,1)) || ...
                                round(point1(3))>length(e1fa(1,1,:,1))
                            stop_list_2(seed_cntr, plane_cntr) = 4;
                            break;
                        end

                        if diffusion_mask(round(point1(1)), round(point1(2)), max([round(point1(3)) 1]))==1

                            % get the eigenvector at point 1
                            e1fa_r_2 = round(point1(1));
                            e1fa_c_2 = round(point1(2));
                            e1fa_s_2 = max([round(point1(3)) 1]);
                            eigvector1_point1 = squeeze(e1fa(e1fa_r_2, e1fa_c_2, e1fa_s_2, 1:3));
                            eigvector1_point1 = check_E1_sign(seed_method, eigvector1_point1, eigvector1_2_old);         %correct, as needed, for antipodal symmetry

                            % store the FA data for termination checking
                            fa_all_2(seed_cntr, plane_cntr,fiber_cntr) = e1fa(e1fa_r_2, e1fa_c_2, e1fa_s_2, 4);

                            % get path direction as the first eigenvector
                            k1 = eigvector1_point1;                                         %initially set the step direction as E1
                            k1 = k1(e1_order);                                              %adjust for frame of reference and row/column indexing
                            k1 = [e1r_sign*k1(1) e1c_sign*k1(2) e1s_sign*k1(3)]';
                            k1(3) = k1(3)/depth_ratio;                                      %account for slice aspect ratio in the z direction.

                            %calculate next point, p2
                            point2 = point1 + (step_incr/2)*k1;

                        else
                            stop_list_2(seed_cntr, plane_cntr)=4;
                            break

                        end %of if point 1 is valid if statement

                        % if mask at point 2 is valid, get the diffusion tensor there and calculate point three
                        if diffusion_mask(round(point2(1)), round(point2(2)), max([round(point2(3)) 1]))==1

                            % get teh eigenvector at point 1
                            e1fa_r_2 = round(point2(1));
                            e1fa_c_2 = round(point2(2));
                            e1fa_s_2 = max([round(point2(3)) 1]);
                            eigvector1_point2 = squeeze(e1fa(e1fa_r_2, e1fa_c_2, e1fa_s_2, 1:3));
                            eigvector1_point2 = check_E1_sign(seed_method, eigvector1_point2, eigvector1_point1);         %correct, as needed, for antipodal symmetry

                            % get the step
                            k2 = eigvector1_point2;                                         %initially set the step direction as E1
                            k2 = k2(e1_order);                                              %adjust for frame of reference and row/column indexing
                            k2 = [e1r_sign*k2(1) e1c_sign*k2(2) e1s_sign*k2(3)]';
                            k2(3) = k2(3)/depth_ratio;                                      %account for slice aspect ratio in the z direction.

                            %calculate next point, p3
                            point3 = point2 + (step_incr/2)*k2;

                            % if mask at point 3 is valid, get the diffusion tensor there and calculate point 4
                            if diffusion_mask(round(point3(1)), round(point3(2)), max([round(point3(3)) 1]))==1

                                % get the eigenvector
                                e1fa_r_2 = round(point3(1));
                                e1fa_c_2 = round(point3(2));
                                e1fa_s_2 = max([round(point3(3)) 1]);
                                eigvector1_point3 = squeeze(e1fa(e1fa_r_2, e1fa_c_2, e1fa_s_2, 1:3));
                                eigvector1_point3 = check_E1_sign(seed_method, eigvector1_point3, eigvector1_point2);         %correct, as needed, for antipodal symmetry

                                % get path direction as the first eigenvector
                                k3 = eigvector1_point3;                                         %initially set the step direction as E1
                                k3 = k3(e1_order);                                              %adjust for frame of reference and row/column indexing
                                k3 = [e1r_sign*k3(1) e1c_sign*k3(2) e1s_sign*k3(3)]';
                                k3(3) = k3(3)/depth_ratio;                                      %account for slice aspect ratio in the z direction.

                                %calculate next point, p4
                                point4 = point3 + step_incr*k3;

                                % if mask at point 4 is valid, get the diffusion tensor there
                                if diffusion_mask(round(point4(1)), round(point4(2)), max([round(point4(3)) 1]))==1
                                    rk_order=4;

                                    % get the eigenvector
                                    e1fa_r_2 = round(point4(1));
                                    e1fa_c_2 = round(point4(2));
                                    e1fa_s_2 = max([round(point4(3)) 1]);
                                    eigvector1_point4 = squeeze(e1fa(e1fa_r_2, e1fa_c_2, e1fa_s_2, 1:3));
                                    eigvector1_point4 = check_E1_sign(seed_method, eigvector1_point4, eigvector1_point3);         %correct, as needed, for antipodal symmetry

                                end                                                                 %of point 4 is valid if statement

                            end                                                                 %of point 3 is valid if statement

                        end                                                                 %of point 2 is valid if statement

                        if rk_order==4

                            %get the step direction from the average of the four eigenvectors.  First, find the mean dyadic tensor
                            dyadic_tensor_1 = eigvector1_point1*eigvector1_point1';
                            dyadic_tensor_2 = eigvector1_point2*eigvector1_point2';
                            dyadic_tensor_3 = eigvector1_point3*eigvector1_point3';
                            dyadic_tensor_4 = eigvector1_point4*eigvector1_point4';
                            dyadic_tensor_mean = (dyadic_tensor_1 + 2*dyadic_tensor_2 + 2*dyadic_tensor_3 + dyadic_tensor_4)/6;
                            [dyadic_tensor_eigvector_m, dyadic_tensor_eigvalue_m] = eig(dyadic_tensor_mean);                     %diagonalize to get e-vectors and e-values

                            %get the derived indices
                            [~, sort_order] = sort(diag(dyadic_tensor_eigvalue_m), 'descend');
                            eigvector1_mean = dyadic_tensor_eigvector_m(:, sort_order(1));
                            eigvector1_2 = check_E1_sign(seed_method, eigvector1_mean, eigvector1_2_old);         %correct, as needed, for antipodal symmetry
                            eigvector1_2_old = eigvector1_2;

                            eigvector1_2_all(fiber_cntr,:) = eigvector1_2;
                            step_dir_2 = eigvector1_2;                                              %initially set the step direction as E1
                            step_dir_2 = step_dir_2(e1_order);                                      %adjust for frame of reference and row/column indexing
                            step_dir_2 = [e1r_sign*step_dir_2(1) e1c_sign*step_dir_2(2) e1s_sign*step_dir_2(3)]';
                            step_dir_2(3) = step_dir_2(3)/depth_ratio;                              %account for the slice aspect ratio

                        elseif rk_order==1

                            eigvector1_2 = eigvector1_point1;
                            eigvector1_2_all(fiber_cntr,:) = eigvector1_2;
                            step_dir_2 = eigvector1_2;                                              %initially set the step direction as E1
                            step_dir_2 = step_dir_2(e1_order);                                      %adjust for frame of reference and row/column indexing
                            step_dir_2 = [e1r_sign*step_dir_2(1) e1c_sign*step_dir_2(2) e1s_sign*step_dir_2(3)]';
                            step_dir_2(3) = step_dir_2(3)/depth_ratio;                              %account for the slice aspect ratio

                        end


                end %of propagation algorithm switch statement


                % 4b) decide whether or not to add the point
                % get indices into tensor matrix from the current point
                e1fa_r_2 = round(next_point_2(1));
                e1fa_c_2 = round(next_point_2(2));
                e1fa_s_2 = max([round(next_point_2(3)) 1]);
                e1fa_s_2 = min([length(e1fa(1,1,:,1,1)) e1fa_s_2]);

                % first criterion - check mask; terminate tract if out of bounds and record as stop criterion = 4
                if diffusion_mask(e1fa_r_2,e1fa_c_2,e1fa_s_2)==0
                    stop_list_2(seed_cntr, plane_cntr) = 4;
                    break;
                end

                %other criteria depend on termination algorithm selected
                switch term_mthd                                                %select between tract termination methods
                    case {'bin1'}                                                %binary OR score, 1 point required for termination
                        %check for out of bounds FA values
                        if fa_all_2(seed_cntr, plane_cntr, fiber_cntr)<fa_low || fa_all_2(seed_cntr, plane_cntr, fiber_cntr)>fa_high
                            stop_list_2(seed_cntr, plane_cntr) = 2;
                            fiber_len_2(seed_cntr, plane_cntr) = fiber_cntr;
                            break;
                        end

                        %check for out of bounds curvature values
                        if fiber_cntr > angle_lookback
                            eigvector1_2_prev = squeeze(eigvector1_2_all((fiber_cntr-angle_lookback),:));
                            step_angle_2 = abs(acosd(dot(eigvector1_2, eigvector1_2_prev)));
                            angle_2_all(fiber_cntr) = step_angle_2;
                            if step_angle_2>angle_thrsh
                                stop_list_2(seed_cntr, plane_cntr) = 3;
                                fiber_len_2(seed_cntr, plane_cntr) = fiber_cntr;
                                break;
                            end
                        end

                    case {'bin2'}                                                %binary OR score, 2 consecutive points required

                        %check for out of bounds FA values
                        if fa_all_2(seed_cntr, plane_cntr, fiber_cntr)<fa_low || fa_all_2(seed_cntr, plane_cntr, fiber_cntr)>fa_high
                            if fa_all_2(seed_cntr, plane_cntr, fiber_cntr-1)<fa_low || fa_all_2(seed_cntr, plane_cntr, fiber_cntr-1)>fa_high
                                stop_list_2(seed_cntr, plane_cntr) = 2;
                                fiber_len_2(seed_cntr, plane_cntr) = fiber_cntr;
                                break;
                            end
                        end

                        %check for out of bounds curvature values
                        if fiber_cntr > (angle_lookback+1)
                            eigvector1_2_prev = squeeze(eigvector1_2_all((fiber_cntr-angle_lookback),:));
                            step_angle_old = abs(acosd(dot(eigvector1_2, eigvector1_2_prev)));

                            eigvector1_last = squeeze(eigvector1_2_all((fiber_cntr-1),:));
                            step_angle_current = abs(acosd(dot(eigvector1_2, eigvector1_last)));
                            angle_2_all(fiber_cntr) = step_angle_current;

                            if step_angle_current>angle_thrsh && step_angle_old>angle_thrsh
                                stop_list_2(seed_cntr, plane_cntr) = 3;
                                fiber_len_2(seed_cntr, plane_cntr) = fiber_cntr;
                                break;
                            end
                        end

                end                                                             %end of termination method switch statement

                % 5) add the point and continue tracking
                fiber_cntr = fiber_cntr + 1;                                     %increment the fiber counter
                fiber_len_2(seed_cntr, plane_cntr) = fiber_cntr;
                next_point_2 = squeeze(fiber_all_2(seed_cntr, plane_cntr,fiber_cntr-1,:)) + ...
                    step_incr*step_dir_2;                                         %multiply step_dir by the step increment before adding the step to the preceding point
                fiber_all_2(seed_cntr, plane_cntr,fiber_cntr,:) = next_point_2;        %add the point
                fa_all_2(seed_cntr, plane_cntr, fiber_cntr) = e1fa(e1fa_r_2, e1fa_c_2, e1fa_s_2);

                roi_flag_2(seed_cntr, plane_cntr)=1;                                 %record as a tracked fiber

            end                                                                 %end of fiber tracking while loop



        end %of seed_points seed for loop

    end %of seed_points plane for loop

    %%%%%%%%%%%%%%%%%%%%% merge datasets %%%%%%%%%%%%%%%%%%%%%

    % initialize variables
    max_fiber_len_1 = max(fiber_len_1);
    max_fiber_len_2 = max(fiber_len_2);
    stop_list = cat(3, stop_list_1, stop_list_2);
    roi_flag = cat(3, roi_flag_1, roi_flag_2);
    fiber_len = cat(3, fiber_len_1, fiber_len_2);

    % assign data
    fiber_all = flip(fiber_all_2(:,:,1:max_fiber_len_2,:), 3);
    fiber_all(:,:,max_fiber_len_2:(max_fiber_len_2+max_fiber_len_1-1),:) = ...
        fiber_all_1(:,:,1:max_fiber_len_1,:);

end % of "other than aponeurosis" tracking methods

%% clean up fiber tracts to account for rounding errors that resulted in inappropriately prolonged tracking loops

% added for v. 2.0

% now work from teh other end of the fiber for non-apo seeding methods
if seed_method(1) ~='a'
% 
%     for row_cntr=1:length(fiber_all(:,1,1,1))
% 
%         for col_cntr=1:length(fiber_all(1,:,1,1))
% 
%             loop_fiber = squeeze(fiber_all(row_cntr, col_cntr, :, :));
%             point1 = find(loop_fiber(:,1), 1);
%             point2 = find(loop_fiber(:,1), 1, 'last');
%             num_points = point2 - point1 + 1;
%             loop_fiber = loop_fiber(point1:point2, :);
% 
%             % examine outer points
%             if num_points>6
% 
%                 % work back from teh end of the fiber tract
%                 end_point = (length(nonzeros(loop_fiber(:,1)))-2);
%                 for fiber_cntr=length(nonzeros(loop_fiber(:,1))):-1:end_point
% 
%                     loop_point = loop_fiber(fiber_cntr, :);
% 
%                     % calculate floor and ceiling values for each point
%                     loop_r = round(loop_point(1));
%                     loop_c = round(loop_point(2));
%                     loop_s = max([1 round(loop_point(3))]);
%                     loop_s = min([loop_s find(sum(sum(diffusion_mask)), 1, 'last')]);
% 
%                     if diffusion_mask(loop_r, loop_c, loop_s)==0
%                         loop_fiber(fiber_cntr, :) = 0;
%                     end
% 
%                     if fiber_cntr==end_point
% 
%                         fiber_all(row_cntr, col_cntr, :,:) = 0;
% 
%                         if sum(sum(loop_fiber)>0)
%                             loop_fiber = [nonzeros(loop_fiber(:,1)) nonzeros(loop_fiber(:,2)) nonzeros(loop_fiber(:,3))];
%                             length_loop_fiber = length(loop_fiber(:,1));
%                             fiber_all(row_cntr, col_cntr, 1:length_loop_fiber, :) = loop_fiber(1:length_loop_fiber, :);
% 
%                         end
% 
%                     end
% 
%                 end % of fiber_cntr loop
% 
%             end % of num_point>6
% 
%         end % of column loop
% 
%     end % of row loop

    % if there is an internal mask, trim points out 
    if sum(sum(sum(internal_mask)))>0

        for row_cntr=1:length(fiber_all(:,1,1,1))

            for col_cntr=1:length(fiber_all(1,:,1,1))

                loop_fiber = squeeze(fiber_all(row_cntr, col_cntr, :,:));
                point1 = find(loop_fiber(:,1), 1);
                point2 = find(loop_fiber(:,1), 1, 'last');
                num_points = point2 - point1 + 1;
                loop_fiber = loop_fiber(point1:point2, :);

                % examine outer points
                if num_points>6

                    for fiber_cntr=1:min([num_points 4])

                        loop_point = loop_fiber(fiber_cntr, :);

                        % calculate floor and ceiling values for each point
                        loop_r = round(loop_point(1));
                        loop_c = round(loop_point(2));
                        loop_s = max([1 round(loop_point(3))]);
                        loop_s = min([loop_s find(sum(sum(diffusion_mask)), 1, 'last')]);

                        if diffusion_mask(loop_r, loop_c, loop_s)==0
                            loop_fiber(fiber_cntr, :) = 0;
                        end

                        if fiber_cntr==3

                            fiber_all(row_cntr, col_cntr, :,:) = 0;

                            if sum(sum(loop_fiber)>0)

                                loop_fiber = [nonzeros(loop_fiber(:,1)) nonzeros(loop_fiber(:,2)) nonzeros(loop_fiber(:,3))];
                                length_loop_fiber = length(loop_fiber(:,1));
                                fiber_all(row_cntr, col_cntr, 1:length_loop_fiber, :) = loop_fiber(1:length_loop_fiber, :);

                            end

                        end

                    end % of num_points>6

                end % of fiber_cntr loop

            end % of column loop

        end % of row loop

    end % of if there is an internal mask with data in it

end %of non-apo seedign methods

%% final stuff

% save outputs in a default file
% save fiber_tract_data fiber_all roi_flag stop_list fiber_len

% plot fiber, mask, and mesh, if desired

if exist('fv_options', 'var') && exist('anat_image', 'var')

    if fv_options.plot_mesh==1 && fv_options.plot_fibers==1
        fiber_visualizer_v11(anat_image, fv_options, roi_mesh, diffusion_mask, fiber_all);
    elseif fv_options.plot_fibers==1
        fiber_visualizer_v11(anat_image, fv_options, [], diffusion_mask, fiber_all);
    end

end

%% end the function
return

%% function to check sign of E1 and adjust

function adjusted_E1 = check_E1_sign(seed_method, E1, direction)

adjusted_E1 = E1;

if size(direction)==1
    direction_sign = sign(direction);

    switch seed_method
        case{'a'}
            if direction_sign*E1(3) < 0                                     %if the x component is <0, reverse it to force tracking in the +X direction.
                adjusted_E1 = -E1;
            end
        case{'c'}
            if direction_sign*E1(1) < 0                                     %if the x component is <0, reverse it to force tracking in the +X direction.
                adjusted_E1 = -E1;
            end
        case{'r'}
            if direction_sign*E1(2) < 0                                     %if the y component is <0, reverse it to force tracking in the +Y direction.
                adjusted_E1 = -E1;
            end
        case{'s'}
            if direction_sign*E1(3) < 0                                     %if the z component is <0, reverse it to force tracking in the +Z direction.
                adjusted_E1 = -E1;
            end
        case{'v'}
            if direction_sign*E1(3) < 0                                      %here, this is just a flag to ensure tracking in opposite directions from the seed point
                adjusted_E1 = -E1;
            end
        case{'e'}
            if direction_sign*E1(3) < 0                                      %here, this is just a flag to ensure tracking in opposite directions from the seed point
                adjusted_E1 = -E1;
            end
        otherwise
            warning('Seed method was either not specified or incorrectly specified; using original eigenvector')
    end

elseif numel(direction)==3                       %if a direction vector was input

    if dot(E1, direction)<0                     %reverse the sign of E1 if the angular difference is >90 and <270
        adjusted_E1 = -E1;
    end

end

return




