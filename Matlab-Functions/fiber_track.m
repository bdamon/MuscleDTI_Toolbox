function [fiber_all, roi_flag, stop_list, fiber_len, fa_all, md_all] = ...
    fiber_track(tensor_m, mask, roi_mesh, ft_options, plot_options, anat_image)
%
%FUNCTION fiber_track
%  [fiber_all, roi_flag, stop_list, fiber_len, fa_all, md_all] = ...
%     fiber_track(tensor_m, mask, roi_mesh, ft_options, plot_options, anat_image);
%
%USAGE
%  The function fiber_track is used to fiber-track a muscle DTI dataset.
%    
%  The required inputs include a 5D matrix hold the diffusion tensor at each
%  voxel and [row column slice] dimensions matching those of the DTMRI data; 
%  the muscle mask, output from define_muscle or other program; the
%  aponeurosis mesh, output from define_roi; and a structure defining the 
%  fiber-tracking options.  This structure allows the user to set options 
%  such as the tracking algorithm, step size, laboratory frame of reference,
%  image orientation, and tract termination method. 
%
%  Fibers are tracked from the aponeurosis mesh according to the selected 
%  propagation algorithm. Each fiber tract is propagated until it reaches the
%  edge of the mask or meets another stop criterion, such as an excessive 
%  inter-segment angle or an out-of-bounds value for fractional anisotropy 
%  (FA).  See the description of the input arguments for additional 
%  information on these criteria. 
%
%  The outputs include the fiber tracts, variables describing the outcomes 
%  of the tracking, and selected data about the tracts.
%
%  The fiber tracts may be viewed using fiber_visualizer, either as part of 
%  the function call to fiber_track or directly from the command line. 
%
%INPUT ARGUMENTS
%  tensor_m: A 5D matrix, with the first-third dimensions matching the 
%  [row column slice] size of the DTI images and the fourth and fifth
%  dimensions holding the 3x3 diffusion tensor at each voxel, calculated 
%  from pre-processing steps
%
%  mask: The mask delimiting the muscle to be fiber-tracked. 
%
%  roi_mesh: The mesh reconstruction of the aponeurosis of muscle fiber 
%    insertion, output from define_roi.
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
%    mesh_dist: The number of pixels to shift the mesh into the muscle,
%      prior to fiber tracking. This can be a (+) or (-) number, depending
%      on the desired direction of the shift.
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
%        -fact: An implementation of the FACT algorithm. FACT is similar to
%          Euler integration, except that the direction is changed as soon 
%          as the tract enters a new voxel.
%
%    step_size: The Euler and 4th-order Runge-Kutta methods require the user
%      to set the fiber-tracking step size, in pixels. A step size of 1 
%      reflects the voxel width
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
%      The FACT algorithm uses its own method for tract termination. Thus,
%      when the propogation algorithm is set to FACT, the user must also
%      create fields called r_crit and num_fact_voxels.  These are
%      used to terminate tracts based on local variability in the first 
%      eigenvector.
%
%    angle_thrsh: A two-element vector containing the angle threshold in
%      degrees and the number of look-back steps.  This is only required when
%      bin1 or bin2 is used.
%
%    fa_thrsh: a two-element vector containing the lower and upper bounds
%      of allowable FA values. This is only required when bin1 or bin2 is
%      used.
%
%    depth_ratio: ratio of slice thickness/in-plane resolution. Note that
%      the function assumes equal in-plane voxel dimensions.
%
%  The following input arguments are optional and are only required if the
%  user wishes toplot the plot tracts from within the fiber_track function.
% 
%  plot_options: If specified, this calls the fiber_visualizer function to
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
%  fa_all: the pointwise FA values on each fiber tract.
%
%  md_all: the pointwise mean diffusivities along each tract
%
%OTHER FUNCTIONS IN THE MUSCLE DTI FIBER-TRACKING TOOLBOX
%  For help visualizing the data, see <a href="matlab: help fiber_visualizer">fiber_visualizer</a>.
%  For help defining the mask, see <a href="matlab: help define_muscle">define_muscle</a>.
%  For help defining the ROI, see <a href="matlab: help define_roi">define_roi</a>.
%  For help smoothing fiber tracts, see <a href="matlab: help fiber_smoother">fiber_smoother</a>.
%  For help quantifying fiber tracts, see <a href="matlab: help fiber_quantifier">fiber_quantifier</a>.
%  For help selecting fiber tracts following their quantification, see <a href="matlab: help fiber_goodness">fiber_goodness</a>.
%
%VERSION INFORMATION
%  v. 0.1
%
%ACKNOWLEDGEMENTS
%  People: Zhaohua Ding, Adam Anderson, Amanda Buck, Anneriet Heemskerk,
%    and Justin Montenegro
%  Grant support: NIH/NIAMS R01 AR050101, NIH/NIAMS R01 AR073831

%% get user options for fiber tracking

%basic tracking options
mesh_dist = ft_options.mesh_dist;
depth_ratio = ft_options.depth_ratio;

%get the frame of reference and image orientation to set conventions for
%propagating points
ref_frame = ft_options.ref_frame;
image_orient = ft_options.image_orient;

switch image_orient
    
    case{'AL'}                                                              %image north is anatomical anterior, image east is anatomical left
        
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
        
    case {'RA'}                                                              %image north is anatomical right, image east is anatomical anterior
        
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
        
    case {'fa'}
        
        term_mthd = 'fact';
        r_crit = ft_options.r_crit;
        num_fact_voxels = ft_options.num_fact_voxels;
        
end


%% fiber track

% initialize matrix to hold the x, y, and z coordinates of the fiber tracts, with units of pixels
fiber_all = zeros(length(roi_mesh(:, 1, 1)), length(roi_mesh(1, :, 1)), 100, 3);    %initialize as roi mesh size by up to 100 points long; can grow, as needed
fa_all = zeros(length(roi_mesh(:, 1, 1)), length(roi_mesh(1, :, 1)), 100);    %initialize as roi mesh size by up to 100 points long; can grow, as needed
md_all = zeros(length(roi_mesh(:, 1, 1)), length(roi_mesh(1, :, 1)), 100);    %initialize as roi mesh size by up to 100 points long; can grow, as needed

% initialize quality-related variables
roi_flag = zeros(size(squeeze(roi_mesh(:, :, 1))));                         %will hold information about tract initiation
stop_list=zeros(length(roi_mesh(:,1,1)),length(roi_mesh(1,:,1)));         %will hold information about the stop criteria
fiber_len=roi_flag;

% fiber track by passing through rows and columns of the mesh
for row_cntr=1:length(roi_mesh(:, 1, 1))                                  %start of the row loop
    
    %progress counter
%     fprintf('%6.2f percent is tracked\n', row_cntr/length(roi_mesh(:,1,1))*100);
    
    for col_cntr=1:length(roi_mesh(1, :, 1))                              %start of the column loop
        
        % initialize some variables
        fiber_cntr = 1;                                                     %counter for the fiber tract points; initialize at 1 for each tract
        angle_all = zeros(100,1);                                           %keep track of all inter-step angles
        eigvector1_all = zeros(100,3);                                      %matrix to hold all of the first eigenvectors (needed for checking inter-point angle stop criterion)
        
        % 1) locate the seed point on the roi_mesh
        seed_point(1:3) = [roi_mesh(row_cntr, col_cntr, 1)+mesh_dist*roi_mesh(row_cntr, col_cntr, 4) ...
            roi_mesh(row_cntr, col_cntr, 2)+mesh_dist*roi_mesh(row_cntr, col_cntr, 5) ...
            roi_mesh(row_cntr, col_cntr, 3)+mesh_dist*roi_mesh(row_cntr, col_cntr, 6)];
        tensor_r = round(seed_point(1));                                    %round the coordiantes to get a row, column, and slice location in the tensor matrix
        tensor_c = round(seed_point(2));
        tensor_s = max([round(seed_point(3)) 1]);                           %just in case the z coordinate rounds to zero, make it a one
        
        % 2) add to the fiber_all matrix
        if mask(tensor_r,tensor_c,tensor_s)==0
            stop_list(row_cntr,col_cntr,1)=4;
            continue;
        end
        
        fiber_all(row_cntr,col_cntr, fiber_cntr,:) = seed_point;
        
        % 3) get the initial direction
        switch prop_algo
            
            case {'eu'}
                
                % get the diffusion tensor
                tensor_indx = [tensor_r tensor_c tensor_s];
                d_tensor = retrieve_tensor(tensor_m, tensor_indx);	%call the function retrieve tensor to the get local diffusion tensor
                [eigvector_m, eigvalue_m] = eig(d_tensor);             %diagonalize to get e-vectors and e-values
                
                %get the derived indices
                [eigvalues, sort_order] = sort(diag(eigvalue_m), 'descend');
                eigvector1 = eigvector_m(:, sort_order(1));
                if eigvector1(3) < 0                                     %if the z component is <0, the eigenvector is upside down, so reverse it.
                    eigvector1 = -eigvector1;
                end
                eigvector1_all(fiber_cntr,:) = eigvector1;
                
                %get the derived indices, FA and E1
                mean_diff = mean(eigvalues(:,1));                                   %calculate the mean diffusivity and FA
                md_all(row_cntr, col_cntr,fiber_cntr) = mean_diff;
                fa_all(row_cntr, col_cntr,fiber_cntr) = sqrt(3/2) * ...
                    sqrt(sum((eigvalues(:,1)-mean_diff).^2)/sum(eigvalues(:,1).^2));
                
                step_dir = eigvector1;                                              %initially set the step direction as E1
                step_dir = step_dir(e1_order);                                      %adjust for frame of reference and row/column indexing
                step_dir = [e1r_sign*step_dir(1) e1c_sign*step_dir(2) e1s_sign*step_dir(3)]';
                step_dir(3) = step_dir(3)/depth_ratio;                              %account for the slice aspect ratio
                    
            case {'rk'} %Runge-Kutta integration (Richard Hamming. 1989 Introduction to Applied Numerical Analysis. Dover Publications, p. 212)
                
                % for point p1, get the diffusion tensor by calling the retrieve_tensor function.
                point1 = squeeze(fiber_all(row_cntr,col_cntr, fiber_cntr,:));
                tensor_indx = [round(point1(1)) round(point1(2)) max([round(point1(3)) 1])];
                d_tensor = retrieve_tensor(tensor_m, tensor_indx);
                [eigvector_m, eigvalue_m] = eig(d_tensor);                     %diagonalize to get e-vectors and e-values
                
                %get the derived indices
                [eigvalues, sort_order] = sort(diag(eigvalue_m), 'descend');
                eigvector1_point1 = eigvector_m(:, sort_order(1));
                if eigvector1_point1(3) < 0                                     %if the z component is <0, the eigenvector is upside down, so reverse it.
                    eigvector1_point1 = -eigvector1_point1;
                end
                
                mean_diff = mean(eigvalues(:,1));                               %calculate the mean diffusivity and FA, but only for the initial point
                md_all(row_cntr, col_cntr,fiber_cntr) = mean_diff;
                fa_all(row_cntr, col_cntr, fiber_cntr) = sqrt(3/2) * ...
                    sqrt(sum((eigvalues(:,1)-mean_diff).^2)/sum(eigvalues(:,1).^2));
                
                % get path direction as the first eigenvector
                k1 = eigvector1_point1;                                         %initially set the step direction as E1
                k1 = k1(e1_order);                                              %adjust for frame of reference and row/column indexing
                k1 = [e1r_sign*k1(1) e1c_sign*k1(2) e1s_sign*k1(3)]';
                k1(3) = k1(3)/depth_ratio;                                      %account for slice aspect ratio in the z direction.
                
                %calculate next point, p2
                point2 = point1 + (step_incr/2)*k1;
                
                % get the diffusion tensor at point p2
                tensor_indx = [round(point2(1)) round(point2(2)) max([round(point2(3)) 1])];
                d_tensor = retrieve_tensor(tensor_m, tensor_indx);
                [eigvector_m, eigvalue_m] = eig(d_tensor);                     %diagonalize to get e-vectors and e-values
                
                %get the derived indices
                [~, sort_order] = sort(diag(eigvalue_m), 'descend');
                eigvector1_point2 = eigvector_m(:, sort_order(1));
                if eigvector1_point2(3) < 0                                     %if the z component is <0, the eigenvector is upside down, so reverse it.
                    eigvector1_point2 = -eigvector1_point2;
                end
                
                % get path direction as the first eigenvector
                k2 = eigvector1_point2;                                         %initially set the step direction as E1
                k2 = k2(e1_order);                                              %adjust for frame of reference and row/column indexing
                k2 = [e1r_sign*k2(1) e1c_sign*k2(2) e1s_sign*k2(3)]';
                k2(3) = k2(3)/depth_ratio;                                      %account for slice aspect ratio in the z direction.
                
                %calculate next point, p3
                point3 = point2 + (step_incr/2)*k2;
                
                % get the diffusion tensor at point p3
                tensor_indx = [round(point3(1)) round(point3(2)) max([round(point3(3)) 1])];
                d_tensor = retrieve_tensor(tensor_m, tensor_indx);
                [eigvector_m, eigvalue_m] = eig(d_tensor);                     %diagonalize to get e-vectors and e-values
                
                %get the derived indices
                [~, sort_order] = sort(diag(eigvalue_m), 'descend');
                eigvector1_point3 = eigvector_m(:, sort_order(1));
                if eigvector1_point3(3) < 0                                     %if the z component is <0, the eigenvector is upside down, so reverse it.
                    eigvector1_point3 = -eigvector1_point3;
                end
                
                % get path direction as the first eigenvector
                k3 = eigvector1_point3;                                         %initially set the step direction as E1
                k3 = k3(e1_order);                                              %adjust for frame of reference and row/column indexing
                k3 = [e1r_sign*k3(1) e1c_sign*k3(2) e1s_sign*k3(3)]';
                k3(3) = k3(3)/depth_ratio;                                      %account for slice aspect ratio in the z direction.
                
                %calculate next point, p4
                point4 = point3 + step_incr*k3;
                
                % get the diffusion tensor at point p4
                tensor_indx = [round(point4(1)) round(point4(2)) max([round(point4(3)) 1])];
                d_tensor = retrieve_tensor(tensor_m, tensor_indx);
                [eigvector_m, eigvalue_m] = eigs(d_tensor);                     %diagonalize to get e-vectors and e-values
                
                %get the derived indices
                [~, sort_order] = sort(diag(eigvalue_m), 'descend');
                eigvector1_point4 = eigvector_m(:, sort_order(1));
                if eigvector1_point4(3) < 0                                     %if the z component is <0, the eigenvector is upside down, so reverse it.
                    eigvector1_point4 = -eigvector1_point4;
                end
                
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
                if eigvector1_mean(3) < 0                                     %if the z component is <0, the eigenvector is upside down, so reverse it.
                    eigvector1_mean = -eigvector1_mean;
                end
                
                eigvector1 = eigvector1_mean;
                eigvector1_all(fiber_cntr,:) = eigvector1;
                step_dir = eigvector1;                                              %initially set the step direction as E1
                step_dir = step_dir(e1_order);                                      %adjust for frame of reference and row/column indexing
                step_dir = [e1r_sign*step_dir(1) e1c_sign*step_dir(2) e1s_sign*step_dir(3)]';
                step_dir(3) = step_dir(3)/depth_ratio;                              %account for the slice aspect ratio
                
            case{'fa'}
                
                % solve diffusion tensor
                tensor_indx = [tensor_r tensor_c tensor_s];
                d_tensor = retrieve_tensor(tensor_m, tensor_indx);	%call the function retrieve tensor to the get local diffusion tensor
                [eigvector_m, eigvalue_m] = eig(d_tensor);                     %diagonalize to get e-vectors and e-values
                
                %get the derived indices
                [eigvalues, sort_order] = sort(diag(eigvalue_m), 'descend');
                eigvector1 = eigvector_m(:, sort_order(1));
                if eigvector1(3) < 0                                     %if the z component is <0, the eigenvector is upside down, so reverse it.
                    eigvector1 = -eigvector1;
                end
                eigvector1_all(fiber_cntr,:) = eigvector1;
                step_dir = eigvector1;                                              %initially set the step direction as E1
                step_dir = step_dir(e1_order);                                      %adjust for frame of reference and row/column indexing
                step_dir = [e1r_sign*step_dir(1) e1c_sign*step_dir(2) e1s_sign*step_dir(3)]';
                step_dir(3) = step_dir(3)/depth_ratio;                              %account for the slice aspect ratio
                
                %get the derived indices, FA and E1
                mean_diff = mean(eigvalues(:,1));                                   %calculate the mean diffusivity and FA
                md_all(row_cntr, col_cntr,fiber_cntr) = mean_diff;
                fa_all(row_cntr, col_cntr,fiber_cntr) = sqrt(3/2) * ...
                    sqrt(sum((eigvalues(:,1)-mean_diff).^2)/sum(eigvalues(:,1).^2));
                
                current_point = squeeze(fiber_all(row_cntr,col_cntr, fiber_cntr,:));
                diff_next_r = diff(floor(current_point(1) + step_dir(1)*(0.01:0.01:1)));
                diff_next_c = diff(floor(current_point(2) + step_dir(2)*(0.01:0.01:1)));
                diff_next_s = diff(floor(current_point(3) + step_dir(3)*(0.01:0.01:1)));
                diff_next_point = sum([diff_next_r' diff_next_c' diff_next_s'], 2);
                switch_next_pixel = min([min(find(diff_next_point))+1 100]);
                step_incr = switch_next_pixel/100;

                
        end
        
        % add the point
        next_point = squeeze(fiber_all(row_cntr,col_cntr, fiber_cntr,:)) + ...
            step_incr(1)*step_dir;                                              %multiply step_dir by the step size before adding the step to the preceding point
        fiber_cntr = fiber_cntr+1;                                              %increment the fiber counter
        
        % get indices into tensor matrix from the current point
        tensor_r = round(next_point(1));
        tensor_c = round(next_point(2));
        tensor_s = max([round(next_point(3)) 1]);                               %take teh maximum with slice 1, in case of low seed point
        
        % check mask; terminate tract if out of bounds. record as stop criterion = 4 and continue to next point
        if mask(tensor_r,tensor_c,tensor_s)==0
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
                    
                    % get the diffusion tensor
                    tensor_indx = [tensor_r tensor_c tensor_s];
                    d_tensor = retrieve_tensor(tensor_m, tensor_indx);                    [eigvector_m, eigvalue_m] = eig(d_tensor);             %diagonalize to get e-vectors and e-values
                    
                    %get the derived indices
                    [eigvalues, sort_order] = sort(diag(eigvalue_m), 'descend');
                    eigvector1 = eigvector_m(:, sort_order(1));
                    if eigvector1(3) < 0                                     %if the z component is <0, the eigenvector is upside down, so reverse it.
                        eigvector1 = -eigvector1;
                    end
                    eigvector1_all(fiber_cntr,:) = eigvector1;
                    
                    %get the derived indices, FA and E1
                    mean_diff = mean(eigvalues(:,1));                                   %calculate the mean diffusivity and FA
                    md_all(row_cntr, col_cntr,fiber_cntr) = mean_diff;
                    fa_all(row_cntr, col_cntr,fiber_cntr) = sqrt(3/2) * ...
                        sqrt(sum((eigvalues(:,1)-mean_diff).^2)/sum(eigvalues(:,1).^2));
                    
                    step_dir = eigvector1;                                              %initially set the step direction as E1
                    step_dir = step_dir(e1_order);                                      %adjust for frame of reference and row/column indexing
                    step_dir = [e1r_sign*step_dir(1) e1c_sign*step_dir(2) e1s_sign*step_dir(3)]';
                    step_dir(3) = step_dir(3)/depth_ratio;                              %account for the slice aspect ratio
                    
                case {'rk'}
                    
                    %create a variable to keep track of the RK order; use order = 1 (Euler's method) unless all four points are inside the mask
                    rk_order=1;
                    
                    % define point p1; if it lies within the mask, get the diffusion tensor by calling the retrieve_tensor function.
                    point1 = squeeze(fiber_all(row_cntr,col_cntr, fiber_cntr,:));
                    
                    if mask(round(point1(1)), round(point1(2)), max([round(point1(3)) 1]))==1
                        
                        tensor_indx = [round(point1(1)) round(point1(2)) max([round(point1(3)) 1])];
                        d_tensor = retrieve_tensor(tensor_m, tensor_indx);
                        [eigvector_m, eigvalue_m] = eig(d_tensor);                     %diagonalize to get e-vectors and e-values
                        
                        %get the derived indices
                        [eigvalues, sort_order] = sort(diag(eigvalue_m), 'descend');
                        eigvector1_point1 = eigvector_m(:, sort_order(1));
                        if eigvector1_point1(3) < 0                                     %if the z component is <0, the eigenvector is upside down, so reverse it.
                            eigvector1_point1 = -eigvector1_point1;
                        end
                        
                        mean_diff = mean(eigvalues(:,1));                               %calculate the mean diffusivity and FA, but only for the initial point
                        md_all(row_cntr, col_cntr,fiber_cntr) = mean_diff;
                        fa_all(row_cntr, col_cntr, fiber_cntr) = sqrt(3/2) * ...
                            sqrt(sum((eigvalues(:,1)-mean_diff).^2)/sum(eigvalues(:,1).^2));
                        
                        % get path direction as the first eigenvector
                        k1 = eigvector1_point1;                                         %initially set the step direction as E1
                        k1 = k1(e1_order);                                              %adjust for frame of reference and row/column indexing
                        k1 = [e1r_sign*k1(1) e1c_sign*k1(2) e1s_sign*k1(3)]';
                        k1(3) = k1(3)/depth_ratio;                                      %account for slice aspect ratio in the z direction.
                        
                        %calculate next point, p2
                        point2 = point1 + (step_incr/2)*k1;
                        
                    else
                        
                        break
                        
                    end
                    
                    % if mask at point 2 is valid, get the diffusion tensor there and calculate point three
                    if mask(round(point2(1)), round(point2(2)), max([round(point2(3)) 1]))==1
                        
                        tensor_indx = [round(point2(1)) round(point2(2)) max([round(point2(3)) 1])];
                        d_tensor = retrieve_tensor(tensor_m, tensor_indx);
                        [eigvector_m, eigvalue_m] = eig(d_tensor);                     %diagonalize to get e-vectors and e-values
                        
                        %get the derived indices
                        [~, sort_order] = sort(diag(eigvalue_m), 'descend');
                        eigvector1_point2 = eigvector_m(:, sort_order(1));
                        if eigvector1_point2(3) < 0                                     %if the z component is <0, the eigenvector is upside down, so reverse it.
                            eigvector1_point2 = -eigvector1_point2;
                        end
                        
                        % get path direction as the first eigenvector
                        k2 = eigvector1_point2;                                         %initially set the step direction as E1
                        k2 = k2(e1_order);                                              %adjust for frame of reference and row/column indexing
                        k2 = [e1r_sign*k2(1) e1c_sign*k2(2) e1s_sign*k2(3)]';
                        k2(3) = k2(3)/depth_ratio;                                      %account for slice aspect ratio in the z direction.
                        
                        %calculate next point, p3
                        point3 = point2 + (step_incr/2)*k2;
                        
                        % if mask at point 3 is valid, get the diffusion tensor there and calculate point 4
                        if mask(round(point3(1)), round(point3(2)), max([round(point3(3)) 1]))==1
                            
                            tensor_indx = [round(point3(1)) round(point3(2)) max([round(point3(3)) 1])];
                            d_tensor = retrieve_tensor(tensor_m, tensor_indx);
                            [eigvector_m, eigvalue_m] = eig(d_tensor);                     %diagonalize to get e-vectors and e-values
                            
                            %get the derived indices
                            [~, sort_order] = sort(diag(eigvalue_m), 'descend');
                            eigvector1_point3 = eigvector_m(:, sort_order(1));
                            if eigvector1_point3(3) < 0                                     %if the z component is <0, the eigenvector is upside down, so reverse it.
                                eigvector1_point3 = -eigvector1_point3;
                            end
                            
                            % get path direction as the first eigenvector
                            k3 = eigvector1_point3;                                         %initially set the step direction as E1
                            k3 = k3(e1_order);                                              %adjust for frame of reference and row/column indexing
                            k3 = [e1r_sign*k3(1) e1c_sign*k3(2) e1s_sign*k3(3)]';
                            k3(3) = k3(3)/depth_ratio;                                      %account for slice aspect ratio in the z direction.
                            
                            %calculate next point, p4
                            point4 = point3 + step_incr*k3;
                            
                            % if mask at point 4 is valid, get the diffusion tensor there
                            if mask(round(point4(1)), round(point4(2)), max([round(point4(3)) 1]))==1
                                rk_order=4;
                                tensor_indx = [round(point4(1)) round(point4(2)) max([round(point4(3)) 1])];
                                d_tensor = retrieve_tensor(tensor_m, tensor_indx);
                                [eigvector_m, eigvalue_m] = eigs(d_tensor);                     %diagonalize to get e-vectors and e-values
                                
                                %get the derived indices
                                [~, sort_order] = sort(diag(eigvalue_m), 'descend');
                                eigvector1_point4 = eigvector_m(:, sort_order(1));
                                if eigvector1_point4(3) < 0                                     %if the z component is <0, the eigenvector is upside down, so reverse it.
                                    eigvector1_point4 = -eigvector1_point4;
                                end
                                
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
                        if eigvector1_mean(3) < 0                                     %if the z component is <0, the eigenvector is upside down, so reverse it.
                            eigvector1_mean = -eigvector1_mean;
                        end
                        
                        eigvector1 = eigvector1_mean;
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
                    
                case {'fa'}
                    
                    % get the diffusion tensor
                    d_tensor = squeeze(tensor_m(tensor_r, tensor_c, tensor_s, :, :));	%call the function retrieve tensor to the get local diffusion tensor
                    [eigvector_m, eigvalue_m] = eigs(d_tensor);             %diagonalize to get e-vectors and e-values
                    
                    %get the derived indices
                    [eigvalues, sort_order] = sort(diag(eigvalue_m), 'descend');
                    eigvector1 = eigvector_m(:, sort_order(1));
                    if eigvector1(3) < 0                                     %if the z component is <0, the eigenvector is upside down, so reverse it.
                        eigvector1 = -eigvector1;
                    end
                    eigvector1_all(fiber_cntr,:) = eigvector1;
                    
                    mean_diff = mean(eigvalues(:,1));                       %calculate the mean diffusivity and FA
                    md_all(row_cntr, col_cntr,fiber_cntr) = mean_diff;
                    fa_all(row_cntr, col_cntr,fiber_cntr) = sqrt(3/2) * ...
                        sqrt(sum((eigvalues(:,1)-mean_diff).^2)/sum(eigvalues(:,1).^2));
                    
                    % get path direction as the first eigenvector
                    step_dir = eigvector1;                                              %initially set the step direction as E1
                    step_dir = step_dir(e1_order);                                      %adjust for frame of reference and row/column indexing
                    step_dir = [e1r_sign*step_dir(1) e1c_sign*step_dir(2) e1s_sign*step_dir(3)]';
                    step_dir(3) = step_dir(3)/depth_ratio;                               %account for the slice aspect ratio
                    
            end                                                             %end of step direction switch statement
            
            % 4b) decide whether or not to add the point
            
            % get indices into tensor matrix from the current point
            tensor_r = round(next_point(1));
            tensor_c = round(next_point(2));
            tensor_s = max([round(next_point(3)) 1]);
            
            % first criterion - check mask; terminate tract if out of bounds and record as stop criterion = 4
            if mask(tensor_r,tensor_c,tensor_s)==0
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
                    
                case {'fact'}
                    
                    %get the data needed to calculate the three nearest voxels
                    dist_mtrx = zeros(3,3,3,3);                             %dist_mtrx will hold the distances from the current point to each surrounding voxel center.
                    tensor_indx_3d = zeros(3,3,3,3);                           %tensor_indx_3d will hold the indices into the tensor matrix
                    
                    for img_r=1:3                                           %loop through the row, column, and slice coordinates of voxel + its nieghbors into the matrix
                        for img_c=1:3
                            for img_s=1:3
                                tensor_indx_3d(img_r, img_c, img_s, :) = ...
                                    [img_r-2 img_c-2 img_s-2] + [tensor_r tensor_c tensor_s];   %indices into the tensor matrix - use below and in the retrive_tensor function call
                                dist_mtrx(img_r, img_c, img_s,:) = squeeze(tensor_indx_3d(img_r, img_c, img_s, :) ...
                                    - fiber_all(row_cntr,col_cntr, fiber_cntr,:));
                            end
                        end
                    end
                    
                    dist_mtrx = dist_mtrx.^2;
                    dist_mtrx = sum(dist_mtrx,4);
                    dist_mtrx = dist_mtrx.^(0.5);
                    
                    % identify the num_fact_voxels voxels nearest to the current point
                    dist_mtrx(2, 2, 2) = 100;   %set the current pixel to distance of 100 - don't want to count it
                    dist_vctr = reshape(dist_mtrx, numel(dist_mtrx), 1);
                    dist_vctr(:,2)=1:length(dist_vctr);
                    dist_vctr=sortrows(dist_vctr, 1);
                    
                    % loop through the nearest neighbor voxels and get the first e-vectors for
                    % the nearest ones; calculate dot products with current e'vector 1; store in the
                    % matrix dot_prod_m
                    dot_prod_m = zeros(1,num_fact_voxels);
                    fact_cntr=1;
                    
                    while fact_cntr<=num_fact_voxels
                        
                        nearest_voxel_idx = find(dist_mtrx==dist_vctr(fact_cntr,1));
                        nearest_voxel_idx = nearest_voxel_idx(1);
                        [nearest_voxel_row, nearest_voxel_col, nearest_voxel_slc] = ind2sub(size(dist_mtrx), nearest_voxel_idx);
                        
                        d_tensor = retrieve_tensor(tensor_m, ...
                            squeeze(tensor_indx_3d(nearest_voxel_row, nearest_voxel_col, nearest_voxel_slc, :)));	%call the function retrieve tensor to the get local diffusion tensor
                        [eigvector_m, eigvalue_m] = eigs(d_tensor);             %diagonalize to get e-vectors and e-values
                        
                        %get the derived indices
                        [~, sort_order] = sort(diag(eigvalue_m), 'descend');
                        eigvector1_fact = eigvector_m(:, sort_order(1));
                        if eigvector1_fact(3) < 0                                     %if the z component is <0, the eigenvector is upside down, so reverse it.
                            eigvector1_fact = -eigvector1_fact;
                        end
                        dot_prod_m(fact_cntr) = dot(eigvector1, eigvector1_fact);
                        
                        %set up for next time through the loop
                        dist_mtrx(nearest_voxel_row, nearest_voxel_col, nearest_voxel_slc)=100;       %don't want to count the same voxel twice
                        fact_cntr=fact_cntr+1;
                        
                    end
                    
                    %apply to stop criterion
                    r_stop=sum(dot_prod_m)/num_fact_voxels;
                    if r_stop < r_crit
                        stop_list(row_cntr, col_cntr) = r_stop;
                        fiber_len(row_cntr, col_cntr) = fiber_cntr;
                        break;
                    end
                    
            end                                                             %end of termination method switch statement
            
            % 4d) in case of FACT, calculate the step
            if prop_algo=='fa'
                
                current_point = squeeze(fiber_all(row_cntr,col_cntr, fiber_cntr,:));        %interpolate at steps of 0.01 pixels, then see when the tract propogates into the next pixel.
                diff_next_r = diff(floor(current_point(1) + step_dir(1)*(0.01:0.01:1)));
                diff_next_c = diff(floor(current_point(2) + step_dir(2)*(0.01:0.01:1)));
                diff_next_s = diff(floor(current_point(3) + step_dir(3)*(0.01:0.01:1)));
                diff_next_point = sum([diff_next_r' diff_next_c' diff_next_s'], 2);
                switch_next_pixel = min([min(find(diff_next_point))+1 100]);
                step_incr = switch_next_pixel/100;                                          %and set that to be the step size
                
            end
            
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


%% plot fiber, mask, and mesh, if desired

if exist('plot_options', 'var') && exist('anat_image', 'var')
    
    fiber_visualizer(anat_image, plot_options, roi_mesh, mask, fiber_all);
    
end

%%

return

