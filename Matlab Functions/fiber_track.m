function [fiber_all, roi_flag, stop_list, fiber_len, fa_all, md_all] = fiber_track(tensor_m, mask, roi_mesh, ft_options, plot_options, anat_image)
%
%FUNCTION fiber_track
%  [fiber_all, roi_flag, stop_list, fiber_len, fa_all, md_all]=...
%     fiber_track(tensor_m, mask, roi_mesh, ft_options, plot_options, anat_image);
%
%USAGE
%    The function fiber_track is used to fiber-track a muscle DTI dataset in
%  the MuscleDTI_Toolbox. 
%    The required inputs include a 5D matrix containing the diffusion tensor; 
%  the mask delimiting the muscle of interest, the mesh reconstruction of the
%  aponeurosis of muscle fiber insertion, and a structure called ft_options.  
%  This structure allows the user to set options such as the tracking algorithm, 
%  step size, frame of reference for the images and diffusion encoding gradients,
%  and tract termination method. 
%    Fibers are tracked from the mesh according to the selected propogation 
%  algorithm until they reach the edge of the mask or meet another stop criterion.
%  Stop criteria are set in the ft_options structure. See the description of 
%  the input arguments for additional information on these variables.
%    The outputs include the fiber tracts, variables describing the outcomes 
%  of the tracking, and selected data about the tracts.
%
%INPUT ARGUMENTS
%  tensor_m: A 5D matrix containing rows, columns, slices, and the 3x3
%    diffusion tensor, calculated from pre-processing steps
%
%  mask: The mask delimiting the muscle to be fiber-tracked. It could be
%    the output of define_muscle or any other image analysis program that
%    creates a binary mask of the same size as the imaging data.   
%
%  roi_mesh: The roi mesh, output from define_roi.  
%
%  ft_options: A structure containing the following fields:
%    ref_frame: The frame of reference in which the diffusion directions
%      are specified. For example, set ft_options.ref_frame='LPS'; if the
%      left, posterior, and superior anatomical positions are (+).
%
%    image_orient: The orientation of the images. Specify the anatomical
%      positions at the north and east edges of the image as A (anterior) or 
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
%          eigenvector. The user must specify the step size in the field
%          ft_options.step_size.
%        -tnsrln: The tensorlines algorithm (Lazar et al, Human Brain Mapping,
%          2003). The tensorlines algorithm combines streamline and tensor
%          deflection components when calculating the propagation direction. 
%          The tensor deflection component can be weighted between deflected
%          and non-deflected terms using the parameter w_punct. w_punct can vary 
%          from 0 to 1, where 1 provides full weighting for the deflected component. 
%          To set w_punct, create a field called ft_options.w_punct. The user 
%          must also specify the value of the largest eigenvalue throughout 
%          the muscle. To do so, create a field called ft_options.eigvalue1_max 
%          and set it to the largest eigenvalue observed in the muscle of interest.
%        -rk4: 4th order Runge-Kutta integration of the first eigenvector.  
%          Note that if the 2nd, 3rd, or 4th order points fall outside of the 
%          mask, the propogation algorithm automatically changes to Euler
%          integration.
%        -fact: The FACT algorithm, as described by Mori et al (Ann Neurol, 
%          1999). FACT is similar to Euler integration, except that the direction 
%          is changed as soon as the tract enters a new voxel.
%
%    step_size: the fiber-tracking step size, in pixels. A step size of 1
%      reflects the voxel width. This is not used for FACT.
%
%    term_mthd: a string variable that specifies the method for determining
%      whether or not to terminate a fiber tract. Any fiber tracking point
%      that falls outside of the image mask will terminate the tract.
%      Other available options include:
%        -bin1: the angle and FA data from the current fiber tracking
%          point are used to decide whether or not to terminate the tract.
%          The angle used is the angle formed by two fiber tracking steps. 
%          The user can decide whether to calculate this angle between
%          a step and its immediate predecessor (1-back) or the step and a
%          step N points ago. The user sets these criteria in the fields
%          term_mthd.angle and term_mthd.FA. If either the angle or FA value
%          is disallowed, then the tract terminates.
%        -bin2: the angle and FA criteria from the current fiber tracking
%          point are combined with those of the preceding point. The
%          user sets the FA criteria as for bin1. If the two consecutive points 
%          have a disallowed FA value, then the tract terminates. For the 
%          angle criterion, the step angle is calculated for teh current point
%          and for one looking back N points (N must be > 1). If the current 
%          step and the preceding step queried have steps that exceed the
%          the angle threshold, then the tract terminates. This option provides 
%          greather tolerance for errors in individual voxels.
%      The FACT algorithm uses its own method for tract termination. Thus,
%      when the propogation algorithm is set to FACT, the user must also
%      create fields called r_crit and num_fact_voxels.  These are
%      the direction used to terminate tracts based on local variability in 
%      the first eigenvector.
%
%    angle_thrsh: A two-element vector containing the angle threshold in
%      degrees and the number of look-back steps (used as described under
%      term_mthd)
%
%    fa_thrsh: a two-element vector containing the lower and upper bounds
%      of allowable FA values (used as described under term_mthd)
%
%    depth_ratio: ratio of slice thickness/in-plane resolution. Note that
%      the function assumes equal in-plane voxel dimensions.
%
%  plot_options: If specified, this calls the fiber_visualizer function to
%    plot the fiber, mask, and roi mesh.
%  
%  anat_image: The structural images, of the same size as the DTI images.
%    These are required only if the user wishes to plot the fiber tracts.
%
%OUTPUT ARGUMENTS
%  fiber_all: the fiber tracts, with units of pixels. The rows and columns
%    correspond to locations on the roi_mesh. Dimension 3 gives point numbers
%    on the tract, and the fourth dimension has row, column, and slice coordinates.
%
%  roi_flag: matrix indicating the presence of fibers that propagated at
%    least 1 point
%
%  stop_list: matrix containing the reason for fiber tract termination
%    (4=mask, 3=curvature, 2=FA, 1=R (for FACT only))
%
%  fiber_len: the length, in points, of each fiber tract. 
%
%  fa_all: the pointwise FA values on each fiber tract.
%
%  md_all: the pointwise mean diffusivities along each tract
%
%OTHER FUNCTIONS IN THE MUSCLE DTI FIBER-TRACKING TOOLBOX
%  For help defining the mask, see <a href="matlab: help define_muscle">define_muscle</a>.
%  For help defining the ROI, see <a href="matlab: help define_roi">define_roi</a>.
%  For help smoothing fiber tracts, see <a href="matlab: help fiber_smoother">fiber_smoother</a>.
%  For help quantifying fiber tracts, see <a href="matlab: help fiber_quantifier">fiber_quantifier</a>.
%  For help selecting fiber tracts following quantification, see <a href="matlab: help fiber_selector">fiber_selector</a>.
%  For help visualizing the data, see <a href="matlab: help fiber_visualizer">fiber_visualizer</a>.
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
            
            case{'LPS'}                                                     %image left, posterior, and superior are (+) directions (DICOM)
                e1_order = [2 1 3];
                e1r_sign = 1;                                               %sign of E1 component in row direction
                e1c_sign = 1;                                               %sign of E1 component in column direction
                e1s_sign = 1;                                               %sign of E1 component in slice direction
                
            case{'LAS'}                                                     %image left, anterior, and superior are (+) directions (NIFTI)
                e1_order = [2 1 3];
                e1r_sign = -1;
                e1c_sign = 1;
                e1s_sign = 1;
                
            case{'RAS'}                                                     %image right, anterior, and superior are (+) directions (NIFTI)
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
        
    case {'tn'}
        
        % tracking options
        step_incr = ft_options.step_size;
        w_punct = ft_options.w_punct;                             %puncture coefficient for tensorlines algorithm
        eigvalue1_max = ft_options.eigvalue1_max;
        
        %check for illegal values
        if w_punct < 0 || w_punct>1
            error('Puncture coefficient is out of bounds')
        end
        
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
    fprintf('%6.2f percent is tracked\n', row_cntr/length(roi_mesh(:,1,1))*100);
    
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
        tensor_s = max([round(seed_point(3)) 1]);                           %just in case teh z coordinate rounds to zero, make it a one
        
        % 2) check to ensure that point is within the mask before tracking; if not, continue to next location on roi_mesh
        if mask(tensor_r,tensor_c,tensor_s)==0
            stop_list(row_cntr,col_cntr,1)=4;
            continue;
        end
        
        % 3) add to the fiber_all matrix
        fiber_all(row_cntr,col_cntr, fiber_cntr,:) = seed_point;
        
        switch prop_algo
            
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
                
            otherwise
                
                % 4) solve diffusion tensor
                d_tensor = retrieve_tensor(tensor_m, [tensor_r tensor_c tensor_s]);	%call the function retrieve tensor to the get local diffusion tensor
                [eigvector_m, eigvalue_m] = eig(d_tensor);                     %diagonalize to get e-vectors and e-values
                
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
                
                if prop_algo == 'tn'                                            %only for the tensorline algorithm, save the incoming E1
                    v_in = eigvector1;
                end
                
                % 5) get initial step direction as the first eigenvector, adjust as
                % needed for image properties
                step_dir = eigvector1;                                              %initially set the step direction as E1
                step_dir = step_dir(e1_order);                                      %adjust for frame of reference and row/column indexing
                step_dir = [e1r_sign*step_dir(1) e1c_sign*step_dir(2) e1s_sign*step_dir(3)]';  
                step_dir(3) = step_dir(3)/depth_ratio;                              %account for the slice aspect ratio
                
                % 6) calculate the next point
                
                if prop_algo=='fa'                                                  %for FACT, calculate the step increment
                    
                    current_point = squeeze(fiber_all(row_cntr,col_cntr, fiber_cntr,:));
                    diff_next_r = diff(floor(current_point(1) + step_dir(1)*(0.01:0.01:1)));
                    diff_next_c = diff(floor(current_point(2) + step_dir(2)*(0.01:0.01:1)));
                    diff_next_s = diff(floor(current_point(3) + step_dir(3)*(0.01:0.01:1)));
                    diff_next_point = sum([diff_next_r' diff_next_c' diff_next_s'], 2);
                    switch_next_pixel = min([min(find(diff_next_point))+1 100]);
                    step_incr = switch_next_pixel/100;
                    
                else                                                                %otherwise, can use the input value
                    
                    step_incr = ft_options.step_size;
                    
                end
                
        end
        
        % add the point
        next_point = squeeze(fiber_all(row_cntr,col_cntr, fiber_cntr,:)) + ...
            step_incr(1)*step_dir;                                             %multiply step_dir by the step size before adding the step to the preceding point
        fiber_cntr = fiber_cntr+1;                                          %increment the fiber counter
        fiber_all(row_cntr,col_cntr, fiber_cntr,:) = next_point;
        
        % 7) begin the fiber tracking loop
        
        while 1
            
            % 7a) get indices into tensor matrix from the current point
            tensor_r = round(next_point(1));
            tensor_c = round(next_point(2));
            tensor_s = max([round(next_point(3)) 1]);
            
            % 7b) check mask; terminate tract if out of bounds. record as
            % stop criterion = 4
            if mask(tensor_r,tensor_c,tensor_s)==0
                stop_list(row_cntr,col_cntr) = 4;
                break;
            end
            
            % 7c) get the step direction
            switch prop_algo
                
                case {'eu'}                                                %if Euler integration is selected
                    
                    % get the diffusion tensor
                    d_tensor = retrieve_tensor(tensor_m, [tensor_r tensor_c tensor_s]);	%call the function retrieve tensor to the get local diffusion tensor
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
                    
                case {'tn'}                                             %if the propagation algorithm is tensor deflection
                    
                    % get the diffusion tensor by calling the retrieve_tensor function.
                    d_tensor = retrieve_tensor(tensor_m, [tensor_r tensor_c tensor_s]);	%call the function retrieve tensor to the get local diffusion tensor
                    [eigvector_m, eigvalue_m] = eig(d_tensor);             %diagonalize to get e-vectors and e-values
                    
                    %get the derived indices
                    [eigvalues, sort_order] = sort(diag(eigvalue_m), 'descend');
                    eigvalue2 = eigvalues(2);
                    eigvalue3 = eigvalues(3);
                    eigvector1 = eigvector_m(:, sort_order(1));
                    if eigvector1(3) < 0                                     %if the z component is <0, the eigenvector is upside down, so reverse it.
                        eigvector1 = -eigvector1;
                    end
                    eigvector1_all(fiber_cntr,:) = eigvector1;
                    v_1 = eigvector1;
                    
                    %get the derived indices, MD and FA
                    mean_diff = mean(eigvalues(:,1));                                   %calculate the mean diffusivity and FA
                    md_all(row_cntr, col_cntr,fiber_cntr) = mean_diff;
                    fa_all(row_cntr, col_cntr,fiber_cntr) = sqrt(3/2) * ...
                        sqrt(sum((eigvalues(:,1)-mean_diff).^2)/sum(eigvalues(:,1).^2));
                    
                    %calculate v_out
                    v_out = d_tensor*v_in;
                    
                    %set coefficient c_l for weighting streamlines vs.
                    %tensor-deflection components
                    c_l = (eigvalue1_max - eigvalue2) / (eigvalue1_max + eigvalue2 + eigvalue3);
                    
                    %calculate step direction
                    v_1 = v_1/norm(v_1);
                    v_in = v_in/norm(v_in);
                    v_out = v_out/norm(v_out);
                    v_prop = c_l*v_1 + (1-c_l)*((1-w_punct)*v_in + w_punct*v_out);    %equation 5 in Weinstein et al
                    
                    step_dir = v_prop;
                    step_dir = step_dir(e1_order);                                      %adjust for frame of reference and row/column indexing
                    step_dir = [e1r_sign*step_dir(1) e1c_sign*step_dir(2) e1s_sign*step_dir(3)]';  
                    step_dir(3) = step_dir(3)/depth_ratio;                              %account for the slice aspect ratio
                    
                    %prepare for next time through the fiber tracking loop
                    v_in = eigvector1;
                    
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
                    d_tensor = retrieve_tensor(tensor_m, [tensor_r tensor_c tensor_s]);	%call the function retrieve tensor to the get local diffusion tensor
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
                    step_dir(3) = step_dir(3)/depth_ratio;                              %account for the slice aspect ratio
                    
            end                                                             %end of step direction switch statement
            
            % 7d) decide whether or not to add the point
            switch term_mthd                                                %sekect between tract termination methods
                
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
                                dist_mtrx(img_r, img_c, img_s,:) = tensor_indx_3d(img_r, img_c, img_s, :) ...
                                    - fiber_all(row_cntr,col_cntr, fiber_cntr,:);
                            end
                        end
                    end
                    
                    dist_mtrx = dist_mtrx.^2;
                    dist_mtrx = sum(dist_mtrx,4);
                    dist_mtrx = dist_mtrx.^(0.5);
                    
                    % identify the num_fact_voxels voxels nearest to the current point
                    dist_mtrx(dist_mtrx==min(min(min(dist_mtrx)))) = 100;   %set the current pixel to distance of 100 - don't want to count it
                    dist_vctr = reshape(dist_mtrx, numel(dist_mtrx), 1);
                    dist_vctr(:,2)=1:length(dist_vctr);
                    dist_vctr=sortrows(dist_vctr);
                    
                    % loop through the nearest neighbor voxels and get the first e-vectors for 
                    % the nearest ones; calculate dot products with current e'vector 1; store in the
                    % matrix dot_prod_m
                    dot_prod_m = zeros(1,num_fact_voxels);
                    fact_cntr=1;
                    
                    while fact_cntr<num_fact_voxels
                        
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
                        fact_cntr=fact_cntr+1;
                        
                    end
                    
                    %apply to stop criterion
                    r_stop=sum(dot_prod_m)/num_fact_voxels;
                    if r_stop > r_crit
                        stop_list(row_cntr, col_cntr) = r_stop;
                        fiber_len(row_cntr, col_cntr) = fiber_cntr;
                        break;
                    end
                    
            end                                                             %end of termination method switch statement
            
            % 8) in case of FACT, calculate the step
            if prop_algo=='fa'
                
                current_point = squeeze(fiber_all(row_cntr,col_cntr, fiber_cntr,:));        %interpolate at steps of 0.01 pixels, then see when the tract propogates into the next pixel. 
                diff_next_r = diff(floor(current_point(1) + step_dir(1)*(0.01:0.01:1)));
                diff_next_c = diff(floor(current_point(2) + step_dir(2)*(0.01:0.01:1)));
                diff_next_s = diff(floor(current_point(3) + step_dir(3)*(0.01:0.01:1)));
                diff_next_point = sum([diff_next_r' diff_next_c' diff_next_s'], 2);
                switch_next_pixel = min([min(find(diff_next_point))+1 100]);
                step_incr = switch_next_pixel/100;                                          %and set that to be the step size
                
            end
            
            % 9) add the point and continue tracking
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

