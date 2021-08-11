function [final_fibers, final_curvature, final_angle, final_distance, qual_mask, num_tracked, mean_fiber_props, mean_apo_props] = ...
    fiber_goodness(fiber_all, angle_list, distance_list, curvature_list, n_points, roi_flag, apo_area, fg_options)
%
%FUNCTION fiber_goodness
%  [final_fibers, final_curvature, final_angle, final_distance, qual_mask, num_tracked, mean_fiber_props, mean_apo_props] = ...
%     fiber_goodness(smoothed_fiber_all, angle_list, distance_list, curvature_list, n_points, roi_flag, apo_area, fg_options);
%
%USAGE
%  The function fiber_goodness is used to assess the goodness of fiber tract  
%  data, and reject outlying fiber tract data, in the MuscleDTI_Toolbox. 
%  
%  The quality algorithm described in Heemskerk et al, 2008 is implemented, 
%  but updated to account for the inclusion of curvature in the architectural
%  computations. Specifically, the fibers are selected for having:
%    1) Monotonically increasing values in the Z direction. In fiber_track, the  
%       Z component of the first eigenvector is forced to be positive; so the  
%       Z components of the unfitted tracts must also increase monotonically.   
%       However, negative steps in the Z direction could result from overfitting   
%       fiber tract smoothing process. Fiber tracts with monotonically increasing 
%       Z values are indicated by entering a value of 1 at the corresponding  
%       [row column] indices in the 1st level of the 3rd dimension of qual_mask.  
%       The number of tracts that meet this criterion are calculated and added
%       to the output argument num_tracked. Tracts that meet this criterion 
%       are advanced to the next level of selection.
%    2) A minimum length, in mm.  This is specified by the user in the min_distance 
%       field of the structure fg_options. Fiber tracts that meet this criterion  
%       are indicated by entering a value of 1 at the corresponding [row column] 
%       indices in the 2nd level of the 3rd dimension of qual_mask. The total 
%       number of tracts that meet this criterion are calculated and added to 
%       to the output argument num_tracked. Tracts that meet this criterion 
%       are advanced to the next level of selection.
%    3) An acceptable pennation angle, in degrees.  This range is defined 
%       by the user in the min_pennation and max_pennation fields of the 
%       structure fg_options. Fiber tracts that meet this criterion are indicated 
%       by entering a value of 1 at the corresponding [row column] indices 
%       in the 3rd level of the 3rd dimension of qual_mask.  The total number of 
%       tracts that meet this criterion are calculated and added to the 
%       to the output argument num_tracked. Tracts that meet this criterion 
%       are advanced to the next level of selection.
%    4) An acceptable curvature value, in m-1.  This range is defined by the 
%       user in the max_curvature field of the structure fg_options. Fiber tracts  
%       that meet this criterion are indicated by entering a value of 1 at the  
%       corresponding [row column] indices in the 4th level of the 3rd dimension  
%       of qual_mask. The total number of tracts that meet this criterion are  
%       calculated and added to to the output argument num_tracked. Tracts that  
%       meet this criterion are advanced to the next level of selection.
%    5) A length, pennation angle, and curvature that lies within the two 
%       standard deviations of the values in the surrounding 24 tracts. Fiber 
%       tracts that meet this criterion are indicated by entering a value of 
%       1 at the corresponding [row column] indices in the 5th level of the  
%       3rd dimension of qual_mask. The total number of tracts that meet this  
%       criterion are calculated and added to the output argument num_tracked. 
%       Tracts that meet this criterion are preserved. 
% The preserved fiber tracts are stored in the matrix final_fibers; their 
% structural properties are stored in the matrices final_curvature, final_angle, 
% and final_distance; and the whole-muscle mean values for length, pennation angle, 
% and curvature are stored in the matrix mean_apo_props.
%
% The use of length, pennation, and curvature criteria require the user to use 
% their knowledge of the expected patterns of muscle geometry to supply values 
% that are reasonable but will not inappropriately bias the results. These 
% selection criteria, as well as the number of tracts that were rejected 
% because of these criteria, should be included in publications.
%
%INPUT ARGUMENTS
% fiber_all: The fiber tracts from which selection will be made.  
%   The matrix could be the output of fiber_track or the smoothed fiber 
%   tracts output from fiber_smoother.
%
% angle_list, distance_list, curvature_list, n_points, apo_area: The outputs 
%   of fiber_quantifier.
%
% roi_flag: A mask indicating fiber tracts that propagated at least one point, 
%   output from fiber_track;
%
% fg_options: a structure containing the following fields:
%   -dwi_res: a three-element vector containing the field of view, matrix
%        size, and slice thickness of the diffusion-weighted images
%   -min_distance: minimum distance for selected tracts, in mm
%   -min_pennation: minimum pennation angle, in degrees
%   -max_pennation: maximum pennation angle, in degrees
%   -max_curvature: maximum curvature, in m^-1
%
%OUTPUT ARGUMENTS
%  final_fibers: the fiber tracts that passed all selection criteria
%
%  final_curvature: pointwise measurements of curvature for the final tracts.
%
%  final_angle: pointwise measurements of pennation angle for the final
%    tracts.
%
%  final_distance: pointwise measurements of cumulative distance for the
%    final tracts
%
%  qual_mask: a 3D matrix of the same row x column size as the roi_mesh, with 6
%    layers corresponding to each stage of the selection process. In each
%    layer, ones correspond to retained fibers and zeros correspond to
%    rejected fibers
%
%  num_tracked: the number of fibers for each of the following steps:
%    1) the number of fiber tracts generated by fiber_track;
%    2-6) the number of fiber tracts that met criteria 1-5 above, respectively.
%
%  mean_fiber_props: A 3D matrix (rows x columns x 5) containing the mean
%    curvature, pennation, and length values along each of the tracts; the
%    amount of aponeurosis area represented by each tract; and the number of
%    points in each tract
%
%  mean_apo_props: A 1 x 3 vector containing the whole-muscle mean values
%    for curvature, pennation, and fiber tract length, weighted by the
%    amount of aponeurosis area represented by each tract.
%
%OTHER FUNCTIONS IN THE MUSCLE DTI FIBER-TRACKING TOOLBOX
%  For help with anisotropic smoothing, see <a href="matlab: help aniso4D_smoothing">aniso4D_smoothing</a>.
%  For help calculating the diffusion tensor, see <a href="matlab: help signal2tensor2">signal2tensor2</a>.
%  For help defining the muscle mask, see <a href="matlab: help define_muscle">define_muscle</a>.
%  For help defining the aponeurosis ROI, see <a href="matlab: help define_roi">define_roi</a>.
%  For help with fiber tracking, see <a href="matlab: help fiber_track">fiber_track</a>.
%  For help smoothing fiber tracts, see <a href="matlab: help fiber_smoother">fiber_smoother</a>.
%  For help quantifying fiber tracts, see <a href="matlab: help fiber_quantifier">fiber_quantifier</a>.
%  For help visualizing fiber tracts and other structures, see <a href="matlab: help fiber_visualizer">fiber_visualizer</a>.
%
%VERSION INFORMATION
%  v. 1.0.0 (initial release), 17 Jan 2021, Bruce Damon
%  v. 1.1.0 (bug fix - corrected reporting of num_tracked fibers), 11 Aug 2021, Bruce Damon
%
%ACKNOWLEDGEMENTS
%  People: Zhaohua Ding, Anneriet Heemskerk
%  Grant support: NIH/NIAMS R01 AR050101, NIH/NIAMS R01 AR073831

%% get options from input structure
min_distance = fg_options.min_distance;
min_pennation = fg_options.min_pennation;
max_pennation = fg_options.max_pennation;
max_curvature = fg_options.max_curvature;

%% intialize output variables

fiber_indices_rows = sum(squeeze(angle_list(:,:,2)), 2);                    %find rows with fiber tracts
fiber_indices_rows = fiber_indices_rows>0;
fiber_indices_cols = sum(squeeze(angle_list(:,:,2)), 1);                    %find columns with fiber tracts
fiber_indices_cols = fiber_indices_cols>0;
first_row = find(fiber_indices_rows, 1);                           %find first and lastrow
last_row = find(fiber_indices_rows, 1, 'last');
first_col = find(fiber_indices_cols, 1);                           %find first and last column
last_col = find(fiber_indices_cols, 1, 'last');

final_fibers=zeros(size(fiber_all));                                        %initialize matrix of fiber fibers; 
final_fibers(first_row:last_row, first_col:last_col, :, :) = ...            %set as tracked fibers; wil prune erroneus results later
    fiber_all(first_row:last_row, first_col:last_col, :, :);
final_angle = angle_list;                                                   %same for geometric measurements of tracts
final_distance = distance_list;
final_curvature = curvature_list;

qual_mask = zeros([size(squeeze(fiber_all(:,:,1,1))) 6]);
qual_mask(:,:,1)=roi_flag;

%% initialize architecture output variables

mean_angle = sum(angle_list, 3)./squeeze(n_points(:,:,2));
mean_angle(isnan(mean_angle)) = 0;

mean_curvature = sum(curvature_list, 3)./squeeze(n_points(:,:,3));

total_distance = squeeze(max(distance_list, [], 3));

%% implement quality checking algorithm

%1) reject fibers that do not monotonically increase in the Z direction

z_positions=squeeze(fiber_all(:,:,:,3));
for row_cntr = first_row:last_row
    
    for col_cntr = first_col:last_col
        
        loop_z=squeeze(z_positions(row_cntr, col_cntr,:));                  % z positions for each fiber tract
        loop_z=nonzeros(loop_z);
        loop_dz=diff(loop_z);                                               %find differences between points
        loop_dz=loop_dz(1:(length(loop_dz)-1));
        
        if length(find(loop_dz<0))>0                                        %look for differences < 0 - indicates down-sloping fiber tracts
            qual_mask(row_cntr,col_cntr,1)=0;                               %write to quality mask
        end
        
    end
    
end
num_tracked(1) = length(find(qual_mask(:,:,1)>0));                          %count number passing through this criterion

%2) reject fibers that are too short
too_short = ones(size(total_distance));                                     %initialize as ones matrix
too_short(total_distance<min_distance) = 0;                                 %then find fiber tracts with length<minimum threshold
qual_mask(:,:,2) = qual_mask(:,:,1).*too_short;                             %write to quality mask
num_tracked(2) = length(find(qual_mask(:,:,2)>0));                          %count number passing through this criterion

%3) reject fibers that out of bounds pennation angle
angles_out_of_range = ones(size(mean_angle));                             	%initialize as ones matrix
angles_out_of_range(mean_angle <= min_pennation) = 0;                       %then find fiber tracts with pennation<minimum threshold
angles_out_of_range(mean_angle >= max_pennation) = 0;                       %or >maximum threshold
qual_mask(:,:,3) = qual_mask(:,:,2).*angles_out_of_range;                   %write to quality mask
num_tracked(3) = length(find(qual_mask(:,:,3)>0));                          %count number passing through this criterion

%4) reject fibers that have excessive curvature
high_curvature = ones(size(mean_curvature));
high_curvature(mean_curvature >= max_curvature) = 0;
qual_mask(:,:,4) = qual_mask(:,:,3).*high_curvature;
num_tracked(4) = length(find(qual_mask(:,:,4)>0));

%6) reject fibers that are very different (>2 SD) from their neighboring pixels in
%   length, curvature, or pennation angle
qual_mask(:,:,5) = qual_mask(:,:,4);                                        %initialize layers 5 and 6 as = index 4
qual_mask(:,:,6) = qual_mask(:,:,4);

for row_cntr = (first_row+2):(last_row-2)
    
    for col_cntr = (first_col+2):(last_col-2)
        
        if qual_mask(row_cntr,col_cntr,5) == 1
            
            row_neighbors = (row_cntr-2):(row_cntr+2);                      %find row indices for 24 nearest neighbors
            col_neighbors = (col_cntr-2):(col_cntr+2);                      %find column indices for 24 nearest neighbors
            local_fibers = qual_mask(row_neighbors,col_neighbors,4);        %get the indices from layer 4 of the matrix
            
            local_angle = mean_angle(row_neighbors, col_neighbors);         %get set of local angles
            local_angle = local_angle.*local_fibers;
            local_angle_non0 = local_angle(local_angle>0);
            mean_local_angle = mean(local_angle_non0);                      %get local mean and SD
            std_local_angle = std(local_angle_non0);
            
            if local_angle(3,3)>(mean_local_angle + 2*std_local_angle) || ...   %set outlier criteria
                    local_angle(3,3)<(mean_local_angle - 2*std_local_angle)
                qual_mask(row_cntr,col_cntr,6) = 0;                         %keep track of this in layer 6 so later results are unaffected
            end
            
            local_curve = mean_curvature(row_neighbors, col_neighbors);     %repeat for curvature
            local_curve = local_curve.*local_fibers;
            local_curve_non0 = local_curve(local_curve>0);
            mean_local_curve = median(local_curve_non0);
            std_local_curve = std(local_curve_non0);
            if local_curve(3,3)>(mean_local_curve + 2*std_local_curve) || ...
                    local_curve(3,3)<(mean_local_curve - 2*std_local_curve)
                qual_mask(row_cntr,col_cntr,6) = 0;
            end
            
            local_length = total_distance(row_neighbors, col_neighbors);    %and for distance
            local_length = local_length.*local_fibers;
            local_length_non0 = local_length(local_length>0);
            mean_local_length = median(local_length_non0);
            std_local_length = std(local_length_non0);
            if local_length(3,3)>(mean_local_length + 2*std_local_length) || ...
                    local_length(3,3)<(mean_local_length - 2*std_local_length)
                qual_mask(row_cntr,col_cntr,6) = 0;
            end
            
        end
        
        % eliminate tracts that failed the tests
        final_fibers(row_cntr,col_cntr,:,1) = final_fibers(row_cntr,col_cntr,:,1)*qual_mask(row_cntr,col_cntr,6);
        final_fibers(row_cntr,col_cntr,:,2) = final_fibers(row_cntr,col_cntr,:,2)*qual_mask(row_cntr,col_cntr,6);
        final_fibers(row_cntr,col_cntr,:,3) = final_fibers(row_cntr,col_cntr,:,3)*qual_mask(row_cntr,col_cntr,6);
        
        %eliminate architecture measurements for failed tracts
        final_curvature(row_cntr,col_cntr,:) = final_curvature(row_cntr,col_cntr,:)*qual_mask(row_cntr,col_cntr,6);
        final_angle(row_cntr,col_cntr,:) = final_angle(row_cntr,col_cntr,:)*qual_mask(row_cntr,col_cntr,6);
        final_distance(row_cntr,col_cntr,:) = final_distance(row_cntr,col_cntr,:)*qual_mask(row_cntr,col_cntr,6);
        
    end
    
end

% calculate output variables, then set index 6 of layer 3 back to zero
num_tracked(5) = length(find(qual_mask(:,:,6)>0));
num_tracked = [length(find(roi_flag)) num_tracked];
qual_final = qual_mask(:,:,6);                                                  %selected fibers
qual_final(isnan(mean_angle)) = 0;                                              %but eliminate NaN values
qual_final(isnan(mean_curvature)) = 0;
qual_final(isnan(total_distance)) = 0;
qual_mask(:,:,5) = qual_mask(:,:,6);                                            %rewrite to index five
qual_mask = qual_mask(:,:,1:5);                                                 %no longer need layer 6                                                    

%mean properties for each tract
mean_curvature = mean_curvature.*qual_mask(:,:,5);                              %get failed fiber tracts out of hte calculation of mean properities
mean_curvature(isnan(mean_curvature)) = 0;                                      %account for division by zero when initially calculating the mean
mean_fiber_props = mean_curvature;                                              %index 1 of third dimension has curvature

mean_angle = mean_angle.*qual_mask(:,:,5);
mean_angle(isnan(mean_angle)) = 0;
mean_fiber_props(:,:,2) = mean_angle;                                         	%index 2 of third dimension has pennation

total_distance = total_distance.*qual_mask(:,:,5);
mean_fiber_props(:,:,3) = total_distance;                                       %index 3 of third dimension has length

mean_fiber_props(:,:,4) = apo_area;                                            	%index 4 of third dimension has aponeurosis area represented by each tract

mean_fiber_props(:,:,5) = n_points(:,:,1).*qual_mask(:,:,5);                    %finally, number of ppints in each tract

%aponeurosis-wide properties - calculate weighted mean (weighted by relative aponeurosis area)
mean_apo_props(1) = sum(sum(mean_curvature.*apo_area.*qual_final))./sum(sum(apo_area.*qual_final));
mean_apo_props(2) = sum(sum(mean_angle.*apo_area.*qual_final))./sum(sum(apo_area.*qual_final));
mean_apo_props(3) = sum(sum(total_distance.*apo_area.*qual_final))./sum(sum(apo_area.*qual_final));

%% end the function

return;

