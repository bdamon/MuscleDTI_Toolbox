function [angle_list, distance_list, curvature_list, fiber_all_mm, n_points, apo_area] = fiber_quantifier(fiber_all, roi_mesh, fq_options)
% 
% FUNCTION fiber_quantifier
%  [angle_list, distance_list, curvature_list, fiber_all_mm, n_points, apo_area] = ...
%    fiber_quantifier(fiber_all, roi_mesh, fq_options)
% 
% USAGE
%    The function fiber_quantifier is used to calculate the muscle architectural
%  parameters pennation angle, fiber tract length, and curvature in the 
%  MuscleDTI_Toolbox. 
%
%  Calculations are only made for fiber tracts having six or more points. 
%  Information about each measurement follows:
%    -Fiber tract length: this is measured by summing the inter-point 
%     distances along the tract.
%    -Pennation: The method for pennation measurements is essentially as 
%     described in Lansdown et al, J Appl Physiol 2007. The approach to 
%     measuring pennation angle traditionally used in ultrasound imaging is 
%     to manually specify the line tangent to a muscle fascicle at the point 
%     of its insertion into the aponeurosis and a second line tangent to the 
%     aponeurosis. The angle formed by these lines is measured. In Lansdown 
%     et al., this concept was extended into 3D space by defining the plane 
%     tangent to the seed point and its normal vector. Also, position vectors 
%     between the seed point and points along the tract were defined. 
%     Pennation angle was defined as the complement to the angle formed by 
%     the normal vector and the position vectors. In the toolbox, this is 
%     modified slightly to calculated the normal vector as the cross product
%     between two tangent lines to the seed point (one lying in the row 
%     direction of the aponeurosis mesh and one lying in the column direction). 
%     This change improves computational efficiency.
% 
%    -Curvature: The method for curvature measurements is described in Damon 
%     et al, Magn Reson Imaging 2012. Briefly, these use a discrete 
%     implementation of the Frenet-Serret equations. Specifically, the 
%     curvature K is defined in
%       dT/ds = K N
%     where T is the tangent line to points along the curve, s is the step 
%     length between points, and N is the normal vector. In fiber_quantifier, 
%     K is calculated by multiplying each side of this equation by the Moore-
%     Penrose pseudoinverse matrix of N.
% 
%     For curvature, the best results are obtained with polynomial-fitted 
%     fiber tracts, calculated using fiber_fitter. Fiber tract length and 
%     pennation angle are unaffected by polynomial fitting.
% 
% INPUT ARGUMENTS
%  fiber_all: A 4D matrix containing the fiber tract points, with units of
%    pixels (in X and Y) or slice number (in Z).  
% 
%  roi_mesh: The mesh reconstruction of the aponeurosis that was used as the
%    seed surface for fiber tracking, output from define_roi.
% 
%  fq_options: A user-defined structure containing the following fields:
%    dwi_res: a three element vector with the FOV, (assumed to be the same for
%      the x and y directions), in-plane matrix size, and the slice thickness
%       of the DTI images. The FOV and slice thickness must be specified in mm.
% 
%    filt_kernel: Pennation angle is calculated as complement to the angle
%      formed by position vectors along the tract and the normal vector to
%      the seed point of the tract. To do so, two lines tangent to the point
%      are calculated. When determining the lines, their slopes are median-
%      filtered over the NxN window specified in filt_kernel; enter this as 
%      a single odd integer, e.g. fq_options.filt_kernel=N.  Note that median 
%      filtering results in a loss of edge information and that the number of 
%      lost rows and columns increases with the size of the filter kernel.
% 
% OUTPUT ARGUMENTS
%  angle_list: The pennation angles, for fiber tracking point numbers
%    2-end. Pennation angles are reported in degrees.
% 
%  distance_list: Distances along the fiber tracts, for fiber tracking point
%    numbers 1-end. Distances are reported in mm.
% 
%  curvature_list: The curvature values, for fiber tracking point numbers
%    starting at 2 and ending 3 points before the tract's end. The latter is to 
%    avoid abrupt changes in curvature.  Curvature values are reported in m-1.
% 
%  fiber_all_mm: The fiber tract points converted to units of mm.
% 
%  n_points: A 3D matrix (rows x columns x 3) containing the number of
%    points used to quantify length, pennation, and curvature in each tract
% 
%  apo_area: A 2D matrix (rows x columns) containing the amount of
%    apoeurosis area represented by each fiber tract
% 
% OTHER FUNCTIONS IN THE MUSCLE DTI FIBER-TRACKING TOOLBOX
%  For help visualizing the data, see <a href="matlab: help fiber_visualizer">fiber_visualizer</a>.
%  For help defining the mask, see <a href="matlab: help define_muscle">define_muscle</a>.
%  For help defining the ROI, see <a href="matlab: help define_roi">define_roi</a>.
%  For help with the fiber tracking program, see <a href="matlab: help fiber_track">fiber_track</a>.
%  For help smoothing fiber tracts, see <a href="matlab: help fiber_smoother">fiber_smoother</a>.
%  For help selecting fiber tracts following their quantification, see <a href="matlab: help fiber_goodness">fiber_goodness</a>.
% 
% VERSION INFORMATION
%  v 0.1
% 
% ACKNOWLEDGEMENTS
%  People: Zhaohua Ding, Adam Anderson, Anneriet Heemskerk
%  Grant support: NIH/NIAMS R01 AR050101, NIH/NIAMS R01 AR073831

%% get options out of structure

dwi_res = fq_options.dwi_res;
filt_kernel = fq_options.filt_kernel;

%% calculate fiber lengths in points

fiber_length = squeeze(fiber_all(:,:,:,1));
fiber_length(fiber_length>0)=1;
fiber_length=sum(fiber_length, 3);

%% initialize output variables

angle_list = zeros(size(squeeze(fiber_all(:,:,:,1))));
distance_list = zeros(size(squeeze(fiber_all(:,:,:,1))));
curvature_list = zeros(size(squeeze(fiber_all(:,:,:,1))));
n_points = zeros([size(squeeze(fiber_all(:,:,1,1))) 3]);
apo_area = zeros(size(squeeze(roi_mesh(:,:,1))));

%% convert spatially dependent matrices from pixels to mm

%account for spatial resolution of the images by adjusting the mesh:
dwi_slicethickness = dwi_res(3);
dwi_fov = dwi_res(1);
dwi_xsize = dwi_res(2);

roi_mesh_mm(:,:,1) = roi_mesh(:,:,1)*(dwi_fov/dwi_xsize);                                   %switch row column to X-Y frame of reference
roi_mesh_mm(:,:,2) = roi_mesh(:,:,2)*(dwi_fov/dwi_xsize);
roi_mesh_mm(:,:,3) = roi_mesh(:,:,3)*dwi_slicethickness;


% convert fiber tracts to mm
fiber_all_mm(:,:,:,1) = squeeze(fiber_all(:,:,:,1))*(dwi_fov/dwi_xsize);                    %convert fiber tracts in pixels to mm
fiber_all_mm(:,:,:,2) = squeeze(fiber_all(:,:,:,2))*(dwi_fov/dwi_xsize);
fiber_all_mm(:,:,:,3) = squeeze(fiber_all(:,:,:,3))*dwi_slicethickness;


%% prepare for pennation measurements

%find the slopes of the tangent lines in the through-plane direction:
dx_dy = zeros(length(roi_mesh_mm(:,1,1))-1, length(roi_mesh_mm(1,:,1))-1);                  %form a zeros matrix to hold dx/dz
dy_dz = zeros(length(roi_mesh_mm(:,1,1))-1, length(roi_mesh_mm(1,:,1))-1);                  %form a zeros matrix to hold dy/dz

for row_cntr = 1:length(roi_mesh_mm(:,1,1))
    inplane_x = squeeze(roi_mesh_mm(row_cntr,:, 2));                                        %find the x positions
    inplane_y = squeeze(roi_mesh_mm(row_cntr,:, 1));                                        %find the y positions
    dx_dy(row_cntr,:) = diff(inplane_x)./diff(inplane_y);                                   %calculate approximate derivative as delta X / delta Y
end

for col_cntr = 1:length(roi_mesh_mm(1,:,1))
    thruplane_y = squeeze(roi_mesh_mm(:, col_cntr, 1));                                     %find the y positions
    thruplane_z = squeeze(roi_mesh_mm(:, col_cntr, 3));                                     %find the z positions
    dy_dz(:,col_cntr) = diff(thruplane_y)./diff(thruplane_z);                               %calculate approximate derivative as delta Y / delta Z
end

%median filter the slopes using user-defined kernel size
dx_dy = medfilt2(dx_dy, [filt_kernel filt_kernel]);
dy_dz = medfilt2(dy_dz, [filt_kernel filt_kernel]);


%% architectural measurements

%median filtering eliminates edge information, so adjust for the kernel
%size when looping through the mesh for architectural measurements
start_row = 1 + floor(filt_kernel/2);
start_col = 1 + floor(filt_kernel/2);
end_row = length(roi_mesh(:,1,1)) - floor(filt_kernel/2);
end_col = length(roi_mesh(1,:,1)) - floor(filt_kernel/2);

% begin the architecture measurement loops
for row_cntr = start_row:end_row
    
    for col_cntr = start_col:end_col
        
        %aponeurosis area calculations
        p0=squeeze(roi_mesh_mm(row_cntr, col_cntr,1:3));                                    %p0-3 are the four points of the quadrilateral on the roi mesh
        p1=squeeze(roi_mesh_mm(row_cntr+1, col_cntr,1:3));                                  %p0-3 are the four points of the quadrilateral on the roi mesh
        p2=squeeze(roi_mesh_mm(row_cntr+1, col_cntr+1,1:3));                                %p0-3 are the four points of the quadrilateral on the roi mesh
        p3=squeeze(roi_mesh_mm(row_cntr, col_cntr+1,1:3));                                  %p0-3 are the four points of the quadrilateral on the roi mesh
        apo_area(row_cntr, col_cntr) = ...                                                  %model the quadrilateral as two irregular triangles with vertices v0-v1-v2 and v2-v3-v0
            norm(cross((p1-p0), (p2-p1)))/2 + norm(cross((p3-p2), (p3-p0)))/2;              %  their area is 1/2 of the magnitude of the indicated cross products
            
        if fiber_length(row_cntr,col_cntr) >= 6                                           	%quantify all tracts of at least 6 points

            %find initial point
            seed_point = squeeze(fiber_all_mm(row_cntr, col_cntr, 1, :));
            x0 = seed_point(2);
            y0 = seed_point(1);
            z0 = seed_point(3);
            
            %get the tangent line in the in-plane direction from the matrix
            %calculated above
            local_dxdy = dx_dy(row_cntr,col_cntr);                                          %local value of the slope, dx/dz
            local_xy_intrcpt = x0 - local_dxdy*y0;                                          %calculate the intercept by rearranging x = b + (dx/dz)*z to solve for the intercept b
            
            %get the tangent line in the through-plane direction:
            local_dydz = dy_dz(row_cntr,col_cntr);                                          %local value of the slope, dy/dz
            local_yz_intrcpt = y0 - local_dydz*z0;                                      	%calculate the intercept by rearranging y = b + (dy/dz)*z to solve for the intercept b
            
            %find two other points in the tangent plane
            y1 = y0 - 1;
            z1 = z0;
            x1 = local_xy_intrcpt + local_dxdy*y1;
            
            x2 = x0;
            z2 = z0 - 1;
            y2 = local_yz_intrcpt + local_dydz*z2;
            
            %calculate the Cartesian equation for the plane and the unit normal vector
            point_0 = seed_point;
            point_1 = [y1 x1 z1]';
            point_2 = [y2 x2 z2]';
            dP1 = point_1-point_0;
            dP2 = point_2-point_0;
            normal_vector = cross(dP1, dP2);
            unit_normal_vector = normal_vector/norm(normal_vector);
            
            %find the fiber length, pennation angle, and curvature as a function of position along the fiber:
            for fiber_cntr = 2:fiber_length(row_cntr,col_cntr)
                
                %fiber length
                delta_p = squeeze(fiber_all_mm(row_cntr, col_cntr, fiber_cntr, :)) - ...
                    squeeze(fiber_all_mm(row_cntr,col_cntr, fiber_cntr-1, :));             	%form position vector from preceding point in the tract
                distance_list(row_cntr,col_cntr,fiber_cntr)  =  ...                         %Euclidean distance from preceding point + the preceding point
                    sqrt(sum(delta_p.^2)) + distance_list(row_cntr,col_cntr,fiber_cntr-1);
                
                %pennation measurements for each point on the tract
                r_vector = squeeze(fiber_all_mm(row_cntr, col_cntr, fiber_cntr, :)) - ...   %form position vectors along the tract
                    seed_point;
                r_unit = r_vector/norm(r_vector);                                          	%r_unit is the unit vector
                theta_degrees = asind(dot(r_unit, unit_normal_vector));                    	%pennation angle is the complement to that formed by the normal vector and the position vector
                if theta_degrees <0 || theta_degrees >90                                    %account for direction of normal vector
                    theta_degrees = asind(dot(r_unit, -unit_normal_vector));
                end
                
                angle_list(row_cntr,col_cntr,fiber_cntr) = theta_degrees;                  	%write to matrix
                if isnan(angle_list(row_cntr,col_cntr,fiber_cntr))                         	%in case of dividing by zero
                    angle_list(row_cntr,col_cntr) = 0;
                end
                
                %keep track of number of points quantified.
                n_points(row_cntr,col_cntr,1) = fiber_cntr;                                 %In 3rd dimension, save #points for distance data in level 1
                n_points(row_cntr,col_cntr,2) = fiber_cntr-1;                               %In 3rd dimension, save #points for pennation data in level 2; 
                                                                                            %use fiber_cntr-1 because measurements start at point #2
                
                if fiber_cntr>2 && fiber_cntr<(fiber_length(row_cntr,col_cntr)-3)           %curvature values blow up at the end of the fiber because the rest of the vector is padded with 0's
                    
                    %curvature measurements use a discrete implementation of the Frenet equations.
                    p1_idx = fiber_cntr-1;                                                  %indices for the three points of interest along the tract - define two pairs of points
                    p2_idx = fiber_cntr;
                    p3_idx = fiber_cntr+1;
%                     loop_fiber_m = squeeze(fiber_all_mm(row_cntr,col_cntr,:,:))/1000;      	%convert from mm to m for curvature measurements
                    
                    delta_p1 = (loop_fiber_m(p1_idx,:))-loop_fiber_m(1,:);                       %three position vectors, one for each point
                    delta_p2 = (loop_fiber_m(p2_idx,:))-loop_fiber_m(1,:);
                    delta_p3 = (loop_fiber_m(p3_idx,:))-loop_fiber_m(1,:);
                    
                    ds21 = sqrt(sum((loop_fiber_m(p2_idx,:)-loop_fiber_m(p1_idx,:)).^2));   %distance between points 1 and 2 and (below) 2 and 3
                    ds32 = sqrt(sum((loop_fiber_m(p3_idx,:)-loop_fiber_m(p2_idx,:)).^2));
                    
                    tangent_vector_2 = (delta_p2-delta_p1)/norm(delta_p2-delta_p1);                           	%normalized tangent lines between the two pairs of points
                    tangent_vector_3 = (delta_p3-delta_p2)/norm(delta_p3-delta_p2);
                    dTds = ((tangent_vector_3-tangent_vector_2)/mean([ds21 ds32]))';      	%dT/ds is the spatial rate of change in tangent lines
                    dTds(isnan(dTds)) = 0;
                    
                    if sum(dTds) ~= 0
                        N_vector = dTds/norm(dTds);                                         %normal to tangent lines
                        curvature_list(row_cntr,col_cntr,fiber_cntr) = pinv(N_vector)*dTds;	%based on dT/ds = curvature * N
                    end
                    
                    n_points(row_cntr,col_cntr,3) = fiber_cntr-1;                           %In 3rd dimension, save # points used for curvature data in level 3
                    
                end                                                                         %of curvature if statement
                
            end                                                                             %of inside fiber loop
            
        end                                                                                 %of outside fiber loop
        
    end                                                                                     %of column loop
    
end                                                                                         %of row loop


%% end the function

return;