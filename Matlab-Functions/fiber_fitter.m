function [fitted_fiber_all, pcoeff_x, pcoeff_y, pcoeff_z] = fiber_fitter(fiber_all, ff_options)
%
%FUNCTION fiber_fitter
%  [fitted_fiber_all, pcoeff_x, pcoeff_y, pcoeff_z] = fiber_fitter(fiber_all, ff_options)
%
%USAGE
%    The function fiber_fitter is used to smooth fiber tracts and increase
%  the spatial resolution of fiber tracts generated using the 
%  MuscleDTI_Toolbox. The x, y, and z positions are separately fitted to Nth order 
%  polynomials as functions of distance along the tract and are uniformly 
%  solved at interpolation distances of interpolation_step. The user selects 
%  the polynomial order.  This procedure is modified from Damon et al, Magn
%  Reson Imaging, 2012 to: 
%    1) fit the tract positions as functions of distance rather than point
%       number. This is required for tracking algorithms that use variable
%       step sizes, such as FACT.
%    2) allow selection of the polynomial order, including different polynomial
%       orders for the X, Y, and X points .  
%  Future development will included automated selection of the most appropriate
%  polynomial order.
%
%INPUT ARGUMENTS 
%  fiber_all: The original fiber tracts, output from fiber_track
%
%  ff_options: A structure containing the following fields:
%    interpolation_step: an interpolation interval for the fitted fiber tract, in
%      units of pixels.  For example, setting interpolation_step to 0.25 would
%      interpolate the fiber tract at intervals of 0.25 pixels.
%
%    p_order: A 3 element vector containing the polynomial orders, [Nx Ny Nz],
%      to use when fitting the tracts. If a single number is entered (i.e., 
%      ff_options.p_order=N), then this order is used to fit the points in the 
%      X, Y and Z directions.
%
%OUTPUT ARGUMENTS
%  fitted_fiber_all: the fiber tracts following Nth order polynomial
%    fitting
%
%  pcoeff_x: a matrix of the Nth order polynomial coefficients for the
%    tracts' x positions 
%
%  pcoeff_y: a matrix of the Nth order polynomial coefficients for the
%    tracts' y positions 
%
%  pcoeff_z: a matrix of the Nth order polynomial coefficients for the
%    tracts' z positions 
%
%OTHER FUNCTIONS IN THE MUSCLE DTI FIBER-TRACKING TOOLBOX
%  For help defining the mask, see <a href="matlab: help define_muscle">define_muscle</a>.
%  For help defining the ROI, see <a href="matlab: help define_roi">define_roi</a>.
%  For help with the fiber tracking program, see <a href="matlab: help fiber_track">fiber_track</a>.
%  For help quantifying fiber tracts, see <a href="matlab: help fiber_quantifier">fiber_quantifier</a>.
%  For help selecting fiber tracts following their quantification, see <a href="matlab: help fiber_selector">fiber_selector</a>.
%  For help visualizing the data, see <a href="matlab: help fiber_visualizer">fiber_visualizer</a>.
%
%VERSION INFORMATION
%  v 0.1
%
%ACKNOWLEDGEMENTS
%  People: Zhaohua Ding, Anneriet Heemskerk
%  Grant support: NIH/NIAMS R01 AR050101, NIH/NIAMS R01 AR073831

%% prepare
interpolation_step=ff_options.interpolation_step;
p_order=ff_options.p_order;
if length(p_order)==1
    p_order=[p_order p_order p_order];
end

%initialize output variable as a zeros matrix
max_length = max(find(squeeze(sum(sum(squeeze(fiber_all(:,:,:,1))))))); %#ok<MXFND>
fitted_fiber_all = ...
    zeros(length(fiber_all(:,1,1,1)), length(fiber_all(1,:,1,1)), (max_length*ceil(1/interpolation_step)), 3);                  %zeros matrix to hold fitted fiber tracts
pcoeff_x = zeros(length(fiber_all(:,1,1,1)), length(fiber_all(1,:,1,1)),(p_order(1)+1));
pcoeff_y = zeros(length(fiber_all(:,1,1,1)), length(fiber_all(1,:,1,1)),(p_order(2)+1));
pcoeff_z = zeros(length(fiber_all(:,1,1,1)), length(fiber_all(1,:,1,1)),(p_order(3)+1));

%% fit each fiber tract

for row_cntr = 1:length(fiber_all(:,1,1,1))
    for col_cntr = 1:length(fiber_all(1,:,1,1))
        
        loop_fiber_length_points = length(find(fiber_all(row_cntr,col_cntr,:,1)));
        
        if loop_fiber_length_points>10
            
            fiber_distance = squeeze(fiber_all(row_cntr,col_cntr,1:loop_fiber_length_points, :));                           %x, y, and z positions
            fiber_distance(2:loop_fiber_length_points,1) = diff(fiber_distance(:,1));                                       %pointwise differences in x positions
            fiber_distance(2:loop_fiber_length_points,2) = diff(fiber_distance(:,2));
            fiber_distance(2:loop_fiber_length_points,3) = diff(fiber_distance(:,3));
            fiber_distance(1,:)=0;                                                                                          %initial point has distance of zed
            fiber_distance = cumsum((sum(fiber_distance.^2, 2)).^0.5);                                                      %calculate distances along fiber tract from initial point

            loop_fiber_x = squeeze(fiber_all(row_cntr,col_cntr,1:loop_fiber_length_points, 1));                         	%get raw tract data in x direction
            x_init=loop_fiber_x(1);
            loop_fiber_x=loop_fiber_x-x_init;
            pcoeff_x(row_cntr,col_cntr,:) = polyfit(fiber_distance, loop_fiber_x, p_order(1));                              %get polynomial coefficients
            loop_fitted_fiber_x = polyval(squeeze(pcoeff_x(row_cntr,col_cntr,:)), ...                                       %fitting in x
                min(fiber_distance):interpolation_step:max(fiber_distance));
            loop_fitted_fiber_x=loop_fitted_fiber_x+x_init;
            pcoeff_x(row_cntr,col_cntr,end)=pcoeff_x(row_cntr,col_cntr,end)+x_init;
            fitted_fiber_all(row_cntr,col_cntr,1:length(loop_fitted_fiber_x),1) = loop_fitted_fiber_x;                   	%copy to output variable
            
            loop_fiber_y = squeeze(fiber_all(row_cntr,col_cntr,1:loop_fiber_length_points, 2));                           	%get raw tract data in y direction
            y_init=loop_fiber_y(1);
            loop_fiber_y=loop_fiber_y-y_init;
            pcoeff_y(row_cntr,col_cntr,:) = polyfit(fiber_distance, loop_fiber_y, p_order(2));                              %get polynomial coefficients
            loop_fitted_fiber_y = polyval(squeeze(pcoeff_y(row_cntr,col_cntr,:)), ...                                       %fitting in y
                min(fiber_distance):interpolation_step:max(fiber_distance));  	
            loop_fitted_fiber_y=loop_fitted_fiber_y + y_init;
            pcoeff_y(row_cntr,col_cntr,end)=pcoeff_y(row_cntr,col_cntr,end) + y_init;
            fitted_fiber_all(row_cntr,col_cntr,1:length(loop_fitted_fiber_y),2) = loop_fitted_fiber_y;                      %copy to output variable
            
            loop_fiber_z = squeeze(fiber_all(row_cntr,col_cntr,1:loop_fiber_length_points, 3));                           	%get raw tract data in z direction
            z_init=loop_fiber_z(1);
            loop_fiber_z=loop_fiber_z - z_init;
            pcoeff_z(row_cntr,col_cntr,:) = polyfit(fiber_distance, loop_fiber_z, p_order(3));                              %get polynomial coefficients
            loop_fitted_fiber_z = polyval(squeeze(pcoeff_z(row_cntr,col_cntr,:)), ...                                       %fitting in z
                min(fiber_distance):interpolation_step:max(fiber_distance));
            loop_fitted_fiber_z = loop_fitted_fiber_z + z_init;
            pcoeff_z(row_cntr,col_cntr,end)=pcoeff_z(row_cntr,col_cntr,end) + z_init;
            fitted_fiber_all(row_cntr,col_cntr,1:length(loop_fitted_fiber_z), 3) = loop_fitted_fiber_z;                  	%copy to output variable
            
        end
    end
end

%% end function

return;
