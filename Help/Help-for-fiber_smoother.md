%
%FUNCTION fiber_smoother
%  [smoothed_fiber_all, pcoeff_x, pcoeff_y, pcoeff_z] = fiber_smoother(fiber_all, smooth_options)
%
%USAGE
%    The function fiber_smoother is used to smooth fiber tracts and increase
%  the spatial resolution of fiber tracts generated using the 
%  MuscleDTI_Toolbox. The x, y, and z positions are fitted to Nth order 
%  polynomials as functions of distance along the tract and are uniformly 
%  solved at interpolation distances of interpolation_step. The user selects 
%  the polynomial order.  This procedure is modified from Damon et al, Magn
%  Reson Imaging, 2012 to: 1) fit the tract positions as functions of distance
%  rather than point number and 2) allow selection of the polynomial order.  
%  The former option is required for tracking algorithms that use variable 
%  step sizes.
%
%INPUT ARGUMENTS 
%  fiber_all: the original fiber tracts, output from fiber_track
%
%  smooth_options: a structure containing the following fields:
%    interpolation_step: an interpolation interval for the fitted fiber tract, in
%      units of pixels.  For example, setting interpolation_step to 0.25 would
%      interpolate the fiber tract at intervals of 0.25 pixels.
%
%    p_order: a 3 element vector containing the polynomial orders, [Nx Ny Nz],
%      to use when fitting the tracts
%
%OUTPUT ARGUMENTS
%  smoothed_fiber_all: the fiber tracts following Nth order polynomial
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
%  In beta-testing mode
%
%ACKNOWLEDGEMENTS
%  People: Zhaohua Ding, Anneriet Heemskerk
%  Grant support: NIH/NIAMS R01 AR050101, NIH/NIAMS R01 AR073831
