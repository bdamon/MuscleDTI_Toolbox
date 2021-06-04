function d_m = signal2tensor2(signal_v, dir_m, b)
%
% FUNCTION signal2tensor2
%  d_m = signal2tensor2(signal_v, dir_m, b);
%
% USAGE
%  signal2tensor2 is used to estimate the diffusion tensor that best fits 
%  the observed signal.
%
%  The user inputs a column vector containing unweighted and diffusion-
%  weighted signals, a matrix describing the diffusion-encoding  
%  directions; and the b-value (assumed to be the same for all directions).
%  The function outputs the calculated tensor, with units of 1/[b]. 
%
% INPUT ARGUMENTS
%  signal_v: A column vector of image signals, the first element of which is
%   the unweighted signal, and the other elements of which correspond to the
%   directions specified by the rows of dir_m 
%
%  dir_m: A matrix containing the X, Y, and Z components of unit-length vectors 
%   describing the diffusion-sensitizing directions; it has dimensions of
%   (Number of Directions x 3)
%
%  b: The diffusion-weighting (b-) value
%
% OUTPUT ARGUMENTS
%  d_m: The diffusion tensor.
%
% OTHER FUNCTIONS IN THE MUSCLE DTI FIBER-TRACKING TOOLBOX
%  For help visualizing the data, see <a href="matlab: help fiber_visualizer">fiber_visualizer</a>.
%  For help defining the mask, see <a href="matlab: help define_muscle">define_muscle</a>.
%  For help defining the aponeurosis ROI, see <a href="matlab: help define_roi">define_roi</a>.
%  For help with the fiber tracking program, see <a href="matlab: help fiber_track">fiber_track</a>.
%  For help fitting fiber tracts, see <a href="matlab: help fiber_fitter">fiber_fitter</a>.
%  For help quantifying fiber tracts, see <a href="matlab: help fiber_quantifier">fiber_quantifier</a>.
%  For help selecting fiber tracts following their quantification, see <a href="matlab: help fiber_goodness">fiber_goodness</a>.
%
% VERSION INFORMATION
%  v. 1.0 2006/12/04 Adam Anderson
%  v. 1.0.1 adds help information; 17 Jan 2021, Bruce Damon



%% start function

% Parameters:
noiseStd = 1;   % Normalized/doesn't matter.

% Get measurement directions:
nDirs = size(dir_m, 1);

% Find design matrix for b=1 (to improve condition number):
design_m = zeros(nDirs+1, 7);   % Fit 7 parameters. 
design_m(1,:) = [1, zeros(1, 6)];
for dir = 1:nDirs
    n_v = dir_m(dir, :);
    b_m = n_v.' * n_v;
    design_m(1+dir, :) = [1, -b_m(1,1), -b_m(2,2), -b_m(3,3), ...
            -2*b_m(1,2), -2*b_m(1,3), -2*b_m(2,3)];
end

% Calculate ln of the signal and its std (column vectors):
ln_v = log(signal_v);
stdLn_v = ( noiseStd ./ signal_v );

% Find weighted least-squares solution:
w_m = diag(1 ./ stdLn_v);   % Matrix of weights.
wDesign_m = w_m * design_m; % Weighted design matrix.

% Find solution to wDesign_m * a_v = w_m * ln_v
a_v = (wDesign_m.' * wDesign_m) \ wDesign_m.' * w_m * ln_v;
lnS0 = a_v(1);  % Not used.

% Form diffusion tensor matrix and scale by actual b (a_v = b*d_v):
d_m = [a_v(2), a_v(5), a_v(6); ...
       a_v(5), a_v(3), a_v(7); ...
       a_v(6), a_v(7), a_v(4)] / b;

return; 

 

