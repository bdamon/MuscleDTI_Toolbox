function d_m = signal2tensor2(signal_v, dir_m, b)
% Usage: d_m = signal2tensor2(signal_v, dir_m, b)
%
%   signal2tensor2 finds the tensor that best fits the observed signal.
% The signal is given by the elements of the column vector signal_v,
% the first element of which is s0, the unweighted signal, and the
% other elements of which correspond to the directions specified by 
% the rows of dir_m (nDirs rows, by 3 columns for x, y, and z). The 
% b factor of all directions is assumed to be b. Note that the 
% units of d_m, the calculated tensor, will be 1/[b].
%   It is assumed that the standard deviation of each signal measurement
% is noiseStd.
%   signal2tensorX is a simplified version of dtiEngineX (useful for
% simulations and demonstrations).
%
% Last modified: 2006/12/04 (AWA)


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

 

