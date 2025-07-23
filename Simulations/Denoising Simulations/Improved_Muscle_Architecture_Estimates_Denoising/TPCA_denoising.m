function [Signal,n_comps] = TPCA_denoising(data, mask, kernel, noise_var)
    %
    % "TPCA_denoising": 4d image denoising by exploiting  data redundancy
    % in the PCA domain using threshold PCA denoising proposed by Henriques
    % et al., Imaging Neuroscience, 2023 
    %
    %  [Signal, n_comps] = TPCA_denoising(data, mask, kernel, sampling, centering, noise_var)
    %       output:
    %           - Signal: [x, y, z, M] denoised data matrix
    %           - n_comps: [x, y, z] number of non-noise pricipal
    %           components identified
    %       input:
    %           - data: [x, y, z, M] 4D DTI volume
    %           - mask:   (optional)  region-of-interest [boolean]
    %           - kernel: (optional)  window size for PCA analyses
    %           - noise_var: [x,y,z] data matrix containing the spatial
    %           noise variance estimates of the data 
    % 
    %  Adapted from:
    % MPPCAdenoising by Jelle Veraart (jelle.veraart@nyumc.org)
    % Copyright (c) 2016 New York Universit and University of Antwerp and
    % PCAdenoising by Rafael Neto Henriques
    % (https://github.com/RafaelNH/PCAdenoising)
    %       
    %      Permission is hereby granted, free of charge, to any non-commercial entity
    %      ('Recipient') obtaining a copy of this software and associated
    %      documentation files (the 'Software'), to the Software solely for
    %      non-commercial research, including the rights to use, copy and modify the
    %      Software, subject to the following conditions: 
    %       
    %        1. The above copyright notice and this permission notice shall be
    %      included by Recipient in all copies or substantial portions of the
    %      Software. 
    %       
    %        2. THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND,
    %      EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIESOF
    %      MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
    %      NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BELIABLE FOR ANY CLAIM,
    %      DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
    %      OTHERWISE, ARISING FROM, OUT OF ORIN CONNECTION WITH THE SOFTWARE OR THE
    %      USE OR OTHER DEALINGS IN THE SOFTWARE. 
    %       
    %        3. In no event shall NYU be liable for direct, indirect, special,
    %      incidental or consequential damages in connection with the Software.
    %      Recipient will defend, indemnify and hold NYU harmless from any claims or
    %      liability resulting from the use of the Software by recipient. 
    %       
    %        4. Neither anything contained herein nor the delivery of the Software to
    %      recipient shall be deemed to grant the Recipient any right or licenses
    %      under any patents or patent application owned by NYU. 
    %       
    %        5. The Software may only be used for non-commercial research and may not
    %      be used for clinical care. 
    %       
    %        6. Any publication by Recipient of research involving the Software shall
    %      cite the references listed below.
    %
    %  Final version written by: Roberto Pineda Guzman
    %  (roberto.guzman@carle.com)
    % 
    % REFERENCES
    %      Veraart, J.; Fieremans, E. & Novikov, D.S. Diffusion MRI noise mapping
    %      using random matrix theory Magn. Res. Med., 2016, early view, doi:
    %      10.1002/mrm.26059
    %      Henriques, R.N.; Ianus, A.; Novello, L.,; Jovicich, J.; Jespersen, S.N.; 
    %      Shemesh, N. Efficient PCA denoising of spatially correlated
    %      redundant MRI data. Imaging Neuroscience, 2023, doi:
    %      10.1162/imag_a_0049
    %      Manjon, J.V.; Coupe, P.; Concha, L.; Buades, A.; Collins, D.L.;
    %      Robles, M. Diffusion weighted image denoising using overcomplete
    %      local PCA. PLoS One, 2013, doi: 10.1371/journal.pone.0073021


 
    if isa(data,'integer') 
        data = single(data);
    end
    [sx, sy, sz, M] = size(data);

       
    if ~exist('mask', 'var') || isempty(mask)
        mask = true([sx, sy, sz]);
    end
    if ~isa(mask,'boolean') 
        mask = mask>0;
    end
  
    if ~exist('kernel', 'var') || isempty(kernel)
        kernel = [5 5 5];
    end
    
    if isscalar(kernel)
        kernel = [kernel, kernel, kernel];
    end
    kernel = kernel + (mod(kernel, 2)-1);   % needs to be odd.
    k = (kernel-1)/2; kx = k(1); ky = k(2); kz = k(3);
    N = prod(kernel);
    
    %warning('image boundaries are not processed.')
    mask(1:k(1), :, :) = 0; % zero mask boundaries where denoising can't be performed  
    mask(sx-k(1)+1:sx, :, :) = 0;
    mask(:, 1:k(2), :) = 0;
    mask(:, sy-k(2)+1:sy, :, :) = 0;           
    mask(:,:,1:k(3)) = 0;
    mask(:,:,sz-k(3)+1:sz) = 0;
    x = []; y = []; z = []; 
    for i = k(3)+1:sz-k(3) % create arrays containing all x,y, and z coordinates of denoised voxels
        [x_mask_slc, y_mask_slc] = find(mask(:,:,i) == 1);
        x = [x; x_mask_slc]; y = [y; y_mask_slc];  z = [z; i*ones(size(y_mask_slc))];
    end 
    x = x(:); y = y(:); z = z(:);

    
    % Declare variables:
    npars = zeros(1, numel(x), 'like', data);
    signal = zeros(M, prod(kernel), numel(x), 'like', data);
    %Signal = zeros(sx, sy, sz, M, 'like', data);
    weight_matrix = zeros(sx, sy, sz, M, 'like', data);
    weighted_signal = zeros(sx, sy, sz, M, 'like', data);
    n_comps = zeros(sx, sy, sz, 'like', data);

    % start denoising
    for nn = 1:numel(x)
        
        % create data matrix 
        X = data(x(nn)-kx:x(nn)+kx, y(nn)-ky:y(nn)+ky, z(nn)-kz:z(nn)+kz, :);
        X = reshape(X, N, M); X = X';

        colmean = mean(X, 1);
        X = X - repmat(colmean, [M, 1]);

        % compute PCA eigenvalues 
        [u, vals, v] = svd(X, 'econ');
        vals = diag(vals).^2 / N; % Scale to eigenvalues of covariance matrix 
        R = min(M, N);

        noise_var_mat = noise_var(x(nn)-kx:x(nn)+kx, y(nn)-ky:y(nn)+ky, z(nn)-kz:z(nn)+kz);
        noise_var_med = median(noise_var_mat,"all");
        
        % TPCA, Section 2.1.6, expression 3) from Henriques et al., Imaging
        % Neuroscience, 2023
        t = find(vals/((1+sqrt(M/N))^2) < noise_var_med, 1);
        
        if isempty(t)
            signal(:, :, nn) = X;  
            t = R+1;
        else
            vals(t:R) = 0;
            s = u*diag(sqrt(N*vals))*v';
            s = s + repmat(colmean, [M, 1]);

        
            signal(:, :, nn) = s;
        end
        npars(nn) = t-1; 
    end

    for nn = 1:numel(x)
        % Signal is weighted based on equation 3 in Manjon et al., 2013
        block_weight = 1/(1+M-npars(nn));
        n_comps(x(nn),y(nn),z(nn)) = npars(nn);
        weighted_signal(x(nn)-k(1):x(nn)+k(1),y(nn)-k(2):y(nn)+k(2),z(nn)-k(3):z(nn)+k(3), :) =...
            weighted_signal(x(nn)-k(1):x(nn)+k(1),y(nn)-k(2):y(nn)+k(2),z(nn)-k(3):z(nn)+k(3), :) + unpatch(signal(:,:,nn), k)*block_weight;
        weight_matrix(x(nn)-k(1):x(nn)+k(1),y(nn)-k(2):y(nn)+k(2),z(nn)-k(3):z(nn)+k(3), :) = ...
            weight_matrix(x(nn)-k(1):x(nn)+k(1),y(nn)-k(2):y(nn)+k(2),z(nn)-k(3):z(nn)+k(3), :) + block_weight*ones(kernel(1),kernel(2),kernel(3),M);
    end

    Signal = weighted_signal./weight_matrix;

end

function data = unpatch(X, k)
    kernel=k+k+1; 
    data = zeros([kernel, size(X, 1)]);
    tmp = zeros(kernel);
    for i = 1:size(X, 1)
        tmp(:) = X(i, :);
        data(:,:,:,i) = tmp;
    end 
end
