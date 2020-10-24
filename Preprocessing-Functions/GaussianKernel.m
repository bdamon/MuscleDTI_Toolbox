function H = GaussianKernel(sz,sigma,sr)
% This function generates a Gaussian kernel
% input: 
%       sz -- [nx, ny, nz] the size of the kernel is 2nx + 1; 2ny + 1...
%       sigma -- the sigma for each dimenstion, note that it could be
%       vector
% output:
%       H -- the kernel

if ~(length(sz) == length(sigma))
    error('the size of sigma and sz does not match');
end
H = ones(2*sz+1);
for i = 1 : length(sz)
    x = [-sz(i):sz(i)]*sr(i);
    kernel_1d = (1/(2*pi*sigma(i)^2)^(1/2))*exp(-x.^2/(2*sigma(i)^2));
    kernel_1d_col = kernel_1d';
    sz_col = (2*sz+1);
    sz_col = circshift(sz_col,[1 i-1]);
    kernel_3d_outorder = repmat(kernel_1d_col,[1 sz_col(2:end)]);
    kernel_3d_inorder = shiftdim(kernel_3d_outorder,length(sz)+1-i);
    H = H.*kernel_3d_inorder;
end
H = H/sum(H(:));
