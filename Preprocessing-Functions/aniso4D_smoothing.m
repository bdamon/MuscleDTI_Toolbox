function img_U_s = aniso4D_smoothing(img_U,sigma,rho,delta_t,sr,type,isnormg,isfasteig)
%
%FUNCTION aniso4D_smoothing
%  img_U_s = aniso4D_smoothing(img_U,sigma,rho,delta_t,sr,type,isnormg,isfasteig);
%
%USAGE
%  The function aniso4D_smoothing is used to smooth diffusion-weighted images.
%
% The user provides the diffusion weighted images and options that control 
% the degree of smoothing. The images are smoothed using the method described 
% by Ding et al. (Ding et al., 2005, Magn Reson Med 53:485-90). This method 
% builds on an anisotropic filtering method described by Weickert et al. 
% that smooths images with the partial differential equation
%   dI/dt=div(T delI)
% Where I is the image intensity, T is the structure tensor that provides 
% smoothing anisotropy, and dt is the iteration time parameter. T is the 
% normalized inverse of the gradient tensor G. G is the convolution of the 
% outer product of delI with a Gaussian kernel having a standard deviation rho.
% The anisotropic smoothing method uses a common definition of G for all 
% diffusion-weighting directions, allowing boundary information that is 
% missing in one weighting direction to be captured from another weighting 
% direction. Thus, the algorithm provides isotropic smoothing inside of 
% structures and anisotropic smoothing at their boundaries. The smoothed 
% images are returned
%
%INPUT ARGUMENTS
%  img_U: The images to be smoothed
%  sigma, rho: The variance of the Gaussian kernel used to smooth the 
%    original images and the structure tensor, respectively
%  delta_t: The time step.
%  sr: The spatial resolution of the images, input as a two-element vector 
%    containing the in-plane resolution and the slice thickness.
%  type: The type of smoothing that is used, input as a string variable. The 
%    options include ‘explicit’, ‘implicit_multisplitting’, ‘implicit_AOS’, 
%    ‘implicit_ADI’, and ‘implicit_AOS_improved’
%  isnormg: determines whether the gradient in the partial differential 
%    equation is normalized; set to either true or false
%  isfasteig: determines whether the fast eigenvector solver is used; set 
%    to either true or false
%
%OUTPUT ARGUMENT
%  img_U_s: The smoothed images 
%
%OTHER FUNCTIONS IN THE MUSCLE DTI FIBER-TRACKING TOOLBOX
%  For help calculating the diffusion tensor, see <a href="matlab: help signal2tensor2">signal2tensor2</a>.
%  For help visualizing the data, see <a href="matlab: help fiber_visualizer">fiber_visualizer</a>.
%  For help defining the mask, see <a href="matlab: help define_muscle">define_muscle</a>.
%  For help defining the ROI, see <a href="matlab: help define_roi">define_roi</a>.
%  For help smoothing fiber tracts, see <a href="matlab: help fiber_smoother">fiber_smoother</a>.
%  For help quantifying fiber tracts, see <a href="matlab: help fiber_quantifier">fiber_quantifier</a>.
%  For help selecting fiber tracts following their quantification, see <a href="matlab: help fiber_goodness">fiber_goodness</a>.
%
%VERSION INFORMATION
%  v. 1.0.0 Zhaohua Ding, 2005
%  v. 1.0.1 Combines all functions into a single m-file; adds help comments. 17 Jan 2021, Bruce Damon
%  v, 1.1.0 Corrects a bug (SolveTriangleTomasBatch function was missing). 10 Aug 2021, Bruce Damon, Zhaohua Ding, and Carly Lockard
%
%ACKNOWLEDGEMENTS
%  People: Adam Anderson, John Gore
%  Grant support: NIH/NIBIB R01 EB000461, NIH/NIBIB R01 EB02777

%%


[img_height, img_width, img_slice, img_type] = size(img_U);

% Smooth image iteratively
Gx = zeros(img_height, img_width, img_slice, img_type);
Gy = zeros(img_height, img_width, img_slice, img_type);
Gz = zeros(img_height, img_width, img_slice, img_type);
J_o = zeros(img_height, img_width, img_slice, img_type, 3, 3);
J_p = zeros(img_height, img_width, img_slice, 3, 3);
NeVal = zeros(3, 1);
D = zeros(img_height, img_width, img_slice, 3, 3);
[x, y, z] = meshgrid((1:img_width)*sr(1), (1:img_height)*sr(2), (1:img_slice)*sr(3));
G_sigma = GaussianKernel([1 1 1],[sigma sigma sigma],sr);
G_rho   = GaussianKernel([1 1 1],[rho rho rho],sr);

% Convolve G_sigma with raw image
for t=1:img_type
    U_sigma(:, :, :, t) = convn(img_U(:,:,:,t), G_sigma, 'same');
    
    % Get gradient info
    [Gx(:, :, :, t), Gy(:, :, :, t), Gz(:, :, :, t)] = gradient(U_sigma(:,:,:,t), sr(1), sr(2), sr(3));
    if isnormg
        G_norm = sqrt(Gx(:,:,:,t).^2 + Gy(:,:,:,t).^2 + Gz(:,:,:,t).^2)+fakezero;
        Gx(:, :, :, t) = Gx(:, :, :, t)./G_norm;
        Gy(:, :, :, t) = Gy(:, :, :, t)./G_norm;
        Gz(:, :, :, t) = Gz(:, :, :, t)./G_norm;
    end
    % Construct structure tensor
    J_o(:, :, :, t, 1, 1) = Gx(:, :, :, t).^2;
    J_o(:, :, :, t, 1, 2) = Gx(:, :, :, t).*Gy(:, :, :, t);
    J_o(:, :, :, t, 1, 3) = Gx(:, :, :, t).*Gz(:, :, :, t);
    J_o(:, :, :, t, 2, 1) = J_o(:, :, :, t, 1, 2);
    J_o(:, :, :, t, 2, 2) = Gy(:, :, :, t).^2;
    J_o(:, :, :, t, 2, 3) = Gy(:, :, :, t).*Gz(:, :, :, t);
    J_o(:, :, :, t, 3, 1) = J_o(:, :, :, t, 1, 3);
    J_o(:, :, :, t, 3, 2) = J_o(:, :, :, t, 2, 3);
    J_o(:, :, :, t, 3, 3) = Gz(:, :, :, t).^2;
end

% Smooth structure tensor with G_rho
for m=1:3
    for n=1:3
        J_p(:, :, :, m, n) =  convn(squeeze(sum(J_o(:, :, :, :, m, n), 4)), G_rho, 'same');
    end
end

if isfasteig
    [eVec,eVal,mask] = eigss_fast(J_p);
    idx = find(mask==0);
    for m = idx'
        [iy,ix,iz] = ind2sub(size(mask),m);
        [eVec(iy,ix,iz,:,:),eValTmp] = eigs(squeeze(J_p(iy, ix, iz, 1:3, 1:3)), 3);
        eVal(iy,ix,iz,:) = diag(eValTmp);
    end
    
    
    mask1 =  eVal(:,:,:,1)<0.1^30;
    mask2 = (eVal(:,:,:,2)<0.1^30)&(~mask1);
    mask3 = (eVal(:,:,:,3)<0.1^30)&(~mask2)&(~mask1);
    mask4 = (~mask1)&(~mask2)&(~mask3);
    idx1 = find(mask1==1);
    idx2 = find(mask2==1);
    idx3 = find(mask3==1);
    idx4 = find(mask4==1);
    eVal = reshape(eVal,[img_height*img_width*img_slice 3]);
    eVal(idx1,:) = 1;
    eVal(idx2,1) = 0;
    eVal(idx2,2:3) = 1.5;
    eVal(idx3,1:2) = 0;
    eVal(idx3,3) = 3;
    eVal(idx4,:) = 1./eVal(idx4,:);
    eVal(idx4,:) = 3*eVal(idx4,:)./repmat(sum(eVal(idx4,:),2),[1 3]);
    eVal = reshape(eVal,[img_height img_width img_slice 3]);
    
    eVal = reshape(eVal,[img_height img_width img_slice 1 3]);
    eVal = repmat(eVal,[1 1 1 3 1]);
    S = eVec.*eVal;
    D(:,:,:,1,1) = squeeze(sum(S(:,:,:,1,:).*eVec(:,:,:,1,:),5));
    D(:,:,:,1,2) = squeeze(sum(S(:,:,:,1,:).*eVec(:,:,:,2,:),5));
    D(:,:,:,1,3) = squeeze(sum(S(:,:,:,1,:).*eVec(:,:,:,3,:),5));
    D(:,:,:,2,2) = squeeze(sum(S(:,:,:,2,:).*eVec(:,:,:,2,:),5));
    D(:,:,:,2,3) = squeeze(sum(S(:,:,:,2,:).*eVec(:,:,:,3,:),5));
    D(:,:,:,3,3) = squeeze(sum(S(:,:,:,3,:).*eVec(:,:,:,3,:),5));
    D(:,:,:,2,1) = D(:,:,:,1,2);
    D(:,:,:,3,1) = D(:,:,:,1,3);
    D(:,:,:,3,2) = D(:,:,:,2,3);
else
    % construct new structure tensors used in the smoothing
    for i=1:img_height
        for j=1:img_width
            for k=1:img_slice
                
                [eVec, eVal] = eigs(squeeze(J_p(i, j, k, 1:3, 1:3)), 3);
                eVal = diag(eVal);
                %[eVec, eVal] = eigss(squeeze(J_p(i, j, k, 1:3, 1:3)));
                if  eVal(1) == 0
                    NeVal(1:3) = [1 1 1];
                    %NeVal(1:3) = [3 3 3];
                elseif  eVal(2) == 0
                    NeVal(1:3) = [0 1.5 1.5];
                    %NeVal(1:3) = [0 3 3];
                elseif  eVal(3) == 0
                    NeVal(1:3) = [0 0 3];
                    %NeVal(1:3) = [0 0 3];
                else
                    NeVal(1:3) = 1./eVal;
                    NeVal = 3*NeVal/sum(NeVal);
                    %NeVal = 3*NeVal/max(NeVal);
                end
                D(i, j, k, 1:3, 1:3) = eVec*[NeVal(1) 0 0; 0 NeVal(2) 0; 0 0 NeVal(3)]*eVec';
            end
        end
    end
    
end

%explicit scheme
switch type
    case 'explicit'
        for t=1:img_type
            img_U(:,:,:,t) = img_U(:,:,:,t) + ...
                delta_t*divergence(x,y,z, ...
                squeeze(D(:,:,:,1,1)).*Gx(:,:,:,t) + squeeze(D(:,:,:,1,2)).*Gy(:,:,:,t) + squeeze(D(:,:,:,1,3)).*Gz(:,:,:,t),...
                squeeze(D(:,:,:,2,1)).*Gx(:,:,:,t) + squeeze(D(:,:,:,2,2)).*Gy(:,:,:,t) + squeeze(D(:,:,:,2,3)).*Gz(:,:,:,t),...
                squeeze(D(:,:,:,3,1)).*Gx(:,:,:,t) + squeeze(D(:,:,:,3,2)).*Gy(:,:,:,t) + squeeze(D(:,:,:,3,3)).*Gz(:,:,:,t));
        end
        
    case 'implicit_multisplitting'
        img_U_tmp = img_U;
        for t=1:img_type
            %explicit
            img_U(:,:,:,t) = img_U(:,:,:,t) + ...
                delta_t*divergence(x,y,z, ...
                squeeze(D(:,:,:,1,2)).*Gy(:,:,:,t) + squeeze(D(:,:,:,1,3)).*Gz(:,:,:,t),...
                squeeze(D(:,:,:,2,1)).*Gx(:,:,:,t)                                      + squeeze(D(:,:,:,2,3)).*Gz(:,:,:,t),...
                squeeze(D(:,:,:,3,1)).*Gx(:,:,:,t) + squeeze(D(:,:,:,3,2)).*Gy(:,:,:,t)                                    );
            %implicit AOS
            img_U(:,:,:,t) = SolveTriangleTomasBatch(squeeze(img_U(:,:,:,t)),squeeze(D(:,:,:,1,1)),delta_t,1,2);
            img_U(:,:,:,t) = SolveTriangleTomasBatch(squeeze(img_U(:,:,:,t)),squeeze(D(:,:,:,2,2)),delta_t,1,1);
            img_U(:,:,:,t) = SolveTriangleTomasBatch(squeeze(img_U(:,:,:,t)),squeeze(D(:,:,:,3,3)),delta_t,1,3);
            %img_U(:,:,:,t) = img_U(:,:,:,t) + (img_Ux + img_Uy + img_Uz);
        end
        
    case 'implicit_AOS'
        img_U_tmp = img_U;
        for t=1:img_type
            %explicit
            img_U(:,:,:,t) = img_U(:,:,:,t) + ...
                delta_t*divergence(x,y,z, ...
                squeeze(D(:,:,:,1,2)).*Gy(:,:,:,t) + squeeze(D(:,:,:,1,3)).*Gz(:,:,:,t),...
                squeeze(D(:,:,:,2,1)).*Gx(:,:,:,t)                                      + squeeze(D(:,:,:,2,3)).*Gz(:,:,:,t),...
                squeeze(D(:,:,:,3,1)).*Gx(:,:,:,t) + squeeze(D(:,:,:,3,2)).*Gy(:,:,:,t)                                    );
            %implicit AOS
            img_Uz = SolveTriangleTomasBatch(squeeze(img_U(:,:,:,t)),squeeze(D(:,:,:,3,3)),3*delta_t,1,3);
            img_Uy = SolveTriangleTomasBatch(squeeze(img_U(:,:,:,t)),squeeze(D(:,:,:,2,2)),3*delta_t,1,1);
            img_Ux = SolveTriangleTomasBatch(squeeze(img_U(:,:,:,t)),squeeze(D(:,:,:,1,1)),3*delta_t,1,2);
            img_U(:,:,:,t) = (img_Ux + img_Uy + img_Uz)/3;
        end
        
    case 'implicit_ADI'
        for t=1:img_type
            img_U(:,:,:,t) = ParobolicEqnSolver3dOnce(squeeze(img_U(:,:,:,t)),D,delta_t,1,1,1);
        end
        
    case 'implicit_AOS_improved'
        for t=1:img_type
            %explicit
            img_U(:,:,:,t) = img_U(:,:,:,t) + ...
                delta_t*divergence(x,y,z, ...
                squeeze(D(:,:,:,1,2)).*Gy(:,:,:,t) + squeeze(D(:,:,:,1,3)).*Gz(:,:,:,t),...
                squeeze(D(:,:,:,2,1)).*Gx(:,:,:,t)                                      + squeeze(D(:,:,:,2,3)).*Gz(:,:,:,t),...
                squeeze(D(:,:,:,3,1)).*Gx(:,:,:,t) + squeeze(D(:,:,:,3,2)).*Gy(:,:,:,t)                                    );
            %implicit AOS
            img_Ux = SolveTriangleTomasBatch(squeeze(img_U(:,:,:,t)),squeeze(D(:,:,:,1,1)).*Gx(:,:,:,t) + squeeze(D(:,:,:,1,2)).*Gy(:,:,:,t) + squeeze(D(:,:,:,1,3)).*Gz(:,:,:,t)./squeeze(Gx(:,:,:,t)),delta_t*3,1,2);
            img_Uy = SolveTriangleTomasBatch(squeeze(img_U(:,:,:,t)),squeeze(D(:,:,:,2,1)).*Gx(:,:,:,t) + squeeze(D(:,:,:,2,2)).*Gy(:,:,:,t) + squeeze(D(:,:,:,2,3)).*Gz(:,:,:,t)./squeeze(Gy(:,:,:,t)),delta_t*3,1,1);
            img_Uz = SolveTriangleTomasBatch(squeeze(img_U(:,:,:,t)),squeeze(D(:,:,:,3,1)).*Gx(:,:,:,t) + squeeze(D(:,:,:,3,2)).*Gy(:,:,:,t) + squeeze(D(:,:,:,3,3)).*Gz(:,:,:,t)./squeeze(Gz(:,:,:,t)),delta_t*3,1,3);
            img_U(:,:,:,t) = (img_Ux + img_Uy + img_Uz)/3;
        end
end
%img_U = img_U/255.0;
% img_U_s = zeros(img_height, img_width, img_slice, img_type+1);
% img_U_s(:,:,:,2:end) = img_U;
% img_U_s(:,:,:,1) = tmpT2;
img_U_s = img_U;

return;


%% eigss_fast

function [evector,evalue,mask]= eigss_fast(D)
%optimization based on Matlab analytical eigenvalue&eigenvector solution "Analytical Computation of the Eigenvalues and Eigenvectors in DT-MRI"
%input: D 3Dx3x3 matrix, if you wanna input a 2D image,please input it as a
%3D array.
%output: evector: the column vectors are eigenvectors
%        evalue:  a row vector, eigenvalue
sz = size(D);
sx = sz(length(sz)-1);
sy = sz(length(sz));
fakezero = 0.1^30;

if (sx==3)
    if (sy==3)
        % I didnot test symmetr here,because I donot wanna introduce any
        % additoanl computation.
        
        %determination of eigenvalue
        I1 = D(:,:,:,1,1)+D(:,:,:,2,2)+D(:,:,:,3,3);
        I2 = (D(:,:,:,1,1).*D(:,:,:,2,2)+D(:,:,:,1,1).*D(:,:,:,3,3)+D(:,:,:,2,2).*D(:,:,:,3,3))-(D(:,:,:,1,2).^2+D(:,:,:,1,3).^2+D(:,:,:,2,3).^2);
        I3 = D(:,:,:,1,1).*D(:,:,:,2,2).*D(:,:,:,3,3)+D(:,:,:,1,2).*D(:,:,:,2,3).*D(:,:,:,3,1)+D(:,:,:,1,3).*D(:,:,:,2,1).*D(:,:,:,3,2)...
            -D(:,:,:,1,3).*D(:,:,:,2,2).*D(:,:,:,3,1)-D(:,:,:,1,2).*D(:,:,:,2,1).*D(:,:,:,3,3)-D(:,:,:,1,1).*D(:,:,:,2,3).*D(:,:,:,3,2);
        
        v = (I1/3).^2-I2/3+fakezero;
        s = (I1/3).^3-I1.*I2/6+I3/2;
        sita = acos(s./v.*sqrt(1./v))/3;
        evalue = zeros([sz(1:3) 3]);
        evalue(:,:,:,1) = I1/3+2*sqrt(v).*cos(sita);
        evalue(:,:,:,2) = I1/3-2*sqrt(v).*cos(pi/3*ones(sz(1:3))+sita);
        evalue(:,:,:,3) = I1-evalue(:,:,:,1)-evalue(:,:,:,2);
        %build the mask for degenerate and non-pos case
        mask = ones(sz(1:3));
        mask = mask.*(I3>fakezero).*(squeeze(D(:,:,:,1,1))>=fakezero).*(squeeze(D(:,:,:,2,2))>=fakezero).*(squeeze(D(:,:,:,3,3))>=fakezero)...
            .*((squeeze(D(:,:,:,1,1)).*squeeze(D(:,:,:,2,2))-squeeze(D(:,:,:,1,2)).^2)>=fakezero)...
            .*((squeeze(D(:,:,:,1,1)).*squeeze(D(:,:,:,3,3))-squeeze(D(:,:,:,1,3)).^2)>=fakezero)...
            .*((squeeze(D(:,:,:,2,2)).*squeeze(D(:,:,:,3,3))-squeeze(D(:,:,:,2,3)).^2)>=fakezero)...
            .*~((abs(squeeze(D(:,:,:,1,3)))<0.1^19).*(abs(squeeze(D(:,:,:,2,3)))<0.1^19))...
            .*~((abs(squeeze(D(:,:,:,1,2)))<0.1^19).*(abs(squeeze(D(:,:,:,1,3)))<0.1^19))...
            .*~((abs(squeeze(D(:,:,:,1,2)))<0.1^19).*(abs(squeeze(D(:,:,:,2,3)))<0.1^19))...
            .*(abs(s./v.*sqrt(1./v))<=1).*(v>0);
        %sum(sum(sum(mask)))
        %sz(1)*sz(2)*sz(3)
        %is it probabatily a degenrate case
        degentype = abs(evalue(:,:,:,1)-evalue(:,:,:,2))<0.1^5;
        % determine the eigenvector
        A = repmat(D(:,:,:,1,1),[1 1 1 2]) - evalue(:,:,:,1:2);
        B = repmat(D(:,:,:,2,2),[1 1 1 2]) - evalue(:,:,:,1:2);
        C = repmat(D(:,:,:,3,3),[1 1 1 2]) - evalue(:,:,:,1:2);
        evector = zeros(sz);
        evector(:,:,:,1,1:2) =((repmat(D(:,:,:,1,2).*D(:,:,:,2,3),[1 1 1 2])-B.*repmat(D(:,:,:,1,3),[1 1 1 2]))...
            .*(repmat(D(:,:,:,1,3).*D(:,:,:,2,3),[1 1 1 2])-C.*repmat(D(:,:,:,1,2),[1 1 1 2])));
        
        evector(:,:,:,2,1:2) =((repmat(D(:,:,:,1,3).*D(:,:,:,2,3),[1 1 1 2])-C.*repmat(D(:,:,:,1,2),[1 1 1 2]))...
            .*(repmat(D(:,:,:,1,3).*D(:,:,:,1,2),[1 1 1 2])-A.*repmat(D(:,:,:,2,3),[1 1 1 2])));
        
        evector(:,:,:,3,1:2) =((repmat(D(:,:,:,1,2).*D(:,:,:,2,3),[1 1 1 2])-B.*repmat(D(:,:,:,1,3),[1 1 1 2]))...
            .*(repmat(D(:,:,:,1,3).*D(:,:,:,1,2),[1 1 1 2])-A.*repmat(D(:,:,:,2,3),[1 1 1 2])));
        %norm of the first 2 eigenvectors
        evector(:,:,:,:,1) = evector(:,:,:,:,1)./(repmat(sqrt(sum(evector(:,:,:,:,1).^2,4)),[1 1 1 3])+fakezero);
        evector(:,:,:,:,2) = evector(:,:,:,:,2)./(repmat(sqrt(sum(evector(:,:,:,:,2).^2,4)),[1 1 1 3])+fakezero);
        %use the tensor product of the first 2 eigenvectors to generate the
        %3rd eigenvector
        evector(:,:,:,1,3) = evector(:,:,:,2,1).*evector(:,:,:,3,2)-evector(:,:,:,2,2).*evector(:,:,:,3,1);
        evector(:,:,:,2,3) = evector(:,:,:,3,1).*evector(:,:,:,1,2)-evector(:,:,:,3,2).*evector(:,:,:,1,1);
        evector(:,:,:,3,3) = evector(:,:,:,1,1).*evector(:,:,:,2,2)-evector(:,:,:,1,2).*evector(:,:,:,2,1);
        
        
        A = repmat(D(:,:,:,1,1),[1 1 1 2]) - evalue(:,:,:,[1 3]);
        B = repmat(D(:,:,:,2,2),[1 1 1 2]) - evalue(:,:,:,[1 3]);
        C = repmat(D(:,:,:,3,3),[1 1 1 2]) - evalue(:,:,:,[1 3]);
        
        evector1 = zeros(sz);
        evector1(:,:,:,1,[1 3]) =((repmat(D(:,:,:,1,2).*D(:,:,:,2,3),[1 1 1 2])-B.*repmat(D(:,:,:,1,3),[1 1 1 2]))...
            .*(repmat(D(:,:,:,1,3).*D(:,:,:,2,3),[1 1 1 2])-C.*repmat(D(:,:,:,1,2),[1 1 1 2])));
        evector1(:,:,:,2,[1 3]) =((repmat(D(:,:,:,1,3).*D(:,:,:,2,3),[1 1 1 2])-C.*repmat(D(:,:,:,1,2),[1 1 1 2]))...
            .*(repmat(D(:,:,:,1,3).*D(:,:,:,1,2),[1 1 1 2])-A.*repmat(D(:,:,:,2,3),[1 1 1 2])));
        evector1(:,:,:,3,[1 3]) =((repmat(D(:,:,:,1,2).*D(:,:,:,2,3),[1 1 1 2])-B.*repmat(D(:,:,:,1,3),[1 1 1 2]))...
            .*(repmat(D(:,:,:,1,3).*D(:,:,:,1,2),[1 1 1 2])-A.*repmat(D(:,:,:,2,3),[1 1 1 2])));
        
        evector1(:,:,:,:,1) = evector1(:,:,:,:,1)./(repmat(sqrt(sum(evector1(:,:,:,:,1).^2,4)),[1 1 1 3])+fakezero);
        evector1(:,:,:,:,3) = evector1(:,:,:,:,3)./(repmat(sqrt(sum(evector1(:,:,:,:,3).^2,4)),[1 1 1 3])+fakezero);
        
        
        evector1(:,:,:,1,2) = evector1(:,:,:,2,1).*evector1(:,:,:,3,3)-evector1(:,:,:,2,3).*evector1(:,:,:,3,1);
        evector1(:,:,:,2,2) = evector1(:,:,:,3,1).*evector1(:,:,:,1,3)-evector1(:,:,:,3,3).*evector1(:,:,:,1,1);
        evector1(:,:,:,3,2) = evector1(:,:,:,1,1).*evector1(:,:,:,2,3)-evector1(:,:,:,1,3).*evector1(:,:,:,2,1);
        
        evector = evector.*repmat(~degentype,[1 1 1 3 3])+evector1.*repmat(degentype,[1 1 1 3 3]);
        
        return;
        
    end
end
return

%% GaussianKernel

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

return


%% SolveTriangleTomasBatch

function U_next = SolveTriangleTomasBatch(U,G,dt,h,dim);

%Do the triangular batching job
sz = size(U);
U = shiftdim(U,dim-1);
G = shiftdim(G,dim-1);
sz = size(U);
U = reshape(U,[sz(1) sz(2)*sz(3)]);
G = reshape(G,[sz(1) sz(2)*sz(3)]);
U = U';
G = G';
for i=1:sz(2)*sz(3)
     U(i,:) = SolveTriangleTomas(squeeze(U(i,:)),squeeze(G(i,:)),dt,h);
%    U(i,:) = SolveTriangleImplicit(squeeze(U(i,:)),squeeze(G(i,:)),dt,h);
end
U = U';
G = G';
U = reshape(U,[sz(1) sz(2) sz(3)]);
G = reshape(G,[sz(1) sz(2) sz(3)]);
U = shiftdim(U,4-dim);
G = shiftdim(G,4-dim);
U_next = U;

return

%% SolveTriangleTomas

function u_next = SolveTriangleTomas(u,g,dt,h)
%Solve (I-dt*A(g))*u_next = u,where u is just a vector with boundary
%condition on two ends

sz = size(u);
d = u;
g = g*dt;
%%%%%%%%%%%slove the linear system using Tomas%%%%%%%%%%%%%%%%%%%%%%
alpha2 = ((g(1)+g(2))/(2*h^2))/(1+(g(1)+g(2))/(2*h^2));
alphan = ((g(end)+g(end-1))/(2*h^2))/(1+(g(end)+g(end-1))/(2*h^2));;
beta2 = d(1)/(1+(g(1)+g(2))/(2*h^2));
betan = d(end)/(1+(g(end)+g(end-1))/(2*h^2));

% %%%%%%%%%%%slove the linear system using Tomas%%%%%%%%%%%%%%%%%%%%%%
% alpha2 = ((g(1)+g(2))/(h^2))/(1+(g(1)+g(2))/(h^2));
% alphan = ((g(end)+g(end-1))/(h^2))/(1+(g(end)+g(end-1))/(h^2));;
% beta2 = d(1)/(1+(g(1)+g(2))/(h^2));
% betan = d(end)/(1+(g(end)+g(end-1))/(h^2));


n = length(d);
%solving the p and q
p = [alpha2];
q = [beta2];
for i=2:n-1
    a = (g(i-1)+g(i))/(2*h^2);
    c = (g(i)+g(i+1))/(2*h^2);
    b = -(a+c);
    a = -a;
    c = -c;
    b = 1-b;
    p_next = -c/(a*p(i-1)+b);
    q_next = (d(i)-a*q(i-1))/(a*p(i-1)+b);
    p = [p p_next];
    q = [q q_next];
end

% sloving the w
w = (alphan*q(end)+betan)/(1-alphan*p(end));

for i=n-1:-1:1
    w_prev = w(1)*p(i)+q(i);
    w = [w_prev w];
end

u_next = w;

return


