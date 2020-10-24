function img_U_s = aniso4D_smoothing_once(img_U,sigma,rho,delta_t,sr,type,isnormg,isfasteig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%take any DWI image img_U and smooth it using anisotropical diffusion%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Input:%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%img_U     --  DWI images, the size should be 4D(y,x,z,diffusion direction)
%sigma,rho --  variance of Gaussian kernel to smooth original image and
%structure tensor respectively.
%delta_t   --  time step
%sr        --  solution for the space variable
%type      --  numerical scheme is employed
%isnormg   --  whether the gradient in the PDE is normalized
%isfasteig --  whether the fast eigenvector solver is used
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%output%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%img_U_s   --  the smoother dwi images, the same size as img_U



%img_U = img_U*255;
% fakezero = 10^-10;
% tmpT2 = img_U(:,:,:,1);
% img_U = img_U(:,:,:,2:end);
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