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
return;