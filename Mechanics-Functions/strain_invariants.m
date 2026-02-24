function [along_fiber_strain,along_fiber_shear,cross_fiber_shear]=...
    strain_invariants(E1_field,C_cauchy_strain,muscle_mask)
%
%FUNCTION strain_invariants
% The function strain_invariants is used to compute the strain invariants 
% that describe the along-fiber strain, along-fiber shear, and cross-fiber
% shear measured in a muscle. These strain invariants can be found on 
% Criscione et al. 2001 (J Mech Phys Solids). See definitions on equations
% 1.3 and 5.3. Variable names follow naming scheme used in Criscione and 
% Blemker et al. 2005 (J Biomech). 
%
%INPUT ARGUMENTS
% E1_field: a 4D array (r x c x s x 3) that contains the 3 components of 
%  the principal diffusion eigenvector in each element of the muscle 3D
%  volume. This eigenvector field can be obtained from the preprocessing
%  scripts available in the MuscleDTI_Toolbox. The eigenvector 
%  components follow MATLAB's image coordinate system convention
% C_cauchy_strain: a 5D array (r x c x s x 3 x 3) that contains the right
%  Cauchy-Green strain tensor in each element of the muscle 3D volume. This
%  array can be obtained from the strain_tensors function. The tensor is
%  MATLAB's image coordinate system convention
% muscle_mask: a binary mask that specifies the muscle ROI, this mask is
%  required if the E1_field has already been masked, as [0,0,0] in the E1_field
%  will result in undefined values for some strain measures. 
%
%OUTPUT ARGUMENTS
% along_fiber_strain: along-fiber strain in each element of the muscle 3D volume
% along_fiber_shear: along-fiber shear in each element of the muscle 3D volume
% cross_fiber_shear: cross-fiber shear in each element of the muscle 3D volume
%
%VERSION INFORMATION
% v1.0.0, Roberto Pineda Guzman, October 7 2024.
%
% ACKNOWLEDGMENTS
% Grant support: NIH/NIAMS R01 AR073831

%%

if nargin<3
    muscle_mask=ones(size(E1_field,1:3)); % Include all volume if mask is not provided
end

along_fiber_stretch=zeros(size(muscle_mask));
along_fiber_shear=zeros(size(muscle_mask));
cross_fiber_shear=zeros(size(muscle_mask));

for s = 1:size(muscle_mask,3)
    for r = 1:size(muscle_mask,1)
        for c = 1:size(muscle_mask,2)
                % estimate strain maps for each voxel in the TA muscle
                if muscle_mask(r,c,s)>0
                    C = squeeze(C_cauchy_strain(r,c,s,:,:)); %Cauchy strain tensor of each muscle element
                    Csqr=C*C; %Square of Cauchy strain tensor
                    fiber_dir = squeeze(E1_field(r,c,s,:)); %Principal muscle fiber orientation in each muscle element
                    fiber_dir=[fiber_dir(1),fiber_dir(2),fiber_dir(3)]; 
                    I1 = trace(C); %1st invariant of Cauchy strain tensor, Criscione 1.3
                    I3 = det(C); %3rd invariant, Criscione 1.3
                    I4 = fiber_dir*C*fiber_dir'; %4th invariant, Criscione 1.3
                    I5 = fiber_dir*Csqr*fiber_dir'; %5th invariant, Criscione 1.3
                    along_fiber_stretch(r,c,s) = sqrt(I4); %along-fiber stretch, Blemker Eq. 3
                    along_fiber_shear(r,c,s) = sqrt((I5/(I4^2))-1); %along-fiber shear, Blemker Eq.3 and Criscione 5.3d
                    cross_fiber_shear(r,c,s) = acosh((I1*I4-I5)/(2*sqrt(I3*I4))); %cross-fiber shear, Criscione 5.3c
                end
        end
    end
end

along_fiber_strain=(along_fiber_stretch-1).*muscle_mask; %engineering strain along muscle fiber direction