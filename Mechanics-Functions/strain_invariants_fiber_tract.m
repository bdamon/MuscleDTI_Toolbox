function [along_fiber_strain,along_fiber_shear,cross_fiber_shear]=...
    strain_invariants_fiber_tract(smooth_fibers,C_cauchy_strain,spacing)
%
%FUNCTION strain_invariants_fiber_tract
% The function strain_invariants_fiber_tract is used to compute the 
% strain invariants that describe the along-fiber strain, along-fiber shear,
% and cross-fibershear measured in the DTI-derived fiber tracts of a muscle. 
% These strain invariants can be found on Criscione et al. 2001 (J Mech Phys Solids). 
% See definitions on equations 1.3 and 5.3. Variable names follow naming scheme
% used in Criscione and Blemker et al. 2005 (J Biomech). 
%
%INPUT ARGUMENTS
% smooth_fibers: a 4D array that contains the row, col, and slc components of 
%  DTI-based fiber tracts from a seed point mesh. This structure is
%  obtained from the fiber_track or fiber_smoother functions of the Muscle 
%  DTI Toolbox. Even though the fiber_tracks output can be used in this
%  function, the smooth tracts obtained from the fiber_smoother function
%  are recommended because strain measurements can be sensitive to noise
%  found in the nonsmoothed fiber-tracts. 
% C_cauchy_strain: a 5D array (r x c x s x 3 x 3) that contains the right
%  Cauchy-Green strain tensor in each element of the muscle 3D volume. This
%  array can be obtained from the strain_tensors function. The tensor is
%  specified in MATLAB's image coordinate system.
% spacing: a 3 x 1 array that includes the physical spacing between the
%  elements of the array in the x,y, and z directions. If the spacing is
%  not specified, an isotropic spacing of 1x1x1 units is assumed
%
%OUTPUT ARGUMENTS
% along_fiber_strain: along-fiber strain in each point of the fiber tracts
% along_fiber_shear: along-fiber shear in each point of the fiber tracts
% cross_fiber_shear: cross-fiber shear in each point of the fiber tracts
%
%VERSION INFORMATION
% v1.0.0, Roberto Pineda Guzman, September 5 2024.
%
% ACKNOWLEDGMENTS
% Grant support: NIH/NIAMS R01 AR073831

%%

if nargin<3
    spacing=[1 1 1]; %Assume isotropic voxels unless previously specified
end

% compute the along-fiber strain of each point in the fiber tracts
start_row = 1;
start_col = 1;
end_row = size(smooth_fibers,1);
end_col = size(smooth_fibers,2);

% calculate fiber lengths in points
fiber_length = squeeze(smooth_fibers(:,:,:,1));
fiber_length(fiber_length>0)=1;
fiber_length=sum(fiber_length, 3);

% initialize variables
along_fiber_strain=zeros(end_row,end_col,size(smooth_fibers,3));
along_fiber_shear=zeros(end_row,end_col,size(smooth_fibers,3));
cross_fiber_shear=zeros(end_row,end_col,size(smooth_fibers,3));

% begin the strain measurement loops
for row_cntr = start_row:end_row
    
    for col_cntr = start_col:end_col
        
        %Measure the tangent vector of each point in the fiber-tracts using discretized versions of the Frenet equations.
        loop_fiber_m = squeeze(smooth_fibers(row_cntr,col_cntr,:,:));
        [~,delta_fiber_position]=gradient(loop_fiber_m); %obtain change in position vector of each point to compute the point's tangent vector
        delta_fiber_position(:,1)=delta_fiber_position(:,1)*spacing(1)/min(spacing); %adjust for anisotropic voxel size
        delta_fiber_position(:,2)=delta_fiber_position(:,2)*spacing(2)/min(spacing);
        delta_fiber_position(:,3)=delta_fiber_position(:,3)*spacing(3)/min(spacing);

        % find the along-fiber strain of each tract
        for fiber_cntr = 1:fiber_length(row_cntr,col_cntr)
            fiber_dir = delta_fiber_position(fiber_cntr,:)/norm(delta_fiber_position(fiber_cntr,:)); %compute tangent vector of each point in the fiber tract
            fiber_dir = [fiber_dir(2),fiber_dir(1),fiber_dir(3)]; % Order vector in x-y-z (from row-col-slc in smooth_fibers array)
            fiber_voxel_position = round(squeeze(smooth_fibers(row_cntr,col_cntr,fiber_cntr,:)));
            C = squeeze(C_cauchy_strain(fiber_voxel_position(1),fiber_voxel_position(2),fiber_voxel_position(3),:,:)); %Cauchy strain of point in the fiber tract
            Csqr = C*C; %Square of Cauchy strain tensor
            I1 = trace(C); %1st invariant of Cauchy strain tensor, Criscionne 1.3
            I3 = det(C); %3rd invariant, Criscionne 1.3
            I4 = fiber_dir*C*fiber_dir'; %4th invariant, Criscionne 1.3
            I5 = fiber_dir*Csqr*fiber_dir'; %5th invariant, Criscionne 1.3
            along_fiber_strain(row_cntr,col_cntr,fiber_cntr) = sqrt(I4)-1; %along-fiber stretch, Blemker Eq. 3
            along_fiber_shear(row_cntr,col_cntr,fiber_cntr) = sqrt((I5/(I4^2))-1); %along-fiber shear, Blemker Eq.3 and Criscionne 5.3d
            cross_fiber_shear(row_cntr,col_cntr,fiber_cntr) = acosh((I1*I4-I5)/(2*sqrt(I3*I4))); %cross-fiber shear, Criscionne 5.3c
        end
    end
end