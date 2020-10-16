function fiber_figure = fiber_visualizer(anat_image, fv_options, roi_mesh, mask, fiber_all)
%
%FUNCTION fiber_visualizer
%  fiber_figure = fiber_visualizer(anat_image, plot_options, roi_mesh, mask, fiber_all)
%
%USAGE
%  The function fiber_visualizer is used to visualize images and the muscle
%  mask, roi mesh, and/or fiber tracts formed using the MuscleDTI_Toolbox.
%
%  The user can call fiber_visualizer from the command line.  The user must  
%  supply theanatomical images, a structure with some plotting options, and  
%  the other variables to be plotted as input arguments. In addition, 
%  define_muscle, define_roi, and fiber_track can be configured to call 
%  fiber_visualizer from within the functions, so that the mask, mesh, and 
%  fiber tracts can be automatically plotted.  Fields of view, matrix sizes, 
%  slice thickness, etc. are appropriately considered so that all structures 
%  are plotted using a consistent measurement scale.
%
%INPUT ARGUMENTS
%  anat_image: The stack of images to be plotted for anatomical reference
%
%  fv_options: A structure containing the following required fields:
%   -anat_dims: A two element vector containing the FOV and the slice
%      thickness of the anatomical images.
%   -anat_slices: A vector containing the slice numbers of the anatomical
%      images to be plotted.
%   -plot_mesh: If set to 1, this field will allow plotting of the aponeurosis 
%       mesh. Otherwise, set to 0.
%   -plot_mask: If set to 1, this field will allow plotting of the mask.
%        Otherwise, set to 0.
%   -plot_fibers:  If set to 1, this field will allow plotting of a single
%      set of fiber tracts. If set to 2, this will allow plotting of two
%      sets of fiber tracts. Otherwise, set to 0.
%
%  Depending on the plot options selected, the following other fields may
%  be required:
%
%  If plot_mesh equals 1, you must also specify:
%   -mesh_size: This two-element vector specifies the in-plane matrix size 
%      of the images used to generate the mesh.
%   -mesh_dims: This two-element vector specifies the FOV and slice thickness 
%      of the images used to create the mesh.
%   -mesh_color: The user creates a color scale using either of the following 
%      two options:
%      *If mesh_color is a 3 element vector of values ranging from 0-1, 
%        the vector is interpreted as RGB levels; the entire mesh is plotted
%        using this color.
%      *If mesh_color is a matrix with size of (#mesh rows) x (#mesh columns) x 3, 
%        and if these values range from 0-1, the matrix will be interpreted  
%        as RGB levels specific to each tract. This could be used to
%        represent the distribution of architectural parameters across the
%        aponeurosis
%   -mesh_dist: If the mesh was shifted for fiber tracking, the user should
%      set this to the value used during fiber tracking.
%
%  If plot_mask equals 1, you must also specify:
%   -mask_size: This two-element vector specifies the in-plane matrix size 
%      of the images used to generate the mask.
%   -mask_dims: This two-element vector specifies the FOV and slice thickness
%      of the images used to create the mask.
%   -mask_color: This three-element vector contains the RGB levels to be used 
%      when plotting the mask 
%
%  If plot_fibers equals 1 or 2, you must also specify:
%    -dti_size: A2-element vector that specifies the matrix size of the images
%       used for fiber tracking.
%    -dti_dims: This two-element vector specifies the FOV and slice thickness
%       of the DT images. The FOV is assumed to be square.
%    -fiber_color: This defines the color of the tracts. Several options are
%       available:
%        *If plot_fibers equals 1 and fiber_color is a 3 element vector of values
%          ranging from 0-1, the vector is interpreted as RGB levels.
%        *If plot_fibers equals 2 and fiber_color is a 2x3 matrix of values
%          ranging from 0-1, each row of the matrix is interpreted as RGB
%          levels for the respective sets of tracts.
%        *If plot_fibers equals 1 and fiber_color is a matrix with size of
%          (#mesh rows) x (#mesh columns) x 3, and if these values range from
%          0-1, the third dimension of the matrix will be interpreted as RGB
%          levels specific to each tract
%        *If plot_fibers equals 2 and fiber_color is a matrix with size of
%          (#mesh rows) x (#mesh columns) x 3 x 2, and if these values range from
%          0-1, the third dimension of the matrix will be interpreted as RGB
%          levels specific to each tract, separately for sets 1 and 2
%    -fiber_skip: Setting fiber_skip to integer values > 1 will skip over
%      fiber tracts when plotting. This may improve visualization and will
%      decrease time for rendering. If not specified, all fibers will be
%      plotted.
%
%  roi_mesh: The roi mesh, defined in define_roi. It is only needed if
%    fv_options.plot_mesh is set to 1. Otherwise, enter an empty matrix
%    ([ ]) as a placeholder.
%
%  mask: A binary mask around the muscle of interest. It could be the
%    output of define_muscle or it could have been defined in another
%    program. It is only needed if fv_options.plot_mask is set to 1. 
%    Otherwise, enter an empty matrix ([ ]) as a placeholder.
%
%  fiber_all: The output of fiber_track (original fiber tracts), fiber_smoother 
%    (smoothed fiber tracts), or fiber_goodness (quality-selected fiber tracts).
%    If plot_fibers equals 1, fiber_all should have size = (#mesh rows) x 
%    (#mesh columns) x (#fiber tract points) x 3. If plot_fibers equals 1, 
%    fiber_all should have size =(#mesh rows) x (#mesh columns) x 
%    (#fiber tract points) x 3 x 2. It is only needed if fv_options.plot_fibers
%    is set to 1. Otherwise, enter an empty matrix ([ ]) as a placeholder.
%
%OUTPUT ARGUMENTS
%  fiber_figure: A Matlab figure structure
%
%OTHER FUNCTIONS IN THE MUSCLE DTI FIBER-TRACKING TOOLBOX
%  For help defining the muscle mask, see <a href="matlab: help define_muscle">define_muscle</a>.
%  For help defining the aponeurosis ROI, see <a href="matlab: help define_roi">define_roi</a>.
%  For help with fiber tracking, see <a href="matlab: help fiber_track">fiber_track</a>.
%  For help smoothing fiber tracts, see <a href="matlab: help fiber_smoother">fiber_smoother</a>.
%  For help quantifying fiber tracts, see <a href="matlab: help fiber_quantifier">fiber_quantifier</a>.
%  For help selecting fiber tracts following their quantification, see <a href="matlab: help fiber_selector">fiber_selector</a>.
%
%VERSION INFORMATION
%  v. 0.1
%
%ACKNOWLEDGEMENTS
%  People: Zhaohua Ding, Hannah Kilpatrick
%  Grant support: NIH/NIAMS R01 AR073831

%% Get basic options from the input arguments

anat_image = double(anat_image);

anat_size = size(anat_image);                                                 %size of the anatomical image
anat_dims=fv_options.anat_dims;
anat_inplane_res = anat_dims(1)/anat_size(1);
anat_slicethick = anat_dims(2);
anat_slices = fv_options.anat_slices;

if fv_options.plot_mesh==1
    mesh_size = fv_options.mesh_size;
    mesh_dims = fv_options.mesh_dims;
    mesh_inplane_res = mesh_dims(1)/mesh_size(1);
    mesh_inplane_ratio = mesh_inplane_res/anat_inplane_res;
    mesh_slicethick = mesh_dims(2);
    mesh_color = fv_options.mesh_color;
    if isfield('fv_options', 'mesh_dist')
        mesh_dist = fv_options.mesh_dist;
    else
        mesh_dist = 0;
    end
end

if fv_options.plot_mask==1
    mask_size = fv_options.mask_size;
    mask_dims = fv_options.mask_dims;
    mask_inplane_res = mask_dims(1)/mask_size(1);
    mask_inplane_ratio = mask_inplane_res/anat_inplane_res;
    mask_slicethick = mask_dims(2);
    mask_color = fv_options.mask_color;
end

if fv_options.plot_fibers>0
    dti_size = fv_options.dti_size;
    dti_dims = fv_options.dti_dims;
    dti_inplane_res = dti_dims(1)/dti_size(1);
    dti_slicethick = dti_dims(2);
    dti_inplane_ratio = dti_inplane_res/anat_inplane_res;
    fiber_color = fv_options.fiber_color;
    if isfield(fv_options, 'fiber_skip')
        fiber_skip = fv_options.fiber_skip;
    else
        fiber_skip = 1;
    end
end

%% as needed, prepare for plotting

if fv_options.plot_mask==1
    
    mask_mesh = zeros(length(find(sum(sum(mask)))), 100, length(find(sum(sum(mask)))));
    mask_mesh1 = zeros(length(find(sum(sum(mask)))), 100, length(find(sum(sum(mask)))));
    mask_mesh2 = zeros(length(find(sum(sum(mask)))), 100, length(find(sum(sum(mask)))));
    n=1;                                                                        %counter to account for absence of mask below msucle of interest
    
    for s=1:length(mask(1,1,:))
        slice_mask1 = mask(:,:,s);                                              %get the mask for each slice
        slice_mask2 = bwmorph(slice_mask1, 'dilate');                           %get a second copy, dilate it once
        if length(find(slice_mask1))>1
            
            mask_boundaries1 = bwboundaries(slice_mask1);                       %gets the boundaries for the mask; in circular order
            rc_indx1 = ind2sub(size(slice_mask1), mask_boundaries1{1});         %convert from indices to row/column
            mask_mesh1(n, 1:99, 1)=imresize(rc_indx1(:,1), [99 1]);             %place x, y, z coordinates into a mesh for the mask; account for resizing factor
            mask_mesh1(n, 1:99, 2)=imresize(rc_indx1(:,2), [99 1]);
            mask_mesh1(n, 1:99, 3)=s;
            
            mask_boundaries2 = bwboundaries(slice_mask2);                       %gets the boundaries for the mask; in circular order
            rc_indx2 = ind2sub(size(slice_mask1), mask_boundaries2{1});         %convert from indices to row/column
            mask_mesh2(n, 1:99, 1)=imresize(rc_indx2(:,1), [99 1]);             %place x, y, z coordinates into a mesh for the mask; account for resizing factor
            mask_mesh2(n, 1:99, 2)=imresize(rc_indx2(:,2), [99 1]);
            mask_mesh2(n, 1:99, 3)=s;
            n=n+1;
        end
    end
    
    mask_mesh1(:, 100,:)=mask_mesh1(:,1,:);                                     %Close the circle - mask 1
    mask_mesh2(:, 100,:)=mask_mesh2(:,1,:);                                     %Close the circle - mask 2
    mask_mesh(:,:,1) = mean(cat(3, mask_mesh1(:,:,1), mask_mesh2(:,:,1)), 3);   %average the x, y, and z positions from the original and dilated masks (moves boundary to center of pixel)
    mask_mesh(:,:,2) = mean(cat(3, mask_mesh1(:,:,2), mask_mesh2(:,:,2)), 3);
    mask_mesh(:,:,3) = mean(cat(3, mask_mesh1(:,:,3), mask_mesh2(:,:,3)), 3);
end

if fv_options.plot_fibers==1                                                  %convert from pixels to desired units, as input by users
    fiber_all(:,:,:,1) = fiber_all(:,:,:,1)*dti_inplane_ratio;
    fiber_all(:,:,:,2) = fiber_all(:,:,:,2)*dti_inplane_ratio;
    fiber_all(:,:,:,3) = fiber_all(:,:,:,3)*dti_slicethick;
elseif fv_options.plot_fibers==2                                                  %convert from pixels to desired units, as input by users
    fiber_all(:,:,:,1,:) = fiber_all(:,:,:,1,:)*dti_inplane_ratio;
    fiber_all(:,:,:,2,:) = fiber_all(:,:,:,2,:)*dti_inplane_ratio;
    fiber_all(:,:,:,3,:) = fiber_all(:,:,:,3,:)*dti_slicethick;
end

%% plot

fiber_figure = figure('Name', 'Plot of Images and Selected Fiber-Tracking Options');

% plot the images
hold on;
for s=1:length(anat_slices)
    loop_image = anat_image(:,:,anat_slices(s));
    loop_image = 240*loop_image/max(max(loop_image));
    surf(ones(anat_size(1), anat_size(2))*anat_slicethick*anat_slices(s), loop_image)
    colormap('gray')
    shading('Interp')
end

if fv_options.plot_mesh==1
    if numel(mesh_color)==3
        mesh_plot=surf((roi_mesh(:, :, 2) + mesh_dist*roi_mesh(:, :, 5))*mesh_inplane_ratio, ...
            (roi_mesh(:, :, 1) + mesh_dist*roi_mesh(:, :, 4))*mesh_inplane_ratio, ...
            (roi_mesh(:, :, 3) + mesh_dist*roi_mesh(:, :, 6))*mesh_slicethick);
        set(mesh_plot, 'FaceColor', mesh_color, 'EdgeColor', mesh_color/2);
    elseif numel(size(mesh_color))==3
        mesh_plot=surf((roi_mesh(:, :, 2) + mesh_dist*roi_mesh(:, :, 5))*mesh_inplane_ratio, ...
            (roi_mesh(:, :, 1) + mesh_dist*roi_mesh(:, :, 4))*mesh_inplane_ratio, ...
            (roi_mesh(:, :, 3) + mesh_dist*roi_mesh(:, :, 6))*mesh_slicethick, ...
            mesh_color);
        set(mesh_plot, 'EdgeAlpha', 0);
    end
end

if fv_options.plot_mask==1
    mask_plot=surf(mask_mesh(:, :, 2)*mask_inplane_ratio, mask_mesh(:, :, 1)*mask_inplane_ratio, mask_mesh(:, :, 3)*mask_slicethick);
    set(mask_plot, 'FaceColor', mask_color, 'EdgeColor', mask_color, 'facealpha', .25, 'edgealpha', 0);
end

if fv_options.plot_fibers==1
    
    if numel(fiber_color)==3
        
        for row_cntr = 1:fiber_skip:(length(roi_mesh(:,1,1)))
            
            for col_cntr = 1:fiber_skip:(length(roi_mesh(1,:,1)))
                
                num_points = length(find(fiber_all(row_cntr,col_cntr,:,1)));
                if num_points>5
                    fiber_plot=plot3(squeeze(fiber_all(row_cntr,col_cntr, 1:num_points, 2)), ...
                        squeeze(fiber_all(row_cntr,col_cntr, 1:num_points, 1)), ...
                        squeeze(fiber_all(row_cntr,col_cntr, 1:num_points, 3)));
                    
                    loop_color = abs(fiber_color + randn(size(fiber_color))*0.1);       %add some random contrast
                    loop_color = loop_color/norm(loop_color);
                    set(fiber_plot, 'color', loop_color)
                end
                
            end
        end
        
    elseif numel(size(fiber_color))==3
        
        for row_cntr = 2:fiber_skip:(length(roi_mesh(:,1,1))-1)
            
            for col_cntr = 2:fiber_skip:(length(roi_mesh(1,:,1))-1)
                
                num_points = length(find(fiber_all(row_cntr,col_cntr,:,1)));
                fiber_plot=plot3(squeeze(fiber_all(row_cntr,col_cntr, 1:num_points, 2)), ...
                    squeeze(fiber_all(row_cntr,col_cntr, 1:num_points, 1)), ...
                    squeeze(fiber_all(row_cntr,col_cntr, 1:num_points, 3)));
                
                loop_color = squeeze(fiber_color(row_cntr, col_cntr,:));
                set(fiber_plot, 'color', loop_color)
                
            end
            
        end
        
    else
        
        error('Incorrect fiber tract color option.')
        
    end
    
elseif fv_options.plot_fibers==2
    
    for f=1:2
        
        fiber_all_f = squeeze(fiber_all(:,:,:,:,f));
        
        if numel(fiber_color)==6
            
            for row_cntr = 1:fiber_skip:(length(roi_mesh(:,1,1)))
                
                for col_cntr = 1:fiber_skip:(length(roi_mesh(1,:,1)))
                    
                    num_points = length(find(fiber_all_f(row_cntr,col_cntr,:,1)));
                    fiber_plot=plot3(squeeze(fiber_all_f(row_cntr,col_cntr, 1:num_points, 2)), ...
                        squeeze(fiber_all_f(row_cntr,col_cntr, 1:num_points, 1)), ...
                        squeeze(fiber_all_f(row_cntr,col_cntr, 1:num_points, 3)));
                    
                    loop_color = fiber_color(f,:);
                    loop_color = abs(loop_color+randn(size(loop_color))*0.1);          %add some contrast to the tracts
                    loop_color = loop_color/norm(loop_color);
                    set(fiber_plot, 'color', loop_color)
                    
                end
            end
            
        elseif numel(size(fiber_color))==4
            
            for row_cntr = 2:fiber_skip:(length(roi_mesh(:,1,1))-1)
                
                for col_cntr = 2:fiber_skip:(length(roi_mesh(1,:,1))-1)
                    
                    num_points = length(find(fiber_all_f(row_cntr,col_cntr,:,1)));
                    fiber_plot=plot3(squeeze(fiber_all_f(row_cntr,col_cntr, 1:num_points, 2)), ...
                        squeeze(fiber_all_f(row_cntr,col_cntr, 1:num_points, 1)), ...
                        squeeze(fiber_all_f(row_cntr,col_cntr, 1:num_points, 3)));
                    
                    loop_color = squeeze(fiber_color(row_cntr, col_cntr, :, f));
                    set(fiber_plot, 'color', loop_color)
                    
                end
                
            end
            
        else
            
            error('Incorrect fiber tract color option.')
            
        end
    end
    
end

set(gca, 'color', 'k')

%% end function

return
