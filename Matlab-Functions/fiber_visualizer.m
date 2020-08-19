function fiber_figure = fiber_visualizer(anat_image, plot_options, roi_mesh, mask, fiber_all)
%
%FUNCTION fiber_visualizer
%  fiber_figure = fiber_visualizer(anat_image, plot_options, roi_mesh, mask, fiber_all)
%
%USAGE
%    The function fiber_visualizer is used to visualize images and the muscle
%  mask, roi mesh, and/or fiber tracts formed using the MuscleDTI_Toolbox.
%  Input the images and a structure containing the plotting options. Depending
%  on the data you wish to view, the roi mesh, mask, and/or fiber tracts must
%  also be input and additional plotting options will be required.
%
%INPUT ARGUMENTS
%  anat_image: The stack of images to be plotted for anatomical reference
%
%  plot_options: A structure containing the following required fields:
%
%    anat_dims: A two element vector containing the FOV and the slice
%      thickness of the anatomical images.
%
%    anat_slices: A vector containing the slice numbers of the anatomical
%      images to be plotted.
%
%    plot_mesh: If set to 1, this field will allow plotting of the roi
%       mesh. Otherwise, set to 0.
%
%    plot_mask: If set to 1, this field will allow plotting of the mask.
%        Otherwise, set to 0.
%
%    plot_fibers:  If set to 1, this field will allow plotting of a single
%      set of fiber tracts. If set to 2, this will allow plotting of two
%      sets of fiber tracts. Otherwise, set to 0.
%
%  Depending on the plot options selected, the following other fields may
%  be required:
%
%  If plot_mesh equals 1, you must also specify:
%    mesh_size: This specifies the in-plane matrix size of the images used to
%      generate the mesh.
%
%    mesh_dims: This specifies the FOV and slice thickness of the images used
%      to create the mesh.
%
%    mesh_color: This is a three element vector containing the RGB levels to
%      be used when plotting the mesh
%
%    mesh_dist: If you shifted the mesh for fiber tracking and want to show
%      add this field and set to the value used during fiber tracking.
%
%  If plot_mask equals 1, you must also specify:
%    mask_size: This specifies the in-plane matrix size of the images used to
%      generate the mask.
%
%    mask_dims: This two-element vector specifies the FOV and slice thickness
%      of the images used to create the mask.
%
%    mask_color: If the mask is to be plotted, create color scale using
%      either of the following two options:
%       -If mesh_color is a 3 element vector of values ranging from 0-1, 
%        the vector is interpreted as RGB levels.
%       -If mesh_color is a matrix with size of (#mesh rows) x (#mesh columns) x 3, 
%        and if these values range from 0-1, the matrix will be interpreted  
%        as RGB levels specific to each tract. This could be used to
%        represent the distribution of architectural parameters across the
%        aponeurosis
%
%  If plot_fibers equals 1 or 2, you must also specify:
%    dti_size: A2-element vector that specifies the matrix size of the images
%      used for fiber tracking.
%
%    dti_dims: This two-element vector specifies the FOV and slice thickness
%      of the DT images. The FOV is assumed to be square.
%
%    fiber_color: Defines the color of the tracts. Several options are available:
%       -If plot_fibers==1 and fiber_color is a 3 element vector of values
%        ranging from 0-1, the vector is interpreted as RGB levels.
%       -If plot_fibers==2 and fiber_color is a 2x3 matrix of values
%        ranging from 0-1, each row of the matrix is interpreted as RGB
%        levels for the respective sets of tracts.
%       -If plot_fibers==1 and fiber_color is a matrix with size of
%        (#mesh rows) x (#mesh columns) x 3, and if these values range from
%        0-1, the third dimension of the matrix will be interpreted as RGB
%        levels specific to each tract
%       -If plot_fibers==2 and fiber_color is a matrix with size of
%        (#mesh rows) x (#mesh columns) x 3 x 2, and if these values range from
%        0-1, the third dimension of the matrix will be interpreted as RGB
%        levels specific to each tract, separately for sets 1 and 2
%
%    fiber_skip: Setting fiber_skip to integer values > 1 will skip over
%      fiber tracts when plotting. This may improve visualization and will
%      decrease time for rendering. If not specified, all fibers will be
%      plotted.
%
%    contrast_multiplier (optional): Used to adjust the brightness of the
%      images being displayed. Set to 1 for default brightness, <1 to darken
%      the images, or >1 to brighten them.
%
%  roi_mesh: The roi mesh, defined in define_roi. It is only needed if
%    plot_options.plot_mesh is set to 1.
%
%  mask: A binary mask around the muscle of interest. It could be the
%    output of define_muscle or it could have been defined in another
%    program. It is only needed if plot_options.plot_mask is set to 1.
%
%  fiber_all: The output of fiber_track (original fiber tracts) or
%    fiber_smoother (smoothed fiber tracts). It is only needed if
%    plot_options.plot_fibers is set to 1. If plot_fibers equals 1, the size
%    should be (#mesh rows) x (#mesh columns) x (#fiber tract points) x 3.
%    If plot_fibers equals 1, the size should be (#mesh rows) x (#mesh columns)
%    x (#fiber tract points) x 3 x 2.
%
%OUTPUT ARGUMENTS
%  fiber_figure: A Matlab figure structure
%
%OTHER FUNCTIONS IN THE MUSCLE DTI FIBER-TRACKING TOOLBOX
%  For help defining the mask, see <a href="matlab: help define_muscle">define_muscle</a>.
%  For help defining the ROI, see <a href="matlab: help define_roi">define_roi</a>.
%  For help with the fiber tracking program, see <a href="matlab: help fiber_track">fiber_track</a>.
%  For help smoothing fiber tracts, see <a href="matlab: help fiber_smoother">fiber_smoother</a>.
%  For help quantifying fiber tracts, see <a href="matlab: help fiber_quantifier">fiber_quantifier</a>.
%  For help selecting fiber tracts following their quantification, see <a href="matlab: help fiber_selector">fiber_selector</a>.
%
%VERSION INFORMATION
%  In beta testing mode
%
%ACKNOWLEDGEMENTS
%  People: Zhaohua Ding, Hannah Kilpatrick
%  Grant support: NIH/NIAMS R01 AR073831

%% Get basic options from the input arguments

anat_image = double(anat_image);

anat_size = size(anat_image);                                                 %size of the anatomical image
anat_dims=plot_options.anat_dims;
anat_inplane_res = anat_dims(1)/anat_size(1);
anat_slicethick = anat_dims(2);
anat_slices = plot_options.anat_slices;

if plot_options.plot_mesh==1
    mesh_size = plot_options.mesh_size;
    mesh_dims = plot_options.mesh_dims;
    mesh_inplane_res = mesh_dims(1)/mesh_size(1);
    mesh_inplane_ratio = mesh_inplane_res/anat_inplane_res;
    mesh_slicethick = mesh_dims(2);
    mesh_color = plot_options.mesh_color;
    if isfield('plot_options', 'mesh_dist')
        mesh_dist = plot_options.mesh_dist;
    else
        mesh_dist = 0;
    end
end

if plot_options.plot_mask==1
    mask_size = plot_options.mask_size;
    mask_dims = plot_options.mask_dims;
    mask_inplane_res = mask_dims(1)/mask_size(1);
    mask_inplane_ratio = mask_inplane_res/anat_inplane_res;
    mask_slicethick = mask_dims(2);
    mask_color = plot_options.mask_color;
end

if plot_options.plot_fibers>0
    dti_size = plot_options.dti_size;
    dti_dims = plot_options.dti_dims;
    dti_inplane_res = dti_dims(1)/dti_size(1);
    dti_slicethick = dti_dims(2);
    dti_inplane_ratio = dti_inplane_res/anat_inplane_res;
    fiber_color = plot_options.fiber_color;
    if isfield(plot_options, 'fiber_skip')
        fiber_skip = plot_options.fiber_skip;
    else
        fiber_skip = 1;
    end
end

if isfield(plot_options, 'contrast_multiplier')
    contrast_multiplier=plot_options.contrast_multiplier;
else
    contrast_multiplier=1;
end

%% as needed, prepare for plotting

if plot_options.plot_mask==1
    
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

if plot_options.plot_fibers==1                                                  %convert from pixels to desired units, as input by users
    fiber_all(:,:,:,1) = fiber_all(:,:,:,1)*dti_inplane_ratio;
    fiber_all(:,:,:,2) = fiber_all(:,:,:,2)*dti_inplane_ratio;
    fiber_all(:,:,:,3) = fiber_all(:,:,:,3)*dti_slicethick;
elseif plot_options.plot_fibers==2                                                  %convert from pixels to desired units, as input by users
    fiber_all(:,:,:,1,:) = fiber_all(:,:,:,1,:)*dti_inplane_ratio;
    fiber_all(:,:,:,2,:) = fiber_all(:,:,:,2,:)*dti_inplane_ratio;
    fiber_all(:,:,:,3,:) = fiber_all(:,:,:,3,:)*dti_slicethick;
end

%% plot

fiber_figure = figure('Name', 'Plot of Images and Selected Fiber-Tracking Options');

%scale the images to max=1
anat_image_norm = contrast_multiplier*255*anat_image/max(max(max(anat_image)));

% plot the images
hold on;
for s=1:length(anat_slices)
    surf(ones(anat_size(1), anat_size(2))*anat_slicethick*anat_slices(s), anat_image_norm(:,:,anat_slices(s)))
    colormap('gray')
    shading('Interp')
end

if plot_options.plot_mesh==1
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

if plot_options.plot_mask==1
    mask_plot=surf(mask_mesh(:, :, 2)*mask_inplane_ratio, mask_mesh(:, :, 1)*mask_inplane_ratio, mask_mesh(:, :, 3)*mask_slicethick);
    set(mask_plot, 'FaceColor', mask_color, 'EdgeColor', mask_color, 'facealpha', .25, 'edgealpha', 0);
end

if plot_options.plot_fibers==1
    
    if numel(fiber_color)==3
        
        for row_cntr = 1:fiber_skip:(length(roi_mesh(:,1,1)))
            
            for col_cntr = 1:fiber_skip:(length(roi_mesh(1,:,1)))
                
                num_points = length(find(fiber_all(row_cntr,col_cntr,:,1)));
                if num_points>5
                    fiber_plot=plot3(squeeze(fiber_all(row_cntr,col_cntr, 1:num_points, 2)), ...
                        squeeze(fiber_all(row_cntr,col_cntr, 1:num_points, 1)), ...
                        squeeze(fiber_all(row_cntr,col_cntr, 1:num_points, 3)));
                    
                    loop_color = abs(fiber_color + randn(size(fiber_color))*0.15);       %add some random contrast
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
    
elseif plot_options.plot_fibers==2
    
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