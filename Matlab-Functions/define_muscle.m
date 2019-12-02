function [mask, alt_mask] = define_muscle(anat_image, slices, alt_mask_size, plot_options)
%
%FUNCTION define_muscle
%  [mask, alt_mask] = define_muscle(anat_image, slices, alt_mask_size, plot_options)
%
%USAGE
%    The function define_muscle is used to define the boundary of a muscle
%  and return the corresponding binary image mask, for use in the MuscleDTI_Toolbox.  
%  A mask defining the muscle boundaries is a required input to define_roi and
%  fiber_track, and it may be used in fiber_visualizer. The mask has the same 
%  dimensions as the input image. An alternatively sized mask may also 
%  be calculated. This is needed if the DTI images and structural images have 
%  different dimensions.
%    After calling the function, three windows are opened: the current slice  
%  (middle window), the preceding slice, and the next slice. Initially, the   
%  zoom feature is enabled so that the user can view the muscle of interest more 
%  closely. Using the roipoly tool, the user defines the ROI in the middle window.   
%  For the main figure window, an interactive tool is opened that allows the user   
%  to adjust the image's window and level settings. After closing the ROI, the  
%  program advances to the next slice. In this and future slices, the level of
%  zoom is automatically set.  This procedure continues until all slices have
%  been defined.
%    A file named mask_file, containing the mask and (if present) the
%  alternatively sized mask, is automatically saved in the working
%  directory.
%
%INPUT ARGUMENTS
%  anat_image: A row x column x slices stack of images
%
%  slices: A two element vector containing the first and last slices to be
%    analyzed
%
%  alt_mask_size: If specified, this is a two element vector containing the row  
%    x column size of a second mask; the same number of slices is assumed.
%
%  plot_options: If included, this calls the fiber_visualizer function to plot 
%    the mask.
%
%OUTPUT ARGUMENTS
%  mask: the binary image mask, with size matching that of the original
%    image
%
%  alt_mask: a second binary image mask, with size matching that of the
%    vector alt_mask_size
%
%OTHER FUNCTIONS IN THE MUSCLE DTI FIBER-TRACKING TOOLBOX
%  For help defining the ROI, see <a href="matlab: help define_roi">define_roi</a>.
%  For help with the fiber tracking program, see <a href="matlab: help fiber_track">fiber_track</a>.
%  For help fitting the fiber tracts, see <a href="matlab: help fiber_fitter">fiber_fitter</a>.
%  For help quantifying fiber tracts, see <a href="matlab: help fiber_quantifier">fiber_quantifier</a>.
%  For help selecting fiber tracts following their quantification, see <a href="matlab: help fiber_selector">fiber_selector</a>.
%  For help visualizing the data, see <a href="matlab: help fiber_visualizer">fiber_visualizer</a>.
%
% VERSION INFORMATION
%  v. 0.1
%
% ACKNOWLEDGMENTS
%  Grant support: NIH/NIAMS R01 AR050101, NIH/NIAMS R01 AR073831

%% display current, preceding, and next slices in three windows; open toolbar to adjust contrast of current slice

% initialize the mask
mask = zeros(size(anat_image));

% create the figure windows
screen_size = get(0,'ScreenSize');
screen_aspect_ratio = screen_size(3) / screen_size(4);
figure('units', 'normalized', 'position', [.05 (1-.3*screen_aspect_ratio) .25 .25*screen_aspect_ratio], 'name', 'Previous Slice')
figure('units', 'normalized', 'position', [.3 .2 .4 .4*screen_aspect_ratio], 'name', 'Current Slice')
figure('units', 'normalized', 'position', [.7 (1-.3*screen_aspect_ratio) .25 .25*screen_aspect_ratio], 'name', 'Next Slice')

% create a loop counter to index values separately from s
n=1;

% ROI selection loop
for s=slices(1):slices(2)
    
    %figure has preceding slice
    figure(1)
    slice1=max([(s-1) 1]);
    imagesc(anat_image(:,:, slice1))
    colormap gray
    title(['Slice #' num2str(slice1)])
    axis off
    if n>1
        hold on
        plot(x_points, y_points, 'c')
        set(gca, 'xlim', [min_col max_col], 'ylim', [min_row max_row]);
    end
    
    %figure 3 has upcoming slice
    figure(3)
    slice3=min([(s+1) length(anat_image(1,1,:))]);
    imagesc(anat_image(:,:,slice3))
    colormap gray
    title(['Slice #' num2str(slice3)])
    axis off
    if n>1
        set(gca, 'xlim', [min_col max_col], 'ylim', [min_row max_row]);
    end
    
    %figure 2 has current slice
    figure(2)
    slice2=s;
    imagesc(anat_image(:,:,slice2))
    colormap gray
    title(['Slice #' num2str(slice2)])
    axis off
    imcontrast                                                          %image contrast tool
    
    if n==1                                                             %first time, user select the level of zoom
        text(5, 5, 'Zoom image and then select Enter to continue', 'color', 'y')
        zoom on
        pause
        
        xlim=get(gca, 'xlim');
        ylim=get(gca, 'ylim');
        min_col = min(xlim);
        max_col = max(xlim);
        min_row = min(ylim);
        max_row = max(ylim);
        
        figure(1)
        set(gca, 'xlim', [min_col max_col], 'ylim', [min_row max_row]);
        figure(3)
        set(gca, 'xlim', [min_col max_col], 'ylim', [min_row max_row]);
        
    else                                                                %after that, set automatically based on previous slice mask
        
        prev_mask=mask(:,:,s-1);
        min_row = find(sum(prev_mask,2), 1 ) - 5;
        max_row = find(sum(prev_mask,2), 1, 'last') + 5;
        min_col = find(sum(prev_mask,1), 1 ) - 5;
        max_col = find(sum(prev_mask,1), 1, 'last') + 5;
        set(gca, 'xlim', [min_col max_col], 'ylim', [min_row max_row]);
        
    end
    
    % get the ROI
    figure(2)
    [mask(:,:,s), x_points, y_points] = roipoly;
    
    %advance the loop counter
    n=n+1;
end

%% as specified by user, create an additional mask

form_alt_mask = ~isempty(alt_mask_size);

if form_alt_mask==1
    
    alt_mask = zeros(alt_mask_size(1), alt_mask_size(2), length(anat_image(1,1,:)));
    
    for s=1:length(anat_image(1,1,:))
        alt_mask(:,:,s) = imresize(squeeze(mask(:,:,s)), alt_mask_size);
    end
    
else                            % form a throwaway variable so that the program doesn't crash
    
    alt_mask=1;
    
end

%% save masks

if form_alt_mask==1
    save mask_file mask alt_mask
else
    save mask_file mask
end

%% plot mask, if desired

plot_mask = isfield(plot_options, 'plot_mask');

if plot_mask==1

    % be sure not to plot unneeded stuff
    plot_options.plot_fibers=0;
    plot_options.plot_mesh=0;
    plot_options.plot_mask=1;
    
    fiber_visualizer(anat_image, plot_options, [], mask, []);
    
end

%% end the function

return
