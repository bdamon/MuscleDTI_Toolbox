function [mask, alt_mask] = define_muscle(anat_image, slices, alt_mask_size, plot_options)
%
%FUNCTION define_muscle
%  [mask, alt_mask] = define_muscle(anat_image, slices, alt_mask_size, plot_options)
%
%USAGE
%  define_muscle is used to define the boundary of a muscle and return its
%  binary image mask. This mask is needed by the functions define_roi and 
%  fiber_track and may be visualized using fiber_visualizer.
%
%  The user provides the anatomical images to be segmented and defines the slice 
%  numbers of interest.  After calling the function, a single figure window 
%  is opened. The initial slice of interest is displayed in the middle pane;
%  the preceding slice (if present) is displayed in the upper left-hand 
%  pane; and the next slice (if present) is displayed in the lower left-
%  hand pane. For the center pane, the zoom tool is enabled; the user clicks
%  and drags the left mouse button to zoom to the muscle of interest. To close
%  the zoom tool, the user selects Enter on their keyboard. All images are
%  then zoomed to this level.
%
%  Using the roipoly tool, the user uses a series of left mouse clicks to 
%  define the muscle as a region of interest (ROI) in the middle pane. 
%  After defining the vertices and closing the ROI, the user can adjust
%  vertex positions.  To complete ROI selection, the user right-clicks the
%  mouse; this brings up bring up a menu, from which the user selects Create
%  Mask to complete ROI selection.  
%
%  Then the program advances to the next slice. In this slice and all
%  subsequent slices, the level of zoom is automatically set as ±20 pixels
%  beyond the previous ROI’s row and column limits.  In the upper left-hand
%  pane, the preceding slice and its ROI are shown. Also shown are gold and
%  red lines depicting the center row and column, respectively, of this ROI.
%  In the lower left-hand pane, the next slice is shown. In the right-center
%  pane, the projections of the image stack and the ROI along the red line
%  are shown.  In the far-right pane, the projections of the image stack
%  and the ROI along the gold line are shown. Examining these windows can
%  help the user maintain consistency in ROI selection. ROI selection
%  continues in this manner until all slices of interest have been defined.
%
%  By default, the mask has the same dimensions as the input image. If the
%  DTI images and structural images have different dimensions from each 
%  other, an alternatively sized mask may also be calculated.  A MATLAB data
%  file named mask_file.mat, containing the mask and the alternatively sized
%  mask (if present), is automatically saved in the working directory. The
%  user is advised to rename this file promptly.
%
%  The mask may be viewed using fiber_visualizer, either as part of the
%  function call to define_muscle or directly from the command line.
%
%INPUT ARGUMENTS
%  anat_image: A row x column x slices stack of images
%
%  slices: A two element vector containing the first and last slices to be
%    analyzed
%
%  alt_mask_size: If specified, this is a two element vector containing the 
%    row x column x slices size of a second mask; the same center position of 
%    the image stack is assumed.
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
%  For help visualizing the data, see <a href="matlab: help fiber_visualizer">fiber_visualizer</a>.
%  For help defining the ROI, see <a href="matlab: help define_roi">define_roi</a>.
%  For help with the fiber tracking program, see <a href="matlab: help fiber_track">fiber_track</a>.
%  For help fitting the fiber tracts, see <a href="matlab: help fiber_smoother">fiber_smoother</a>.
%  For help quantifying fiber tracts, see <a href="matlab: help fiber_quantifier">fiber_quantifier</a>.
%  For help selecting fiber tracts following their quantification, see <a href="matlab: help fiber_goodness">fiber_goodness</a>.
%
% VERSION INFORMATION
%  v. 0.2
%
% ACKNOWLEDGMENTS
%  Grant support: NIH/NIAMS R01 AR050101, NIH/NIAMS R01 AR073831

%% display current, preceding, and next slices in three windows; open toolbar to adjust contrast of current slice

% initialize the mask
mask = zeros(size(anat_image));

% create the figure windows
screen_size = get(0,'ScreenSize');
figure(101)
set(gcf, 'position', [.05*screen_size(3) .1*screen_size(4) .8*screen_size(3) .6*screen_size(4)])

% create a loop counter to index values separately from s
n=1;

% ROI selection loop
for s=slices(1):slices(2)
    
    figure(101)
    
    %upper left has the preceding slice
    figure(101)
    subplot(2,5,1)
    slice_prev=max([(s-1) 1]);                                                  %make sure it doesn't crash if user starts on slice #1 of the stack
    image_prev = anat_image(:,:,slice_prev);
    image_prev = image_prev/max(max(image_prev));
    imagesc(image_prev)
    colormap gray
    title(['Previous Slice (#' num2str(slice_prev) ')'])
    axis image, axis off
    
    [counts, img_values] = imhist(image_prev);
    set(gca, 'clim', [0 img_values(find(cumsum(counts)>(0.995*sum(counts)), 1))])      %adjust windowing                                               %image contrast tool
    
    if n>1
        hold on
        plot(x_points, y_points, 'c-', 'linewidth', 2)
        plot([mean_col mean_col], [min_row max_row], 'r--', 'linewidth', 2)
        plot([min_col max_col], [mean_row mean_row], 'color', [1 .72 0], 'linewidth', 2)
        set(gca, 'xlim', [min_col max_col], 'ylim', [min_row max_row]);
       
    end
    
    %lower left has upcoming slice
    subplot(2,5,6)
    slice_next = min([(s+1) length(anat_image(1,1,:))]);
    image_next = anat_image(:,:,slice_next);
    image_next = image_next/max(max(image_next));
    imagesc(image_next)
    colormap gray
    title(['Upcoming Slice (#' num2str(slice_next) ')'])
    axis image, axis off
    
    [counts, img_values] = imhist(image_next);
    set(gca, 'clim', [0 img_values(find(cumsum(counts)>(0.995*sum(counts)), 1))])      %adjust windowing                                               %image contrast tool
    
    if n>1
        set(gca, 'xlim', [min_col max_col], 'ylim', [min_row max_row]);
    end
    
    %middle pane has current slice
    subplot(2,5, [2 3 7 8])
    slice_current=s;
    image_current = anat_image(:,:,slice_current);
    image_current = image_current/max(max(image_current));
    imagesc(image_current)
    colormap gray
    title(['Current Slice (#' num2str(slice_current) ')'])
    axis image, axis off
    
    [counts, img_values] = imhist(image_current);
    set(gca, 'clim', [0 img_values(find(cumsum(counts)>(0.995*sum(counts)), 1))])      %adjust windowing                                               %image contrast tool
    
    
    if n==1                                                             %first time, user select the level of zoom
        subplot(2,5, [2 3 7 8])
        text(5, 5, 'Zoom image and then select Enter to continue', 'color', 'y')
        zoom on
        pause
        
        xlim=get(gca, 'xlim');
        ylim=get(gca, 'ylim');
        min_col = min(xlim);
        max_col = max(xlim);
        mean_col = round(mean([min_col max_col]));
        
        min_row = min(ylim);
        max_row = max(ylim);
        mean_row = round(mean([min_row max_row]));
        
        subplot(2 ,5, 1)
        set(gca, 'xlim', [min_col max_col], 'ylim', [min_row max_row]);
        subplot(2, 5, 6)
        set(gca, 'xlim', [min_col max_col], 'ylim', [min_row max_row]);
        
    else                                                                %after that, set automatically based on previous slice mask
        
%         subplot(2,5, [2 3 7 8])
        
        prev_mask=mask(:,:,s-1);
        min_row = find(sum(prev_mask,2), 1 ) - 20;                           %zoom level is +/-20 pixels from previous
        min_row = max([1 min_row]);                                         %but you can't have an axis limit <1
        max_row = find(sum(prev_mask,2), 1, 'last') + 20;
        max_row = min([length(mask(:,1,1)) max_row]);                   %but you can't have an axis limit bigger than the image
        mean_row = round(mean([min_row max_row]));
        
        min_col = find(sum(prev_mask,1), 1 ) - 20;
        min_col = max([1 min_col]);                                         %but you can't have an axis limit <1
        max_col = find(sum(prev_mask,1), 1, 'last') + 20;
        max_col = min([length(mask(1,:,1)) max_col]);                 %but you can't have an axis limit bigger than the image
        mean_col = round(mean([min_col max_col]));
        set(gca, 'xlim', [min_col max_col], 'ylim', [min_row max_row]);
        
    end
    
    % get the ROI
    subplot(2, 5, [2 3 7 8])
    [mask(:,:,s), x_points, y_points] = roipoly;
    loop_mask_c = flipud(squeeze(mask(:, mean_col, :))');
    
    %view in sagital and coronal planes
    show_image_c = flipud(squeeze(anat_image(:, mean_col, :))');
    show_image_c = show_image_c/max(max(show_image_c));
    show_image_c(:,:,2) = show_image_c(:,:,1);
    show_image_c(:,:,3) = show_image_c(:,:,1);
    show_image_c = 0.75*show_image_c;
    show_image_c(:,:,1) = show_image_c(:,:,1) + loop_mask_c*.25;
    
    subplot(2, 5, [4 9])
    imagesc(show_image_c)
    title('Projection Along Red (- -) Line')
    axis off
    
    loop_mask_r = flipud(squeeze(mask(mean_row, :, :))');
    show_image_r = flipud(squeeze(anat_image(mean_row, :, :))');
    show_image_r = show_image_r/max(max(show_image_r));
    show_image_r(:,:,2) = show_image_r(:,:,1);
    show_image_r(:,:,3) = show_image_r(:,:,1);
    show_image_r = 0.75*show_image_r;
    show_image_r(:,:,1) = show_image_r(:,:,1) + loop_mask_r*.25;
    show_image_r(:,:,2) = show_image_r(:,:,2) + loop_mask_r*.18;
    
    subplot(2, 5, [5 10])
    imagesc(show_image_r)
    title('Projection Along Gold Line')
    axis off
    
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
