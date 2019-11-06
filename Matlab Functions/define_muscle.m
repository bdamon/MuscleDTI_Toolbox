function [mask, alt_mask] = define_muscle(anat_image, slices, alt_mask_size, plot_options)
%
%FUNCTION define_muscle
%  [mask, alt_mask] = define_muscle(anat_image, slices, alt_mask_size, plot_options)
%
%USAGE
%    The function define_muscle is used to define the boundary of a muscle
%  and return the corresponding binary image mask, for use in the MuscleDTI_Toolbox.  
%  Three images are displayed: the current slice (middle panel), the preceding 
%  slice with its ROI (if present), and the next slice. Using the roipoly 
%  tool, the user defines the ROI in the middle slice. For the main figure window,
%  an interactive tool is opened that allows the user to adjust the image's
%  window and level settings. After closing the ROI, the program advances
%  to the next slice until all slices have been defined.
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
%  plot_options: If specified, this calls the fiber_visualizer function to
%    plot the mask.
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
%  For help smoothing fiber tracts, see <a href="matlab: help fiber_fitter">fiber_fitter</a>.
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

mask = zeros(size(anat_image));

figure('position', [100 200 500 500], 'name', 'Previous Slice')
figure('position', [400 100 700 700], 'name', 'Current Slice')
figure('position', [700 200 500 500], 'name', 'Next Slice')
n=1;
for s=slices(1):slices(2)
    
    if s>1
        figure(1)
        imagesc(anat_image(:,:,s-1))
        colormap gray
        title(['Slice #' num2str(s-1)])
        axis off
        if n>1
            hold on
            plot(x_points, y_points, 'c')
        end
    end
    
    if s<=length(anat_image(1,1,:))
        figure(3)
        imagesc(anat_image(:,:,s+1))
        colormap gray
        title(['Slice #' num2str(s+1)])
        axis off
    end
    
    figure(2)
    imagesc(anat_image(:,:,s))
    colormap gray
    title(['Slice #' num2str(s)])
    axis off
    imcontrast
    [mask(:,:,s), x_points, y_points] = roipoly;
    n=n+1;
end

%% as specified by user, create an additional mask

form_alt_mask = ~isempty(alt_mask_size);

if form_alt_mask==1
    
    alt_mask = zeros(alt_mask_size(1), alt_mask_size(2), length(anat_image(1,1,:)));
    
    for s=1:length(anat_image(1,1,:))
        alt_mask(:,:,s) = imresize(squeeze(mask(:,:,s)), alt_mask_size);
    end
    
end

%% save masks

if form_alt_mask==1
    save mask_file mask alt_mask
else
    save mask_file mask
end

%% plot mask, if desired

if exist('plot_options', 'var')
    
    % be sure not to plot unneeded stuff
    plot_options.plot_fibers=0;
    plot_options.plot_mesh=0;
    
    fiber_visualizer(anat_image, plot_options, [], mask, []);
    
end

%% end the function

return

