function roi_mesh=define_roi(anat_image, mask, defroi_options, plot_options)
%
%FUNCTION define_roi
%  roi_mesh=define_roi(anat_image, mask, defroi_options, plot_options)
%
%USAGE
%    The function define_roi is used to digitize the aponeurosis of muscle fiber
%  insertion in the MuscleDTI_Toolbox.  The digitized points are used to
%  reconstruct a mesh; the mesh is used as the seed surface for fiber tracking.
%  There are two options for defining the aponeurosis:
%    -Manual: The user is prompted initially to select two points, which define
%     the level of zoom to be used throughout the entire process. Then the user
%     advances through the slices to select the aponeurosis. The mask is provided
%     as a guide; the aponeurosis must fall within its boundaries. The selected
%     points can form a line or close to form a polygon. At each slice, the
%     user is given the option of repeating the procedure in case of error.
%     For each figure window, an interactive tool is opened that allows the user
%     to adjust the image's window and level settings.
%    -Automatic: The aponeurosis is automatically segmented from within the
%     region of the image represented by the muscle mask. TWo segmentation
%     methods (edge detection and k-means clustering) are used, and the
%     segmented region is defined as the consensus of the two results. The
%     region is displayed and the user is allowed to correct misassignments.
%     The boundaries of the segmented region are smoothed using a Savitsky-
%     Golay filter and used to form the mesh.
%  The mesh is initially formed with resolution n_row x n_col.  To smooth
%  the mesh, it is then downsampled by a size factor of four. Finally, the smoothed
%  mesh is used to create a high resolution mesh at the desired size. A file
%  called roi_mesh_file.mat is automatically saved in the working directory.
%    If the input argument plot_options is included, the mesh and mask are
%  plotted using the function fiber_visualizer.
%
%INPUT ARGUMENTS
%  anat_image: The imaging data.  If input as a structure, then the imaging
%    data are assumed to exist in a field called anat_image.Data.  If
%    specified as a matrix, the data are used directly.
%
%  mask: The mask, as defined by the function define_mask
%
%  defroi_options: A structure containing the following fields:
%    slices: A two-element vector containing the first and last slices that the
%      user wishes to digitize.
%
%    dti_size: The size of the DTI image dataset (rows x columns x slices),
%      input as a three element vector.
%
%    mesh_size: A two-element vector containing the numbers of rows (n_row) and
%     columns (n_col) desired in the output mesh.
%
%    method: a string variable set either to 'manual' or 'auto'. The manual
%      and automatic options are described above.
%
%  plot_options: If specified, this calls the fiber_visualizer function to
%    plot the mask and roi mesh.
%
%OUTPUT ARGUMENTS
%  roi_mesh: a 3D matrix containing the reconstructed mesh with size rows x
%    columns x 6. In the 3rd dimension, levels 1-3 hold the row-column-slice
%    data and levels 4-6 hold the normal vector to the roi surface at the
%    point {row column slice}.
%
%OTHER FUNCTIONS IN THE MUSCLE DTI FIBER-TRACKING TOOLBOX
%  For help defining the mask, see <a href="matlab: help define_muscle">define_muscle</a>.
%  For help with the fiber tracking program, see <a href="matlab: help fiber_track">fiber_track</a>.
%  For help smoothing fiber tracts, see <a href="matlab: help fiber_smoother">fiber_smoother</a>.
%  For help quantifying fiber tracts, see <a href="matlab: help fiber_quantifier">fiber_quantifier</a>.
%  For help selecting fiber tracts following their quantification, see <a href="matlab: help fiber_selector">fiber_selector</a>.
%  For help visualizing the data, see <a href="matlab: help fiber_visualizer">fiber_visualizer</a>.
%
%VERSION INFORMATION
%  v 0.1
%
%ACKNOWLEDGMENTS
%  People: Zhaohua Ding
%  Grant support: NIH/NIAMS R01 AR050101, NIH/NIAMS R01 AR073831

%% preliminary stuff
dti_size = defroi_options.dti_size;

frst_slice = defroi_options.slices(1);
last_slice = defroi_options.slices(2);

n_row = defroi_options.mesh_size(1);
n_col = defroi_options.mesh_size(2);

select_method = defroi_options.method(1:2);

if isstruct(anat_image)
    working_image=anat_image.Data;
else
    working_image=anat_image;
end
sz_work_imag=size(working_image);


%% digitization loops

switch select_method
    
    case{'ma'}                                          %manual definition
        
        k=frst_slice;
        while k<=last_slice
            
            %get the roi data
            curr_slice=k;
            if k>1
                figure(1)
                
                prev_slice_rgb(:,:,1) = squeeze(working_image(:,:,curr_slice-1,:));
                prev_slice_rgb(:,:,2) = squeeze(working_image(:,:,curr_slice-1,:));
                prev_slice_rgb(:,:,3) = squeeze(working_image(:,:,curr_slice-1,:));
                prev_slice_rgb = prev_slice_rgb/max(max(max(prev_slice_rgb)));
                prev_slice_rgb = prev_slice_rgb*0.75;
                
                imagesc(prev_slice_rgb)
                title(['Slice #' num2str(curr_slice-1)])
                axis off
                axis image
            end
            
            if k<=length(working_image(1,1,:))
                figure(3)
                
                next_slice_rgb(:,:,1) = squeeze(working_image(:,:,curr_slice+1,:));
                next_slice_rgb(:,:,2) = squeeze(working_image(:,:,curr_slice+1,:));
                next_slice_rgb(:,:,3) = squeeze(working_image(:,:,curr_slice+1,:));
                next_slice_rgb = next_slice_rgb/max(max(max(next_slice_rgb)));
                next_slice_rgb = next_slice_rgb*0.75;
                
                imagesc(next_slice_rgb)
                title(['Slice #' num2str(curr_slice+1)])
                axis off
                axis image
            end
            
            figure(2)
            
            curr_slice_rgb(:,:,1) = squeeze(working_image(:,:,curr_slice,:));
            curr_slice_rgb(:,:,2) = squeeze(working_image(:,:,curr_slice,:));
            curr_slice_rgb(:,:,3) = squeeze(working_image(:,:,curr_slice,:));
            curr_slice_rgb = curr_slice_rgb/max(max(max(curr_slice_rgb)));
            curr_slice_rgb = curr_slice_rgb*0.75;
            curr_slice_rgb(:,:,1) = curr_slice_rgb(:,:,1) + mask(:,:,curr_slice)*0.25;
            
            imagesc(squeeze(curr_slice_rgb))
            title(['Slice #' num2str(curr_slice)])
            axis off
            axis image
            
            if k==frst_slice
                text_hndl1=text(10, 20, 'Click two points to zoom image');              %zoom using two ginputs - give instruction to user
                set(text_hndl1, 'fontsize', 10, 'color', 'w');
                [x_tmp y_tmp]=ginput(2);
                set(text_hndl1, 'visible', 'off');
            end
            
            axis([min(x_tmp) max(x_tmp) min(y_tmp) max(y_tmp)]);
            text_hndl2=text(min(x_tmp)+2, min(y_tmp)+2, ...                       %define the ROI.  Must be done in an A to P manner so that points will be ordered appropriately in subsequent steps
                strvcat('Define the ROI using left mouse clicks beginning anteriorly', '  and continuing posteriorly.  Use rightmouse button to terminate the ROI'));
            set(text_hndl2, 'fontsize', 10, 'color', 'w');
            roi_curx = [];
            roi_cury = [];
            roi_curz = [];
            
            while 1
                [col_fix, row_fix, b] = ginput(1);
                
                if  b==3
                    break;
                end
                
                if  b==1
                    hold on
                    plot(col_fix, row_fix, 'w.');
                    roi_curx = [roi_curx row_fix];
                    roi_cury = [roi_cury col_fix];
                    roi_curz = [roi_curz k];
                end
                
            end
            
            set(text_hndl2, 'visible', 'off');
            
            text_hndl3=text(min(x_tmp)+2, min(y_tmp)+2, ...                             %wait for a click before continuing
                'Click left mouse button to continue or right mouse button to repeat this ROI');
            set(text_hndl3, 'fontsize', 10, 'color', 'y');
            [~, ~, b]=ginput(1);
            if b==1
                roi_surfx(curr_slice, :) = interp1(roi_curx, 1:(length(roi_curx)-1)/50:length(roi_curx), 'pchip');
                roi_surfy(curr_slice, :) = interp1(roi_cury, 1:(length(roi_cury)-1)/50:length(roi_cury), 'pchip');
                roi_surfz(curr_slice, :) = interp1(roi_curz, 1:(length(roi_curz)-1)/50:length(roi_curz), 'pchip');
                k=k+1;
            elseif b==3
                k=k;
            end
            
            set(text_hndl3, 'visible', 'off');
            
            clear text_hndl1 text_hndl2 text_hndl3 junk1 junk 2                         %get ready for the next time through the loop
            x_tmp=get(gca, 'XLim');
            y_tmp=get(gca, 'YLim');
            
        end
        
    case{'au'}                                  %automated definition
        
        apo_segmented_initial=zeros(size(working_image));
        apo_segmented_final=zeros(size(working_image));
        
        
        for k=frst_slice:last_slice
            
            curr_slice=k;
            
            %set level of zoom based on current mask size
            loop_mask = mask(:,:,curr_slice);
            mask_rows=[min(find(sum(loop_mask,2)))-6 max(find(sum(loop_mask,2)))+6];
            mask_cols=[min(find(sum(loop_mask,1)))-6 max(find(sum(loop_mask,1)))+6];
            
            loop_img = double(working_image(:,:,k));
            if k>frst_slice
                prev_segmented=apo_segmented;
            end
            
            loop_roi=mask(:,:,k);                                        %draw roi arund muscle of itnerest
            roi_edge=bwboundaries(loop_roi);
            roi_edge_xy=roi_edge{1};
            roi_erode=bwmorph(loop_roi, 'erode', 2);                  %erode the roi
            roi_dilate=bwmorph(loop_roi, 'dilate', 2);                  %dilate the roi
            
            %find edges within the mask
            edge_segmented=edge(loop_img.*roi_dilate);                    %find edges of and within dilated ROI
            edge_segmented=edge_segmented.*roi_erode;                           %clean out the edge of mask, though
            
            %also do a k-means clustering:
            [kmeans_segmented, kmeans_centers] = imsegkmeans(int16(loop_img.*roi_erode), 3);
            kmeans_segmented = double(kmeans_segmented).*double(roi_erode);
            kmeans_segmented(loop_img>(median(kmeans_centers)))=0;
            kmeans_segmented(kmeans_segmented>0)=1;
            
            % if there is a previous slice, add that in there too:
            if k>frst_slice
                apo_segmented = bwmorph(kmeans_segmented, 'dilate') + edge_segmented + prev_segmented;
            else
                apo_segmented = bwmorph(kmeans_segmented, 'dilate') + edge_segmented;%
            end
            apo_segmented(apo_segmented<2)=0;
            apo_segmented(apo_segmented>0)=1;
            
            apo_segmented=bwmorph(apo_segmented, 'clean');              %remove isolated pixels
            apo_segmented=bwmorph(apo_segmented, 'close');
            apo_segmented_initial(:,:,curr_slice)=apo_segmented;
            
            %visualize/correct
            loop_img_rgb=cat(3, loop_img, loop_img, loop_img);
            loop_img_rgb=loop_img_rgb/max(max(max(loop_img_rgb)));
            
            %visualize points in Panel 1
            figure(100)
            if k==frst_slice
                set(gcf, 'position', [0 50 1535 430]);
            else
                set(gcf, 'position', figure_position);
            end
            
            subplot(1,3,1), imagesc(loop_img_rgb), axis image
            [pixel_rows, pixel_cols] = ind2sub(size(apo_segmented), find(apo_segmented));
            hold on
            plot(pixel_cols, pixel_rows, 'm.', 'markersize', 5)
            title('Current Voxels')
            set(gca, 'xlim', mask_cols, 'ylim', mask_rows)
            
            %Visualize region, interact with panel 2
            loop_img_region=loop_img_rgb*0.75;
            loop_img_region(:,:,1)=loop_img_region(:,:,1)+.25*apo_segmented;
            subplot(1,3,2), imagesc(loop_img_region), axis image
            hold on
            plot(roi_edge_xy(:,2), roi_edge_xy(:,1), 'c');
            set(gcf, 'name', 'Follow the Prompts in Command Window');
            title('Interactive Window')
            xlabel(['Slice ' num2str(curr_slice) ' of Slices ' num2str(frst_slice) ' to ' num2str(last_slice)])
            set(gca, 'xlim', mask_cols, 'ylim', mask_rows)
            
            subplot(1,3,3)
            imagesc(loop_img_rgb)
            hold on
            axis image
            set(gca, 'xlim', mask_cols, 'ylim', mask_rows)
            if sum(sum(apo_segmented))>0
                edge_pixels = bwboundaries(apo_segmented);
                edge_xy=edge_pixels{1};
                clear edge_xy_smooth
                edge_xy_smooth(:,2)=smooth(edge_xy(:,2), 'sgolay');
                edge_xy_smooth(:,1)=smooth(edge_xy(:,1), 'sgolay');
                edge_xy_smooth(end,:)=edge_xy_smooth(1,:);
                plot(edge_xy_smooth(:,2), edge_xy_smooth(:,1), 'm.', 'markersize', 5);
            end
            
            correct_yn='c';
            initialize_correction=1;
            while correct_yn(1)=='c' || correct_yn(1)=='C'
                
                subplot(1,3,1), set(gca, 'xlim', mask_cols, 'ylim', mask_rows);
                subplot(1,3,2), set(gca, 'xlim', mask_cols, 'ylim', mask_rows);
                
                clc
                if initialize_correction==0
                    correct_yn=input('Select a to accept this result or c to correct pixels: ', 's');
                    if isempty(correct_yn)
                        correct_yn='c';
                    end
                else
                    initialize_correction=0;
                end
                
                set(gcf, 'name', 'Follow Prompt in Command Window');
                if correct_yn(1)=='a' || correct_yn(1)=='A'
                    break
                end
                
                disp('Use the left mouse button to select voxels for removal; use the right mouse button to select voxels to add;')
                disp('press ENTER when finished.')
                [col_fix,row_fix,button] = ginput;
                row_fix=round(row_fix);
                col_fix=round(col_fix);
                for b=1:length(button)
                    if button(b)==1
                        apo_segmented(row_fix(b), col_fix(b))=0;
                    elseif button(b)==3
                        apo_segmented(row_fix(b), col_fix(b))=1;
                    end
                end
                
                %visualize result
                figure(100)
                
                subplot(1,3,1)
                imagesc(loop_img_rgb), axis image
                [pixel_rows, pixel_cols] = ind2sub(size(apo_segmented), find(apo_segmented));
                set(gca, 'xlim', mask_cols, 'ylim', mask_rows)
                hold on
                plot(roi_edge_xy(:,2), roi_edge_xy(:,1), 'c');
                plot(pixel_cols, pixel_rows, 'm.', 'markersize', 5)
                
                loop_img_region=loop_img_rgb*0.75;
                loop_img_region(:,:,1)=loop_img_region(:,:,1)+.25*apo_segmented;
                subplot(1,3,2), imagesc(loop_img_region), axis image
                hold on
                plot(roi_edge_xy(:,2), roi_edge_xy(:,1), 'c');
                title('Interactive Window')
                xlabel(['Slice ' num2str(curr_slice) ' of Slices ' num2str(frst_slice) ' to ' num2str(last_slice)])
                set(gca, 'xlim', mask_cols, 'ylim', mask_rows)
                
                edge_pixels = bwboundaries(apo_segmented);
                edge_xy=edge_pixels{1};
                clear edge_xy_smooth
                edge_xy_smooth(:,2)=smooth(edge_xy(:,2), 'sgolay');
                edge_xy_smooth(:,1)=smooth(edge_xy(:,1), 'sgolay');
                edge_xy_smooth(end,:)=edge_xy_smooth(1,:);
                
                subplot(1,3,3)
                imagesc(loop_img_rgb)
                hold on
                plot(edge_xy_smooth(:,2), edge_xy_smooth(:,1), 'm.', 'markersize', 5);
                axis image
                set(gca, 'xlim', mask_cols, 'ylim', mask_rows);
                
            end
            
            clf
            
            %add to roi_surf matrices
            roi_curx = edge_xy_smooth(:,1);
            roi_cury = edge_xy_smooth(:,2);
            roi_curz = ones(size(roi_curx))*curr_slice;
            
            roi_surfx(curr_slice, :) = interp1(roi_curx, 1:(length(roi_curx)-1)/100:length(roi_curx), 'pchip');
            roi_surfy(curr_slice, :) = interp1(roi_cury, 1:(length(roi_cury)-1)/100:length(roi_cury), 'pchip');
            roi_surfz(curr_slice, :) = interp1(roi_curz, 1:(length(roi_curz)-1)/100:length(roi_curz), 'pchip');
            
            %add to matrix of final segmentations
            apo_segmented_final(:,:,curr_slice)=apo_segmented;
            
            figure_position = get(gcf, 'position');
            
        end
        
end


%% prompt user to correct any slices

if select_method=='au'
    clc
    correct_slices=input('Do you want to correct any slices (y/n)? ', 's');
    while correct_slices(1)=='y' || correct_slices(1)=='Y'
        
        curr_slice=input('What slice do you want to correct? ');
        if curr_slice>last_slice || curr_slice<frst_slice
            curr_slice=input(['Value entered is out of range.  Please enter a number from ' num2str(frst_slice) ' to ' num2str(last_slice) ':']);
        end
        
        %set level of zoom based on current mask size
        loop_mask = mask(:,:,curr_slice);
        mask_rows=[min(find(sum(loop_mask,2)))-6 max(find(sum(loop_mask,2)))+6];
        mask_cols=[min(find(sum(loop_mask,1)))-6 max(find(sum(loop_mask,1)))+6];
        
        %set loop-specific variables:
        loop_roi=mask(:,:,k);                                        %draw roi arund muscle of itnerest
        roi_edge=bwboundaries(loop_roi);
        roi_edge_xy=roi_edge{1};
        loop_img = double(working_image(:,:,curr_slice));
        
        %get the last result for this slice
        apo_segmented=apo_segmented_final(:,:,curr_slice);
        
        %visualize/correct
        loop_img_rgb=cat(3, loop_img, loop_img, loop_img);
        loop_img_rgb=loop_img_rgb/max(max(max(loop_img_rgb)));
        
        %visualize points in Panel 1
        figure(100)
        set(gcf, 'position', figure_position);
        
        subplot(1,3,1), imagesc(loop_img_rgb), axis image
        [pixel_rows, pixel_cols] = ind2sub(size(apo_segmented), find(apo_segmented));
        hold on
        plot(pixel_cols, pixel_rows, 'm.', 'markersize', 5)
        title('Current Voxels')
        
        %Visualize region, interact with panel 2
        loop_img_region=loop_img_rgb*0.75;
        loop_img_region(:,:,1)=loop_img_region(:,:,1)+.25*apo_segmented;
        subplot(1,3,2), imagesc(loop_img_region), axis image
        hold on
        plot(roi_edge_xy(:,2), roi_edge_xy(:,1), 'c');
        set(gcf, 'name', 'Follow the Prompts in Command Window');
        title('Interactive Window')
        xlabel(['Slice ' num2str(curr_slice) ' of Slices ' num2str(frst_slice) ' to ' num2str(last_slice)])
        zoom on
        
        %look at initial boundary result in panel 3
        edge_pixels = bwboundaries(apo_segmented);
        edge_xy=edge_pixels{1};
        clear edge_xy_smooth
        edge_xy_smooth(:,2)=smooth(edge_xy(:,2), 'sgolay');
        edge_xy_smooth(:,1)=smooth(edge_xy(:,1), 'sgolay');
        edge_xy_smooth(end,:)=edge_xy_smooth(1,:);
        
        subplot(1,3,3)
        imagesc(loop_img_rgb)
        hold on
        plot(edge_xy_smooth(:,2), edge_xy_smooth(:,1), 'm.', 'markersize', 5);
        axis image
        set(gca, 'xlim', mask_cols, 'ylim', mask_rows);
        
        correct_yn='c';
        initialize_correction=1;
        while correct_yn(1)=='c' || correct_yn(1)=='C'
            
            subplot(1,3,1), set(gca, 'xlim', mask_cols, 'ylim', mask_rows);
            subplot(1,3,2), set(gca, 'xlim', mask_cols, 'ylim', mask_rows);
            
            clc
            if initialize_correction==0
                correct_yn=input('Select a to accept this result or c to correct pixels: ', 's');
                if isempty(correct_yn)
                    correct_yn='c';
                end
            else
                initialize_correction=0;
            end
            
            set(gcf, 'name', 'Follow Prompt in Command Window');
            if correct_yn(1)=='a' || correct_yn(1)=='A'
                break
            end
            
            disp('Use the left mouse button to select voxels for removal; use the right mouse button to select voxels to add;')
            disp('press ENTER when finished.')
            [col_fix,row_fix,button] = ginput;
            row_fix=round(row_fix);
            col_fix=round(col_fix);
            for b=1:length(button)
                if button(b)==1
                    apo_segmented(row_fix(b), col_fix(b))=0;
                elseif button(b)==3
                    apo_segmented(row_fix(b), col_fix(b))=1;
                end
            end
            
            %visualize result
            figure(100)
            
            subplot(1,3,1)
            imagesc(loop_img_rgb), axis image
            [pixel_rows, pixel_cols] = ind2sub(size(apo_segmented), find(apo_segmented));
            set(gca, 'xlim', mask_cols, 'ylim', mask_rows);
            hold on
            plot(roi_edge_xy(:,2), roi_edge_xy(:,1), 'c');
            plot(pixel_cols, pixel_rows, 'm.', 'markersize', 5)
            
            loop_img_region=loop_img_rgb*0.75;
            loop_img_region(:,:,1)=loop_img_region(:,:,1)+.25*apo_segmented;
            subplot(1,3,2), imagesc(loop_img_region), axis image
            hold on
            plot(roi_edge_xy(:,2), roi_edge_xy(:,1), 'c');
            title('Interactive Window')
            xlabel(['Slice ' num2str(curr_slice) ' of Slices ' num2str(frst_slice) ' to ' num2str(last_slice)])
            set(gca, 'xlim', mask_cols, 'ylim', mask_rows);
            
            edge_pixels = bwboundaries(apo_segmented);
            edge_xy=edge_pixels{1};
            clear edge_xy_smooth
            edge_xy_smooth(:,2)=smooth(edge_xy(:,2), 'sgolay');
            edge_xy_smooth(:,1)=smooth(edge_xy(:,1), 'sgolay');
            edge_xy_smooth(end,:)=edge_xy_smooth(1,:);
            
            subplot(1,3,3)
            imagesc(loop_img_rgb)
            hold on
            plot(edge_xy_smooth(:,2), edge_xy_smooth(:,1), 'm.', 'markersize', 5);
            axis image
            set(gca, 'xlim', mask_cols, 'ylim', mask_rows);
            
        end
        
        % add to apo_segmented_final
        apo_segmented_final(:,:,curr_slice)=apo_segmented;
        
        %add to roi_surf matrices
        roi_curx = edge_xy_smooth(:,1);
        roi_cury = edge_xy_smooth(:,2);
        roi_curz = ones(size(roi_curx))*curr_slice;
        
        roi_surfx(curr_slice, :) = interp1(roi_curx, 1:(length(roi_curx)-1)/100:length(roi_curx), 'pchip');
        roi_surfy(curr_slice, :) = interp1(roi_cury, 1:(length(roi_cury)-1)/100:length(roi_cury), 'pchip');
        roi_surfz(curr_slice, :) = interp1(roi_curz, 1:(length(roi_curz)-1)/100:length(roi_curz), 'pchip');
        
        correct_slices=input('Do you want to correct any more slices (y/n)? ', 's');
        
        figure_position = get(gcf, 'position');
        clf
    end
    
end


%% create the mesh

%convert to dimensions of DTI image
roi_surfx=roi_surfx*dti_size(1)/sz_work_imag(1);
roi_surfy=roi_surfy*dti_size(2)/sz_work_imag(2);
roi_surfz=roi_surfz*dti_size(3)/sz_work_imag(3);

%resample to desired size:
roi_mesh = zeros(n_row, n_col, 6);
roi_mesh(:,:,1)=imresize(roi_surfx(frst_slice:last_slice,:), [n_row n_col]);
roi_mesh(:,:,2)=imresize(roi_surfy(frst_slice:last_slice,:), [n_row n_col]);
roi_mesh(:,:,3)=imresize(roi_surfz(frst_slice:last_slice,:), [n_row n_col]);

% find normal to mesh at each point:
mesh_row_vec = circshift(roi_mesh(:, :, 1:3), [0 -1 0]) - roi_mesh(:, :, 1:3);
mesh_col_vec = circshift(roi_mesh(:, :, 1:3), [-1 0 0]) - roi_mesh(:, :, 1:3);
mesh_row_vec(:, n_col, :) = mesh_row_vec(:, n_col-1, :);
mesh_col_vec(n_row, :, :) = mesh_col_vec(n_row-1, :, :);

roi_mesh(:, :, 4:6) = cross(mesh_row_vec, mesh_col_vec, 3);
roi_mesh(:, :, 4:6) = smooth3(roi_mesh(:, :, 4:6));
roi_norm = (roi_mesh(:, :, 4).^2 + roi_mesh(:, :, 5).^2 + roi_mesh(:, :,6).^2).^0.5;
roi_mesh(:, :, 4) = roi_mesh(:, :, 4)./roi_norm;
roi_mesh(:, :, 5) = roi_mesh(:, :, 5)./roi_norm;
roi_mesh(:, :, 6) = roi_mesh(:, :, 6)./roi_norm;
roi_mesh(:, :, 4:6) = smooth3(roi_mesh(:, :, 4:6));
roi_norm = (roi_mesh(:, :, 4).^2 + roi_mesh(:, :, 5).^2 + roi_mesh(:, :,6).^2).^0.5;
roi_mesh(:, :, 4) = roi_mesh(:, :, 4)./roi_norm;
roi_mesh(:, :, 5) = roi_mesh(:, :, 5)./roi_norm;
roi_mesh(:, :, 6) = roi_mesh(:, :, 6)./roi_norm;


%% save data files
save roi_mesh_file roi_mesh
save apo_segmentation_results apo_segmented_final apo_segmented_initial

%% plot mesh, if desired
close all

if exist('plot_options', 'var')
    
    % be sure not to plot unneeded stuff
    plot_options.plot_fibers=0;
    plot_options.plot_mask=0;
    
    fiber_visualizer(.75*(double(anat_image)/max(max(max(double(anat_image))))), plot_options, roi_mesh, [], []);
    
end


%% end the function
return;

