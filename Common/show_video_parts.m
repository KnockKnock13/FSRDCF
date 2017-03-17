function update_visualization_func = show_video(img_files, video_path, resize_image)
%SHOW_VIDEO
%   Visualizes a tracker in an interactive figure, given a cell array of
%   image file names, their path, and whether to resize the images to
%   half size or not.
%
%   This function returns an UPDATE_VISUALIZATION function handle, that
%   can be called with a frame number and a bounding box [x, y, width,
%   height], as soon as the results for a new frame have been calculated.
%   This way, your results are shown in real-time, but they are also
%   remembered so you can navigate and inspect the video afterwards.
%   Press 'Esc' to send a stop signal (returned by UPDATE_VISUALIZATION).
%
%   Joao F. Henriques, 2014
%   http://www.isr.uc.pt/~henriques/


	%store one instance per frame
	num_frames = numel(img_files);
	boxes = cell(num_frames,1);
    points = cell(num_frames,1);

	%create window
	[fig_h, axes_h, unused, scroll] = videofig(num_frames, @redraw, [], [], @on_key_press);  %#ok, unused outputs
	set(fig_h, 'IntegerHandle','off','Name', ['Tracker - ' video_path])%Number
	axis off;
	
	%image and rectangle handles start empty, they are initialized later
	im_h = [];
	rect_h = [];
    pointsC_h = [];
    pointsL_h = [];
    pointsR_h = [];
    pointsU_h = [];
    pointsD_h = [];
    pointsM_h = [];
    pointsT_h = [];
    
	
	update_visualization_func = @update_visualization;
	stop_tracker = false;
	

	function stop = update_visualization(frame, box, point)
		%store the tracker instance for one frame, and show it. returns
		%true if processing should stop (user pressed 'Esc').
		boxes{frame} = box;
        points{frame} = point;
		scroll(frame);        
		stop = stop_tracker;
	end

	function redraw(frame)
		%render main image
		im = imread([video_path img_files{frame}]);
% 		if size(im,3) > 1,
% 			im = rgb2gray(im);
% 		end
		if resize_image,
			im = imresize(im, 0.5);
		end
		
		if isempty(im_h),  %create image
			im_h = imshow(im, 'Border','tight', 'InitialMag',200, 'Parent',axes_h);
		else  %just update it
			set(im_h, 'CData', im)
		end
		
		%render target bounding box for this frame
		if isempty(rect_h),  %create it for the first time
			rect_h = rectangle('Position',[0,0,1,1], 'EdgeColor','g', 'Parent',axes_h);%
		end
		if ~isempty(boxes{frame}),
			set(rect_h, 'Visible', 'on', 'Position', boxes{frame});%
		else
			set(rect_h, 'Visible', 'off');
        end
%         
%         %
        hold on
        if(isempty(pointsC_h))  %create it for the first time
			pointsC_h = plot(0,0,'r+','Parent',axes_h);
        end       
         if(~isempty(pointsC_h))  %create it for the first time
			set(pointsC_h,'XData',points{frame}{1}(2),'YData',points{frame}{1}(1));
         else
             set(pointsC_h, 'Visible', 'off');
         end
%          
         if(isempty(pointsL_h))
             pointsL_h = scatter(0,0,'g<','Parent',axes_h);
         end
        if(~isempty(pointsL_h))
            set(pointsL_h,'XData',points{frame}{2}(2),'YData',points{frame}{2}(1));
        else
            set(pointsL_h, 'Visible', 'off');
        end
%         
        if(isempty(pointsR_h))
             pointsR_h = scatter(0,0,'b>','Parent',axes_h);
         end
        if(~isempty(pointsR_h))
            set(pointsR_h,'XData',points{frame}{3}(2),'YData',points{frame}{3}(1));
        else
            set(pointsR_h, 'Visible', 'off');
        end
%         
        if(isempty(pointsU_h))
             pointsU_h = scatter(0,0,'y^','Parent',axes_h);
         end
        if(~isempty(pointsU_h))
            set(pointsU_h,'XData',points{frame}{4}(2),'YData',points{frame}{4}(1));
        else
            set(pointsU_h, 'Visible', 'off');
        end
%         
         if(isempty(pointsD_h))
             pointsD_h = scatter(0,0,'cv','Parent',axes_h);
         end
        if(~isempty(pointsD_h))
            set(pointsD_h,'XData',points{frame}{5}(2),'YData',points{frame}{5}(1));
        else
            set(pointsD_h, 'Visible', 'off');
        end
%         
        if(isempty(pointsM_h))
             pointsM_h = scatter(0,0,'r*','Parent',axes_h);
         end
        if(~isempty(pointsM_h))
            set(pointsM_h,'XData',points{frame}{6}(2),'YData',points{frame}{6}(1));
        else
            set(pointsM_h, 'Visible', 'off');
        end
%         
        if(isempty(pointsT_h))
             pointsT_h = scatter(0,0,'md','Parent',axes_h);
         end
        if(~isempty(pointsT_h))
            set(pointsT_h,'XData',points{frame}{7}(2),'YData',points{frame}{7}(1));
        else
            set(pointsT_h, 'Visible', 'off');
        end
        hold off
	end

	function on_key_press(key)
		if strcmp(key, 'escape'),  %stop on 'Esc'
			stop_tracker = true;
		end
	end

end

