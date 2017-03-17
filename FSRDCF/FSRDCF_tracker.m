function [ results ] =FSRDCF_tracker( params )

start_up =params.start_up;
if(start_up)
    tic();
    start_up_time = 0;
end
search_area_scale = params.search_area_scale;
output_sigma_factor = params.output_sigma_factor;
lambda = params.lambda;
learning_rate = params.learning_rate;
refinement_iterations = params.refinement_iterations;
filter_max_area = params.filter_max_area;
nScales = params.number_of_scales;
scale_step = params.scale_step;
interpolate_response = params.interpolate_response;
num_GS_iter = params.num_GS_iter;

features = params.t_features;

s_frames = params.s_frames;
pos = floor(params.init_pos);
target_sz = floor(params.wsize);

num_frames = numel(s_frames);

init_target_sz = target_sz;
% label = params.label;
offset = 1;
debug = params.debug;
%
global WW Y 
global update_visualization_func

%set the feature ratio to the feature-cell size
featureRatio = params.t_global.cell_size;

search_area = prod(init_target_sz / featureRatio * search_area_scale);

% when the number of cells are small, choose a smaller cell size
if isfield(params.t_global, 'cell_selection_thresh')
    if search_area < params.t_global.cell_selection_thresh * filter_max_area
        params.t_global.cell_size = min(featureRatio, max(1, ceil(sqrt(prod(init_target_sz * search_area_scale)/(params.t_global.cell_selection_thresh * filter_max_area)))));
        
        featureRatio = params.t_global.cell_size;
        search_area = prod(init_target_sz / featureRatio * search_area_scale);
    end
end

global_feat_params = params.t_global;

if search_area > filter_max_area
    currentScaleFactor = sqrt(search_area / filter_max_area);
else
    currentScaleFactor = 1.0;
end

% target size at the initial scale
base_target_sz = target_sz / currentScaleFactor;

%window size, taking padding into account
switch params.search_area_shape
    case 'proportional'
        sz = floor( base_target_sz * search_area_scale);     % proportional area, same aspect ratio as the target
    case 'square'
        sz = repmat(sqrt(prod(base_target_sz * search_area_scale)), 1, 2); % square area, ignores the target aspect ratio
    case 'fix_padding'
        sz = base_target_sz + sqrt(prod(base_target_sz * search_area_scale) + (base_target_sz(1) - base_target_sz(2))/4) - sum(base_target_sz)/2; % const padding
    otherwise
        error('Unknown "params.search_area_shape". Must be ''proportional'', ''square'' or ''fix_padding''');
end


% set the size to exactly match the cell size
sz = round(sz / featureRatio) * featureRatio;
use_sz = floor(sz/featureRatio);
if(mod(use_sz(1),2)==0)
    use_sz(1) = use_sz(1) - 1; 
    sz(1) = sz(1) - featureRatio;
end
if(mod(use_sz(2),2)==0)
    use_sz(2) = use_sz(2) - 1; 
    sz(2) = sz(2) - featureRatio;
end
%
shiftL = shift_m(use_sz(1));
shiftR = shift_m(use_sz(2));

% construct the label function
output_sigma = sqrt(prod(floor(base_target_sz/featureRatio))) * output_sigma_factor;
rg = circshift(-floor((use_sz(1)-1)/2):ceil((use_sz(1)-1)/2), [0 -floor((use_sz(1)-1)/2)]);
cg = circshift(-floor((use_sz(2)-1)/2):ceil((use_sz(2)-1)/2), [0 -floor((use_sz(2)-1)/2)]);
[rs, cs] = ndgrid( rg,cg);
y = exp(-0.5 * (((rs.^2 + cs.^2) / output_sigma^2)));
yf = fft2(y);
yf = shiftL*yf*shiftR;
Y = yf(:);

if interpolate_response == 1
    interp_sz = use_sz * featureRatio;
else
    interp_sz = use_sz;
end

% construct cosine window
cos_window = single(hann(use_sz(1))*hann(use_sz(2))');

% the search area size
support_sz = prod(use_sz);

% Calculate feature dimension
im = imread([params.path s_frames{1}]);
if size(im,3) == 3
    if all(all(im(:,:,1) == im(:,:,2)))
        colorImage = false;
    else
        colorImage = true;
    end
else
    colorImage = false;
end
% compute feature dimensionality
feature_dim = 0;
for n = 1:length(features)
    
    if ~isfield(features{n}.fparams,'useForColor')
        features{n}.fparams.useForColor = true;
    end;
    
    if ~isfield(features{n}.fparams,'useForGray')
        features{n}.fparams.useForGray = true;
    end;
    
    if (features{n}.fparams.useForColor && colorImage) || (features{n}.fparams.useForGray && ~colorImage)
        feature_dim = feature_dim + features{n}.fparams.nDim;
    end;
end;
% regulation window
[reg_win ,reg] = reg_window(use_sz, base_target_sz/featureRatio);
reg_win = sparse(reg_win);
WW = eval(['blkdiag(reg_win' repmat(',reg_win', 1, feature_dim-1) ');']);
%Scale handing
multires_pixel_template = zeros(sz(1), sz(2), size(im,3), nScales, 'uint8');
if nScales > 0
    scale_exp = (-floor((nScales-1)/2):ceil((nScales-1)/2));
    
    scaleFactors = scale_step .^ scale_exp;
    
    %force reasonable scale changes
    min_scale_factor = scale_step ^ ceil(log(max(5 ./ sz)) / log(scale_step));
    max_scale_factor = scale_step ^ floor(log(min([size(im,1) size(im,2)] ./ base_target_sz)) / log(scale_step));
end

if interpolate_response >= 3
    % Pre-computes the grid that is used for socre optimization
    ky = circshift(-floor((use_sz(1) - 1)/2) : ceil((use_sz(1) - 1)/2), [1, -floor((use_sz(1) - 1)/2)]);
    kx = circshift(-floor((use_sz(2) - 1)/2) : ceil((use_sz(2) - 1)/2), [1, -floor((use_sz(2) - 1)/2)])';
    newton_iterations = params.newton_iterations;
end
xlf_old_b = zeros(support_sz, feature_dim);
if(debug)
    fig = figure(2);
    colormap hsv;
end
rect_position = zeros(num_frames, 4);
if(~start_up)
    time = 0;
end
for frame = 1:num_frames,
    %load image
    im = imread([params.path s_frames{frame}]);
    if size(im,3) > 1 && colorImage == false
        im = im(:,:,1);
        multires_pixel_template = zeros(sz(1), sz(2), size(im,3), nScales, 'uint8');
    end
    
    if(~start_up)
       tic();
    end
    
    if(frame>1)
        old_pos = inf(size(pos));
        iter = 1;
        %translation search
        while iter <= refinement_iterations && any(old_pos ~= pos)
            % Get multi-resolution image
            for scale_ind = 1:nScales
                multires_pixel_template(:,:,:,scale_ind) = ...
                    get_pixels(im, pos-offset, round(sz*currentScaleFactor*scaleFactors(scale_ind)), sz);
            end
            
            xt = bsxfun(@times,get_features(multires_pixel_template,features,global_feat_params),cos_window);
            
            xtf = fft2(xt);%/unitary
            
            responsef = permute(sum(bsxfun(@times, hf, xtf), 3), [1 2 4 3]);
            
            % if we undersampled features, we want to interpolate the
            % response so it has the same size as the image patch
            if interpolate_response == 2
                % use dynamic interp size
                interp_sz = floor(size(y) * featureRatio * currentScaleFactor);
            end
            responsef_padded = resizeDFT2(responsef, interp_sz);
            
            % response
            response = ifft2(responsef_padded, 'symmetric');%*unitary
            
            % find maximum
            if interpolate_response == 3
                error('Invalid parameter value for interpolate_response');
            elseif interpolate_response == 4
                [disp_row, disp_col, sind] = resp_newton(response, responsef_padded, newton_iterations, ky, kx, use_sz);
            else
                [row, col, sind] = ind2sub(size(response), find(response == max(response(:)), 1));
                disp_row = mod(row - 1 + floor((interp_sz(1)-1)/2), interp_sz(1)) - floor((interp_sz(1)-1)/2);
                disp_col = mod(col - 1 + floor((interp_sz(2)-1)/2), interp_sz(2)) - floor((interp_sz(2)-1)/2);
            end
            
            % calculate translation
            switch interpolate_response
                case 0
                    translation_vec = round([disp_row, disp_col] * featureRatio * currentScaleFactor * scaleFactors(sind));
                case 1
                    translation_vec = round([disp_row, disp_col] * currentScaleFactor * scaleFactors(sind));
                case 2
                    translation_vec = round([disp_row, disp_col] * scaleFactors(sind));
                case 3
                    translation_vec = round([disp_row, disp_col] * featureRatio * currentScaleFactor * scaleFactors(sind));
                case 4
                    translation_vec = round([disp_row, disp_col] * featureRatio * currentScaleFactor * scaleFactors(sind));
            end
            
            % set the scale
            currentScaleFactor = currentScaleFactor * scaleFactors(sind);
            % adjust to make sure we are not to large or to small
            if currentScaleFactor < min_scale_factor
                currentScaleFactor = min_scale_factor;
            elseif currentScaleFactor > max_scale_factor
                currentScaleFactor = max_scale_factor;
            end
            
            % update position
            old_pos = pos;
            pos = pos + translation_vec;
            
            iter = iter + 1;
        end
    end
    if(params.visualization) 
        visual = pos-floor(target_sz/2);
        update_visualization_func(frame,[visual(2),visual(1),base_target_sz(2)*currentScaleFactor,base_target_sz(1)*currentScaleFactor]);
    end
    % extract training sample image region
    pixels = get_pixels(im,pos-offset,round(sz*currentScaleFactor),sz);
    
    % extract features and do windowing
    xl = bsxfun(@times,get_features(pixels,features,global_feat_params),cos_window);
    xlf = fft2(xl);%/unitary
    
    xlf_reshaped = reshape(xlf, [support_sz, feature_dim]);
    
    if(frame>1)
        autoc = sparse(double(conj(xlf_reshaped(:,:)).*xlf_reshaped(:,:)));
        xlf_old_a = (1-learning_rate)*xlf_old_a + learning_rate*autoc;
        xlf_old_b = (1-learning_rate)*xlf_old_b + learning_rate*xlf_reshaped;
    end
    if(frame == 1)
        if(params.visualization)
            visual = pos-floor(target_sz/2);
            update_visualization_func(frame,[visual(2),visual(1),base_target_sz(2)*currentScaleFactor,base_target_sz(1)*currentScaleFactor]);
        end
        xxf = sparse(double(conj(xlf_reshaped(:,:)).*xlf_reshaped(:,:)));
        A = diag(xxf(:)) + WW;
        b = bsxfun(@times, Y, conj(xlf_reshaped(:,:)));
        b = double(b(:));
        w = A\b;
        xlf_old_b = xlf_reshaped;
        xlf_old_a = xxf;
    end
     w = gauss_seidel( xlf_old_a,xlf_old_b, w,num_GS_iter);
     hf = reshape(single(reshape(full(w), [support_sz, feature_dim])), [use_sz, feature_dim]);
%      for i=1: feature_dim
%          hf(:,:,i) = shiftL*hf(:,:,i)*shiftR;
%      end
    target_sz = floor(base_target_sz * currentScaleFactor);
    
    %save position and calculate FPS
    rect_position(frame,:) = [pos([2,1]) - floor(target_sz([2,1])/2), target_sz([2,1])];
    if(debug)
     h = ifft2(hf);%
     figure(fig);
     surf(real(h(:,:,1)))
    end
    if(start_up)
     start_up_time = start_up_time + toc();
    else
      time = time + toc();
    end
    if(start_up)
        results.res = [];
        results.fps = start_up_time;
        return ;
    end
    if(frame ==8)
        paperim = 1;
    end
end

fps = numel(s_frames) / time;
results.type = 'rect';
% results.fps = fps;
results.res = rect_position;
results.fps = fps;
end

