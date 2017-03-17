function [ patch ] = crop( I, rect,XWorldLimits, YWorldLimits )
%crop an image path 
%   
    if(nargin<4)
        patch = imcrop(I, [rect(1:2) rect(3:4)-1]);
        s = size(patch);
        error = rect([3,4])-s([2,1]);
        if(s(1)==0 && s(2)==0)
            patch = zeros(rect([4,3]));
        else
            if(any(error>0))
                for i=1:error(1)
                    patch = [patch patch(:,s(2))];
                end
                for i=1:error(2)
                    patch = [patch;patch(s(1),:)];
                end
            end
        end
    else
        patch = imcrop(XWorldLimits, YWorldLimits, I, [rect(1:2) rect(3:4)-1]);
        s = size(patch);
        error = rect([3,4])-s([2,1]);
        if(any(error~=0))
            patch = imcrop(XWorldLimits, YWorldLimits, I, [rect(1:2) rect(3:4)-1+error]);
        end
    end
end

