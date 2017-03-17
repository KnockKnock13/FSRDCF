function [ imageFile, label ] = loadImgLabel( num )
%load image imformation according to the parameter num
%   num --> the index of dataset in the seqs or the name of the dataset
%   imageFile --> struct(name, path, startFrame, endFrame, nz, ext, init_rect, len, frames)
%                 frames --> cell(len, 1) are the paths of every image
%   label --> matrix(len,4) are the groundtruth
%             [a b c d] (a,b) is the coordinate of the TOP-LEFT of the
%             rectangel,(c,d) are the width and height of the rectangel

    global seqs
    
    if(~isa(num,'double'))
        for i=1:length(seqs)
            if(strcmp(seqs{i}.name,num))
                num=i;
                break;
            end
        end
    end
    
    imageFile = seqs{num};
    nz	= strcat('%0',num2str(imageFile.nz),'d');
    imageFile.len = imageFile.endFrame - imageFile.startFrame + 1;
    imageFile.frames = cell(imageFile.len,1);
    %% image name list to read later
    for i=1:imageFile.len
            image_no = imageFile.startFrame + (i-1);
            id = sprintf(nz,image_no);
            imageFile.frames{i} = strcat(id,'.',imageFile.ext);
    end
    %% ground truth files
    label = load(strcat(imageFile.path,'..\','groundtruth_rect.txt'));
    %label(:,[1,2]) = label(:,[1,2]) + floor(label(:,[3,4])/2);
end

