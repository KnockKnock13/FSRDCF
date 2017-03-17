clear all

addpath('../Common');
save = 0;
path = 'C:\Users\KnockKnock\Pictures\data_seq_row';%the folder that includes the sub-folders OTB50 and OTB100
%%
global seqs
seqs = configSeqs(path);
global update_visualization_func
for i = 1:1
[imgFile,label_R] = loadImgLabel(i);
imgFile.init_rect = label_R(1,:);
update_visualization_func = show_video(imgFile.frames, imgFile.path, false);

% Run FSRDCF
 results = run_FSRDCF(imgFile,label_R);
end