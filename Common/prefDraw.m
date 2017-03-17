clear all
addpath('../../expriment analysis/SRDCF');
addpath('./rstEval');
%%
global seqs
seqs = configSeqs;
anno = [];

thresholdSetOverlap = 0:0.05:1;
thresholdSetError = 0:50;
idxNum = 1;

successNumOverlap = zeros(idxNum,length(thresholdSetOverlap));
successNumErr = zeros(idxNum,length(thresholdSetError));

for i = 1:100
    [~,label_R] = loadImgLabel(i);
    anno = [anno;label_R];
end
numImg = size(anno,1);

 rect = load('label.txt');
 rect = rect(end-numImg+1:end,:);
 overlap=calOverlap(rect, anno);
 center = [rect(:,1)+(rect(:,3)-1)/2 rect(:,2)+(rect(:,4)-1)/2];
%  center = rect(:,[1,2]);%[rect(:,1)+(rect(:,3)-1)/2 rect(:,2)+(rect(:,4)-1)/2];
%  rect(:,[1,2]) = [rect(:,1)-(rect(:,3)-1)/2 rect(:,2)-(rect(:,4)-1)/2];
 centerG = [anno(:,1)+(anno(:,3)-1)/2 anno(:,2)+(anno(:,4)-1)/2];
 errCenter = sqrt(sum(((center - centerG).^2),2));
 errCoverage = calcRectInt(rect,anno);
 
 aveErrCoverage = sum(errCoverage)/numImg;

aveErrCenter = sum(errCenter)/numImg;
for idx = 1:idxNum
     for tIdx=1:length(thresholdSetOverlap)
        successNumOverlap(idx,tIdx) = sum(errCoverage >thresholdSetOverlap(tIdx));
     end
     for tIdx=1:length(thresholdSetError)
        successNumErr(idx,tIdx) = sum(errCenter <= thresholdSetError(tIdx));
     end
end
aveSuccessRatePlot = successNumOverlap/(numImg+eps);
rankO = mean(aveSuccessRatePlot);
aveSuccessRatePlotErr = successNumErr/(numImg+eps);
rankC = aveSuccessRatePlotErr(20)
precision = load('precision.txt');
precision20 = sum(precision(:,21))/size(precision,1)