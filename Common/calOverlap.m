function overlap=calOverlap(A,B)
% each row of rt and gt is a rectangle [left,top,width,height]
% the function calculate the overlap score between these two rectangles

leftA = A(:,1);
topA = A(:,2);
rightA = leftA + A(:,3) - 1;
bottomA = topA + A(:,4) - 1;

leftB = B(:,1);
topB = B(:,2);
rightB = leftB + B(:,3) - 1;
bottomB = topB + B(:,4) - 1;

tmp = (max(0, min(rightA, rightB) - max(leftA, leftB)+1 )) .* (max(0, min(bottomA, bottomB) - max(topA, topB)+1 ));
areaA = A(:,3) .* A(:,4);
areaB = B(:,3) .* B(:,4);
overlap = tmp./(areaA+areaB-tmp);