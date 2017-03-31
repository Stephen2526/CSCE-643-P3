function [ norPts, T ] = normalise2dpts( pts )
%NORMALISE Summary of this function goes here
%   Detailed explanation goes here
% verify input matrix size
if size(pts,1) ~= 3
    error('input matrix must be 3*N');
end
% scale the thirt coord to 1
pts(1,:) = pts(1,:)./pts(3,:);
pts(2,:) = pts(2,:)./pts(3,:);
pts(3,:) = 1;

% find centroid of the points
c = mean(pts(1:2,:),2);
%c
% translate points to make centriod the center
trnPts(1,:) = pts(1,:)-c(1);
trnPts(2,:) = pts(2,:)-c(2);

% mean dist to centroid
dists = sqrt(trnPts(1,:).^2+trnPts(2,:).^2);
meandist = mean(dists(:));
%meandist

scale = meandist/sqrt(2);

T = [1/scale    0     -c(1)/scale;
     0       1/scale  -c(2)/scale;
     0          0         1      ];
%T
norPts = T*pts;
end

