function [P0] = DLT(pts_wld, pts_img)
%DLT Summary of this function goes here
%   DLT algorithm

% if size(arrS,2) < 4
%     error('points is not enough');
% end
%construct A matrix
A = [];
for i = 1:size(pts_wld,2)
    % For each point, build matrix for its two equations
    current = [pts_wld(1,i) pts_wld(2,i) pts_wld(3,i) 1 0 0 0 0 -(pts_wld(1,i) * pts_img(1,i)) -(pts_wld(2,i) * pts_img(1,i)) -(pts_wld(3,i) * pts_img(1,i)) -pts_img(1,i);
              0 0 0 0 -pts_wld(1,i) -pts_wld(2,i) -pts_wld(3,i) -1 (pts_wld(1,i) * pts_img(2,i)) (pts_wld(2,i) * pts_img(2,i)) (pts_wld(3,i) * pts_img(2,i)) pts_img(2,i)];
    
    % Concatenate the equations to one matrix
    A = [A; current];
end

%find null space of A to solve h
[U,S,V] = svd(A);
p = V(:,end);
%display('A*h = '); A*h

% Respape h to 3 by 3 homography matrix
P0 = reshape(p, 4, 3)';

end

