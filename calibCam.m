% main part of calibration
% Using golden stardard algorithm to estimate P matrix
%
%clear all;  
%% load images
imgs = loadImg;

%% obtain points data
% points on the image plane

% strcut array to store parameters
para.wintx = 5;
para.winty = 5;
para.n_sq_x = 6;
para.n_sq_y = 4;
para.dX = 30; %mm
para.dY = 30; %mm
para.num_plane = 2; %number of planes

% extract points on image and 3D world
[x,~,X] = extractCorner(para, imgs);
%extractCorner;

for i=1:size(x,2)
    %% points data normalzation
    %  add homo coor
    x = [x{i};ones(1,size(x{i},2))];
    X = [X{i};ones(1,size(X{i},2))];
    [nor_pts2d, T_2d] = normalise2dpts(x);
    [not_pts3d, T_3d] = normalise3dpts(X);
    %% obtain an initial solution using DLT
    P0 = DLT(not_pts3d, nor_pts2d)
    %% Minimize geometric error
    nor_P = MLE(P0, not_pts3d, nor_pts2d)
    %% Denormalization
    P = T_2d\nor_P*T_3d
end
