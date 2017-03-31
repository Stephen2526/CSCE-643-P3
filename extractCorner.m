function [ x, xx, X ] = extractCorner( para, imgs)
%extractCorner Summary of this function goes here
%   extract point positions on image plane
%% set up parameters
wintx = para.wintx;
winty = para.winty;
n_sq_x = para.n_sq_x;
n_sq_y = para.n_sq_y; 
dX = para.dX; %mm
dY = para.dY; %mm
num_plane = para.num_plane;

%% initialize x
x = cell(1:size(imgs,2));
X = cell(1:size(imgs,2));
xx = cell(1:size(imgs,2));
if num_plane == 1
    x(:) = {zeros(2,(n_sq_x+1)*(n_sq_y+1))};
    xx(:) = {zeros(2,(n_sq_x+1)*(n_sq_y+1))};
    X(:) = {zeros(3,(n_sq_x+1)*(n_sq_y+1))};
else %==2
    x(:) = {zeros(2,2*(n_sq_x+1)*(n_sq_y+1))};
    xx(:) = {zeros(2,2*(n_sq_x+1)*(n_sq_y+1))};
    X(:) = {zeros(3,2*(n_sq_x+1)*(n_sq_y+1))};
end

%% compute x and X
for i = 1:size(imgs,2)
    if num_plane == 1
        [x{i},xx{i},X{i}] = cornExtract(para, imgs{i});
    else %==2
        [x_l,xx_l,~] = cornExtract(para, imgs{i});
        [x_r,xx_r,~] = cornExtract(para, imgs{i});
        drawMarker(para,x_l,x_r,imgs{i});
        x{i} = [x_l x_r];
        xx{i} = [xx_l xx_r];
        
        % compute X of orthogonal pattern
        % left board
        Np = (n_sq_x+1)*(n_sq_y+1);
        X_l = reshape(((0:n_sq_x)*dX)'*ones(1,n_sq_y+1),Np,1)';
        Y_l = reshape(ones(n_sq_x+1,1)*(0:n_sq_y)*dY,Np,1)';
        Z_l = zeros(1,Np);
        % right board
        X_r = reshape(((0:n_sq_x)*dX)'*ones(1,n_sq_y+1),Np,1)';
        Y_r = ones(1,Np)*(n_sq_y+1)*dY;
        Z_r = reshape(ones(n_sq_x+1,1)*(1:n_sq_y+1)*dY,Np,1)';
        X{i} = [X_l X_r;Y_l Y_r;Z_l Z_r];
    end
end

end

