function [ meanOfd2, sumOfd2, varErr ] = testModel( para, img, P, xc )
%testModel Summary of this function goes here
%   Detailed explanation goes here
dX = para.dX;
dY = para.dY;
n_sq_x = para.n_sq_x;
n_sq_y = para.n_sq_y;

%% define testing world points
% middle seven
mid_pts = [(0:n_sq_x)*dX; ones(1,(n_sq_x+1))*(n_sq_y+1)*dY; zeros(1,(n_sq_x+1))];
% left seven
left_pts = [(0:n_sq_x)*dX; ones(1,(n_sq_x+1))*(-1)*dY; zeros(1,(n_sq_x+1))];
% right seven
right_pts = [(0:n_sq_x)*dX; ones(1,(n_sq_x+1))*(n_sq_y+1)*dY; ones(1,(n_sq_x+1))*(n_sq_y+2)*dY];
test_wldpts = [left_pts mid_pts right_pts];

test_wldpts = [test_wldpts; ones(1,size(test_wldpts,2))];
%% mapping
test_imgpts = P*test_wldpts;
% dehomo
test_imgpts = test_imgpts./repmat(test_imgpts(3,:),3,1);

figure(2);
imshow(img);
axis on; hold on
plot(test_imgpts(1,:),test_imgpts(2,:),'r+');
plot(xc(1,:), xc(2,:),'b+');
drawnow;
hold off;

%% compute error
diff_pts = xc-test_imgpts(1:2,:);
sumOfd2 = sum(sum(diff_pts.^2));
meanOfd2 = mean(sum(diff_pts.^2));
varErr = var(sum(diff_pts.^2));
end

