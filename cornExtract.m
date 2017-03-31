function [ x, xx, X ] = cornExtract( para, img )
%cornExtract Summary of this function goes here
%   Detailed explanation goes here
% extract corner in the specific area
% INPUT: an image
% OUTPUT: x-corner coordinates on image plane, X-corner coordinates in 3D
% world
%% set up parameters
% wintx = 5;
% winty = 5;
% n_sq_x = 5;
% n_sq_y = 7; 
% dX = 30; %mm
% dY = 30; %mm

wintx = para.wintx;
winty = para.winty;
n_sq_x = para.n_sq_x;
n_sq_y = para.n_sq_y; 
dX = para.dX; %mm
dY = para.dY; %mm


fprintf(1,'Using (wintx,winty)=(%d,%d) - Window size = %dx%d - square size = (%d,%d)\n',wintx,winty,2*wintx+1,2*winty+1,dX,dY);

%% manually select four boundary corners
figure(1);
imshow(img);
axis on;


title('Click on the four extreme corners of the rectangular pattern (first corner = origin)...');
disp('Click on the four extreme corners of the rectangular complete pattern (the first clicked corner is the origin)...');

x= [];y = [];
figure(1); hold on;
for count = 1:4
    [xi,yi] = ginput4(1);
    [xxi] = cornerfinder([xi;yi],img,winty,wintx);
    xi = xxi(1);
    yi = xxi(2);
    figure(1);
    plot(xi,yi,'+','color',[ 1.000 0.0 0.0 ],'linewidth',2);
    plot(xi + [wintx+.5 -(wintx+.5) -(wintx+.5) wintx+.5 wintx+.5],yi + [winty+.5 winty+.5 -(winty+.5) -(winty+.5)  winty+.5],'-','color',[ 1.000 0.0 0.0 ],'linewidth',2);
    x = [x;xi];
    y = [y;yi];
    plot(x,y,'-','color',[ 1.000 0.0 0.0 ],'linewidth',2);
    drawnow;
end;
plot([x;x(1)],[y;y(1)],'-','color',[ 1.000 0.0 0.0 ],'linewidth',2);
drawnow;
hold off;
%w = waitforbuttonpress;
%% refine four corners
[Xc,good,bad,type] = cornerfinder([x';y'],img,winty,wintx); % the four corners
x = Xc(1,:)';
y = Xc(2,:)';

% figure(1);hold on;
% plot([x;x(1)],[y;y(1)],'-','color',[ 0.0 0.0 1.0 ],'linewidth',2);
% drawnow;
% hold off;

% Sort the corners:
x_mean = mean(x);
y_mean = mean(y);
x_v = x - x_mean;
y_v = y - y_mean;

theta = atan2(-y_v,x_v);
%[junk,ind] = sort(theta);

[junk,ind] = sort(mod(theta-theta(1),2*pi));

%ind = ind([2 3 4 1]);

ind = ind([4 3 2 1]); %-> New: the Z axis is pointing uppward

x = x(ind);
y = y(ind);
% x1= x(1); x2 = x(2); x3 = x(3); x4 = x(4);
% y1= y(1); y2 = y(2); y3 = y(3); y4 = y(4);
% 
% 
% % Find center:
% p_center = cross(cross([x1;y1;1],[x3;y3;1]),cross([x2;y2;1],[x4;y4;1]));
% x5 = p_center(1)/p_center(3);
% y5 = p_center(2)/p_center(3);
% 
% % center on the X axis:
% x6 = (x3 + x4)/2;
% y6 = (y3 + y4)/2;
% 
% % center on the Y axis:
% x7 = (x1 + x4)/2;
% y7 = (y1 + y4)/2;
% 
% % Direction of displacement for the X axis:
% vX = [x6-x5;y6-y5];
% vX = vX / norm(vX);
% 
% % Direction of displacement for the X axis:
% vY = [x7-x5;y7-y5];
% vY = vY / norm(vY);
% 
% % Direction of diagonal:
% vO = [x4 - x5; y4 - y5];
% vO = vO / norm(vO);
% 
% delta = 30;
% 
% 
% figure(1); 
% imshow(img);
% axis on;
% hold on;
% plot([x;x(1)],[y;y(1)],'g-');
% plot(x,y,'og');
% hx=text(x6 + delta * vX(1) ,y6 + delta*vX(2),'X');
% set(hx,'color','g','Fontsize',14);
% hy=text(x7 + delta*vY(1), y7 + delta*vY(2),'Y');
% set(hy,'color','g','Fontsize',14);
% hO=text(x4 + delta * vO(1) ,y4 + delta*vO(2),'O','color','g','Fontsize',14);
% hold off;

%% extract the cornor points within the region
%  compute initial positions of corners
% Compute the inside points through computation of the planar homography (collineation)

% a00 = [x(1);y(1);1];
% a10 = [x(2);y(2);1];
% a11 = [x(3);y(3);1];
% a01 = [x(4);y(4);1];

a00 = [x(4);y(4);1];
a10 = [x(3);y(3);1];
a11 = [x(2);y(2);1];
a01 = [x(1);y(1);1];
% Compute the planar collineation: (return the normalization matrix as well)

[Homo,Hnorm,inv_Hnorm] = compute_homography([a00 a10 a11 a01],[0 1 1 0;0 0 1 1;1 1 1 1]);


% Build the grid using the planar collineation:

x_l = ((0:n_sq_x)'*ones(1,n_sq_y+1))/n_sq_x;
y_l = (ones(n_sq_x+1,1)*(0:n_sq_y))/n_sq_y;
pts = [x_l(:) y_l(:) ones((n_sq_x+1)*(n_sq_y+1),1)]';

XX = Homo*pts;
XX = XX(1:2,:) ./ (ones(2,1)*XX(3,:));

figure(1);imshow(img);hold on;axis on;
plot(XX(1,:),XX(2,:),'r+');
title('The initial guessed postion of corners');
hold off;
% points on the straignt lines
xx = XX;

%  refine the position of corners
disp('Refine the position of corners...');
grid_pts = cornerfinder(XX,img,winty,wintx);

%  mark the corners on the image
figure(1);hold on;
plot(grid_pts(1,:),grid_pts(2,:),'b+');
title('The measured and corrected postion of corners');
hold off;
%savefig('before_undisto.fig');

%  save position of corners in image space
%grid_pts = grid_pts - 1; % subtract 1 to bring the origin to (0,0)
x = grid_pts;

%% calculate corner position in 3D space
Np = (n_sq_x+1)*(n_sq_y+1);
Xi = reshape(((0:n_sq_x)*dX)'*ones(1,n_sq_y+1),Np,1)';
Yi = reshape(ones(n_sq_x+1,1)*(0:n_sq_y)*dY,Np,1)';
Zi = zeros(1,Np);

Xgrid = [Xi;Yi;Zi];
X = Xgrid;
disp('All the points coordinates saved (x-on the image, X-3D)');
end

