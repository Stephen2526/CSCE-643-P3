function [ ] = drawMarker( para, x_l, x_r, img )
%drawMarker Summary of this function goes here
%   Detailed explanation goes here
n_sq_x = para.n_sq_x;
n_sq_y = para.n_sq_y;

%4-----1
%|     |
%|     |
%3-----2
x1= x_l(1,end-n_sq_x); x2 = x_l(1,end); x3 = x_l(1,n_sq_x+1); x4 = x_l(1,1);
y1= x_l(2,end-n_sq_x); y2 = x_l(2,end); y3 = x_l(2,n_sq_x+1); y4 = x_l(2,1);


% Find center:
p_center = cross(cross([x1;y1;1],[x3;y3;1]),cross([x2;y2;1],[x4;y4;1]));
x5 = p_center(1)/p_center(3);
y5 = p_center(2)/p_center(3);

% center on the X axis:
x6 = (x3 + x4)/2;
y6 = (y3 + y4)/2;

% center on the Y axis:
x7 = (x1 + x4)/2;
y7 = (y1 + y4)/2;

% Direction of displacement for the X axis:
vX = [x6-x5;y6-y5];
vX = vX / norm(vX);

% Direction of displacement for the X axis:
vY = [x7-x5;y7-y5];
vY = vY / norm(vY);

% Direction of diagonal:
vO = [x4 - x5; y4 - y5];
vO = vO / norm(vO);

delta = 30;


figure(1); 
imshow(img);
axis on;
hold on;
%draw axis
plot([x1;x4;x3],[y1;y4;y3],'g-');
plot(x4,y4,'og');
hx=text(x6 + delta * vX(1) ,y6 + delta*vX(2),'X');
set(hx,'color','g','Fontsize',14);
hy=text(x7 + delta*vY(1), y7 + delta*vY(2),'Y');
set(hy,'color','g','Fontsize',14);
hO=text(x4 + delta * vO(1) ,y4 + delta*vO(2),'O','color','g','Fontsize',14);

%draw corners
plot(x_l(1,:),x_l(2,:),'r+');
plot(x_r(1,:),x_r(2,:),'r+');
title('The detected postion of corners');
hold off;

end

