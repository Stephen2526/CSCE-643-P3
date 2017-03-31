function [ xc ] = pickNPts( para, img, N )
%pickNPts Summary of this function goes here
%   manually pick N corner points

% load parameters
wintx = para.wintx;
winty = para.winty;
% initial x
x = zeros(2,N);

fprintf(1,'Using (wintx,winty)=(%d,%d) - Window size = %dx%d \n',wintx,winty,2*wintx+1,2*winty+1);

figure(1);
imshow(img);
axis on; hold on;
for count = 1:N
    [xi,yi] = ginput4(1);
    [xxi] = cornerfinder([xi;yi],img,winty,wintx);
    xi = xxi(1);
    yi = xxi(2);
    %figure(1);
    plot(xi,yi,'+','color',[ 1.000 0.0 0.0 ],'linewidth',2);
    plot(xi + [wintx+.5 -(wintx+.5) -(wintx+.5) wintx+.5 wintx+.5],yi + [winty+.5 winty+.5 -(winty+.5) -(winty+.5)  winty+.5],'-','color',[ 1.000 0.0 0.0 ],'linewidth',2);
    x(1,count) = xi;
    x(2,count) = yi;
    drawnow;
end
hold off;

% refine corners
[xc,~,~,~] = cornerfinder([x(1,:);x(2,:)],img,winty,wintx);

end

