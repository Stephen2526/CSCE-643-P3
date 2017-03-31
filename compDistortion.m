function [para_opt] = compDistortion(para)
% compute the radial distortion coefficients and remove the distortion

% %% Extract points for using
% bond_pts_x = [];
% bond_pts_x = [bond_pts_x x(1, 1:n_sq_x)];% left |
% bond_pts_x = [bond_pts_x x(1, n_sq_x+1 : n_sq_x+1 : (n_sq_x+1)*(n_sq_y+1)-1)];% bottom _
% bond_pts_x = [bond_pts_x x(1, (n_sq_x+1)*(n_sq_y+1) : -1 : (n_sq_x+1)*(n_sq_y+1)-(n_sq_x-1))];% right |
% bond_pts_x = [bond_pts_x x(1, (n_sq_x+1)*(n_sq_y+1)-(n_sq_x) : -(n_sq_x+1) : (n_sq_x+2))];% bottom _
% %bond_pts_x = bond_pts_x + 1; %move origin to (1 1)
% 
% bond_pts_y = [];
% bond_pts_y = [bond_pts_y x(2, 1:n_sq_x)];% left |
% bond_pts_y = [bond_pts_y x(2, n_sq_x+1 : n_sq_x+1 : (n_sq_x+1)*(n_sq_y+1)-1)];% bottom _
% bond_pts_y = [bond_pts_y x(2, (n_sq_x+1)*(n_sq_y+1) : -1 : (n_sq_x+1)*(n_sq_y+1)-(n_sq_x-1))];% right |
% bond_pts_y = [bond_pts_y x(2, (n_sq_x+1)*(n_sq_y+1)-(n_sq_x) : -(n_sq_x+1) : (n_sq_x+2))];% bottom _
% %bond_pts_y = bond_pts_y + 1; %move origin to (1 1)


% % superpose points on the image
% cor_ind = [1 n_sq_x+1 n_sq_x+1+n_sq_y n_sq_x+1+n_sq_y+n_sq_x];

% figure(2);
% imshow(img);
% axis on; hold on;
% plot(bond_pts_x(cor_ind), bond_pts_y(cor_ind),'og');
% plot(bond_pts_x([cor_ind 1]), bond_pts_y([cor_ind 1]),'g-');
% plot(bond_pts_x,bond_pts_y,'b+');
% title('Image with radial distortion');
% hold off;

%% load image
imgs = loadImg();
img = imgs{1};

%% extract points
[x,xx,~] = extractCorner(para,imgs);
x = x{1}
xx = xx{1}
%% optimization process 
%  parameter set: k1 k2 k3 k4 xc yc
%  give an initial value of parameters
para_ini = [1.0e-5 1.00e-6 0.0000 0.0000 0.0000 0.0000 floor(size(img,2)/2) floor(size(img,1)/2)];
para_ini
%  minimization
options = optimoptions(@lsqnonlin,'MaxFunctionEvaluations',5000,'Display','iter-detailed');
options.Algorithm = 'levenberg-marquardt';
options.StepTolerance = 1.000000e-12;
options.FunctionTolerance = 1.000000e-8;
%para_opt = lsqnonlin(@funcDist1,para_ini,[],[],options,cor_ind,n_sq_x,n_sq_y,bond_pts_x,bond_pts_y);
para_opt = lsqnonlin(@funcDist2,para_ini,[],[],options,x,xx);
para_opt


%% undistort iamge using optimaized parameters
f = 1;
p1 = 0;
p2 = 0;
ppx = para_opt(end-1);
ppy = para_opt(end);
k1 = para_opt(1);
k2 = para_opt(2);
k3 = para_opt(3);
k4 = para_opt(4);
k5 = para_opt(5);
k6 = para_opt(6);
undisto_img = undistortimage(img, f, ppx, ppy, k1, k2, k3, k4, k5, k6, p1, p2);
figure(3);
imshow(undisto_img);

end
%% cost function definition (consider the boundary)
function [dists_sqrt] = funcDist1(para_opt, cor_ind, n_sq_x, n_sq_y, bond_pts_x, bond_pts_y)

%cor_pts = [bond_pts_x(cor_ind);bond_pts_y(cor_ind);ones(1,4)]; % four corner points(00,10,11,01)
cor_pts = [bond_pts_x(cor_ind);bond_pts_y(cor_ind)]; % four corner points(00,10,11,01)

%  straignt lines in homogeneous form
%str_lines=[cross(cor_pts(:,1), cor_pts(:,2)), cross(cor_pts(:,2), cor_pts(:,3)), cross(cor_pts(:,3), cor_pts(:,4)), cross(cor_pts(:,4), cor_pts(:,1))];

num_pts = size(bond_pts_x,2);
dists_sqrt = zeros(num_pts,1);

r2 = (bond_pts_x-para_opt(5)).^2 + (bond_pts_y-para_opt(6)).^2;
r = sqrt(r2); %radial distance

ref_pts_x = bond_pts_x + (para_opt(1)*r+para_opt(2)*(r.^2)+para_opt(3)*(r.^3)+para_opt(4)*(r.^4)).*(bond_pts_x-para_opt(5));
ref_pts_y = bond_pts_y + (para_opt(1)*r+para_opt(2)*(r.^2)+para_opt(3)*(r.^3)+para_opt(4)*(r.^4)).*(bond_pts_y-para_opt(6));
ref_pts = [ref_pts_x;ref_pts_y];

for i = 1:num_pts
    if i>1 && i<n_sq_x+1
        denomr_l = sqrt((cor_pts(1,2)-cor_pts(1,1))^2 + (cor_pts(2,2)-cor_pts(2,1))^2);
        dists_sqrt(i) = sqrt(abs((cor_pts(1,1)-ref_pts(1,i))*(cor_pts(2,2)-cor_pts(2,1)) - (cor_pts(2,1)-ref_pts(2,i))*(cor_pts(1,2)-cor_pts(1,1)))/denomr_l);
    elseif i>n_sq_x+1 && i<n_sq_x+1+n_sq_y
        denomr_b = sqrt((cor_pts(1,3)-cor_pts(1,2))^2 + (cor_pts(2,3)-cor_pts(2,2))^2);
        dists_sqrt(i) = sqrt(abs((cor_pts(1,2)-ref_pts(1,i))*(cor_pts(2,3)-cor_pts(2,2)) - (cor_pts(2,2)-ref_pts(2,i))*(cor_pts(1,3)-cor_pts(1,2)))/denomr_b);
    elseif i>n_sq_x+1+n_sq_y && i<n_sq_x+1+n_sq_y+n_sq_x
        denomr_r = sqrt((cor_pts(1,4)-cor_pts(1,3))^2 + (cor_pts(2,4)-cor_pts(2,3))^2);
        dists_sqrt(i) = sqrt(abs((cor_pts(1,3)-ref_pts(1,i))*(cor_pts(2,4)-cor_pts(2,3)) - (cor_pts(2,3)-ref_pts(2,i))*(cor_pts(1,4)-cor_pts(1,3)))/denomr_r);
    elseif i>n_sq_x+1+n_sq_y+n_sq_x
        denomr_u = sqrt((cor_pts(1,1)-cor_pts(1,4))^2 + (cor_pts(2,1)-cor_pts(2,4))^2);
        dists_sqrt(i) = sqrt(abs((cor_pts(1,4)-ref_pts(1,i))*(cor_pts(2,1)-cor_pts(2,4)) - (cor_pts(2,4)-ref_pts(2,i))*(cor_pts(1,1)-cor_pts(1,4)))/denomr_u);
    else
        
    end        
end
end

%% cost function two (consider all the corners L(r)=k1*r + k2*r^2 +k3*r^3 +k4*r^4 +k5*r^5 +k6*r^6)
function [dists_all_sqrt] = funcDist2(para_opt, x, xx)
    x = x+1; % move origin to (1 1)
    %num_pts = size(x,2);
    %dists_all_sqrt = zeros(num_pts,1);
    
    %transform points
    r2 = (x(1,:)-para_opt(5)).^2 + (x(2,:)-para_opt(6)).^2;
    r = sqrt(r2); %radial distance
    ref_pts_x = x(1,:) + (para_opt(1)*r+para_opt(2)*(r.^2)+para_opt(3)*(r.^3)+para_opt(4)*(r.^4)+para_opt(5)*(r.^5)+para_opt(6)*(r.^6)).*(x(1,:)-para_opt(5));
    ref_pts_y = x(2,:) + (para_opt(1)*r+para_opt(2)*(r.^2)+para_opt(3)*(r.^3)+para_opt(4)*(r.^4)+para_opt(5)*(r.^5)+para_opt(6)*(r.^6)).*(x(2,:)-para_opt(6));
    ref_pts = [ref_pts_x;ref_pts_y];
    
    dists_all_sqrt = sqrt(sum((ref_pts-xx).^2,1));
end

%% cost function two (consider all the corners L(r)=K1*r + k2*r^2 +k4*r^4 +k6*r^6)
function [dists_all_sqrt] = funcDist3(para_opt, x, xx)
    x = x+1; % move origin to (1 1)
    %num_pts = size(x,2);
    %dists_all_sqrt = zeros(num_pts,1);
    
    %transform points
    r2 = (x(1,:)-para_opt(5)).^2 + (x(2,:)-para_opt(6)).^2;
    r = sqrt(r2); %radial distance
    ref_pts_x = x(1,:) + (para_opt(1)*r+para_opt(2)*(r.^2)+para_opt(3)*(r.^4)+para_opt(4)*(r.^6)).*(x(1,:)-para_opt(5));
    ref_pts_y = x(2,:) + (para_opt(1)*r+para_opt(2)*(r.^2)+para_opt(3)*(r.^4)+para_opt(4)*(r.^6)).*(x(2,:)-para_opt(6));
    ref_pts = [ref_pts_x;ref_pts_y];
    
    dists_all_sqrt = sqrt(sum((ref_pts-xx).^2,1));
end