function [ minP ] = MLE( P0, pts_wld, pts_img )
%MLE Summary of this function goes here
%   MLE algorithm
%% Initial value
P0_tr = P0';
p0 = P0_tr(:)

%% Geometric minimization of Sampson error

% mininization process
options = optimoptions('lsqnonlin','FiniteDifferenceType','central','Display','iter-detailed','algorithm','Levenberg-Marquardt');
%options.StepTolerance = 1.000000e-08;
minp = lsqnonlin(@funGeome,p0,[],[],options, pts_wld, pts_img);
minP = reshape(minp, 4, 3)';
end

% function [ err ] = funSamp(p, ptsA, ptsB)
% % This is the user-defined function to express the objective function (Sampson error)
% %Arguemnt
% % denote elements in h
% h11 = p(1);h12 = p(2);h13 = p(3);
% h21 = p(4);h22 = p(5);h23 = p(6);
% h31 = p(7);h32 = p(8);h33 = p(9);
% 
% err = zeros(size(ptsA,2),1);
% for i = 1:size(ptsA,2)
%    eps = [-ptsB(3,i)*(ptsA(1,i)*h21 + ptsA(2,i)*h22 + h23) + ptsB(2,i)*(ptsA(1,i)*h31 + ptsA(2,i)*h32 + h33);
%           ptsB(3,i)*(ptsA(1,i)*h11 + ptsA(2,i)*h12 + h13) - ptsB(1,i)*(ptsA(1,i)*h31 + ptsA(2,i)*h32 + h33)];
%    J = [-ptsB(3,i)*h21+ptsB(2,i)*h31 , -ptsB(3,i)*h22+ptsB(2,i)*h32 , 0 , ptsA(1,i)*h31+ptsA(2,i)*h32+h33;
%         ptsB(3,i)*h11-ptsB(1,i)*h31 , -ptsB(3,i)*h12-ptsB(1,i)*h32 , -ptsA(1,i)*h31-ptsA(2,i)*h32-h33 , 0];
%    err(i) = sqrt(eps'*inv(J*J')*eps);
% end
% % display('err:');
% % err'*err
% end

function [err] = funGeome(minp, pts_wld, pts_img)
% cost function for minimize geometric error
P_tmp = reshape(minp, 4, 3)';
pts_trs = P_tmp*pts_wld;
%make third coor 1
coor_thd = repmat(pts_trs(3,:),3,1);
pts_trs = pts_trs./coor_thd;
dis_diff = pts_img - pts_trs;
err = sqrt(dis_diff(1,:).^2 + dis_diff(2,:).^2);
end