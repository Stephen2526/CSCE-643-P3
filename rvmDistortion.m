%
%remove radial distortion from images
%
%  load images
imgs = loadImg();

% remove distortion using optimized parameters
f = 1;
p1 = 0;
p2 = 0;
ppx = para_opt(5);
ppy = para_opt(6);
k1 = para_opt(1);
k2 = para_opt(2);
k3 = para_opt(3);
k4 = para_opt(4);
undis_img = cell(size(imgs,2));
for i=1 : size(imgs,2)
    undis_img{i} = undistortimage(imgs{i}, f, ppx, ppy, k1, k2, k3, k4, p1, p2);
    dir_tmp = strcat('imgs/undisImg',int2str(i),'_new.jpg');
    imwrite(undis_img{i},dir_tmp);
end
