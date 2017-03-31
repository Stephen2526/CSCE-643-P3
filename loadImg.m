function [imgs] = loadImg()
%loadImg.m
%
%load all images in a folder to the memory
%User enters the folder dir, image basename and image format
%


%% input dir, basename and format
dir_name = input('The folder where images exist (relative dir): ', 's');
base_img = input('The basename (without number) of the image: ', 's');
format_img = input('The format of the image (without point): ', 's');

%% open the folder and read images
curr_dir = pwd;
full_dirName = strcat(curr_dir,'\',dir_name,'\', base_img,'*.',format_img);
img_files = dir(full_dirName);
imgs = cell(1,length(img_files));
for i = 1 : length(img_files)
    file_name = strcat(curr_dir,'\',dir_name,'\',img_files(i).name);
    fprintf(1, 'loading image %s \n', file_name);
    imgs{i} = imread(file_name);
end
end
