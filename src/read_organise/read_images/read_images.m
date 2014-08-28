function [images] = read_images(folder,show)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Read and show the images from
%% a given folder.
%% NOTE: make sure that your folder only contains the images you
%%        you want to read and show.
%%
%% Inputs: 1. folder -> name of folder containing files
%%         3. show ->  0/1 = show/not show the images
%%
%% Output: 1. images -> matrix containing all the images
%%
%% Execute:
%% - images = read_images('resources/images/',1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dcm = dir(folder);

for i=3:size(dcm,1)
        
    filename = strcat(folder,dcm(i).name);
    im = imread(filename);

    images{i-2} = im;

    if show == 1
        h = imshow(images{i-2});title(['Image: ',dcm(i).name]);
        pause;
        close;
    end

end

if show == 1
    close all;
end