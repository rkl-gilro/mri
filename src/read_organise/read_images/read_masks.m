function [images,alpha_ch] = read_masks(folder,show)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Read and show the images from
%% a given folder containing the segmentation masks
%% NOTE: make sure that the mask images in your folder have the form
%%      '% mask.png'. If not, you can change 'line 22'.
%%
%% Inputs: 1. folder -> name of folder containing files
%%         3. show ->  0/1 = show/not show the images
%%
%% Output: 1. images -> matrix containing all the images (labels)
%%         2. alpha_ch -> binary matrix (masks)
%%
%% Execute:
%% - [images,alpha_ch] = read_masks('resources/images/masks/',1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dcm = dir(folder);

for i=3:size(dcm,1)
    
    if strcmp(dcm(i).name(end-7:end),'mask.png')
        
        filename = strcat(folder,dcm(i).name);
        [im map alpha] = imread(filename);

        images{i-2} = im;
        alpha_ch{i-2} = alpha;
        
        if show == 1
            h = imshow(images{i-2});title(['Image: ',dcm(i).name]);
            set(h,'AlphaData',alpha_ch{i-2});
            pause;
            close;
        end
    end
end

if show == 1
    close all;
end