function [output] = read_dcm(folder,show)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Read and show the images from
%% a given folder with dicom files
%% NOTE: make sure that your folder only contains the images you
%%        you want to read and show.
%%
%% Inputs: 1. folder -> name of folder containing files
%%         3. show ->  0/1 = show/not show the images
%%
%% Output: 1. output -> struct matrix containing all the images
%%
%% Execute:
%% - images = read_images('resources/images/',1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dcm = dir(folder);

for i=3:size(dcm,1)
        
    filename = strcat(folder,dcm(i).name)
    
    im = dicomread(filename);
    info = dicominfo(filename);
    
    output.image(:,:,i-2)  = im;
    output.image_info{i-2} = info;
 

    if show == 1
        h = imshow(output.image(:,:,i-2),[]);title(['Image: ',dcm(i).name]);
        pause;
        close;
    end

end

if show == 1
    close all;
end