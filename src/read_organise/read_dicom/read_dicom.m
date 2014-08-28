function [images, dicom_inf, series_des, name_series] = read_dicom(folder,resize_factor,show)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Read and show the DICOM images from
%% a given folder.
%% Inputs: 1. folder -> name of folder containing dcm files
%%         2. resize_factor -> resize images from dcm files
%%         3. show ->  0/1 = show/not show the images
%%
%% Output: 1. images -> cell of images obtain from dcm files and resize
%%                      given the resize_factor
%%         2. dicom_inf -> cell of dicom information
%%         3. series_des ->
%%         4. name_series ->
%%
%% Execute:
%% - [im, im_info] = read_dicom('resources/patients/patient1/',.5,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<2
    resize_factor = 1;
    show = false;
end
if nargin<3
    show = false;
end

dcm = dir(folder);
num_series = 1;
series_des = {};
name_series = {};

h = waitbar(0,'Please wait...');

for i=3:size(dcm,1)

    [im, info] = dicom(strcat(folder,'/', dcm(i).name), 0);
    im = imresize(im, resize_factor,'bicubic');
    images{i-2}    = im;
    dicom_inf{i-2} = info;

    %% Let us compare the series description with the rest of the already read
    Fx = @(x) strcmp(x, info.SeriesDescription);
    inds = cellfun(Fx, name_series);
    
    %% If it is a new series description
    if isempty(find(inds, 1))
        
        name_series{num_series}          = info.SeriesDescription;
        series_des{num_series}.images{1} = im;  
        series_des{num_series}.info{1}   = info; 
        
        num_series = num_series + 1;
    else
        series_des{inds}.images{end + 1} = im;  
        series_des{inds}.info{end + 1}   = info; 
    end
    
    if show == 1
        imshow(im,[]);title(['Image: ',dcm(i).name]);
        pause;
    end
    
    waitbar((i-2) / (size(dcm,1)-2))
    
end

close(h)

if show == 1
    close all;
end