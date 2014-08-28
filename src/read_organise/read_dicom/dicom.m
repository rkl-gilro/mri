function [im, file] = dicom(filename, show)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Given a dicome file, it returns the image and its dicom info
%% NOTE: you can return only the tags that you need by uncommenting from 
%%       'line 20', and returning 'imval' instead of file
%%
%% Inputs:  1. filename -> string with the name of the file
%%          2. show -> 1 if show image, 0 otherwise
%% Outputs: 1. im -> image
%%          2. file -> struct containing all the dicom info
%%
%% Execute:
%% - [im,im_info] = dicom('resources/patients/patient1/exp0000',1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<2
    show = false;
end

file = dicominfo(filename);
im   = dicomread(filename);

%% This part can be changed to whatever tags you need, if you don't want 
%% to return all the dicom information
% imval.Width = file.Width;
% imval.Height = file.Height;
% imval.BitDepth = file.BitDepth;
% imval.PixelSpacing = file.PixelSpacing;
% imval.BitsAllocated = file.BitsAllocated;
% imval.BitsStored = file.BitsStored;
% imval.HighBit = file.HighBit;
% imval.PatientName = file.PatientName;

if show
    figure;
    imshow(im,[]);
end