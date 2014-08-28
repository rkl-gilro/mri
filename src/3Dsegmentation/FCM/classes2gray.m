function output = classes2gray(class)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% From the output of BCFCM3D, 3D data with 3 different channels
%% for each class, get the 3D data in grayscale.
%%
%% Inputs:  1. class -> segmented 3D data (dbt images)
%%                      using fcm3D (3 classes)
%%
%% Outputs: 1. output -> 3D segmented data (grayscale)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:size(class,3)
    output(:,:,i) = class(:,:,i,1).*0 + class(:,:,i,2).*.5 +class(:,:,i,3).*1; 
end
