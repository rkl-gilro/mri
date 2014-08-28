%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Main program to compute Otsu thresholding method
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; close all;

%% Select the view (axial, sagittal or coronal)
D = views.axial;
%D = views.sagittal;
%D = views.coronal;
%% Preprocessing step
D = aniso2D(D);

%% Set the number of labels
n = 5;

%% Compute thresholding using Otsu method
%% Otsu returns a labeled image between 1-n
for i=1:size(D,3)
    seg_D(:,:,i) = otsu(D(:,:,i),n);
end

show_results(seg_D);