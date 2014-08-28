%% Main program for watershed using entropy as the input
clc; close all;

%% Select the view (axial, sagittal or coronal)
D = views.axial;
%D = views.sagittal;
%D = views.coronal;
%% Preprocessing step
D = aniso2D(D);

%% Compute the entropy and the watershed algorithm
for i=1:size(D,3)
    seg_D(:,:,i) = entropyfilt(uint16(D(:,:,i)));
    water(:,:,i) = watershed(seg_D(:,:,i));
end

%% Show results
show_results(water);
