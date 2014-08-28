%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Main co-occurrence matrix
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
%% Select input image
%im = views.axial(:,:,10);
im = views.sagittal(:,:,12);
%im = views.coronal(:,:,12);

%% Preprocessing step
im = anisodiff2D(im, 10,1/7,50,2);

%im = convert2u8(im);
%im = im2double(im);
%% Set co-occurrence matrix parameters
num_of_levels = 8;
offset = 1;

% Directions
% 0     [  0  D ]
% 45    [ -D  D ]
% 90    [ -D  0 ]
% 135   [ -D -D ]

angle = offset.*[0 1];
neigh_range = [55 55];

tic;

gl = GLCM(im, num_of_levels, angle, neigh_range);

clear angle im neigh_range num_of_levels offset

toc
%% Show results
figure;
subplot(221);imshow(gl(:,:,1),[]);title('Contrast');
subplot(222);imshow(gl(:,:,2),[]);title('Correlation');
subplot(223);imshow(gl(:,:,3),[]);title('Energy');
subplot(224);imshow(gl(:,:,4),[]);title('Homogeneity');