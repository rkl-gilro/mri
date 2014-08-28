%% This is a demo showing how to use this toolbox

%   Copyright by Quan Wang, 2012/04/25
%   Please cite: Quan Wang. HMRF-EM-image: Implementation of the 
%   Hidden Markov Random Field Model and its Expectation-Maximization 
%   Algorithm. arXiv:1207.3510 [cs.CV], 2012.

%clear;clc;close all;



%I=imread('Beijing World Park 8.JPG');
%Y=rgb2gray(I);
%Y = patient_or.sagittal(:,:,10);
%Y = patient_or.axial(:,:,10);
%Y = views.coronal(:,:,12);

Y=double(ax);
%Y=gaussianBlur(Y,3);
Y = anisodiff2D(Y, 10, 1/7, 30, 2);
Z = edge(uint16(Y),'canny');
%Z = entropyfilt(uint16(Y));

k=4;
EM_iter=7; % max num of iterations
MAP_iter=10; % max num of iterations

tic
fprintf('Performing Fuzzy C-means segmentation\n');
%[X mu sigma]=image_kmeans(Y,k);
[X mu sigma] = fuzzy_cmeans(Y,k);
%imwrite(uint8(X*120),'initial labels.png');

[out mu sigma]=HMRF_EM(X,Y,Z,mu,sigma,k,EM_iter,MAP_iter);
%imwrite(uint8(X*120),'final labels.png');

figure;
subplot(221);imshow(Y,[]);
subplot(222);imshow(Z,[]);
subplot(223);imshow(X,[]);
subplot(224);imshow(out,[]);

toc