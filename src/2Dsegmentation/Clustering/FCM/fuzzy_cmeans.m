function [im_out mean_cl_ord std_cl_ord] = fuzzy_cmeans(im,c)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Compute Fuzzy C-means of a given image im
%%
%% Inputs: 1. im -> input image
%%         2. c -> number of clusters
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[rows cols ch] = size(im);

if ch > 1
    im = im(:,:,1);
end
%im = double(im);
%im = imadjust(im);

data = [im(:) im(:)]; % data array

% Options for FCM
opt = [NaN 100 NaN 0];

[center,U,obj_fcn] = fcm(data, c, opt); % Fuzzy C-means 
         
% Finding the pixels for each class
% Calculate the mean for its class in order
% to obtain 1 for dense tissue, and 0 for background
maxU = max(U);

mean_cl = [];
std_cl = [];
for i=1:c
    index{i} = find(U(i,:) == maxU);
    mean_cl = [mean_cl mean2(im(index{i}))];
    std_cl = [std_cl std2(im(index{i}))];
end
[mean_cl_ord ind]= sort(mean_cl);
std_cl_ord = std_cl(ind);

for i=1:c
    index_ord{i} = index{ind(i)};
end

% Assigning pixel to each class by giving them a specific value
fcmImage(1:length(data))=0; 
rate = 1/(c-1);
for i=1:c
    fcmImage(index_ord{i})= (i-1);%*rate;
end

% Reshape the array to a image
im_out = reshape(fcmImage,rows,cols);

%figure;
% subplot(121);imagesc(im);
% subplot(122);imagesc(im_out);
% subplot(121);imshow(im,[]);
%imshow(im_out,[]);