function statistics = Laws_filter(im,mask)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Compute the mean and standard deviation after applying a Laws mask
%%
%% Inputs:  1. im -> input image
%%          2. mask -> Laws filter mask
%% Outputs: 1. statistics -> a struct containing the mean, absmean and
%%                           stdev
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
  disp('An image and a mask are needed')
  return 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[rg_rows rg_cols] = size(mask);
[rows cols] = size(im);

Law_im = imfilter(im,mask);

for i=1+rg_rows:rows-rg_rows
    for j=1+rg_cols:cols-rg_cols

        statistics.mean(i,j) = mean2(Law_im(i-rg_rows:i+rg_rows,j-rg_cols:j+rg_cols));
        statistics.absmean(i,j) = abs(statistics.mean(i,j));
        statistics.stdev(i,j) = std2(Law_im(i-rg_rows:i+rg_rows,j-rg_cols:j+rg_cols));
        
    end
end