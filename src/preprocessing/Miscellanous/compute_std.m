function new_im = compute_std(im,range)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Compute the standard deviation of a given image
%%
%% Inputs:  1. im -> input image
%%          2. range -> an array [r c] of range size
%% Outputs: 1. new_im -> output image
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2
  range = [3 3];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[rows cols] = size(im);
new_im = zeros(rows,cols);

im_r = zeros(rows + 2*range(1),cols + 2*range(2));
im_r(1+range(1):end - range(1),1+range(2):end - range(2)) = im;


for i=1+range(1):rows+range(1)
    for j=1+range(2):cols+range(2)

        i_n = i-range(1);
        j_n = j-range(2);
        new_im(i_n,j_n) = std2(im_r(i-range(1):i+range(1),j-range(2):j+range(2)));
        
    end
end

