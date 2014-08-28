function neig = block(pos,size_neig,im)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% It returns a neighborhood of a pixel, 
%% given a specific size in an image
%%
%% Inputs:  1. pos -> 1x2 vector (i,j) position of the pixel
%%          2. size_neig -> 1x2 vector size of the neighborhood, i.e. [3 3], it should be odd
%%          3. im -> mxn image 
%% Outputs: 1. neig -> size_neig matrix 
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

im_pad = zeros(size(im) + (size_neig-1));

pad_r = (size_neig(1)-1)/2;
pad_c = (size_neig(2)-1)/2;

im_pad(1+pad_r:end-pad_r,1+pad_c:end-pad_c) = im;

i = pos(1) + pad_r;
j = pos(2) + pad_c;

neig = im_pad(i-pad_r:i+pad_r,j-pad_c:j+pad_c);