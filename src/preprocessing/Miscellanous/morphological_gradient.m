function out = morphological_gradient(im,se)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  Compute morphological gradient for a given image and structuring 
%%  element
%%
%%  Inputs:  1. im -> input image
%%           2. se -> structuring element
%%  Outputs: 1. out -> output image 
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out = imdilate(im, se) - imerode(im, se);