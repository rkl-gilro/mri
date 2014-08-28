function new_im = convert2u8(im)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Convert image to uint8
%%
%% Inputs:  1. im -> input image
%% Outputs: 1. new_im -> output image
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
im = double(im);
min_v_n = 0;
max_v_n = 255;

min_v = min(min(im));
max_v = max(max(im));

new_im = min_v_n + (im - min_v).*((max_v_n - min_v_n)/(max_v - min_v));

new_im = uint8(new_im);