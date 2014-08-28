function val = bilinear_interpolation(i,j,neig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Bilinear interpolation given (i,j) pixel, and the four values used for it
%%
%% Inputs:  1. i -> row coordinate of the pixel (real number)
%%          2. j -> col coordinate of the pixel (real number)
%%          3. neig -> 2 x 2 matrix containing the gray values of the four
%%                     neighborhood pixels 
%%                     [I(0,0) I(0,1);I(1,0) I(1,1)]
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i = i - floor(i);
j = j - floor(j);

val = [1-i i] * neig * [1-j j]';


