function [vol_wrap, position] = wrap_volume(vol, show)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  Crop a volume
%%
%% Inputs:  1. vol  -> mxnxs volume
%%          2. show ->  (optional) 0/1 show or not the resulting volume
%% Outputs: 1. vol_wrap -> the resulting cropped volume
%%          2. position -> the (x,y) coordinates of the left up corner &
%%                             (x,y) coordinates of the right bottom corner 
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check the inputs
if nargin < 2
    show = 0;
end

hf = figure;
imshow(vol(:,:,round(size(vol,3)/2)),[]);title('Double click to close')

h = imrect;
position = wait(h);

close all;

%% Calculate the coordinates of the rectangle
position(3) = position(3) + position(1);
position(4) = position(4) + position(2);

%% Check the boundaries
position(2) = max(1, position(2));
position(1) = max(1, position(1));

position(4) = min(size(vol,2), position(4));
position(3) = min(size(vol,2), position(3));


%% Crop the volume
vol_wrap = vol(position(2):position(4), position(1):position(3),:);

if show
    show_results(vol_wrap)
end