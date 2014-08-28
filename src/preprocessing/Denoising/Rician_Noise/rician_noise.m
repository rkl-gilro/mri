function [out, backgr, h] = rician_noise(rima, parameters, show )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Denoising Filter for Rician noise based on 
%% Optimized Blockwise Non Local Means Denoising Filter ( ONLM ) (ACM folder)
%% 
%% Inputs:  1. rima       -> mxnxs volume data
%%          2. parameters -> (optional) array of | length 2 ( [M alpha] ) or 
%%                                               | length 3 ([M alpha h])
%%          3. show       -> (optional) 0/1 - not show/show the results
%% 
%% Outputs: 1. out        -> the denoised volume
%%          2. backgr     -> in case parameter h is estimated, the selected
%%                           background for its computation
%%          3. h          -> the estimated (or not) the parameter h
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check the input arguments
if nargin < 3
    show  = 0;
end
if nargin < 2
    alpha = 1;
    M     = 3;
    parameters = [M, alpha];
end

rima = double(rima);
s = size(rima);
backgr = [];

%% The input parameters is an array of 2 elements [M alpha]
if length(parameters) == 2
    
    M = parameters(1);
    alpha = parameters(2);
    
    %% Estimate 'h' value by selecting a region in the middle slice
    figure; imshow(rima(:,:,round(s(3)/2)),[]);
    handle = impoly;
    position = getPosition(handle);
    
    r_min = min(position(:,2));
    r_max = max(position(:,2));
    c_min = min(position(:,1));
    c_max = max(position(:,1));
    
    backgr = rima(r_min:r_max,c_min:c_max,:);
    
    close;
    %% Calculate the standard deviation
    
    %% Calculate the second order moment
    fun2 = @(x) sum( x(:).^2 ) / (length(x(:)) - 1);
    second_mean = zeros(size(backgr));
    
    for i=1:size(backgr,3)
        
        second_mean(:,:,i) = nlfilter(backgr(:,:,i),[7 7],fun2);
        
    end
    %% Calculate the mode of the distribution of second-order moment
    second_mean_array = second_mean(:);
    indices    =  find(diff([second_mean(:); realmax]) > 0); % indices where repeated values change
    [~ ,i]     =  max (diff([0; indices]));        % longest persistence length of repeated values
    mode       =  second_mean_array(indices(i));
    h = ceil(sqrt(mode / 2)) * 2
    
    %% Optimized Blockwise Non Local Means Denoising Filter
    out = ornlm(rima, M+2, alpha, h);
    
end

%% The input parameters is an array of 3 elements [M alpha h]
if length(parameters) == 3
    
    M = parameters(1);
    alpha = parameters(2);
    h = parameters(3);
    
    %%  Optimized Blockwise Non Local Means Denoising Filter
    out = ornlm(rima, M+2, alpha, h);
    
end

%% Display the results if show = 1
if show
    figure;
    subplot(121);imshow(rima(:,:,round(s(3)/2)),[]);title('Original Image');
    subplot(122);imshow(out(:,:,round(s(3)/2)),[]);title('Denoised Image');
end


