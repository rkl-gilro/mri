function GLCM_im = GLCM(im, num_of_levels, offset, range)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Compute the co-occurrence matrix for the given image
%%
%% Inputs:  1. im -> input image
%%          2. num_of_levels -> number of gray levels
%%          3. offset -> an array for the angle and offset (e.g [0 1] )
%%          4. range -> array for the neighborhood range (e.g [3 3] )
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2
    num_of_levels = 8;
    offset = [0 1];
    range = [3 3];
end
if nargin < 3
    offset = [0 1];
    range = [3 3];
end
if nargin < 4
    range = [3 3];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[rows cols] = size(im);
GLCM_im = zeros(rows,cols,4);

range = [floor((range(1)-1)/2) floor((range(2)-1)/2)]
im_r = zeros(rows + 2*range(1),cols + 2*range(2));
im_r(1+range(1):end - range(1),1+range(2):end - range(2)) = im;

[rows_n cols_n] = size(im_r);

for i=1+range(1):rows_n-range(1)
    
    for j=1+range(2):cols_n - range(2)
        
        i_n = i-range(1);
        j_n = j-range(2);
        
        neighbour = im_r(i-range(1):i+range(1),j-range(2):j+range(2));
        
        % Compute the co-occurrence matrix for the given neighborhood
        GLCM_local = graycomatrix(neighbour,'NumLevels',num_of_levels,'Offset',offset,'GrayLimits' ,[min(min(im)) max(max(im))]);
        stats = graycoprops(GLCM_local);
        
        % Contrast
        GLCM_im(i_n,j_n,1) = mean(stats.Contrast);
        
        % Correlation
        GLCM_im(i_n,j_n,2) = mean(stats.Correlation);
        
        % Energy
        GLCM_im(i_n,j_n,3) = mean(stats.Energy);
        
        % Homo
        GLCM_im(i_n,j_n,4) = stats.Homogeneity;
        
    end
    
end