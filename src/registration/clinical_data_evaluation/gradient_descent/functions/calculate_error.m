function [diff, diff_st, diff1, diff2] = calculate_error(dim1, dim2, dim3, var_cell, vol1, vol2, M1, M2, show)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Calculate the intensity differences between the intersection of two given 
%% volumes
%%
%% Inputs:  1. dim1     -> (scalar) # of slices of the first volume
%%          2. dim2     -> (scalar) # of slices of the second volume
%%          3. dim3     -> (scalar) # of points (discretization) per intersection
%%          4. var_cell -> (cell) containing the 3D intersection points
%%          5. vol1     -> (matrix) N1 x M1 x L1 first volume 
%%          6. vol2     -> (matrix) N2 x M2 x L2 second volume 
%%          7. M1       -> (matrix) 3x4 transform 3D point to 2D pixel  (vol1)
%%          8. M2       -> (matrix) 3x4 transform 3D point to 2D pixel  (vol2)
%%
%% Outputs: 1. diff      -> (array) of size N, number of intersections, with the 
%%                     difference between intensity values in the intersections
%%          2. diff_st   -> (struct) statistic values of the 'diff' 
%%          4. diff1     -> (array) of size N, number of intersections, with the 
%%                          intensity values of the intersections in view 1
%%          5. diff2     -> (array) of size N, number of intersections, with the 
%%                          intensity values of the intersections in view 2
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check input 
if nargin < 9
    show = 0;
end

%% Define the array dimension
diff  = zeros(dim1 , dim2 , dim3);
diff1 = zeros(dim1 , dim2 , dim3);
diff2 = zeros(dim1 , dim2 , dim3);

% ind = 1;

for i = 1:dim1
    for j = 1:dim2
        for k = 1:dim3
            

            %% First plane %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [i1, j1, i2, j2, real_v]  = compute_coord(M1{i}, [var_cell{i,j,k} 1], size(vol1,1), size(vol1,2));
            
            neig = [vol1(i1, j1, i)   vol1(i1, j2, i);...
                    vol1(i2, j1, i)   vol1(i2, j2, i)];
            
            diff1(i,j,k) = bilinear_interpolation(real_v(2),real_v(1),double(neig));
            
            %% Second plane %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [i1, j1, i2, j2, real_v]  = compute_coord(M2{j}, [var_cell{i,j,k} 1], size(vol2,1), size(vol2,2));
            
            neig = [vol2(i1, j1, j)   vol2(i1, j2, j);...
                    vol2(i2, j1, j)   vol2(i2, j2, j)];

            diff2(i,j,k) = bilinear_interpolation(real_v(2),real_v(1),double(neig));
            
            %% Compute the difference between the intensity values
            diff(i,j,k) = abs(diff1(i,j,k) - diff2(i,j,k));
%             ind = ind + 1;
            
        end
    end
end

diff_st.mean = mean(diff(:));
diff_st.std  = std(diff(:));
diff_st.max  = max(diff(:));

%% Plot the results
if show
    
    figure;
    plot(0:dim3-1,squeeze(diff1(20,9,:)), 'r-');hold on
    plot(0:dim3-1,squeeze(diff2(20,9,:)), 'b-');hold on
    
end
