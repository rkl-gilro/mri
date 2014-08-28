function show_error(diff1, diff2, slice1, slice2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  Displays the error in the intensity values in the intersection of 
%%  two slices from different views
%%
%%  Inputs:  1. diff1  -> (array) 
%%           2. diff2  ->
%%           3. slice1 ->
%%           4. slice2 ->
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~, ~, dim3] = size(diff1);

figure;
plot(0:dim3-1,squeeze(diff1(slice1, slice2, :)), 'r-');hold on
plot(0:dim3-1,squeeze(diff2(slice1, slice2, :)), 'b-');hold on
    
