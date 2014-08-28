function [diff, statistics] = compare_grid(grid1, grid2, show)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  Calculate the distances between two given grids, the source vs the target,
%%  already computed for the registration problem.
%% 
%%  Inputs:  1.
%%           2.
%%           3. show
%%  Outputs: 1. diff  -> array containing the distances between points
%%           2. statistics -> struct with statistics, e.g mean std max ...
%%                            of the array 'diff'
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check the input arguments
if size(grid1.X,1) ~= size(grid2.X,1)
    disp('Error the size of the grids is not the same \n');
    return;
else
    num_points = size(grid1.X,1);
end

if nargin < 3
    show = 0;
end

diff = zeros(1, num_points);

%% Compute the distances (l2) between the 2 grids
for i=1:num_points
    
    diff(i) = sqrt( sum((grid1.X(i,:) - grid2.X(i,:)).^2) ); %% Euclidean distance
    
end

%% Compute some statistics, add more if needed
statistics.mean = mean(diff);
statistics.min  = min(diff);
statistics.max  = max(diff);
statistics.std  = std(diff);

%% Plot the results 
if show
    
    figure;
    tetramesh(grid1, 'EdgeColor', 'r', 'FaceColor', 'none'); hold on
    tetramesh(grid2, 'EdgeColor', 'b', 'FaceColor', 'none'); hold on
    plot3(grid1.X(:,1), grid1.X(:,2), grid1.X(:,3), 'r+');hold on %% first grid
    plot3(grid2.X(:,1), grid2.X(:,2), grid2.X(:,3), 'b*');hold on %% second grid
    
    for i=1:num_points
        plot3([grid1.X(i,1);grid2.X(i,1)], [grid1.X(i,2);grid2.X(i,2)], [grid1.X(i,3);grid2.X(i,3)],'g-');hold on
    end
    
end

