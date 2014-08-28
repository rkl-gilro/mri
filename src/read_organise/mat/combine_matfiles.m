function [faces vertices] = combine_matfiles(option, mat_cell, info_view, smooth, combine_ply, show)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  Surface Reconstruction:
%%  Given the 3 mat files of the contours of the 3 different views 
%%  combine all the vertex and calculate the new mesh. They should correspond
%%  to the same organ or region.
%%
%%  NOTE: the coordinate points in the mat files are 2D space (i,j) coordinates
%%
%%  Inputs:  1. mat_cell -> (cell) mat files containing the vertex 
%%                        of the segmented organ in different views
%%           2. infor_view -> (cell) view info for computing the transform M
%%                           (i,j) --> (x,y,z)
%%           3. smooth  -> 0/1 :=  no smooth/smooth
%%           4. combine_ply(optional) -> (string) path of the ply file to save the 
%%                                   vertex and faces information
%%           5. show    -> if 1 shows the calculate surface, 0 (default) otherwise  
%% Outputs:  1. faces -> Nfx3 matrix, where each row is a face and the columns
%%                              correspond to the vertex of this face
%%           2. vertices -> Nv x 3 matrix with rows number of contour points, 
%%                    and columns (x,y,z) coordinates
%%
%% EXECUTE
%%
%% files = {'resources/manual_seg/LASSALAS/lassalas_uterurs_cor2.mat','resources/manual_seg/LASSALAS/lassalas_uterurs_ax2.mat',...
%%          'resources/manual_seg/LASSALAS/lassalas_uterurs_sag2.mat'};
%% info  = {views.coronal_info, views.axial_info, views.sagittal_info};
%% combine_matfiles(files, info);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save_file = 1;

%% Check input arguments
if length(mat_cell) ~= length(info_view)
    
    disp('Error: number of files is not equal to number of info views')
    
end

N = length(mat_cell);

%% Compute the (x,y,z) contour points in RCS
vert = [];


for i=1:N
    
    a = load(mat_cell{i});
    
    %% struct of a is a.vertex{j}(N,1:3)
    for j=1:length(a.vertex)
        if ~isempty(a.vertex{j})
            %% Compute the M transformations for each info_view
            [M{i}{j}, ~, ~] = compute_M_M1(info_view{i}{j}, 0);
            %% Get the points in RCS
            p = M{i}{j} * [a.vertex{j}(:,1) a.vertex{j}(:,2) ones(size(a.vertex{j},1),1)]';
            vert = [vert; p(1:3)];
        end
    end
    
end

%% Compute the surface from points 'vert'
[faces, vertices] = surface_reconstruction( vert, smooth );

%% Display results
%% Plot of the output triangulation
if show

    figure;trimesh(faces, vertices(:,1),vertices(:,2),vertices(:,3),'facecolor','c','edgecolor','b','edgealpha',0,'facelighting','flat');camlight
    title('Flow Smooth Triangulation','fontsize',14)
    axis equal
    hold on

end

%% Save the results in combine_ply 
if save_file
    save_ply(vertices, faces, combine_ply);
end
