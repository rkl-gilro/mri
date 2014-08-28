function [faces, vertices] = combine_plyfiles(ply_cell, smooth, combine_ply, show)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  Surface Reconstruction:
%%  Given the 3 ply files of the contours of the 3 different views 
%%  combine all the vertex and calculate the new mesh. They should correspond
%%  to the same organ or region.
%%
%%  NOTE: the coordinate points in the ply files (x,y,z) are in the 3D space
%%        i.e, in the RCS !
%%
%%  Inputs:  1. ply_cell -> (cell) ply files containing the vertex 
%%                        of the segmented organ in different views
%%           2. smooth  -> 0/1 :=  no smooth/smooth
%%           2. combine_ply(optional) -> (string) path of the ply file to save the 
%%                                   vertex and faces information
%%           3. show    -> if 1 shows the calculate surface, 0 (default) otherwise  
%% Outputs:  1. faces -> Nfx3 matrix, where each row is a face and the columns
%%                              correspond to the vertex of this face
%%           2. vertices -> Nv x 3 matrix with rows number of contour points, 
%%                    and columns (x,y,z) coordinates
%%
%% EXECUTE
%%
%% contour2mesh('resources/manual_seg/LASSALAS/lassalas_uterurs_cor2.mat',views.coronal_info,'lassalas_cor.ply');
%% contour2mesh('resources/manual_seg/LASSALAS/lassalas_uterurs_ax2.mat',views.axial_info,'lassalas_ax.ply');
%% contour2mesh('resources/manual_seg/LASSALAS/lassalas_uterurs_sag2.mat',views.sagittal_info,'lassalas_sag.ply');
%% combine_plyfiles({'lassalas_ax.ply','lassalas_sag.ply','lassalas_cor.ply'});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save_file = 1;

var = [];

%% Check the number of inputs 
if nargin == 4
    
    for i = 1:length(ply_cell)
        [va, ~] = read_ply(ply_cell{i});
        var = [var; va];
    
    end
end

if nargin == 3
    
    for i = 1:length(ply_cell)
        
        [va, ~] = read_ply(ply_cell{i});
        var = [var; va];
    
    end
    
    show = 0;
end

if nargin == 2
    save_file = 0;
    show = 0;
    
    for i = 1:length(ply_cell)
        
        [va, ~] = read_ply(ply_cell{i});
        var = [var; va];
    
    end

end

if nargin == 1
    save_file = 0;
    show = 0;
    smooth = 1;
    
    for i = 1:length(ply_cell)
        
        [va, ~] = read_ply(ply_cell{i});
        var = [var; va];
    
    end

end

%% Compute the surface from points 'var'
[faces, vertices] = surface_reconstruction( var, smooth );

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