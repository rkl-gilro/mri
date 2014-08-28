function [faces, vertices] = surface_reconstruction( points, smooth )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Surface Reconstruction using Crust algorithm from a given set of points
%%
%% Inputs:  1. points   -> Nx3 matrix of points
%%          2. smooth   -> 0/1 no smoothing - smooth surface
%% Outputs: 1. faces    -> Nf x 3 matrix containing the faces of the surface
%%          2. vertices -> Nv x 3 matrix containing the vertices 
%%                        of the surface
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,I,~] = unique(points,'first','rows');
I = sort(I);
points = points(I,:);

%% Compute the surface reconstruction 
[t] = Crust(points);

FV.vertices = points;
FV.faces = double(t);

%% Smooth the resulting surface 
% FV2 = smoothpatch(FV,1,5);
% out = meannorm_trismooth(FV.vertices, FV.faces );
out2 = lpflow_trismooth( FV.vertices, FV.faces ); 

if smooth 
    
    faces = FV.faces;
    vertices = out2;
        
else
    
    faces    = FV.faces;
    vertices = FV.vertices;
    
end
