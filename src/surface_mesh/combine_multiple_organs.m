function st = combine_multiple_organs( ply_cell )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Combines different '.ply' files for different organs in a single struct
%% NOTE: The '.ply' files contain the vertices and faces of an organ
%% 
%% Inputs:  1. ply_cell -> cell containing the 'N' ply files to combine
%%
%% Outputs: 1. st  -> struct with vertices & faces, same length as the number
%%                    of ply files 'N'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i = 1:length(ply_cell)
    
    [faces, vertices] = read_ply(ply_cell{i});
    st(i).vertices = vertices;
    st(i).faces    = faces;
    
end





