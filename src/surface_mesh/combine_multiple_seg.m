function [faces, vertices] = compute_faces_vertex(st, views_info, view, save_ply, show)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Number of contours
num_contours = length(st.vertex{1}.contour);

for i = 1:num_contours
    
    for j = 1:length(st.vertex)
        cont{j} = st.vertex{j}.contour{i};
    end
    ply_file = strcat(save_ply, num2str(i), '.ply');
    [faces, vertices] = contour2mesh(cont, views_info, view, ply_file, show );

end