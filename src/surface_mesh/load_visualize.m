function load_visualize(file_name, views, show)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Main program for visualizing the slices and the manual segmentation
%% in RCS
%%
%% Inputs: 1. file_name -> '.mat' or '.ply' file containing the vertices and 
%%                         faces of the segmented organs
%%         2. views     -> 
%%         3. show      -> (optional) 0/1 not show/show the results
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check the inputs
if nargin < 3
    show = 0;
end

if iscell(file_name)
    
    for i=1:length(file_name)
        
        [vertices, faces] = read_ply(file_name{i});
        st(i).vertices = vertices;
        st(i).faces    = faces;
        
    end
    
    if show
        view3d(views, st);
    end
    
else
    %% Check the file extension %%
    if strcmp(file_name(end-2:end),'mat')
        
        st = load(file_name);
        if show
            view3d(views, st.struct_view);
        end
        
    end
    
    if strcmp(file_name(end-2:end),'ply')
        
        [vertices, faces] = read_ply(file_name);
        st.vertices = vertices;
        st.faces    = faces;
        if show
            view3d(views, st);
        end
        
    end
end




