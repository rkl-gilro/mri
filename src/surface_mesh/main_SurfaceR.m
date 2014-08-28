%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Main program for visualizing the slices and the manual segmentation
%% in RCS
%% 
%% The .mat files have the following struct: 'st.vertex{s,1}.contour{n,1}'
%% where 's' is the number of slices, and 'n' is the number of contours in 
%% the current slice
%% NOTE: This file contains all the different contours we wanted to segment
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

show = 1;

%% Name of the files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Name for the combined ply file (with vertices and faces) in case we want
%% to save it
combine_file = 'combine_u';

%% Axial file
load_file_ax = 'resources/manual_seg/bourasset/mat/myoma_ax.mat';
save_plyax = 'uterus_ax_';


%% Sagittal file
load_file_sag = 'resources/manual_seg/bourasset/mat/myoma_sag.mat';
save_plysag = 'uterus_sag_';

%% Coronal file
load_file_cor = 'resources/manual_seg/bourasset/mat/myoma_cor.mat';
save_plycor = 'uterus_cor_';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Registration %%
disp('-------------------- Registration -------------------------------------------');
% main_optimization;

%% Store the new data %%
% images.axial = new_axial;
% images.axial_info = axial_m;
% 
% images.sagittal = new_sagittal;
% images.sagittal_info = sag_m;
% 
% images.coronal = new_coronal;
% images.coronal_info = cor_m;

images.axial = views.axial;
images.axial_info = views.axial_info;

images.sagittal = views.sagittal;
images.sagittal_info = views.sagittal_info;

images.coronal = views.coronal;
images.coronal_info = views.coronal_info;

%% The manual segmentation 
disp('-------------------- Load Axial segmentation -------------------------------------------');

ax_st = load(load_file_ax); %% It is a struct: ax_st.vertex{s,1}.contour{n,1}

[faces_ax, vertices_ax] = compute_faces_vertex(ax_st, views.axial_info, 1, save_plyax, show);

disp('-------------------- Load Sagittal segmentation -------------------------------------------');

sag_st = load(load_file_sag); %% st.vertex{s,1}.contour{n,1}

[faces_sag, vertices_sag] = compute_faces_vertex(sag_st, views.sagittal_info, 2, save_plysag, show);

disp('-------------------- Load Coronal segmentation -------------------------------------------');

cor_st = load(load_file_cor);%% st.vertex{s,1}.contour{n,1}

[faces_cor, vertices_cor] = compute_faces_vertex(cor_st, views.coronal_info, 3, save_plycor, show);

%% Plot of the output triangulation
if show

    figure;
    trisurf(faces_cor,vertices_cor(:,1),vertices_cor(:,2),vertices_cor(:,3),'facecolor','c','edgecolor','b')%plot della superficie trattata
    hold on
    trisurf(faces_ax, vertices_ax(:,1), vertices_ax(:,2), vertices_ax(:,3), 'facecolor','r','edgecolor','b')%plot della superficie trattata
    hold on
    trisurf(faces_sag,vertices_sag(:,1),vertices_sag(:,2),vertices_sag(:,3),'facecolor','g','edgecolor','b')%plot della superficie trattata
    hold on
    
end

%% Combine the ply files for the final surface reconstruction
num_contours = length(ax_st.vertex{1}.contour); %% the number of contours is the same in the 3 views

for i = 1:num_contours
    
    ply_file_ax  = strcat(save_plyax, num2str(i), '.ply');
    ply_file_sag = strcat(save_plysag,num2str(i), '.ply');
    ply_file_cor = strcat(save_plycor,num2str(i), '.ply');
    
    combine_ply = strcat(combine_file,num2str(i), '.ply');
    
    [faces, vertices] = combine_plyfiles({ply_file_sag, ply_file_cor},1,combine_ply,1);
    st(i).vertices = vertices;
    st(i).faces    = faces;
    
end

%% Show the results of the reconstruction and the corresponding MRI
if show
    view3d(images, st);
end


