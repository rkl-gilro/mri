%% Example for input .mat file with the previous struct .vertex{}
% file_ax  = 'resources/manual_seg/LASSALAS/mat/uterus/lassalas_uterus_ax.mat';
% file_sag = 'resources/manual_seg/LASSALAS/mat/uterus/lassalas_uterus_sag.mat';
% file_cor = 'resources/manual_seg/LASSALAS/mat/uterus/lassalas_uterus_cor.mat';
% 
% contour2mesh(file_cor,cor_m,  'lassalas_cor.ply');
% contour2mesh(file_ax, axial_m,'lassalas_ax.ply');
% contour2mesh(file_sag,sag_m,  'lassalas_sag.ply');
% 
% [f_ut v_ut] = combine_plyfiles('lassalas_ax.ply','lassalas_sag.ply','lassalas_cor.ply');
% 
% file_ax  = 'resources/manual_seg/LASSALAS/mat/myoma/axial1.mat';
% file_sag = 'resources/manual_seg/LASSALAS/mat/myoma/sag1.mat';
% file_cor = 'resources/manual_seg/LASSALAS/mat/myoma/cor1.mat';
% 
% contour2mesh(file_cor, cor_m,  'lassalas_myoma_cor.ply');
% contour2mesh(file_ax,  axial_m,'lassalas_myoma_ax.ply');
% contour2mesh(file_sag, sag_m,  'lassalas_myoma_sag.ply');
% 
% [f_my v_my] = combine_plyfiles('lassalas_myoma_ax.ply','lassalas_myoma_sag.ply','lassalas_myoma_cor.ply');

% filename_my = 'lassalas_myoma.ply';
% filename_ut = 'lassalas_uterus.ply';
% [v_ut, f_ut] = read_ply(filename_ut);
% [v_my, f_my] = read_ply(filename_my);

% cont(1).vertices = v_ut;
% cont(1).faces    = f_ut;
% 
% cont(2).vertices = v_my;
% cont(2).faces    = f_my;

%% Directly from the mat file with all the organs segmented
% st = load('resources/manual_seg/LASSALAS/mat/lassalas_uterus_myoma.mat');
% 
% view3d(views, st.struct_view);