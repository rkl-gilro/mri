
global lava_flex_ax
global lava_flex_sag
global lava_flex_cor

global vol_ax_eval1
global vol_sag_eval1
global vol_cor_eval1

global deformation_type


R = makeresampler('cubic', 'replicate');
TDIMS_A = [1 2 3];
TDIMS_B = [1 2 3];

TMAP_B = [];
F = 0;


if deformation_type.rigid
    
    disp('--------- Apply deformations axial -----')
    %% Calculate the new volume for axial 
    % rotation matrix
    rot = compute_rotm(deformation_type.angles_ax);
    
    % translate to the origin first, then apply rigid transformation
%     T1 = [eye(3) -(size(lava_flex_ax)' + 1)/2;0 0 0 1];
    T1 = [1 0 0 0
    0 1 0 0
    0 0 1 0
    -(size(lava_flex_ax) + 1)/2 1];

    rig = [rot zeros(3,1); 0 0 0 1];
%     T2 = [eye(3) (size(lava_flex_ax)' + 1)/2;0 0 0 1];
    T2 = [1 0 0 0
    0 1 0 0
    0 0 1 0
    (size(lava_flex_ax) + 1)/2 1];

    Transl = [eye(3) deformation_type.translation_ax';0 0 0 1];
    
    T_transf = T1 * rig' * T2 * Transl';
    
    tform = maketform('affine', T_transf);
    TSIZE_B = size(lava_flex_ax);
    
    deformation_type.T_ax = T_transf;
    
%     vol_ax_eval1 = affine3d(lava_flex_ax, inv(T_transf), 1:size(lava_flex_ax,1), 1:size(lava_flex_ax,2), 1:size(lava_flex_ax,3), 'spline');
    vol_ax_eval1 = tformarray(lava_flex_ax, tform, R, TDIMS_A, TDIMS_B, TSIZE_B, TMAP_B, F);
%     vol_ax_eval1 = affine_transform(lava_flex_ax, T_transf, 0);
    disp('--------- Apply deformations sagittal -----')
    %% Calculate the new volume for sagittal 
    rot = compute_rotm(deformation_type.angles_sag);
    
%     T1 = [eye(3) -(size(lava_flex_sag)' + 1)/2;0 0 0 1];
    T1 = [1 0 0 0
    0 1 0 0
    0 0 1 0
    -(size(lava_flex_sag) + 1)/2 1];

    rig = [rot zeros(3,1); 0 0 0 1];
%     T2 = [eye(3) (size(lava_flex_ax)' + 1)/2;0 0 0 1];
    T2 = [1 0 0 0
    0 1 0 0
    0 0 1 0
    (size(lava_flex_sag) + 1)/2 1];

%     rig = [rot zeros(3,1); 0 0 0 1];
%     T2 = [eye(3) (size(lava_flex_sag)' + 1)/2;0 0 0 1];
    
    Transl = [eye(3) deformation_type.translation_sag';0 0 0 1];
    
    T_transf = T1 * rig' * T2 * Transl';
    tform = maketform('affine', T_transf);
    TSIZE_B = size(lava_flex_sag);
    
    deformation_type.T_sag = T_transf;
    
%     vol_sag_eval = affine3d(lava_flex_sag, inv(T_transf), 1:size(lava_flex_sag,1), 1:size(lava_flex_sag,2), 1:size(lava_flex_sag,3), 'spline');
    vol_sag_eval1 = tformarray(lava_flex_sag, tform, R, TDIMS_A, TDIMS_B, TSIZE_B, TMAP_B, F);
%     vol_sag_eval1 = affine_transform(lava_flex_sag, T_transf, 0);
    disp('--------- Apply deformations coronal -----')
    %% Calculate the new volume for coronal 
    rot = compute_rotm(deformation_type.angles_cor);
    
%     T1  = [eye(3) -(size(lava_flex_cor)' + 1)/2;0 0 0 1];
    T1 = [1 0 0 0
    0 1 0 0
    0 0 1 0
    -(size(lava_flex_cor) + 1)/2 1];

    rig = [rot zeros(3,1); 0 0 0 1];
%     T2 = [eye(3) (size(lava_flex_ax)' + 1)/2;0 0 0 1];
    T2 = [1 0 0 0
    0 1 0 0
    0 0 1 0
    (size(lava_flex_cor) + 1)/2 1];
%     rig = [rot zeros(3,1); 0 0 0 1];
%     T2  = [eye(3) (size(lava_flex_cor)' + 1)/2;0 0 0 1];
    
    Transl = [eye(3) deformation_type.translation_cor';0 0 0 1];
     
    T_transf = T1 * rig' * T2 * Transl';
    tform = maketform('affine', T_transf);
    TSIZE_B = size(lava_flex_cor);
    
    deformation_type.T_cor = T_transf;
    
%     vol_cor_eval1 = affine3d(lava_flex_cor, inv(T_transf), 1:size(lava_flex_cor,1), 1:size(lava_flex_cor,2), 1:size(lava_flex_cor,3), 'spline');
    vol_cor_eval1 = tformarray(lava_flex_cor, tform, R, TDIMS_A, TDIMS_B, TSIZE_B, TMAP_B, F);
%     vol_cor_eval1 = affine_transform(lava_flex_cor, T_transf, 0);
else
    %% Non rigid deformation based on b-splines
    
    percent = .2;
    
    [vol_ax_eval1,  grid_def] = random_deformation_bspline(  lava_flex_ax, deformation_type.current_mag_def, percent );
    deformation_type.T_ax = grid_def;
    
    [vol_sag_eval1, grid_def] = random_deformation_bspline( lava_flex_sag, 0, percent ); % deformation_type.current_mag_def
    deformation_type.T_sag = grid_def;
    
    [vol_cor_eval1, grid_def] = random_deformation_bspline( lava_flex_cor, 0, percent ); % deformation_type.current_mag_def
    deformation_type.T_cor = grid_def;
    
end

