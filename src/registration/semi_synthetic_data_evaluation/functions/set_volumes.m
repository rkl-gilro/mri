global crop_volume
global crop_rectangle

global lava_flex

if crop_volume
    [lava_flex_n, crop_rectangle] = wrap_volume(lava_flex);
else
    lava_flex_n = lava_flex;
end
global size1 size2 size3

[size1, size2, size3] = size(lava_flex_n);

global lava_flex_ax
global lava_flex_sag
global lava_flex_cor

%% Rician noise

M     = 3;
alpha = 1;
show  = 0;
h   = 26;

disp('--------- Image denoising -----')
disp('-------------Rician Noise ---------')

% lava_flex = rician_noise(lava_flex_n, [M, alpha, h], show);
lava_flex = lava_flex_n;

disp('--------- Apply deformations axial -----')
    %% Calculate the new volume for axial 
%     R = makeresampler('cubic', 'replicate');
% TDIMS_A = [1 2 3];
% TDIMS_B = [1 2 3];
% 
% TMAP_B = [];
% F = 0;
%     % rotation matrix
%     rot = compute_rotm([0 0 .2]);
%     
%     % translate to the origin first, then apply rigid transformation
% %     T1 = [eye(3) -(size(lava_flex_ax)' + 1)/2;0 0 0 1];
%     T1 = [1 0 0 0
%     0 1 0 0
%     0 0 1 0
%     -(size(lava_flex) + 1)/2 1];
% 
%     rig = [rot zeros(3,1); 0 0 0 1];
% %     T2 = [eye(3) (size(lava_flex_ax)' + 1)/2;0 0 0 1];
%     T2 = [1 0 0 0
%     0 1 0 0
%     0 0 1 0
%     (size(lava_flex) + 1)/2 1];
% 
%     
%     T_transf = T1 * rig' * T2;
%     
%     tform = maketform('affine', T_transf);
%     TSIZE_B = size(lava_flex);
%     
% %     vol_ax_eval1 = affine3d(lava_flex_ax, inv(T_transf), 1:size(lava_flex_ax,1), 1:size(lava_flex_ax,2), 1:size(lava_flex_ax,3), 'spline');
%     lava_flex = tformarray(lava_flex, tform, R, TDIMS_A, TDIMS_B, TSIZE_B, TMAP_B, F);
% disp('-------------Anisotropic diffusion ---------')
% 
% for i=1:size3
%     lava_flex(:,:,i) = anisodiff2D(lava_flex_n(:,:,i), 20, 1/7, 30, 1);
% end

disp('--------- Define each direction -----')
lava_flex_ax  = lava_flex;

for k = 1:size3
    for i = 1:size1
        lava_flex_sag(k,i,:) = lava_flex(i,:,k); 
    end
end


for k = 1:size3
    for j = 1:size2
        lava_flex_cor(k,j,:) = lava_flex(:,j,k); 
    end
end

disp('--------- Calculate transformations ----------------------------------------------------')
global lava_axM
global lava_axM_1
global lava_sagM
global lava_sagM_1
global lava_corM
global lava_corM_1

global X_ax_v
global Y_ax_v
global Z_ax_v

global X_sag_v
global Y_sag_v
global Z_sag_v

global X_cor_v
global Y_cor_v
global Z_cor_v

global ortho
ortho = 1;


if crop_volume
    
    disp('--------- Calculate M & M^(-1) and plane eq. for each slice in the first direction -----')
    %% Axial
    [lava_axM, lava_axM_1, X_ax_v, Y_ax_v, Z_ax_v, ~]      = calculate_transforms(lava_flex_info, size(lava_flex_ax), [crop_rectangle(2)-1 crop_rectangle(4)-1], [crop_rectangle(1)-1 crop_rectangle(3)-1], ortho,  1);
    
    disp('--------- Calculate  M & M^(-1) and plane eq. for each slice in the second direction -----')
    %% Sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [lava_sagM, lava_sagM_1, X_sag_v, Y_sag_v, Z_sag_v, ~] = calculate_transforms_synt_wrap(lava_flex_info{1}, size(lava_flex_sag), [0 size(lava_flex_sag,1)-1], [crop_rectangle(2)-1 crop_rectangle(4)-1], [crop_rectangle(1)-1 crop_rectangle(3)-1], ortho, lava_axM{1}, 2);
    
    disp('--------- Calculate  M & M^(-1) and plane eq. for each slice in the third direction -----')
    %% Coronal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [lava_corM, lava_corM_1, X_cor_v, Y_cor_v, Z_cor_v, ~] = calculate_transforms_synt_wrap(lava_flex_info{1}, size(lava_flex_cor), [0 size(lava_flex_cor,1)-1], [crop_rectangle(1)-1 crop_rectangle(3)-1], [crop_rectangle(2)-1 crop_rectangle(4)-1], ortho, lava_axM{1}, 3);
    
else
    
    disp('--------- Calculate M & M^(-1) and plane eq. for each slice in the first direction -----')
    %% Axial    
    [lava_axM, lava_axM_1, X_ax_v, Y_ax_v, Z_ax_v, ~]      = calculate_transforms(lava_flex_info, size(lava_flex_ax), [0 size(lava_flex_ax,1)-1], [0 size(lava_flex_ax,2)-1], ortho,  1);
    
    disp('--------- Calculate  M & M^(-1) and plane eq. for each slice in the second direction -----')
    %% Sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [lava_sagM, lava_sagM_1, X_sag_v, Y_sag_v, Z_sag_v, ~] = calculate_transforms_synt(lava_flex_info{1}, size(lava_flex_sag), [0 size(lava_flex_sag,1)-1], [0 size(lava_flex_sag,2)-1], ortho, lava_axM{1}, 2);
    
    disp('--------- Calculate  M & M^(-1) and plane eq. for each slice in the third direction -----')
    %% Coronal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [lava_corM, lava_corM_1, X_cor_v, Y_cor_v, Z_cor_v, ~] = calculate_transforms_synt(lava_flex_info{1}, size(lava_flex_cor), [0 size(lava_flex_cor,1)-1], [0 size(lava_flex_cor,2)-1], ortho, lava_axM{1}, 3);
end
