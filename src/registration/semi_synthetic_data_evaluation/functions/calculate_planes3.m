disp('--------- Calculate transformations ----------------------------------------------------')
global crop_volume

global lava_flex_ax
global lava_flex_sag
global lava_flex_cor

global lava_flex_info

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