function measurements( ) 

global measure_errors

global vol_ax
global vol_sag
global vol_cor

global new_axial
global new_sagittal
global new_coronal

global axial_m1
global sag_m1
global cor_m1

global var_cell1
global var_cell2
global var_cell3

global optimizer

show = 0;

%% First intersection
[diff1_a, diff_st1_a, diff_ax1_a, diff_sag1_a] = calculate_error(size(vol_ax,3), size(vol_sag,3), optimizer.t, var_cell1, new_axial, new_sagittal, axial_m1, sag_m1, show);
[diff1_b, diff_st1_b, diff_ax1_b, diff_sag1_b] = calculate_error(size(vol_ax,3), size(vol_sag,3), optimizer.t, var_cell1,    vol_ax,      vol_sag, axial_m1, sag_m1, show);

%% Second intersection
[diff2_a, diff_st2_a, diff_ax2_a, diff_cor2_a] = calculate_error(size(vol_ax,3), size(vol_cor,3), optimizer.t, var_cell2, new_axial, new_coronal, axial_m1, cor_m1, show);
[diff2_b, diff_st2_b, diff_ax2_b, diff_cor2_b] = calculate_error(size(vol_ax,3), size(vol_cor,3), optimizer.t, var_cell2,    vol_ax,     vol_cor, axial_m1, cor_m1, show);

%% Third intersection
[diff3_a, diff_st3_a, diff_cor3_a, diff_sag3_a] = calculate_error(size(vol_cor,3), size(vol_sag,3), optimizer.t, var_cell3, new_coronal, new_sagittal, cor_m1, sag_m1, show);
[diff3_b, diff_st3_b, diff_cor3_b, diff_sag3_b] = calculate_error(size(vol_cor,3), size(vol_sag,3), optimizer.t, var_cell3,     vol_cor,      vol_sag, cor_m1, sag_m1, show);


measure_errors.diff1_a     = diff1_a;
measure_errors.diff1_b     = diff1_b;
measure_errors.diff_st1_a  = diff_st1_a;
measure_errors.diff_st1_b  = diff_st1_b;
measure_errors.diff_ax1_a  = diff_ax1_a;
measure_errors.diff_ax1_b  = diff_ax1_b;
measure_errors.diff_sag1_a = diff_sag1_a;
measure_errors.diff_sag1_b = diff_sag1_b;

measure_errors.diff2_a     = diff2_a;
measure_errors.diff2_b     = diff2_b;
measure_errors.diff_st2_a  = diff_st2_a;
measure_errors.diff_st2_b  = diff_st2_b;
measure_errors.diff_ax2_a  = diff_ax2_a;
measure_errors.diff_ax2_b  = diff_ax2_b;
measure_errors.diff_cor2_a = diff_cor2_a;
measure_errors.diff_cor2_b = diff_cor2_b;

measure_errors.diff3_a     = diff3_a;
measure_errors.diff3_b     = diff3_b;
measure_errors.diff_st3_a  = diff_st3_a;
measure_errors.diff_st3_b  = diff_st3_b;
measure_errors.diff_cor3_a = diff_cor3_a;
measure_errors.diff_cor3_b = diff_cor3_b;
measure_errors.diff_sag3_a = diff_sag3_a;
measure_errors.diff_sag3_b = diff_sag3_b;

measure_errors.total_diff_a = [diff1_a(:); diff3_a(:); diff3_a(:)];
measure_errors.total_diff_b = [diff1_b(:); diff3_b(:); diff3_b(:)];