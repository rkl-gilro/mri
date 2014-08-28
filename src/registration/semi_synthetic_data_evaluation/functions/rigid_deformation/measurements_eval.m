function measurements_eval( ) 

global measure_errors_eval

global vol_ax_eval
global vol_sag_eval
global vol_cor_eval

global new_axial
global new_sagittal
global new_coronal

global axial_M1
global sag_M1
global cor_M1

global var_cell1_v
global var_cell2_v
global var_cell3_v

global optimizer

show = 0;

%% First intersection
[diff1_a, diff_st1_a, diff_ax1_a, diff_sag1_a] = calculate_error(size(vol_ax_eval,3), size(vol_sag_eval,3), optimizer.t, var_cell1_v,   new_axial, new_sagittal, axial_M1, sag_M1, show);
[diff1_b, diff_st1_b, diff_ax1_b, diff_sag1_b] = calculate_error(size(vol_ax_eval,3), size(vol_sag_eval,3), optimizer.t, var_cell1_v, vol_ax_eval, vol_sag_eval, axial_M1, sag_M1, show);

%% Second intersection
[diff2_a, diff_st2_a, diff_ax2_a, diff_cor2_a] = calculate_error(size(vol_ax_eval,3), size(vol_cor_eval,3), optimizer.t, var_cell2_v,   new_axial,  new_coronal, axial_M1, cor_M1, show);
[diff2_b, diff_st2_b, diff_ax2_b, diff_cor2_b] = calculate_error(size(vol_ax_eval,3), size(vol_cor_eval,3), optimizer.t, var_cell2_v, vol_ax_eval, vol_cor_eval, axial_M1, cor_M1, show);

%% Third intersection
[diff3_a, diff_st3_a, diff_cor3_a, diff_sag3_a] = calculate_error(size(vol_cor_eval,3), size(vol_sag_eval,3), optimizer.t, var_cell3_v,  new_coronal, new_sagittal, cor_M1, sag_M1, show);
[diff3_b, diff_st3_b, diff_cor3_b, diff_sag3_b] = calculate_error(size(vol_cor_eval,3), size(vol_sag_eval,3), optimizer.t, var_cell3_v, vol_cor_eval, vol_sag_eval, cor_M1, sag_M1, show);


measure_errors_eval.diff1_a     = diff1_a;
measure_errors_eval.diff1_b     = diff1_b;
measure_errors_eval.diff_st1_a  = diff_st1_a;
measure_errors_eval.diff_st1_b  = diff_st1_b;
measure_errors_eval.diff_ax1_a  = diff_ax1_a;
measure_errors_eval.diff_ax1_b  = diff_ax1_b;
measure_errors_eval.diff_sag1_a = diff_sag1_a;
measure_errors_eval.diff_sag1_b = diff_sag1_b;

measure_errors_eval.diff2_a     = diff2_a;
measure_errors_eval.diff2_b     = diff2_b;
measure_errors_eval.diff_st2_a  = diff_st2_a;
measure_errors_eval.diff_st2_b  = diff_st2_b;
measure_errors_eval.diff_ax2_a  = diff_ax2_a;
measure_errors_eval.diff_ax2_b  = diff_ax2_b;
measure_errors_eval.diff_cor2_a = diff_cor2_a;
measure_errors_eval.diff_cor2_b = diff_cor2_b;

measure_errors_eval.diff3_a     = diff3_a;
measure_errors_eval.diff3_b     = diff3_b;
measure_errors_eval.diff_st3_a  = diff_st3_a;
measure_errors_eval.diff_st3_b  = diff_st3_b;
measure_errors_eval.diff_cor3_a = diff_cor3_a;
measure_errors_eval.diff_cor3_b = diff_cor3_b;
measure_errors_eval.diff_sag3_a = diff_sag3_a;
measure_errors_eval.diff_sag3_b = diff_sag3_b;

measure_errors_eval.total_diff_a = [diff1_a(:); diff3_a(:); diff3_a(:)];
measure_errors_eval.total_diff_b = [diff1_b(:); diff3_b(:); diff3_b(:)];