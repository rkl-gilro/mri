function F2 = Intersection2(n2, b2c_ncoord1_int2, axial_m1, sub_3_int2, vol_ax, b2c_ncoord2_int2, cor_m1, sub_2_int2, vol_cor)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F2 = zeros(1, n2);

[r_ax c_ax, ~] = size(vol_ax);
[r_cor, c_cor, ~] = size(vol_cor);

for i = 1: n2
    
    %% axial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tmp_v1_ax = axial_m1{sub_3_int2(i)} * [b2c_ncoord1_int2(i,:) 1]'; % 3D point to 2D point in the frame coordinates axial_m1_par2{i}
    
    fl  = floor(tmp_v1_ax(1) + 1);
    fl2 = floor(tmp_v1_ax(2) + 1);
    cl  = ceil(tmp_v1_ax(1) + 1);
    cl2 = ceil(tmp_v1_ax(2) + 1);
    
    min_max_r1ax = min(max(fl2,1),r_ax);
    min_max_r2ax = min(max(cl2,1),r_ax);
    min_max_c1ax = min(max(fl,1), c_ax);
    min_max_c2ax = min(max(cl,1), c_ax);
    
    neig = [vol_ax(min_max_r1ax, min_max_c1ax, sub_3_int2(i))   vol_ax(min_max_r1ax, min_max_c2ax, sub_3_int2(i));...
            vol_ax(min_max_r2ax, min_max_c1ax, sub_3_int2(i))   vol_ax(min_max_r2ax, min_max_c2ax, sub_3_int2(i))];
    
    new_im_ax = bilinear_interpolation(tmp_v1_ax(2),tmp_v1_ax(1),double(neig));
    
    %% coronal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tmp_v1_cr =  cor_m1{sub_2_int2(i)} * [b2c_ncoord2_int2(i,:) 1]'; % 3D point to 2D point in the frame coordinates var_cell{i,j,k}

    fl  = floor(tmp_v1_cr(1) + 1);
    fl2 = floor(tmp_v1_cr(2) + 1);
    cl  = ceil(tmp_v1_cr(1) + 1);
    cl2 = ceil(tmp_v1_cr(2) + 1);

    min_max_r1cr = min(max(fl2,1),r_cor);
    min_max_r2cr = min(max(cl2,1),r_cor);
    min_max_c1cr = min(max(fl,1), c_cor);
    min_max_c2cr = min(max(cl,1), c_cor);
    
    neig = [vol_cor(min_max_r1cr, min_max_c1cr, sub_2_int2(i))   vol_cor(min_max_r1cr, min_max_c2cr, sub_2_int2(i));...
            vol_cor(min_max_r2cr, min_max_c1cr, sub_2_int2(i))   vol_cor(min_max_r2cr, min_max_c2cr, sub_2_int2(i))];
    
    new_im_cor = bilinear_interpolation(tmp_v1_cr(2),tmp_v1_cr(1),double(neig));
    
    %% Function
    F2(i) = new_im_ax - new_im_cor;
    
end