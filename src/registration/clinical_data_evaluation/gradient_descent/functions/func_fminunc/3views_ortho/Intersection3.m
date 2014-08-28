function F3 = Intersection3(n3, b2c_ncoord1_int3, sag_m1, sub_2_int3, vol_sag, b2c_ncoord2_int3, cor_m1, sub_3_int3, vol_cor)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F3 = zeros(1,n3);

[r_cor, c_cor, ~] = size(vol_cor);
[r_sag, c_sag, ~] = size(vol_sag);

for i = 1: n3
    
    %% sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tmp_v1_sg = sag_m1{sub_2_int3(i)} * [b2c_ncoord1_int3(i,:) 1]'; % 3D point to 2D point in the frame coordinates sag_m1_par2{i}
    
    fl  = floor(tmp_v1_sg(1) + 1);
    fl2 = floor(tmp_v1_sg(2) + 1);
    cl  = ceil(tmp_v1_sg(1)  + 1);
    cl2 = ceil(tmp_v1_sg(2)  + 1);
    
    min_max_r1sg = min(max(fl2,1),r_sag);
    min_max_r2sg = min(max(cl2,1),r_sag);
    min_max_c1sg = min(max(fl,1), c_sag);
    min_max_c2sg = min(max(cl,1), c_sag);
    
    neig = [vol_sag(min_max_r1sg, min_max_c1sg, sub_2_int3(i))   vol_sag(min_max_r1sg, min_max_c2sg, sub_2_int3(i));...
            vol_sag(min_max_r2sg, min_max_c1sg, sub_2_int3(i))   vol_sag(min_max_r2sg, min_max_c2sg, sub_2_int3(i))];
    
    new_im_sag = bilinear_interpolation(tmp_v1_sg(2), tmp_v1_sg(1), double(neig));
    
    %% coronal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tmp_v1_cr = cor_m1{sub_3_int3(i)} * [b2c_ncoord2_int3(i,:) 1]'; % 3D point to 2D point in the frame coordinates var_cell{i,j,k}cor_m1_par2{i}

    fl  = floor(tmp_v1_cr(1) + 1);
    fl2 = floor(tmp_v1_cr(2) + 1);
    cl  = ceil(tmp_v1_cr(1)  + 1);
    cl2 = ceil(tmp_v1_cr(2)  + 1);

    min_max_r1cr = min(max(fl2,1),r_cor);
    min_max_r2cr = min(max(cl2,1),r_cor);
    min_max_c1cr = min(max(fl,1), c_cor);
    min_max_c2cr = min(max(cl,1), c_cor);
    
    neig = [vol_cor(min_max_r1cr, min_max_c1cr, sub_3_int3(i))   vol_cor(min_max_r1cr, min_max_c2cr, sub_3_int3(i));...
            vol_cor(min_max_r2cr, min_max_c1cr, sub_3_int3(i))   vol_cor(min_max_r2cr, min_max_c2cr, sub_3_int3(i))];
    
    new_im_cor = bilinear_interpolation(tmp_v1_cr(2), tmp_v1_cr(1), double(neig));

    %% Function
    F3(i) = new_im_sag - new_im_cor;
    
end