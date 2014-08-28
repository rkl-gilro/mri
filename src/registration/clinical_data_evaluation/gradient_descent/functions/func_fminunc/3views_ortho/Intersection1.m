function F1 = Intersection1(n1, b2c_ncoord1_int1, axial_m1, sub_3_int1, vol_ax, b2c_ncoord2_int1, sag_m1, sub_2_int1, vol_sag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F1 = zeros(1, n1);

[r_ax c_ax, ~] = size(vol_ax);
[r_sag, c_sag, ~] = size(vol_sag);

for i = 1: n1
    
    %% axial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tmp_v1_ax = axial_m1{sub_3_int1(i)} * [b2c_ncoord1_int1(i,:) 1]'; % 3D point to 2D point in the frame coordinates axial_m1_par1{i}
    
    fl  = floor(tmp_v1_ax(1) + 1);
    fl2 = floor(tmp_v1_ax(2) + 1);
    cl  = ceil(tmp_v1_ax(1) + 1);
    cl2 = ceil(tmp_v1_ax(2) + 1);
    
    min_max_r1ax = min(max(fl2,1),r_ax);
    min_max_r2ax = min(max(cl2,1),r_ax);
    min_max_c1ax = min(max(fl,1), c_ax);
    min_max_c2ax = min(max(cl,1), c_ax);
    
    neig = [vol_ax(min_max_r1ax, min_max_c1ax, sub_3_int1(i))   vol_ax(min_max_r1ax, min_max_c2ax, sub_3_int1(i));...
            vol_ax(min_max_r2ax, min_max_c1ax, sub_3_int1(i))   vol_ax(min_max_r2ax, min_max_c2ax, sub_3_int1(i))];
    
    new_im_ax = bilinear_interpolation(tmp_v1_ax(2), tmp_v1_ax(1), double(neig));
    
    %% sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tmp_v1_sg = sag_m1{sub_2_int1(i)} * [b2c_ncoord2_int1(i,:) 1]'; % 3D point to 2D point in the frame coordinates var_cell{i,j,k} sag_m1_par1{i}

    fl  = floor(tmp_v1_sg(1) + 1);
    fl2 = floor(tmp_v1_sg(2) + 1);
    cl  = ceil(tmp_v1_sg(1) + 1);
    cl2 = ceil(tmp_v1_sg(2) + 1);

    min_max_r1sg = min(max(fl2,1),r_sag);
    min_max_r2sg = min(max(cl2,1),r_sag);
    min_max_c1sg = min(max(fl,1), c_sag);
    min_max_c2sg = min(max(cl,1), c_sag);
    
    neig = [vol_sag(min_max_r1sg, min_max_c1sg, sub_2_int1(i))   vol_sag(min_max_r1sg, min_max_c2sg, sub_2_int1(i));...
            vol_sag(min_max_r2sg, min_max_c1sg, sub_2_int1(i))   vol_sag(min_max_r2sg, min_max_c2sg, sub_2_int1(i))];
    
    new_im_sag = bilinear_interpolation(tmp_v1_sg(2), tmp_v1_sg(1), double(neig));
     
    %% Function
    F1(i) = new_im_ax - new_im_sag;
    
end