function F_obj = myfun_unc2_t(mesh0_X, S_struct)


global source_tri_v

global var_array1_v
global var_array2_v
global var_array3_v

global vol_ax_eval
global vol_sag_eval
global vol_cor_eval

global sub_2_int1_v
global sub_3_int1_v

global sub_2_int2_v
global sub_3_int2_v

global sub_2_int3_v
global sub_3_int3_v

global current_tr_int1_v
global c2b_coord_int1_v
global current_tr_int2_v
global c2b_coord_int2_v
global current_tr_int3_v
global c2b_coord_int3_v
global list_edges_v

global axial_M1
global sag_M1
global cor_M1

size_sour = size( source_tri_v.X,1);

for q = 1:size(mesh0_X,1)
    

    
    mesh0_X_m = reshape(mesh0_X(q,:)',3 * size_sour,3);
    mesh0_X_mat = mesh0_X_m + [source_tri_v.X; source_tri_v.X;source_tri_v.X];

    % Define the new mesh
%     mesh0_X_mat(1:size_sour,:) = [mesh0_X_m(1:size_sour,:) zeros(size_sour,1)] + source_tri_v.X;
%     mesh0_X_mat(size_sour+1:2*size_sour,:) = [zeros(size_sour,1) mesh0_X_m(size_sour+1:2*size_sour,:)] + source_tri_v.X;
%     mesh0_X_mat(2*size_sour+1:3*size_sour,:) = [mesh0_X_m(2*size_sour+1:end,1) zeros(size_sour,1) mesh0_X_m(2*size_sour+1:end,2)] + source_tri_v.X;
    
    
    tri1 = TriRep( source_tri_v.Triangulation,mesh0_X_mat(1:size_sour,:)); % define the new mesh for the axial
    tri2 = TriRep( source_tri_v.Triangulation,mesh0_X_mat(size_sour+1:2*size_sour,:)); % define the new mesh for the sagittal
    tri3 = TriRep( source_tri_v.Triangulation,mesh0_X_mat(2*size_sour+1:end,:));
    
    % First intersection
    b2c_ncoord1_int1 = baryToCart(tri1,  current_tr_int1_v,  c2b_coord_int1_v);
    b2c_ncoord2_int1 = baryToCart(tri2,  current_tr_int1_v,  c2b_coord_int1_v);
    
    % Second intersection
    b2c_ncoord1_int2 = baryToCart(tri1,  current_tr_int2_v,  c2b_coord_int2_v);
    b2c_ncoord2_int2 = baryToCart(tri3,  current_tr_int2_v,  c2b_coord_int2_v);
    
    % Third intersection
    b2c_ncoord1_int3 = baryToCart(tri2,  current_tr_int3_v,  c2b_coord_int3_v);
    b2c_ncoord2_int3 = baryToCart(tri3,  current_tr_int3_v,  c2b_coord_int3_v);
    
    
    lambda = .01;
    sqrt_lambda = sqrt(lambda);
    
    n1 = size(var_array1_v,1);
    n2 = size(var_array2_v,1);
    n3 = size(var_array3_v,1);
    
    n = n1 + n2 + n3;
    
    incr =  n + size_sour;
    
    F = zeros(size(var_array1_v,1) + size(var_array2_v,1) + size(var_array3_v,1) + 3*size_sour,1);
    
    rows_ax = size(vol_ax_eval,1);
    cols_ax = size(vol_ax_eval,2);
    
    rows_sag = size(vol_sag_eval,1);
    cols_sag = size(vol_sag_eval,2);
    
    rows_cor = size(vol_cor_eval,1);
    cols_cor = size(vol_cor_eval,2);
    
    %% First intersection %%
    for i = 1: n1
        
        %% axial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        tmp_v1 = axial_M1{sub_3_int1_v(i)} * [b2c_ncoord1_int1(i,:) 1]'; % 3D point to 2D point in the frame coordinates
        
        fl  = floor(tmp_v1(1) + 1);
        fl2 = floor(tmp_v1(2) + 1);
        cl  = ceil(tmp_v1(1) + 1);
        cl2 = ceil(tmp_v1(2) + 1);
        
        min_max_r  = min(max(fl2,1),rows_ax);
        min_max_r2 = min(max(cl2,1),rows_ax);
        min_max_c  = min(max(fl,1), cols_ax);
        min_max_c2 = min(max(cl,1), cols_ax);
        
        neig = [vol_ax_eval(min_max_r, min_max_c, sub_3_int1_v(i))   vol_ax_eval(min_max_r, min_max_c2, sub_3_int1_v(i));...
                vol_ax_eval(min_max_r2,min_max_c, sub_3_int1_v(i))   vol_ax_eval(min_max_r2,min_max_c2, sub_3_int1_v(i))];
        
        new_im_ax = bilinear_interpolation(tmp_v1(2),tmp_v1(1),double(neig));
        
        
        
        %% sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        tmp_v1 = sag_M1{sub_2_int1_v(i)} * [b2c_ncoord2_int1(i,:) 1]'; % 3D point to 2D point in the frame coordinates var_cell{i,j,k}
        
        fl  = floor(tmp_v1(1) + 1);
        fl2 = floor(tmp_v1(2) + 1);
        cl  = ceil(tmp_v1(1) + 1);
        cl2 = ceil(tmp_v1(2) + 1);
        
        min_max_r  = min(max(fl2,1),rows_sag);
        min_max_r2 = min(max(cl2,1),rows_sag);
        min_max_c  = min(max(fl,1), cols_sag);
        min_max_c2 = min(max(cl,1), cols_sag);
        
        neig = [vol_sag_eval(min_max_r, min_max_c,sub_2_int1_v(i))   vol_sag_eval(min_max_r, min_max_c2,sub_2_int1_v(i));...
            vol_sag_eval(min_max_r2,min_max_c,sub_2_int1_v(i))   vol_sag_eval(min_max_r2,min_max_c2,sub_2_int1_v(i))];
        
        new_im_sag = bilinear_interpolation(tmp_v1(2),tmp_v1(1),double(neig));
        
        
        %% Function
        F(i) = new_im_ax - new_im_sag;
        
    end
    
    %% Second intersection %%
    for i = 1: n2
        
        %% axial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        tmp_v1 = axial_M1{sub_3_int2_v(i)} * [b2c_ncoord1_int2(i,:) 1]'; % 3D point to 2D point in the frame coordinates axial_M1_par2{i}
        
        fl  = floor(tmp_v1(1) + 1);
        fl2 = floor(tmp_v1(2) + 1);
        cl  = ceil(tmp_v1(1) + 1);
        cl2 = ceil(tmp_v1(2) + 1);
        
        min_max_r  = min(max(fl2,1),rows_ax);
        min_max_r2 = min(max(cl2,1),rows_ax);
        min_max_c  = min(max(fl,1), cols_ax);
        min_max_c2 = min(max(cl,1), cols_ax);
        
        neig = [vol_ax_eval(min_max_r, min_max_c, sub_3_int2_v(i))   vol_ax_eval(min_max_r, min_max_c2, sub_3_int2_v(i));...
                vol_ax_eval(min_max_r2,min_max_c, sub_3_int2_v(i))   vol_ax_eval(min_max_r2,min_max_c2, sub_3_int2_v(i))];
        
        new_im_ax = bilinear_interpolation(tmp_v1(2),tmp_v1(1),double(neig));
        
        
        
        %% coronal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        tmp_v1 =  cor_M1{sub_2_int2_v(i)} * [b2c_ncoord2_int2(i,:) 1]'; % 3D point to 2D point in the frame coordinates var_cell{i,j,k}
        
        fl  = floor(tmp_v1(1) + 1);
        fl2 = floor(tmp_v1(2) + 1);
        cl  = ceil(tmp_v1(1) + 1);
        cl2 = ceil(tmp_v1(2) + 1);
        
        min_max_r  = min(max(fl2,1),rows_cor);
        min_max_r2 = min(max(cl2,1),rows_cor);
        min_max_c  = min(max(fl,1), cols_cor);
        min_max_c2 = min(max(cl,1), cols_cor);
        
        neig = [vol_cor_eval(min_max_r, min_max_c,sub_2_int2_v(i))   vol_cor_eval(min_max_r, min_max_c2,sub_2_int2_v(i));...
                vol_cor_eval(min_max_r2,min_max_c,sub_2_int2_v(i))   vol_cor_eval(min_max_r2,min_max_c2,sub_2_int2_v(i))];
        
        new_im_cor = bilinear_interpolation(tmp_v1(2),tmp_v1(1),double(neig));
        
        
        %% Function
        F(i+n1) = new_im_ax - new_im_cor;
        
        
    end
    
    %% Third intersection %%
    for i = 1: n3
        
        %% sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        tmp_v1 = sag_M1{sub_2_int3_v(i)} * [b2c_ncoord1_int3(i,:) 1]'; % 3D point to 2D point in the frame coordinates sag_M1_par2{i}
        
        fl  = floor(tmp_v1(1) + 1);
        fl2 = floor(tmp_v1(2) + 1);
        cl  = ceil(tmp_v1(1) + 1);
        cl2 = ceil(tmp_v1(2) + 1);
        
        min_max_r  = min(max(fl2,1),rows_sag);
        min_max_r2 = min(max(cl2,1),rows_sag);
        min_max_c  = min(max(fl,1), cols_sag);
        min_max_c2 = min(max(cl,1), cols_sag);
        
        neig = [vol_sag_eval(min_max_r, min_max_c, sub_2_int3_v(i))   vol_sag_eval(min_max_r, min_max_c2, sub_2_int3_v(i));...
                vol_sag_eval(min_max_r2,min_max_c, sub_2_int3_v(i))   vol_sag_eval(min_max_r2,min_max_c2, sub_2_int3_v(i))];
        
        new_im_sag = bilinear_interpolation(tmp_v1(2),tmp_v1(1),double(neig));
        
        
        
        %% coronal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        tmp_v1 = cor_M1{sub_3_int3_v(i)} * [b2c_ncoord2_int3(i,:) 1]'; % 3D point to 2D point in the frame coordinates var_cell{i,j,k}cor_M1_par2{i}
        
        fl  = floor(tmp_v1(1) + 1);
        fl2 = floor(tmp_v1(2) + 1);
        cl  = ceil(tmp_v1(1) + 1);
        cl2 = ceil(tmp_v1(2) + 1);
        
        min_max_r  = min(max(fl2,1),rows_cor);
        min_max_r2 = min(max(cl2,1),rows_cor);
        min_max_c  = min(max(fl,1), cols_cor);
        min_max_c2 = min(max(cl,1), cols_cor);
        
        neig = [vol_cor_eval(min_max_r, min_max_c,sub_3_int3_v(i))   vol_cor_eval(min_max_r, min_max_c2,sub_3_int3_v(i));...
                vol_cor_eval(min_max_r2,min_max_c,sub_3_int3_v(i))   vol_cor_eval(min_max_r2,min_max_c2,sub_3_int3_v(i))];
        
        new_im_cor = bilinear_interpolation(tmp_v1(2),tmp_v1(1),double(neig));
        
        
        %% Function
        F(i + n1 + n2) = new_im_sag - new_im_cor;
        
        
    end
    for i = 1:size( source_tri_v.X,1)
        
        F(i +  n)               =  sqrt_lambda * (sum( (tri1.X(i,:) - mean(tri1.X( list_edges_v{i},:))).^2 ));
        F(i + incr)             =  sqrt_lambda * (sum( (tri2.X(i,:) - mean(tri2.X( list_edges_v{i},:))).^2 ));
        F(i + n + 2*size_sour)  =  sqrt_lambda * (sum( (tri3.X(i,:) - mean(tri3.X( list_edges_v{i},:))).^2 ));
    end
    
    F_tmp(q) = sum(F.^2);
    
end

F_obj = F_tmp';