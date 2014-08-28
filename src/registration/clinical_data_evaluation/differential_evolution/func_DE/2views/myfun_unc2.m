function F_obj = myfun_unc2(mesh0_X, S_struct)


global source_tri
global var_array1
global sag_m1
global axial_m1
global vol_ax
global vol_sag
global vol_cor
global sub_1
global sub_2
global sub_3
global current_tr
global c2b_coord
global list_edges
global neig_ax
global neig_sag
global neig_cor
global gradx_ax
global grady_ax
global gradx_sag
global grady_sag
global gradx_cor
global grady_cor
global axial_m1_par
global sag_m1_par

for q = 1:size(mesh0_X,1)
    
    mesh0_X_mat = reshape(mesh0_X(q,:)',2 * size( source_tri.X,1),3);
    mesh0_X_mat = mesh0_X_mat + [source_tri.X;source_tri.X];
    
    tri1 = TriRep( source_tri.Triangulation,mesh0_X_mat(1:size( source_tri.X,1),:)); % define the new mesh for the axial
    tri2 = TriRep( source_tri.Triangulation,mesh0_X_mat(size( source_tri.X,1)+1:end,:)); % define the new mesh for the sagittal
    
    
    b2c_ncoord1 = baryToCart(tri1,  current_tr,  c2b_coord);
    b2c_ncoord2 = baryToCart(tri2,  current_tr,  c2b_coord);
    
    lambda = .01;
    sqrt_lambda = sqrt(lambda);
    
    n = size(var_array1,1);
    
    incr =  n + size(source_tri.X,1);
    
    F = zeros(size(var_array1,1) + size(source_tri.X,1) + size(source_tri.X,1),1);
    
    rows = size(vol_ax,1);
    cols = size(vol_ax,2);
    
    for i = 1: n
        
        %% axial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        tmp_v1 = axial_m1_par{i} * [b2c_ncoord1(i,:) 1]'; % 3D point to 2D point in the frame coordinates
        
        fl  = floor(tmp_v1(1) + 1);
        fl2 = floor(tmp_v1(2) + 1);
        cl  = ceil(tmp_v1(1) + 1);
        cl2 = ceil(tmp_v1(2) + 1);
        
        min_max_r  = min(max(fl2,1),rows);
        min_max_r2 = min(max(cl2,1),rows);
        min_max_c  = min(max(fl,1),cols);
        min_max_c2 = min(max(cl,1),cols);
        
        neig = [vol_ax(min_max_r, min_max_c, sub_3(i))   vol_ax(min_max_r, min_max_c2, sub_3(i));...
            vol_ax(min_max_r2,min_max_c, sub_3(i))   vol_ax(min_max_r2,min_max_c2, sub_3(i))];
        
        new_im_ax = bilinear_interpolation(tmp_v1(2),tmp_v1(1),double(neig));
        
        
        
        %% sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        tmp_v1 = sag_m1_par{i} * [b2c_ncoord2(i,:) 1]'; % 3D point to 2D point in the frame coordinates var_cell{i,j,k}
        
        fl  = floor(tmp_v1(1) + 1);
        fl2 = floor(tmp_v1(2) + 1);
        cl  = ceil(tmp_v1(1) + 1);
        cl2 = ceil(tmp_v1(2) + 1);
        
        min_max_r  = min(max(fl2,1),rows);
        min_max_r2 = min(max(cl2,1),rows);
        min_max_c  = min(max(fl,1),cols);
        min_max_c2 = min(max(cl,1),cols);
        
        neig = [vol_sag(min_max_r, min_max_c,sub_2(i))   vol_sag(min_max_r, min_max_c2,sub_2(i));...
                vol_sag(min_max_r2,min_max_c,sub_2(i))   vol_sag(min_max_r2,min_max_c2,sub_2(i))];
        
        new_im_sag = bilinear_interpolation(tmp_v1(2),tmp_v1(1),double(neig));
        
        
        %% Function
        F(i) = new_im_ax - new_im_sag;
        
    end
    
    
    for i = 1:size( source_tri.X,1)
        
        F(i +  n)   =  sqrt_lambda*(sum( (tri1.X(i,:) - mean(tri1.X( list_edges{i},:))).^2 ));
        F(i + incr) =  sqrt_lambda*(sum( (tri2.X(i,:) - mean(tri2.X( list_edges{i},:))).^2 ));
    end
    
    F_tmp(q) = sum(F.^2);
    
end
F_obj = F_tmp';