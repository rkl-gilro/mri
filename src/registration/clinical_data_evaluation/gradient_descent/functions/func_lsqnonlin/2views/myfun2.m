function F = myfun2(mesh0_X,in,out)

tri1 = TriRep(in.source_tri.Triangulation,mesh0_X(1:size(in.source_tri.X,1),:)); % define the new mesh for the axial
tri2 = TriRep(in.source_tri.Triangulation,mesh0_X(size(in.source_tri.X,1)+1:end,:)); % define the new mesh for the sagittal


b2c_ncoord1 = baryToCart(tri1, in.current_tr, in.c2b_coord);
b2c_ncoord2 = baryToCart(tri2, in.current_tr, in.c2b_coord);

incr = in.n + in.size_source;

for i = 1:in.n
    
    %% axial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tmp_v1 = in.axial_m1{in.sub_3(i)} * [b2c_ncoord1(i,:) 1]'; % 3D point to 2D point in the frame coordinates
    
    fl  = floor(tmp_v1(1) + 1);
    fl2 = floor(tmp_v1(2) + 1);
    cl  = ceil(tmp_v1(1) + 1);
    cl2 = ceil(tmp_v1(2) + 1);
    
    min_max_r = min(max(fl2,1),in.rows);
    min_max_r2 = min(max(cl2,1),in.rows);
    min_max_c = min(max(fl,1),in.cols);
    min_max_c2 = min(max(cl,1),in.cols);
    
    neig = [in.vol_ax(min_max_r,min_max_c,in.sub_3(i))   in.vol_ax(min_max_r,min_max_c2,in.sub_3(i));...
            in.vol_ax(min_max_r2,min_max_c,in.sub_3(i))  in.vol_ax(min_max_r2,min_max_c2,in.sub_3(i))];
    new_im_ax = bilinear_interpolation(tmp_v1(2),tmp_v1(1),double(neig));
    
    
    
    %% sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tmp_v1 = in.sag_m1{in.sub_2(i)} * [b2c_ncoord2(i,:) 1]'; % 3D point to 2D point in the frame coordinates var_cell{i,j,k}

    fl  = floor(tmp_v1(1) + 1);
    fl2 = floor(tmp_v1(2) + 1);
    cl  = ceil(tmp_v1(1) + 1);
    cl2 = ceil(tmp_v1(2) + 1);

    min_max_r  = min(max(fl2,1),in.rows);
    min_max_r2 = min(max(cl2,1),in.rows);
    min_max_c  = min(max(fl,1),in.cols);
    min_max_c2 = min(max(cl,1),in.cols);
    
    neig = [in.vol_sag(min_max_r,min_max_c,in.sub_2(i))   in.vol_sag(min_max_r,min_max_c2,in.sub_2(i));...
            in.vol_sag(min_max_r2,min_max_c,in.sub_2(i)) in.vol_sag(min_max_r2, min_max_c2,in.sub_2(i))];
    new_im_sag = bilinear_interpolation(tmp_v1(2),tmp_v1(1),double(neig));
    
    
    %% Function
    out.F(i) = new_im_ax - new_im_sag;
    
end


for i = 1:size(in.source_tri.X,1)
    
    out.F(i + in.n) = in.sqrt_lambda*(norm( (tri1.X(i,:) - mean(tri1.X(in.list_edges{i},:))) ));
    out.F(i + incr) = in.sqrt_lambda*(norm( (tri2.X(i,:) - mean(tri2.X(in.list_edges{i},:))) ));
end

F = out.F;

