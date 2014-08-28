function F_obj = myfun_meta(mesh0_X, S_struct)


global source_tri
global var_array1
global var_array2
global var_array3

global vol_ax
global vol_sag
global vol_cor

global sub_2_int1
global sub_3_int1

global sub_2_int2
global sub_3_int2

global sub_2_int3
global sub_3_int3
global current_tr_int1
global c2b_coord_int1
global current_tr_int2
global c2b_coord_int2
global current_tr_int3
global c2b_coord_int3
global list_edges

global axial_m1
global sag_m1
global cor_m1

size_sour = size( source_tri.X,1);

for q = 1:size(mesh0_X,1)
    
    mesh0_X_m = reshape(mesh0_X(q,:)',3 * size_sour,2);
%     mesh0_X_mat = mesh0_X_m + [source_tri.X(:,1:2); source_tri.X(:,2:3);source_tri.X(:,1) source_tri.X(:,3)];

    mesh1 = mesh0_X_m(1:size_sour,:);
    mesh2 = mesh0_X_m(size_sour + 1 :2 * size_sour,: );
    mesh3 = mesh0_X_m(2*size_sour + 1 :end,: );
    
%     mesh11 = reshape(mesh1',2,size_sour/6);
%     mesh21 = reshape(mesh2',2,size_sour/6);
%     mesh31 = reshape(mesh3',2,size_sour/6);
    
    
    tri1 = TriRep(source_tri.Triangulation, [mesh1+ source_tri.X(:,1:2) source_tri.X(:,3)]);  
    tri2 = TriRep(source_tri.Triangulation, [source_tri.X(:,1) mesh2+ source_tri.X(:,2:3)]);  
    tri3 = TriRep(source_tri.Triangulation, [mesh3(:,1) + source_tri.X(:,1) source_tri.X(:,2) mesh3(:,2) + source_tri.X(:,3)]);  
    
    % First intersection
    b2c_ncoord1_int1 = baryToCart(tri1,  current_tr_int1,  c2b_coord_int1);
    b2c_ncoord2_int1 = baryToCart(tri2,  current_tr_int1,  c2b_coord_int1);
    
    % Second intersection
    b2c_ncoord1_int2 = baryToCart(tri1,  current_tr_int2,  c2b_coord_int2);
    b2c_ncoord2_int2 = baryToCart(tri3,  current_tr_int2,  c2b_coord_int2);
    
    % Third intersection
    b2c_ncoord1_int3 = baryToCart(tri2,  current_tr_int3,  c2b_coord_int3);
    b2c_ncoord2_int3 = baryToCart(tri3,  current_tr_int3,  c2b_coord_int3);
    
    lambda = .1;
    sqrt_lambda = sqrt(lambda);
    
    n1 = size(var_array1,1);
    n2 = size(var_array2,1);
    n3 = size(var_array3,1);
    
    incr = size(source_tri.X,1);
    total_var = 6 * incr; % 9
    xyz = 2; % 3
    % incr_row = in.n + in.size_source;
    
    F1  = zeros(size(var_array1,1),1); % + size(source_tri.X,1) + size(source_tri.X,1),1);
    F2  = zeros(size(var_array2,1),1);
    F3  = zeros(size(var_array3,1),1);
    F_s = zeros(3*size(mesh0_X,1),1);
    
    r_ax  = size(vol_ax, 1);
    c_ax  = size(vol_ax, 2);
    r_sag = size(vol_sag,1);
    c_sag = size(vol_sag,2);
    r_cor = size(vol_sag,1);
    c_cor = size(vol_sag,2);
    
    %% First intersection %%
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
    
    %% Second intersection %%
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
        F2(i) = new_im_cor - new_im_ax;
        
    end
    %% Third intersection %%
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
    
    lapl_tri1 = zeros(size( source_tri.X,1), xyz);
    lapl_tri2 = zeros(size( source_tri.X,1), xyz);
    lapl_tri3 = zeros(size( source_tri.X,1), xyz);
    
    % Smooth term %%
    for i = 1:size( source_tri.X,1)
        
        mean_tmp1 = mean(tri1.X(list_edges{i},1:2)); % axial (x,y)
        mean_tmp2 = mean(tri2.X(list_edges{i},2:3)); % sagittal (y,z)
        mean_tmp3 = mean([tri3.X(list_edges{i},1) tri3.X(list_edges{i},3)]); % coronal (x,z)
        
        %% L1 norm
        %     lapl_tri1 = abs(tri1.X(i,1:2) - mean_tmp1);
        %     lapl_tri2 = abs(tri2.X(i,2:3) - mean_tmp2);
        %     lapl_tri3 = abs([tri3.X(i,1) tri3.X(i,3)] - mean_tmp3);
        
        %% L2 norm
        lapl_tri1(i,:) = tri1.X(i,1:2) - mean_tmp1;
        lapl_tri2(i,:) = tri2.X(i,2:3) - mean_tmp2;
        lapl_tri3(i,:) = [tri3.X(i,1) tri3.X(i,3)] - mean_tmp3;
        
        %% Functional %%
        F_s( i )         =  lambda * sum( lapl_tri1(i,:).^2 ); % lambda * lapl_tri1; %
        F_s(i + incr)    =  lambda * sum( lapl_tri2(i,:).^2 ); % lambda * lapl_tri2; %
        F_s(i + 2*incr)  =  lambda * sum( lapl_tri3(i,:).^2 ); % lambda * lapl_tri3; %
        
    end
    
    F(q) = (.5/(n1 + n2 + n3 + 3*size( source_tri.X,1)) ) .* (sum(F1.^2) + sum(F2.^2) + sum(F3.^2) + sum(F_s));
    
    
end

F_obj = F';