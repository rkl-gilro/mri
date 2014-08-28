function [F J] = myfun_nonortho(mesh0_X, in)

mesh1 = mesh0_X(1:length(mesh0_X)/3);
mesh2 = mesh0_X(length(mesh0_X)/3 + 1 :2*length(mesh0_X)/3);
mesh3 = mesh0_X(2*length(mesh0_X)/3 + 1 :end );

mesh11 = reshape(mesh1',3,length(mesh0_X)/9);
mesh21 = reshape(mesh2',3,length(mesh0_X)/9);
mesh31 = reshape(mesh3',3,length(mesh0_X)/9);

tri1 = TriRep(in.source_tri.Triangulation, mesh11'); % define the new mesh for the axial      mesh0_X(:,1:in.size_source)  mesh0_X(:,1:size(mesh0_X,2)/2)' mesh11'
tri2 = TriRep(in.source_tri.Triangulation, mesh21'); % define the new mesh for the sagittal  mesh21' mesh0_X(:,in.size_source+1:end)  mesh0_X(:,size(mesh0_X,2)/2+1:end)' mesh21'
tri3 = TriRep(in.source_tri.Triangulation, mesh31');

% First intersection
b2c_ncoord1 = baryToCart(tri1, in.current_tr, in.c2b_coord);
b2c_ncoord2 = baryToCart(tri2, in.current_tr, in.c2b_coord);

% Second intersection
b2c_ncoord1_int2 = baryToCart(tri1,  in.current_tr2,  in.c2b_coord2);
b2c_ncoord2_int2 = baryToCart(tri3,  in.current_tr2,  in.c2b_coord2);

% Third intersection
b2c_ncoord1_int3 = baryToCart(tri2,  in.current_tr3,  in.c2b_coord3);
b2c_ncoord2_int3 = baryToCart(tri3,  in.current_tr3,  in.c2b_coord3);

incr_row = in.n + in.size_source;
incr_col = 3 * in.size_source;

incr = size(in.source_tri.X,1);
total_var = 9 * incr; % 9
xyz = 3; 

sub_3 = in.sub_3;
sub_2 = in.sub_2;

vol_ax  = in.vol_ax;
vol_sag = in.vol_sag;
vol_cor = in.vol_cor;

r_ax  = size(vol_ax, 1);
c_ax  = size(vol_ax, 2);
r_sag = size(vol_sag,1);
c_sag = size(vol_sag,2);
r_cor = size(vol_cor,1);
c_cor = size(vol_cor,2);

F1 = zeros(in.n1, 1);
F2 = zeros(in.n2, 1);
F3 = zeros(in.n3, 1);


% F = zeros(in.n + 2 * in.size_source, 1);
% Js = zeros(in.n + 2 * in.size_source, 6 * in.size_source);
%  F = zeros(in.n, 1);
 
current_tr  = in.current_tr;
c2b_coord   = in.c2b_coord;
current_tr2 = in.current_tr2;
c2b_coord2  = in.c2b_coord2;
current_tr3 = in.current_tr3;
c2b_coord3  = in.c2b_coord3;

axi_m1_par = in.axial_m1_par;
sag_m1_par = in.sag_m1_par;
cor_m1_par1 = in.cor_m1_par1;
axi_m1_par2 = in.axial_m1_par2;
sag_m1_par2 = in.sag_m1_par2;
cor_m1_par2 = in.cor_m1_par2;


for i = 1 : in.n
    
    %% axial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    tmp_v1_ax = axi_m1_par{i} * [b2c_ncoord1(i,:) 1]'; % 3D point to 2D point in the frame coordinates axial_m1{sub_3(i)}
    
    %% Make sure the indexes are correct and use bilinear interpolation
    
    fl1ax = floor(tmp_v1_ax(1) + 1);
    fl2ax = floor(tmp_v1_ax(2) + 1);
    cl1ax = ceil(tmp_v1_ax(1) + 1);
    cl2ax = ceil(tmp_v1_ax(2) + 1);
    
    min_max_r1ax = min(max(fl2ax,1),r_ax);
    min_max_r2ax = min(max(cl2ax,1),r_ax);
    min_max_c1ax = min(max(fl1ax,1),c_ax);
    min_max_c2ax = min(max(cl1ax,1),c_ax);
    
    neigax = [vol_ax(min_max_r1ax, min_max_c1ax, sub_3(i)) vol_ax(min_max_r1ax, min_max_c2ax, sub_3(i));...
              vol_ax(min_max_r2ax, min_max_c1ax, sub_3(i)) vol_ax(min_max_r2ax, min_max_c2ax, sub_3(i))];
    
    new_im_ax = bilinear_interpolation(tmp_v1_ax(2), tmp_v1_ax(1), double(neigax));
       
    %% sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tmp_v1_sg = sag_m1_par{i} * [b2c_ncoord2(i,:) 1]'; % 3D point to 2D point in the frame coordinates var_cell{i,j,k} sag_m1{sub_2(i)}
    
    %% Make sure the indexes are correct
    
    fl1sg = floor(tmp_v1_sg(1) + 1);
    fl2sg = floor(tmp_v1_sg(2) + 1);
    cl1sg = ceil(tmp_v1_sg(1) + 1);
    cl2sg = ceil(tmp_v1_sg(2) + 1);
    
    min_max_r1sg = min(max(fl2sg,1),r_sag);
    min_max_r2sg = min(max(cl2sg,1),r_sag);
    min_max_c1sg = min(max(fl1sg,1),c_sag);
    min_max_c2sg = min(max(cl1sg,1),c_sag);
    
    neigsg = [vol_sag(min_max_r1sg, min_max_c1sg, sub_2(i)) vol_sag(min_max_r1sg, min_max_c2sg,sub_2(i));...
              vol_sag(min_max_r2sg, min_max_c1sg, sub_2(i)) vol_sag(min_max_r2sg, min_max_c2sg,sub_2(i))];
    
    new_im_sag = bilinear_interpolation(tmp_v1_sg(2), tmp_v1_sg(1), double(neigsg));

    %% Function
    
    F1(i) = new_im_ax - new_im_sag;

end

%% Second intersection %%
for i = 1: in.n2
    
    %% axial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tmp_v1_ax = axi_m1_par2{i} * [b2c_ncoord1_int2(i,:) 1]'; % 3D point to 2D point in the frame coordinates axial_m1_par2{i}
    
    fl  = floor(tmp_v1_ax(1) + 1);
    fl2 = floor(tmp_v1_ax(2) + 1);
    cl  = ceil(tmp_v1_ax(1) + 1);
    cl2 = ceil(tmp_v1_ax(2) + 1);
    
    min_max_r1ax = min(max(fl2,1),r_ax);
    min_max_r2ax = min(max(cl2,1),r_ax);
    min_max_c1ax = min(max(fl,1), c_ax);
    min_max_c2ax = min(max(cl,1), c_ax);
    
    neig = [vol_ax(min_max_r1ax, min_max_c1ax, in.sub_3_int2(i))   vol_ax(min_max_r1ax, min_max_c2ax, in.sub_3_int2(i));...
            vol_ax(min_max_r2ax, min_max_c1ax, in.sub_3_int2(i))   vol_ax(min_max_r2ax, min_max_c2ax, in.sub_3_int2(i))];
    
    new_im_ax = bilinear_interpolation(tmp_v1_ax(2),tmp_v1_ax(1),double(neig));
    
    %% coronal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tmp_v1_cr =  cor_m1_par1{i} * [b2c_ncoord2_int2(i,:) 1]'; % 3D point to 2D point in the frame coordinates var_cell{i,j,k}

    fl  = floor(tmp_v1_cr(1) + 1);
    fl2 = floor(tmp_v1_cr(2) + 1);
    cl  = ceil(tmp_v1_cr(1) + 1);
    cl2 = ceil(tmp_v1_cr(2) + 1);

    min_max_r1cr = min(max(fl2,1),r_cor);
    min_max_r2cr = min(max(cl2,1),r_cor);
    min_max_c1cr = min(max(fl,1), c_cor);
    min_max_c2cr = min(max(cl,1), c_cor);
    
    neig = [vol_cor(min_max_r1cr, min_max_c1cr, in.sub_2_int2(i))   vol_cor(min_max_r1cr, min_max_c2cr, in.sub_2_int2(i));...
            vol_cor(min_max_r2cr, min_max_c1cr, in.sub_2_int2(i))   vol_cor(min_max_r2cr, min_max_c2cr, in.sub_2_int2(i))];
    
    new_im_cor = bilinear_interpolation(tmp_v1_cr(2),tmp_v1_cr(1),double(neig));
    
    %% Function
    F2(i) = new_im_ax - new_im_cor;
    
end

%% Third intersection %%
for i = 1: in.n3
    
    %% sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tmp_v1_sg = sag_m1_par2{i} * [b2c_ncoord1_int3(i,:) 1]'; % 3D point to 2D point in the frame coordinates sag_m1_par2{i}
    
    fl  = floor(tmp_v1_sg(1) + 1);
    fl2 = floor(tmp_v1_sg(2) + 1);
    cl  = ceil(tmp_v1_sg(1)  + 1);
    cl2 = ceil(tmp_v1_sg(2)  + 1);
    
    min_max_r1sg = min(max(fl2,1),r_sag);
    min_max_r2sg = min(max(cl2,1),r_sag);
    min_max_c1sg = min(max(fl,1), c_sag);
    min_max_c2sg = min(max(cl,1), c_sag);
    
    neig = [vol_sag(min_max_r1sg, min_max_c1sg, in.sub_2_int3(i))   vol_sag(min_max_r1sg, min_max_c2sg, in.sub_2_int3(i));...
            vol_sag(min_max_r2sg, min_max_c1sg, in.sub_2_int3(i))   vol_sag(min_max_r2sg, min_max_c2sg, in.sub_2_int3(i))];
    
    new_im_sag = bilinear_interpolation(tmp_v1_sg(2), tmp_v1_sg(1), double(neig));
    
    %% coronal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tmp_v1_cr = cor_m1_par2{i} * [b2c_ncoord2_int3(i,:) 1]'; % 3D point to 2D point in the frame coordinates var_cell{i,j,k}cor_m1_par2{i}

    fl  = floor(tmp_v1_cr(1) + 1);
    fl2 = floor(tmp_v1_cr(2) + 1);
    cl  = ceil(tmp_v1_cr(1)  + 1);
    cl2 = ceil(tmp_v1_cr(2)  + 1);

    min_max_r1cr = min(max(fl2,1),r_cor);
    min_max_r2cr = min(max(cl2,1),r_cor);
    min_max_c1cr = min(max(fl,1), c_cor);
    min_max_c2cr = min(max(cl,1), c_cor);
    
    neig = [vol_cor(min_max_r1cr, min_max_c1cr, in.sub_3_int3(i))   vol_cor(min_max_r1cr, min_max_c2cr, in.sub_3_int3(i));...
            vol_cor(min_max_r2cr, min_max_c1cr, in.sub_3_int3(i))   vol_cor(min_max_r2cr, min_max_c2cr, in.sub_3_int3(i))];
    
    new_im_cor = bilinear_interpolation(tmp_v1_cr(2), tmp_v1_cr(1), double(neig));

    %% Function
    F3(i) = new_im_sag - new_im_cor;
    
end

% for i = 1:in.size_source
%     
%     mean_tmp1 = mean(tri1.X(in.list_edges{i},:));
%     mean_tmp2 = mean(tri2.X(in.list_edges{i},:));
%     
%     lapl_tri1 = abs(tri1.X(i,:) - mean_tmp1);
%     lapl_tri2 = abs(tri2.X(i,:) - mean_tmp2);
%     
%     norm2_lapl_tri1 = sum( lapl_tri1 );
%     norm2_lapl_tri2 = sum( lapl_tri2 );
%     
%     F(i + in.n)     = in.sqrt_lambda * norm2_lapl_tri1;
%     F(i + incr_row) = in.sqrt_lambda * norm2_lapl_tri2;
% 
% end

F = [F1;F2;F3];


if nargout > 1
    
%     J = zeros(in.n + 2*in.size_source, 6 * in.size_source);
    
    for i = 1 : in.n1
        
        %% axial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        tmp_v1_ax = axi_m1_par{i} * [b2c_ncoord1(i,:) 1]'; % 3D point to 2D point in the frame coordinates axial_m1{sub_3(i)}
        
        %% Make sure the indexes are correct and use bilinear interpolation
        
        fl1ax = floor(tmp_v1_ax(1) + 1);
        fl2ax = floor(tmp_v1_ax(2) + 1);
        cl1ax = ceil(tmp_v1_ax(1) + 1);
        cl2ax = ceil(tmp_v1_ax(2) + 1);
        
        min_max_r1ax = min(max(fl2ax,1),r_ax);
        min_max_r2ax = min(max(cl2ax,1),r_ax);
        min_max_c1ax = min(max(fl1ax,1),c_ax);
        min_max_c2ax = min(max(cl1ax,1),c_ax);

        neigax_gradx = [in.gradx_ax(min_max_r1ax, min_max_c1ax, sub_3(i)) in.gradx_ax(min_max_r1ax, min_max_c2ax,sub_3(i));...
                        in.gradx_ax(min_max_r2ax, min_max_c1ax, sub_3(i)) in.gradx_ax(min_max_r2ax, min_max_c2ax,sub_3(i))];
        
        neigax_grady = [in.grady_ax(min_max_r1ax, min_max_c1ax, sub_3(i)) in.grady_ax(min_max_r1ax, min_max_c2ax,sub_3(i));...
                        in.grady_ax(min_max_r2ax, min_max_c1ax, sub_3(i)) in.grady_ax(min_max_r2ax, min_max_c2ax,sub_3(i))];
        
        grad_ax = [-bilinear_interpolation(tmp_v1_ax(2), tmp_v1_ax(1), double(neigax_gradx))  -bilinear_interpolation(tmp_v1_ax(2), tmp_v1_ax(1), double(neigax_grady))];
        
        %% sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        tmp_v1_sg = sag_m1_par{i} * [b2c_ncoord2(i,:) 1]'; % 3D point to 2D point in the frame coordinates var_cell{i,j,k} sag_m1{sub_2(i)}
        
        %% Make sure the indexes are correct
        
        fl1sg = floor(tmp_v1_sg(1) + 1);
        fl2sg = floor(tmp_v1_sg(2) + 1);
        cl1sg = ceil(tmp_v1_sg(1) + 1);
        cl2sg = ceil(tmp_v1_sg(2) + 1);
        
        min_max_r1sg = min(max(fl2sg,1),r_sag);
        min_max_r2sg = min(max(cl2sg,1),r_sag);
        min_max_c1sg = min(max(fl1sg,1),c_sag);
        min_max_c2sg = min(max(cl1sg,1),c_sag);
        
        neigsg_gradx = [in.gradx_sag(min_max_r1sg, min_max_c1sg, sub_2(i)) in.gradx_sag(min_max_r1sg, min_max_c2sg,sub_2(i));...
                        in.gradx_sag(min_max_r2sg, min_max_c1sg, sub_2(i)) in.gradx_sag(min_max_r2sg, min_max_c2sg,sub_2(i))];
        
        neigsg_grady = [in.grady_sag(min_max_r1sg, min_max_c1sg, sub_2(i)) in.grady_sag(min_max_r1sg, min_max_c2sg,sub_2(i));...
                        in.grady_sag(min_max_r2sg, min_max_c1sg, sub_2(i)) in.grady_sag(min_max_r2sg, min_max_c2sg,sub_2(i))];
        
        grad_sag = [-bilinear_interpolation(tmp_v1_sg(2), tmp_v1_sg(1), double(neigsg_gradx))  -bilinear_interpolation(tmp_v1_sg(2), tmp_v1_sg(1), double(neigsg_grady))];

        
        %% Jacobian of warp function qi = w(pi,X), pi is a given data point, and X is the mesh
        
        w  = zeros(3, total_var);
        w2 = zeros(3, total_var);
%         w  = zeros(3, 2 * incr_col);
%         w2 = zeros(3, 2 * incr_col);
        
        %% Get the index of the points in the mesh (tetrahedron) that contains point i
        %% We multiply by 3 because we have 3 coordinates per each point
        ind_w1 = (tri1.Triangulation(current_tr(i),:)-1) * 3 + 1;
        
        w(1,ind_w1 ) = c2b_coord(i,:);
        w2(1,:)      = circshift(w(1,:)', incr_col);
        
        w(2,:)  = circshift(w(1,:)', 1);
        w(3,:)  = circshift(w(1,:)', 2);

        w2(2,:) = circshift(w2(1,:)',1);
        w2(3,:) = circshift(w2(1,:)',2);
        

        %% Jacobian
                
        J1(i,:) =  (grad_ax * axi_m1_par{i}(1:2,1:3) * w) - (grad_sag * sag_m1_par{i}(1:2,1:3) * w2);  

    end
    %% Second intersection %%
    for i = 1: in.n2
        
        %% axial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        tmp_v1_ax = axi_m1_par2{i} * [b2c_ncoord1_int2(i,:) 1]'; % 3D point to 2D point in the frame coordinates axial_m1_par2{i}
        
        fl  = floor(tmp_v1_ax(1) + 1);
        fl2 = floor(tmp_v1_ax(2) + 1);
        cl  = ceil(tmp_v1_ax(1) + 1);
        cl2 = ceil(tmp_v1_ax(2) + 1);
        
        min_max_r1ax = min(max(fl2,1),r_ax);
        min_max_r2ax = min(max(cl2,1),r_ax);
        min_max_c1ax = min(max(fl,1), c_ax);
        min_max_c2ax = min(max(cl,1), c_ax);
        
        %% Gradient in axial view
        
        neigax_gradx = [in.gradx_ax(min_max_r1ax, min_max_c1ax, in.sub_3_int2(i)) in.gradx_ax(min_max_r1ax, min_max_c2ax, in.sub_3_int2(i));...
                        in.gradx_ax(min_max_r2ax, min_max_c1ax, in.sub_3_int2(i)) in.gradx_ax(min_max_r2ax, min_max_c2ax, in.sub_3_int2(i))];
        
        neigax_grady = [in.grady_ax(min_max_r1ax, min_max_c1ax, in.sub_3_int2(i)) in.grady_ax(min_max_r1ax, min_max_c2ax, in.sub_3_int2(i));...
                        in.grady_ax(min_max_r2ax, min_max_c1ax, in.sub_3_int2(i)) in.grady_ax(min_max_r2ax, min_max_c2ax, in.sub_3_int2(i))];
        
        grad_ax = [-bilinear_interpolation(tmp_v1_ax(2), tmp_v1_ax(1), double(neigax_gradx))  -bilinear_interpolation(tmp_v1_ax(2), tmp_v1_ax(1), double(neigax_grady))];
        
        %% coronal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        tmp_v1_cr =  cor_m1_par1{i} * [b2c_ncoord2_int2(i,:) 1]'; % 3D point to 2D point in the frame coordinates var_cell{i,j,k}
        
        fl  = floor(tmp_v1_cr(1) + 1);
        fl2 = floor(tmp_v1_cr(2) + 1);
        cl  = ceil(tmp_v1_cr(1) + 1);
        cl2 = ceil(tmp_v1_cr(2) + 1);
        
        min_max_r1cr = min(max(fl2,1),r_cor);
        min_max_r2cr = min(max(cl2,1),r_cor);
        min_max_c1cr = min(max(fl,1), c_cor);
        min_max_c2cr = min(max(cl,1), c_cor);
        
        %% Gradient in coronal view
        
        neigcr_gradx = [in.gradx_cor(min_max_r1cr, min_max_c1cr, in.sub_2_int2(i)) in.gradx_cor(min_max_r1cr, min_max_c2cr, in.sub_2_int2(i));...
                        in.gradx_cor(min_max_r2cr, min_max_c1cr, in.sub_2_int2(i)) in.gradx_cor(min_max_r2cr, min_max_c2cr, in.sub_2_int2(i))];
        
        neigcr_grady = [in.grady_cor(min_max_r1cr, min_max_c1cr, in.sub_2_int2(i)) in.grady_cor(min_max_r1cr, min_max_c2cr, in.sub_2_int2(i));...
                        in.grady_cor(min_max_r2cr, min_max_c1cr, in.sub_2_int2(i)) in.grady_cor(min_max_r2cr, min_max_c2cr, in.sub_2_int2(i))];
        
        grad_cor = [-bilinear_interpolation(tmp_v1_cr(2), tmp_v1_cr(1), double(neigcr_gradx))  -bilinear_interpolation(tmp_v1_cr(2), tmp_v1_cr(1), double(neigcr_grady))];
                
        %% Jacobian
        w  = zeros(3, total_var);
        w2 = zeros(3, total_var);
        
        ind_w = (tri1.Triangulation(current_tr2(i),:)-1) * xyz + 1;
        
        w( 1,ind_w ) = c2b_coord2(i,:);
        w( 2,:) = circshift(w(1,:)', 1);
        w( 3,:) = circshift(w(1,:)', 2);

        w2(1,:) = circshift(w(1,:)', 2 * xyz * incr );
        w2(2,:) = circshift(w2(1,:)',1);
        w2(3,:) = circshift(w2(1,:)',2); 
        
        J2(i,:) = (grad_ax * axi_m1_par2{i}(1:2,1:3) * w) - (grad_cor * cor_m1_par1{i}(1:2,1:3) * w2);
        
    end
    
    %% Third intersection %%
    for i = 1: in.n3
        
        %% sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        tmp_v1_sg = sag_m1_par2{i} * [b2c_ncoord1_int3(i,:) 1]'; % 3D point to 2D point in the frame coordinates sag_m1_par2{i}
        
        fl  = floor(tmp_v1_sg(1) + 1);
        fl2 = floor(tmp_v1_sg(2) + 1);
        cl  = ceil(tmp_v1_sg(1)  + 1);
        cl2 = ceil(tmp_v1_sg(2)  + 1);
        
        min_max_r1sg = min(max(fl2,1),r_sag);
        min_max_r2sg = min(max(cl2,1),r_sag);
        min_max_c1sg = min(max(fl,1), c_sag);
        min_max_c2sg = min(max(cl,1), c_sag);
        
        %% Gradient in sagittal view
        
        neigsg_gradx = [in.gradx_sag(min_max_r1sg, min_max_c1sg, in.sub_2_int3(i)) in.gradx_sag(min_max_r1sg, min_max_c2sg, in.sub_2_int3(i));...
                        in.gradx_sag(min_max_r2sg, min_max_c1sg, in.sub_2_int3(i)) in.gradx_sag(min_max_r2sg, min_max_c2sg, in.sub_2_int3(i))];
        
        neigsg_grady = [in.grady_sag(min_max_r1sg, min_max_c1sg, in.sub_2_int3(i)) in.grady_sag(min_max_r1sg, min_max_c2sg, in.sub_2_int3(i));...
                        in.grady_sag(min_max_r2sg, min_max_c1sg, in.sub_2_int3(i)) in.grady_sag(min_max_r2sg, min_max_c2sg, in.sub_2_int3(i))];
        
        grad_sag = [-bilinear_interpolation(tmp_v1_sg(2), tmp_v1_sg(1), double(neigsg_gradx))  -bilinear_interpolation(tmp_v1_sg(2), tmp_v1_sg(1), double(neigsg_grady))];
        
        %% coronal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        tmp_v1_cr = cor_m1_par2{i} * [b2c_ncoord2_int3(i,:) 1]'; % 3D point to 2D point in the frame coordinates var_cell{i,j,k}cor_m1_par2{i}
        
        fl  = floor(tmp_v1_cr(1) + 1);
        fl2 = floor(tmp_v1_cr(2) + 1);
        cl  = ceil(tmp_v1_cr(1)  + 1);
        cl2 = ceil(tmp_v1_cr(2)  + 1);
        
        min_max_r1cr = min(max(fl2,1),r_cor);
        min_max_r2cr = min(max(cl2,1),r_cor);
        min_max_c1cr = min(max(fl,1), c_cor);
        min_max_c2cr = min(max(cl,1), c_cor);
        
        %% Gradient in coronal view
        
        neigcr_gradx = [in.gradx_cor(min_max_r1cr, min_max_c1cr, in.sub_3_int3(i)) in.gradx_cor(min_max_r1cr, min_max_c2cr, in.sub_3_int3(i));...
                        in.gradx_cor(min_max_r2cr, min_max_c1cr, in.sub_3_int3(i)) in.gradx_cor(min_max_r2cr, min_max_c2cr, in.sub_3_int3(i))];
        
        neigcr_grady = [in.grady_cor(min_max_r1cr, min_max_c1cr, in.sub_3_int3(i)) in.grady_cor(min_max_r1cr, min_max_c2cr, in.sub_3_int3(i));...
                        in.grady_cor(min_max_r2cr, min_max_c1cr, in.sub_3_int3(i)) in.grady_cor(min_max_r2cr, min_max_c2cr, in.sub_3_int3(i))];
        
        grad_cor = [-bilinear_interpolation(tmp_v1_cr(2), tmp_v1_cr(1), double(neigcr_gradx))  -bilinear_interpolation(tmp_v1_cr(2), tmp_v1_cr(1), double(neigcr_grady))];

        %% Jacobian
        w_tmp  = zeros(1, total_var);
        w  = zeros(3, total_var);
        w2 = zeros(3, total_var);
        
        ind_w = (tri3.Triangulation(current_tr3(i),:)-1) * xyz + 1;
        
        w_tmp( 1,ind_w ) = c2b_coord3(i,:);
        
        w(1,:)   = circshift(w_tmp(1,:)', xyz * incr);
        w(2,:) = circshift(w(1,:)',  1);
        w(3,:) = circshift(w(1,:)',  2); 
        
        w2(1,:)  = circshift(w_tmp(1,:)', 2 * xyz * incr);
        w2(2,:) = circshift(w2(1,:)',1);
        w2(3,:) = circshift(w2(1,:)',2); 
        
        J3(i,:) = (grad_cor * cor_m1_par2{i}(1:2,1:3) * w2) - (grad_sag * sag_m1_par2{i}(1:2,1:3) * w); %
        
    end
    
    
%     for i = 1:in.size_source
%         
%         jsi1 = (i - 1) * 3 + 1;
%         jsi2 = jsi1 + incr_col;
%         J(i + in.n    , jsi1:jsi1+2) = in.sqrt_lambda;%.* lapl_tri1; %(in.sqrt_lambda * ( 1 / norm_lapl_tri1 ) ) .* lapl_tri1;
%         J(i + incr_row, jsi2:jsi2+2) = in.sqrt_lambda;%.* lapl_tri2; %(in.sqrt_lambda * ( 1 / norm_lapl_tri2 ) ) .* lapl_tri2;
%         
%         for j = 1:length(in.list_edges{i})
%             
%             edgein1 = (in.list_edges{i}(j)-1) * 3 + 1;
%             edgein2 = edgein1 + incr_col;
%             J(i + in.n,     edgein1:edgein1+2) = J(i + in.n,     jsi1:jsi1+2).*[in.scalar(i) in.scalar(i) in.scalar(i)];
%             J(i + incr_row, edgein2:edgein2+2) = J(i + incr_row, jsi2:jsi2+2).*[in.scalar(i) in.scalar(i) in.scalar(i)];
%             
%         end
%         
%     end

     J = [J1;J2;J3];
end




