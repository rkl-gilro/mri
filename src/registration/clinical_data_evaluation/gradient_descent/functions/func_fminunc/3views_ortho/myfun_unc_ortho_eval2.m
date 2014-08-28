function [F_obj, J] = myfun_unc_ortho_eval2(mesh0_X)


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

global sub_2_int1_v2
global sub_3_int1_v2

global sub_2_int2_v2
global sub_3_int2_v2

global sub_2_int3_v2
global sub_3_int3_v2
global current_tr_int1_v2
global c2b_coord_int1_v2
global current_tr_int2_v2
global c2b_coord_int2_v2
global current_tr_int3_v2
global c2b_coord_int3_v2

global tetra2
% global list_edges_v
global list_edges_v2

global scalar_v2

global axial_M1
global sag_M1
global cor_M1

global grady_ax
global gradx_ax

global grady_sag
global gradx_sag

global grady_cor
global gradx_cor

global crop_volume
global crop_rectangle

global out_ind

mesh1 = mesh0_X(1                       : length(mesh0_X)/3);
mesh2 = mesh0_X(length(mesh0_X)/3 + 1   : 2 * length(mesh0_X)/3 );
mesh3 = mesh0_X(2*length(mesh0_X)/3 + 1 : end );

mesh11 = reshape(mesh1',2,length(mesh0_X)/6);
mesh21 = reshape(mesh2',2,length(mesh0_X)/6);
mesh31 = reshape(mesh3',2,length(mesh0_X)/6);

%% Define a new mesh of points first intersection
tmp_points = zeros(size(source_tri_v.X));
tmp_points(out_ind==1,:) = source_tri_v.X(out_ind==1,:);
tmp_points(out_ind==0,:) = [mesh11' source_tri_v.X(out_ind==0,3)]; 

tri1 = TriRep( source_tri_v.Triangulation,[tmp_points(:,1) tmp_points(:,2) tmp_points(:,3)]); % define the new mesh for the axial
tri12 =  TriRep( tetra2, [tmp_points(out_ind==0,1) tmp_points(out_ind==0,2) tmp_points(out_ind==0,3)]); 

%% Define a new mesh of points second intersection
tmp_points = zeros(size(source_tri_v.X));
tmp_points(out_ind==1,:) = source_tri_v.X(out_ind==1,:);
tmp_points(out_ind==0,:) = [source_tri_v.X(out_ind==0,1) mesh21']; 

tri2 = TriRep( source_tri_v.Triangulation,[tmp_points(:,1) tmp_points(:,2) tmp_points(:,3)]); % define the new mesh for the sagittal
tri22 =  TriRep( tetra2, [tmp_points(out_ind==0,1) tmp_points(out_ind==0,2) tmp_points(out_ind==0,3)]); 

%% Define a new mesh of points third intersection
tmp_points = zeros(size(source_tri_v.X));
tmp_points(out_ind==1,:) = source_tri_v.X(out_ind==1,:);
tmp_points(out_ind==0,:) = [mesh31(1,:)' source_tri_v.X(out_ind==0,2) mesh31(2,:)']; 

tri3 = TriRep( source_tri_v.Triangulation, [tmp_points(:,1) tmp_points(:,2) tmp_points(:,3)]); % define the new mesh for the coronal
tri32 =  TriRep( tetra2, [tmp_points(out_ind==0,1) tmp_points(out_ind==0,2) tmp_points(out_ind==0,3)]); 

% First intersection
b2c_ncoord1_int1 = baryToCart(tri12,  current_tr_int1_v2,  c2b_coord_int1_v2);
b2c_ncoord2_int1 = baryToCart(tri22,  current_tr_int1_v2,  c2b_coord_int1_v2);

% Second intersection
b2c_ncoord1_int2 = baryToCart(tri12,  current_tr_int2_v2,  c2b_coord_int2_v2);
b2c_ncoord2_int2 = baryToCart(tri32,  current_tr_int2_v2,  c2b_coord_int2_v2);

% Third intersection
b2c_ncoord1_int3 = baryToCart(tri22,  current_tr_int3_v2,  c2b_coord_int3_v2);
b2c_ncoord2_int3 = baryToCart(tri32,  current_tr_int3_v2,  c2b_coord_int3_v2);

lambda = .01;
sqrt_lambda = sqrt(lambda);

n1 = size(var_array1_v,1);
n2 = size(var_array2_v,1);
n3 = size(var_array3_v,1);

incr = length( find(out_ind == 0) );%size(source_tri_v.X,1);
total_var = 6 * incr; % 9
xyz = 2; % 3

F1  = zeros(size(var_array1_v,1),1); % + size(source_tri_v.X,1) + size(source_tri_v.X,1),1);
F2  = zeros(size(var_array2_v,1),1); 
F3  = zeros(size(var_array3_v,1),1); 
F_s = zeros(3*size(mesh0_X,1),1);

r_ax = size(vol_ax_eval,1);
c_ax = size(vol_ax_eval,2);

r_sag = size(vol_sag_eval,1);
c_sag = size(vol_sag_eval,2);

r_cor = size(vol_cor_eval,1);
c_cor = size(vol_cor_eval,2);

max_val = 1;%max([max(vol_ax_eval(:)) max(vol_sag_eval(:)) max(vol_cor_eval(:))]);

%% First intersection %%
for i = 1: n1
    
    %% axial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if var_array1_v(i,1) ~= -Inf
        
        tmp_v1_ax = axial_M1{sub_3_int1_v(i)} * [b2c_ncoord1_int1(i,:) 1]'; % 3D point to 2D point in the frame coordinates axial_M1_par1{i}
        
        
        if crop_volume
            fl1 = floor(tmp_v1_ax(1) - crop_rectangle(1) + 1); % Cols
            fl2 = floor(tmp_v1_ax(2) - crop_rectangle(2) + 1); % Rows
            cl1 = ceil(tmp_v1_ax(1)  - crop_rectangle(1) + 1); % Cols
            cl2 = ceil(tmp_v1_ax(2)  - crop_rectangle(2) + 1); % Rows
        else
            fl1 = floor(tmp_v1_ax(1) + 1); % Cols
            fl2 = floor(tmp_v1_ax(2) + 1); % Rows
            cl1 = ceil(tmp_v1_ax(1)  + 1); % Cols
            cl2 = ceil(tmp_v1_ax(2)  + 1); % Rows
        end
    
        min_max_r1ax = min(max(fl2,1),r_ax);
        min_max_r2ax = min(max(cl2,1),r_ax);
        min_max_c1ax = min(max(fl1,1), c_ax);
        min_max_c2ax = min(max(cl1,1), c_ax);
        
        neig = [vol_ax_eval(min_max_r1ax, min_max_c1ax, sub_3_int1_v(i))   vol_ax_eval(min_max_r1ax, min_max_c2ax, sub_3_int1_v(i));...
                vol_ax_eval(min_max_r2ax, min_max_c1ax, sub_3_int1_v(i))   vol_ax_eval(min_max_r2ax, min_max_c2ax, sub_3_int1_v(i))];
        
        new_im_ax = bilinear_interpolation(tmp_v1_ax(2),tmp_v1_ax(1),double(neig));
        
           
        %% sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        tmp_v1_sg = sag_M1{sub_2_int1_v(i)} * [b2c_ncoord2_int1(i,:) 1]'; % 3D point to 2D point in the frame coordinates var_cell{i,j,k} sag_M1_par1{i}
   
        if crop_volume
            fl1 = floor(tmp_v1_sg(1,:) - crop_rectangle(2) + 1);
            fl2 = floor(tmp_v1_sg(2,:) + 1);
            cl1 = ceil(tmp_v1_sg(1,:)  - crop_rectangle(2) + 1);
            cl2 = ceil(tmp_v1_sg(2,:)  + 1);
        else
            fl1 = floor(tmp_v1_sg(1,:) + 1);
            fl2 = floor(tmp_v1_sg(2,:) + 1);
            cl1 = ceil(tmp_v1_sg(1,:)  + 1);
            cl2 = ceil(tmp_v1_sg(2,:)  + 1);
        end
    
    
        min_max_r1sg = min(max(fl2,1), r_sag);
        min_max_r2sg = min(max(cl2,1), r_sag);
        min_max_c1sg = min(max(fl1,1), c_sag);
        min_max_c2sg = min(max(cl1,1), c_sag);
        
        neig = [vol_sag_eval(min_max_r1sg, min_max_c1sg,sub_2_int1_v(i))   vol_sag_eval(min_max_r1sg, min_max_c2sg,sub_2_int1_v(i));...
                vol_sag_eval(min_max_r2sg, min_max_c1sg,sub_2_int1_v(i))   vol_sag_eval(min_max_r2sg, min_max_c2sg,sub_2_int1_v(i))];
        
        new_im_sag = bilinear_interpolation(tmp_v1_sg(2),tmp_v1_sg(1),double(neig));
        
        
        %% Function
        F1(i) = (new_im_ax - new_im_sag)/max_val;
    else
        F1(i) = 0
    end
    
end

%% Second intersection %%
for i = 1: n2
    
    %% axial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if var_array2_v(i,1) ~= -Inf
        
        tmp_v1_ax = axial_M1{sub_3_int2_v(i)} * [b2c_ncoord1_int2(i,:) 1]'; % 3D point to 2D point in the frame coordinates axial_M1_par2{i}
  
        if crop_volume
            fl1 = floor(tmp_v1_ax(1) - crop_rectangle(1) + 1); % Cols
            fl2 = floor(tmp_v1_ax(2) - crop_rectangle(2) + 1); % Rows
            cl1 = ceil(tmp_v1_ax(1)  - crop_rectangle(1) + 1); % Cols
            cl2 = ceil(tmp_v1_ax(2)  - crop_rectangle(2) + 1); % Rows
        else
            fl1 = floor(tmp_v1_ax(1) + 1); % Cols
            fl2 = floor(tmp_v1_ax(2) + 1); % Rows
            cl1 = ceil(tmp_v1_ax(1)  + 1); % Cols
            cl2 = ceil(tmp_v1_ax(2)  + 1); % Rows
        end
        
        min_max_r1ax = min(max(fl2,1), r_ax);
        min_max_r2ax = min(max(cl2,1), r_ax);
        min_max_c1ax = min(max(fl1,1), c_ax);
        min_max_c2ax = min(max(cl1,1), c_ax);
        
        neig = [vol_ax_eval(min_max_r1ax, min_max_c1ax, sub_3_int2_v(i))   vol_ax_eval(min_max_r1ax, min_max_c2ax, sub_3_int2_v(i));...
                vol_ax_eval(min_max_r2ax, min_max_c1ax, sub_3_int2_v(i))   vol_ax_eval(min_max_r2ax, min_max_c2ax, sub_3_int2_v(i))];
        
        new_im_ax = bilinear_interpolation(tmp_v1_ax(2),tmp_v1_ax(1),double(neig));
        
        
        
        %% coronal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        tmp_v1_cr =  cor_M1{sub_2_int2_v(i)} * [b2c_ncoord2_int2(i,:) 1]'; % 3D point to 2D point in the frame coordinates var_cell{i,j,k}
  
        if crop_volume
            fl1 = floor(tmp_v1_cr(1,:) - crop_rectangle(1) + 1);
            fl2 = floor(tmp_v1_cr(2,:) + 1);
            cl1 = ceil(tmp_v1_cr(1,:)  - crop_rectangle(1) + 1);
            cl2 = ceil(tmp_v1_cr(2,:)  + 1);
        else
            fl1 = floor(tmp_v1_cr(1,:) + 1);
            fl2 = floor(tmp_v1_cr(2,:) + 1);
            cl1 = ceil(tmp_v1_cr(1,:)  + 1);
            cl2 = ceil(tmp_v1_cr(2,:)  + 1);
        end
    
        min_max_r1cr = min(max(fl2,1),r_cor);
        min_max_r2cr = min(max(cl2,1),r_cor);
        min_max_c1cr = min(max(fl1,1), c_cor);
        min_max_c2cr = min(max(cl1,1), c_cor);
        
        neig = [vol_cor_eval(min_max_r1cr, min_max_c1cr,sub_2_int2_v(i))   vol_cor_eval(min_max_r1cr, min_max_c2cr,sub_2_int2_v(i));...
                vol_cor_eval(min_max_r2cr, min_max_c1cr,sub_2_int2_v(i))   vol_cor_eval(min_max_r2cr, min_max_c2cr,sub_2_int2_v(i))];
        
        new_im_cor = bilinear_interpolation(tmp_v1_cr(2),tmp_v1_cr(1),double(neig));
        
        
        %% Function
        F2(i) = (new_im_ax - new_im_cor)/max_val;
    else
        F2(i) = 0
    end
    
end
%% Third intersection %%
for i = 1: n3
    
    %% sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if var_array3_v(i,1) ~= -Inf
        
        tmp_v1_sg = sag_M1{sub_2_int3_v(i)} * [b2c_ncoord1_int3(i,:) 1]'; % 3D point to 2D point in the frame coordinates sag_M1_par2{i}

        if crop_volume
            fl1 = floor(tmp_v1_sg(1,:) - crop_rectangle(2) + 1);
            fl2 = floor(tmp_v1_sg(2,:) + 1);
            cl1 = ceil(tmp_v1_sg(1,:)  - crop_rectangle(2) + 1);
            cl2 = ceil(tmp_v1_sg(2,:)  + 1);
        else
            fl1 = floor(tmp_v1_sg(1,:) + 1);
            fl2 = floor(tmp_v1_sg(2,:) + 1);
            cl1 = ceil(tmp_v1_sg(1,:)  + 1);
            cl2 = ceil(tmp_v1_sg(2,:)  + 1);
        end
        
        min_max_r1sg = min(max(fl2,1), r_sag);
        min_max_r2sg = min(max(cl2,1), r_sag);
        min_max_c1sg = min(max(fl1,1), c_sag);
        min_max_c2sg = min(max(cl1,1), c_sag);
        
        neig = [vol_sag_eval(min_max_r1sg, min_max_c1sg, sub_2_int3_v(i))   vol_sag_eval(min_max_r1sg, min_max_c2sg, sub_2_int3_v(i));...
                vol_sag_eval(min_max_r2sg, min_max_c1sg, sub_2_int3_v(i))   vol_sag_eval(min_max_r2sg, min_max_c2sg, sub_2_int3_v(i))];
        
        new_im_sag = bilinear_interpolation(tmp_v1_sg(2),tmp_v1_sg(1),double(neig));
        
        
        
        %% coronal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        tmp_v1_cr = cor_M1{sub_3_int3_v(i)} * [b2c_ncoord2_int3(i,:) 1]'; % 3D point to 2D point in the frame coordinates var_cell{i,j,k}cor_M1_par2{i}

        if crop_volume
            fl1 = floor(tmp_v1_cr(1,:) - crop_rectangle(1) + 1);
            fl2 = floor(tmp_v1_cr(2,:) + 1);
            cl1 = ceil(tmp_v1_cr(1,:)  - crop_rectangle(1) + 1);
            cl2 = ceil(tmp_v1_cr(2,:)  + 1);
        else
            fl1 = floor(tmp_v1_cr(1,:) + 1);
            fl2 = floor(tmp_v1_cr(2,:) + 1);
            cl1 = ceil(tmp_v1_cr(1,:)  + 1);
            cl2 = ceil(tmp_v1_cr(2,:)  + 1);
        end
        
        min_max_r1cr = min(max(fl2,1),r_cor);
        min_max_r2cr = min(max(cl2,1),r_cor);
        min_max_c1cr = min(max(fl1,1), c_cor);
        min_max_c2cr = min(max(cl1,1), c_cor);
        
        neig = [vol_cor_eval(min_max_r1cr, min_max_c1cr,sub_3_int3_v(i))   vol_cor_eval(min_max_r1cr, min_max_c2cr,sub_3_int3_v(i));...
                vol_cor_eval(min_max_r2cr, min_max_c1cr,sub_3_int3_v(i))   vol_cor_eval(min_max_r2cr, min_max_c2cr,sub_3_int3_v(i))];
        
        new_im_cor = bilinear_interpolation(tmp_v1_cr(2),tmp_v1_cr(1),double(neig));
        
        
        %% Function
        F3(i) = (new_im_cor - new_im_sag)/max_val;
    else
        F3(i) = 0
    end
    
end

lapl_tri1 = zeros(incr, xyz);
lapl_tri2 = zeros(incr, xyz);
lapl_tri3 = zeros(incr, xyz);

aux_array = find(out_ind == 0) ;

tri1_in = tri1.X(out_ind == 0,:);
tri2_in = tri3.X(out_ind == 0,:);
tri3_in = tri3.X(out_ind == 0,:);

%% Smooth term %%
for i = 1:incr % size( source_tri_v.X,1)
    
%     F_s( i )         =  sqrt_lambda * (sum( (tri1.X(i,:) - mean(tri1.X( list_edges_v{i},:))).^2 ));
%     F_s(i + incr)    =  sqrt_lambda * (sum( (tri2.X(i,:) - mean(tri2.X( list_edges_v{i},:))).^2 ));
%     F_s(i + 2*incr)  =  sqrt_lambda * (sum( (tri3.X(i,:) - mean(tri3.X( list_edges_v{i},:))).^2 ));
    
    mean_tmp1 = mean(tri12.X(list_edges_v2{i},1:2)); %mean(tri1.X(list_edges_v2{i},1:2)); % axial (x,y)
    mean_tmp2 = mean(tri22.X(list_edges_v2{i},2:3)); % mean(tri2.X(list_edges_v2{i},2:3)); % sagittal (y,z)
    mean_tmp3 = mean([tri32.X(list_edges_v2{i},1) tri32.X(list_edges_v2{i},3)]); %mean([tri3.X(list_edges_v2{i},1) tri3.X(list_edges_v2{i},3)]); % coronal (x,z)
    
    %% L1 norm
%     lapl_tri1 = abs(tri1.X(i,1:2) - mean_tmp1);
%     lapl_tri2 = abs(tri2.X(i,2:3) - mean_tmp2);
%     lapl_tri3 = abs([tri3.X(i,1) tri3.X(i,3)] - mean_tmp3);

    %% L2 norm
    lapl_tri1(i,:) = tri12.X(i,1:2) - mean_tmp1;%tri1.X(aux_array(i),1:2) - mean_tmp1;
    lapl_tri2(i,:) = tri22.X(i,2:3) - mean_tmp2;%tri2.X(aux_array(i),2:3) - mean_tmp2;
    lapl_tri3(i,:) = [tri32.X(i,1) tri32.X(i,3)] - mean_tmp3;%[tri3.X(aux_array(i),1) tri3.X(aux_array(i),3)] - mean_tmp3;
    
    %% Functional %%
    F_s( i )         =  lambda * sum( lapl_tri1(i,:).^2 ); % lambda * lapl_tri1; % 
    F_s(i + incr)    =  lambda * sum( lapl_tri2(i,:).^2 ); % lambda * lapl_tri2; % 
    F_s(i + 2*incr)  =  lambda * sum( lapl_tri3(i,:).^2 ); % lambda * lapl_tri3; % 
    
end


% F_obj = sum(F1.^2) + sum(F2.^2) + sum(F3.^2) + sum(F_s.^2);
% F_obj = (.5/( n1 + n2 + n3 ) ).* (sum(F1.^2) + sum(F2.^2) + sum(F3.^2) );

F_obj = (.5/( n1 + n2 + n3 + 3*incr ) ).* (sum(F1.^2) + sum(F2.^2) + sum(F3.^2) + sum(F_s));

if nargout > 1
    
    Js_int1 = zeros(n1, total_var); % size(var_array1,1)
    Js_int2 = zeros(n2, total_var);
    Js_int3 = zeros(n3, total_var);
    Js_s    = zeros(3*incr, total_var);

    incr_col = xyz * incr;
    
    %% First intersection %%
    for i = 1: n1
        
        %% axial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         inters = slice_int1(i);
        
        tmp_v1_ax = axial_M1{sub_3_int1_v(i)} * [b2c_ncoord1_int1(i,:) 1]'; % 3D point to 2D point in the frame coordinates axial_m1_par1{i}
        
        if crop_volume
            fl1 = floor(tmp_v1_ax(1) - crop_rectangle(1) + 1); % Cols
            fl2 = floor(tmp_v1_ax(2) - crop_rectangle(2) + 1); % Rows
            cl1 = ceil(tmp_v1_ax(1)  - crop_rectangle(1) + 1); % Cols
            cl2 = ceil(tmp_v1_ax(2)  - crop_rectangle(2) + 1); % Rows
        else
            fl1 = floor(tmp_v1_ax(1) + 1); % Cols
            fl2 = floor(tmp_v1_ax(2) + 1); % Rows
            cl1 = ceil(tmp_v1_ax(1)  + 1); % Cols
            cl2 = ceil(tmp_v1_ax(2)  + 1); % Rows
        end
        
        min_max_r1ax = min(max(fl2,1), r_ax);
        min_max_r2ax = min(max(cl2,1), r_ax);
        min_max_c1ax = min(max(fl1,1), c_ax);
        min_max_c2ax = min(max(cl1,1), c_ax);
        
        %% Gradient in axial view
        
        neigax_gradx = [gradx_ax(min_max_r1ax, min_max_c1ax, sub_3_int1_v(i)) gradx_ax(min_max_r1ax, min_max_c2ax, sub_3_int1_v(i));...
                        gradx_ax(min_max_r2ax, min_max_c1ax, sub_3_int1_v(i)) gradx_ax(min_max_r2ax, min_max_c2ax, sub_3_int1_v(i))];
        
        neigax_grady = [grady_ax(min_max_r1ax, min_max_c1ax, sub_3_int1_v(i)) grady_ax(min_max_r1ax, min_max_c2ax, sub_3_int1_v(i));...
                        grady_ax(min_max_r2ax, min_max_c1ax, sub_3_int1_v(i)) grady_ax(min_max_r2ax, min_max_c2ax, sub_3_int1_v(i))];
        
        grad_ax = [-bilinear_interpolation(tmp_v1_ax(2), tmp_v1_ax(1), double(neigax_gradx))  -bilinear_interpolation(tmp_v1_ax(2), tmp_v1_ax(1), double(neigax_grady))];
        
        %% sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        tmp_v1_sg = sag_M1{sub_2_int1_v(i)} * [b2c_ncoord2_int1(i,:) 1]'; % 3D point to 2D point in the frame coordinates var_cell{i,j,k} sag_m1_par1{i}
        
        if crop_volume
            fl1 = floor(tmp_v1_sg(1,:) - crop_rectangle(2) + 1);
            fl2 = floor(tmp_v1_sg(2,:) + 1);
            cl1 = ceil(tmp_v1_sg(1,:)  - crop_rectangle(2) + 1);
            cl2 = ceil(tmp_v1_sg(2,:)  + 1);
        else
            fl1 = floor(tmp_v1_sg(1,:) + 1);
            fl2 = floor(tmp_v1_sg(2,:) + 1);
            cl1 = ceil(tmp_v1_sg(1,:)  + 1);
            cl2 = ceil(tmp_v1_sg(2,:)  + 1);
        end
        
        min_max_r1sg = min(max(fl2,1), r_sag);
        min_max_r2sg = min(max(cl2,1), r_sag);
        min_max_c1sg = min(max(fl1,1), c_sag);
        min_max_c2sg = min(max(cl1,1), c_sag);
        
        %% Gradient in sagittal view
        neigsg_gradx = [gradx_sag(min_max_r1sg, min_max_c1sg, sub_2_int1_v(i)) gradx_sag(min_max_r1sg, min_max_c2sg, sub_2_int1_v(i));...
                        gradx_sag(min_max_r2sg, min_max_c1sg, sub_2_int1_v(i)) gradx_sag(min_max_r2sg, min_max_c2sg, sub_2_int1_v(i))];
        
        neigsg_grady = [grady_sag(min_max_r1sg, min_max_c1sg, sub_2_int1_v(i)) grady_sag(min_max_r1sg, min_max_c2sg, sub_2_int1_v(i));...
                        grady_sag(min_max_r2sg, min_max_c1sg, sub_2_int1_v(i)) grady_sag(min_max_r2sg, min_max_c2sg, sub_2_int1_v(i))];
        
        grad_sag = [-bilinear_interpolation(tmp_v1_sg(2), tmp_v1_sg(1), double(neigsg_gradx))  -bilinear_interpolation(tmp_v1_sg(2), tmp_v1_sg(1), double(neigsg_grady))];
        
        %% Jacobian
        w  = zeros(3, total_var);
        w2 = zeros(3, total_var);
        
        ind_w = (tri22.Triangulation(current_tr_int1_v2(i),:)-1)*xyz + 1;
        
        % axial (x,y)
        w(1,ind_w ) = c2b_coord_int1_v2(i,:); 
        w(2,:)  = circshift(w(1,:)', 1);
        % sagittal (y,z)
        w2(2,:) = circshift(w(1,:)', xyz * incr);
        w2(3,:) = circshift(w2(2,:)',1);
            
        Js_int1(i,:) = (F1(i)) .*  ( grad_ax * axial_M1{sub_3_int1_v(i)}(1:2,1:3) * w  - grad_sag * sag_M1{sub_2_int1_v(i)}(1:2,1:3) * w2 ); 
    end
    
    %% Second intersection %%
    for i = 1: n2
        
        %% axial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         inters = slice_int2(i);
        inters = i;
        tmp_v1_ax = axial_M1{sub_3_int2_v(inters)} * [b2c_ncoord1_int2(inters,:) 1]'; % 3D point to 2D point in the frame coordinates axial_m1_par2{i}

        if crop_volume
            fl1 = floor(tmp_v1_ax(1) - crop_rectangle(1) + 1); % Cols
            fl2 = floor(tmp_v1_ax(2) - crop_rectangle(2) + 1); % Rows
            cl1 = ceil(tmp_v1_ax(1)  - crop_rectangle(1) + 1); % Cols
            cl2 = ceil(tmp_v1_ax(2)  - crop_rectangle(2) + 1); % Rows
        else
            fl1 = floor(tmp_v1_ax(1) + 1); % Cols
            fl2 = floor(tmp_v1_ax(2) + 1); % Rows
            cl1 = ceil(tmp_v1_ax(1)  + 1); % Cols
            cl2 = ceil(tmp_v1_ax(2)  + 1); % Rows
        end
        
        min_max_r1ax = min(max(fl2,1), r_ax);
        min_max_r2ax = min(max(cl2,1), r_ax);
        min_max_c1ax = min(max(fl1,1), c_ax);
        min_max_c2ax = min(max(cl1,1), c_ax);
        
        %% Gradient in axial view
        
        neigax_gradx = [gradx_ax(min_max_r1ax, min_max_c1ax, sub_3_int2_v(inters)) gradx_ax(min_max_r1ax, min_max_c2ax, sub_3_int2_v(inters));...
                        gradx_ax(min_max_r2ax, min_max_c1ax, sub_3_int2_v(inters)) gradx_ax(min_max_r2ax, min_max_c2ax, sub_3_int2_v(inters))];
        
        neigax_grady = [grady_ax(min_max_r1ax, min_max_c1ax, sub_3_int2_v(inters)) grady_ax(min_max_r1ax, min_max_c2ax, sub_3_int2_v(inters));...
                        grady_ax(min_max_r2ax, min_max_c1ax, sub_3_int2_v(inters)) grady_ax(min_max_r2ax, min_max_c2ax, sub_3_int2_v(inters))];
        
        grad_ax = [-bilinear_interpolation(tmp_v1_ax(2), tmp_v1_ax(1), double(neigax_gradx))  -bilinear_interpolation(tmp_v1_ax(2), tmp_v1_ax(1), double(neigax_grady))];
        
        %% coronal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        tmp_v1_cr =  cor_M1{sub_2_int2_v(inters)} * [b2c_ncoord2_int2(inters,:) 1]'; % 3D point to 2D point in the frame coordinates var_cell{i,j,k}

        if crop_volume
            fl1 = floor(tmp_v1_cr(1,:) - crop_rectangle(1) + 1);
            fl2 = floor(tmp_v1_cr(2,:) + 1);
            cl1 = ceil(tmp_v1_cr(1,:)  - crop_rectangle(1) + 1);
            cl2 = ceil(tmp_v1_cr(2,:)  + 1);
        else
            fl1 = floor(tmp_v1_cr(1,:) + 1);
            fl2 = floor(tmp_v1_cr(2,:) + 1);
            cl1 = ceil(tmp_v1_cr(1,:)  + 1);
            cl2 = ceil(tmp_v1_cr(2,:)  + 1);
        end
        
        min_max_r1cr = min(max(fl2,1), r_cor);
        min_max_r2cr = min(max(cl2,1), r_cor);
        min_max_c1cr = min(max(fl1,1), c_cor);
        min_max_c2cr = min(max(cl1,1), c_cor);
        
        %% Gradient in coronal view
        
        neigcr_gradx = [gradx_cor(min_max_r1cr, min_max_c1cr, sub_2_int2_v(inters)) gradx_cor(min_max_r1cr, min_max_c2cr, sub_2_int2_v(inters));...
                        gradx_cor(min_max_r2cr, min_max_c1cr, sub_2_int2_v(inters)) gradx_cor(min_max_r2cr, min_max_c2cr, sub_2_int2_v(inters))];
        
        neigcr_grady = [grady_cor(min_max_r1cr, min_max_c1cr, sub_2_int2_v(inters)) grady_cor(min_max_r1cr, min_max_c2cr, sub_2_int2_v(inters));...
                        grady_cor(min_max_r2cr, min_max_c1cr, sub_2_int2_v(inters)) grady_cor(min_max_r2cr, min_max_c2cr, sub_2_int2_v(inters))];
        
        grad_cor = [-bilinear_interpolation(tmp_v1_cr(2), tmp_v1_cr(1), double(neigcr_gradx))  -bilinear_interpolation(tmp_v1_cr(2), tmp_v1_cr(1), double(neigcr_grady))];
                
        %% Jacobian
        w  = zeros(3, total_var);
        w2 = zeros(3, total_var);
        
        ind_w = (tri12.Triangulation(current_tr_int2_v2(inters),:)-1) * xyz + 1;
        
        % axial (x,y)
        w( 1,ind_w ) = c2b_coord_int2_v2(inters,:);
        w( 2,:) = circshift(w(1,:)', 1);
        % coronal (x,z)
        w2(1,:) = circshift(w(1,:)', 2 * xyz * incr );
        w2(3,:) = circshift(w2(1,:)',1);
      
        Js_int2(i,:) = (F2(i)) .*  ( - grad_cor * cor_M1{sub_2_int2_v(inters)}(1:2,1:3) * w2 + grad_ax * axial_M1{sub_3_int2_v(inters)}(1:2,1:3) * w );
        
    end
    %% Third intersection %%
    for i = 1: n3
        
        %% sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         inters = slice_int3(i);
        inters = i;
        tmp_v1_sg = sag_M1{sub_2_int3_v(inters)} * [b2c_ncoord1_int3(inters,:) 1]'; % 3D point to 2D point in the frame coordinates sag_m1_par2{i}
        
        if crop_volume
            fl1 = floor(tmp_v1_sg(1,:) - crop_rectangle(2) + 1);
            fl2 = floor(tmp_v1_sg(2,:) + 1);
            cl1 = ceil(tmp_v1_sg(1,:)  - crop_rectangle(2) + 1);
            cl2 = ceil(tmp_v1_sg(2,:)  + 1);
        else
            fl1 = floor(tmp_v1_sg(1,:) + 1);
            fl2 = floor(tmp_v1_sg(2,:) + 1);
            cl1 = ceil(tmp_v1_sg(1,:)  + 1);
            cl2 = ceil(tmp_v1_sg(2,:)  + 1);
        end
        
        min_max_r1sg = min(max(fl2,1), r_sag);
        min_max_r2sg = min(max(cl2,1), r_sag);
        min_max_c1sg = min(max(fl1,1), c_sag);
        min_max_c2sg = min(max(cl1,1), c_sag);
        
        %% Gradient in sagittal view
        
        neigsg_gradx = [gradx_sag(min_max_r1sg, min_max_c1sg, sub_2_int3_v(inters)) gradx_sag(min_max_r1sg, min_max_c2sg, sub_2_int3_v(inters));...
                        gradx_sag(min_max_r2sg, min_max_c1sg, sub_2_int3_v(inters)) gradx_sag(min_max_r2sg, min_max_c2sg, sub_2_int3_v(inters))];
        
        neigsg_grady = [grady_sag(min_max_r1sg, min_max_c1sg, sub_2_int3_v(inters)) grady_sag(min_max_r1sg, min_max_c2sg, sub_2_int3_v(inters));...
                        grady_sag(min_max_r2sg, min_max_c1sg, sub_2_int3_v(inters)) grady_sag(min_max_r2sg, min_max_c2sg, sub_2_int3_v(inters))];
        
        grad_sag = [-bilinear_interpolation(tmp_v1_sg(2), tmp_v1_sg(1), double(neigsg_gradx))  -bilinear_interpolation(tmp_v1_sg(2), tmp_v1_sg(1), double(neigsg_grady))];
        
        %% coronal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        tmp_v1_cr = cor_M1{sub_3_int3_v(inters)} * [b2c_ncoord2_int3(inters,:) 1]'; % 3D point to 2D point in the frame coordinates var_cell{i,j,k}cor_m1_par2{i}
        
        if crop_volume
            fl1 = floor(tmp_v1_cr(1,:) - crop_rectangle(1) + 1);
            fl2 = floor(tmp_v1_cr(2,:) + 1);
            cl1 = ceil(tmp_v1_cr(1,:)  - crop_rectangle(1) + 1);
            cl2 = ceil(tmp_v1_cr(2,:)  + 1);
        else
            fl1 = floor(tmp_v1_cr(1,:) + 1);
            fl2 = floor(tmp_v1_cr(2,:) + 1);
            cl1 = ceil(tmp_v1_cr(1,:)  + 1);
            cl2 = ceil(tmp_v1_cr(2,:)  + 1);
        end
        
        min_max_r1cr = min(max(fl2,1), r_cor);
        min_max_r2cr = min(max(cl2,1), r_cor);
        min_max_c1cr = min(max(fl1,1), c_cor);
        min_max_c2cr = min(max(cl1,1), c_cor);
        
        %% Gradient in coronal view
        
        neigcr_gradx = [gradx_cor(min_max_r1cr, min_max_c1cr, sub_3_int3_v(inters)) gradx_cor(min_max_r1cr, min_max_c2cr, sub_3_int3_v(inters));...
                        gradx_cor(min_max_r2cr, min_max_c1cr, sub_3_int3_v(inters)) gradx_cor(min_max_r2cr, min_max_c2cr, sub_3_int3_v(inters))];
        
        neigcr_grady = [grady_cor(min_max_r1cr, min_max_c1cr, sub_3_int3_v(inters)) grady_cor(min_max_r1cr, min_max_c2cr, sub_3_int3_v(inters));...
                        grady_cor(min_max_r2cr, min_max_c1cr, sub_3_int3_v(inters)) grady_cor(min_max_r2cr, min_max_c2cr, sub_3_int3_v(inters))];
        
        grad_cor = [-bilinear_interpolation(tmp_v1_cr(2), tmp_v1_cr(1), double(neigcr_gradx))  -bilinear_interpolation(tmp_v1_cr(2), tmp_v1_cr(1), double(neigcr_grady))];

        %% Jacobian
        w_tmp  = zeros(1, total_var);
        w  = zeros(3, total_var);
        w2 = zeros(3, total_var);
        
        ind_w = (tri32.Triangulation(current_tr_int3_v2(inters),:)-1) * xyz + 1;
        
        w_tmp( 1,ind_w ) = c2b_coord_int3_v2(inters,:);
        
        % sagittal (y,z)
        w(2,:) = circshift(w_tmp(1,:)', xyz * incr);
        w(3,:) = circshift(w(2,:)',  1);
        % coronal (x,z) 
        w2(1,:)  = circshift(w_tmp(1,:)', 2 * xyz * incr);
        w2(3,:) = circshift(w2(1,:)',1);
        
        Js_int3(i,:) = (F3(i)) .* ( - grad_sag * sag_M1{sub_2_int3_v(inters)}(1:2,1:3) * w + grad_cor * cor_M1{sub_3_int3_v(inters)}(1:2,1:3) * w2 ); %
        
    end
    
    %% Smooth term %%
    for i = 1:incr
        
        %% Gradient %%
        jsi1 = (i - 1) * xyz + 1;
        jsi2 = jsi1 + incr_col; % incr_col = xyz * incr, where incr = size( source_tri.X,1) & xyz # of variables
        jsi3 = jsi2 + incr_col;
        
        Js_s(i    ,        jsi1:jsi1+1) =   (lambda) .* lapl_tri1(i,:); %( sqrt_lambda / (sqrt ( sum(lapl_tri1(i,:).^2)) + .5) ) .* lapl_tri1(i,:); %(lambda) .* lapl_tri1(i,:); %sqrt_lambda; 
        Js_s(i + incr,     jsi2:jsi2+1) =   (lambda) .* lapl_tri2(i,:); %( sqrt_lambda / (sqrt ( sum(lapl_tri2(i,:).^2)) + .5) ) .* lapl_tri2(i,:); %(lambda) .* lapl_tri2(i,:); %sqrt_lambda; 
        Js_s(i + 2 * incr, jsi3:jsi3+1) =   (lambda) .* lapl_tri3(i,:); %( sqrt_lambda / (sqrt ( sum(lapl_tri3(i,:).^2)) + .5) ) .* lapl_tri3(i,:); %(lambda) .* lapl_tri3(i,:); %sqrt_lambda;  
       
        
        for j = 1:length(list_edges_v2{i})
            
%             if list_edges_v{aux_array(i)}(j) > min(aux_array) && list_edges_v{aux_array(i)}(j) < max(aux_array)
                
                edgein1 = (list_edges_v2{i}(j) - 1) * xyz + 1;

                edgein2 = edgein1 + incr_col;
                edgein3 = edgein2 + incr_col;
                
                Js_s(i   ,         edgein1:edgein1+1) = Js_s(i ,         jsi1:jsi1+1).*[scalar_v2(i) scalar_v2(i)];
                Js_s(i + incr,     edgein2:edgein2+1) = Js_s(i + incr,   jsi2:jsi2+1).*[scalar_v2(i) scalar_v2(i)];
                Js_s(i + 2 * incr, edgein3:edgein3+1) = Js_s(i + 2*incr, jsi3:jsi3+1).*[scalar_v2(i) scalar_v2(i)];
                
                %             Js_s(i   ,         edgein1:edgein1+2) = Js_s(i ,         jsi1:jsi1+2).*[scalar(i) scalar(i) scalar(i)];
                %             Js_s(i + incr,     edgein2:edgein2+2) = Js_s(i + incr,   jsi2:jsi2+2).*[scalar(i) scalar(i) scalar(i)];
                %             Js_s(i + 2 * incr, edgein3:edgein3+2) = Js_s(i + 2*incr, jsi3:jsi3+2).*[scalar(i) scalar(i) scalar(i)];
%             end
            
        end
        
    end

    J =  (1/(n1 + n2 + n3 + 3*incr) ).*(sum(Js_int1) + sum(Js_int2) + sum(Js_int3) + sum(Js_s));

end