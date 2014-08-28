function [F_obj J] = myfun_unc(mesh0_X)


global source_tri
global var_array1
global sag_m1
global axial_m1
global vol_ax
global vol_sag
global vol_cor
global sub_1
global sub_2_int1
global sub_3_int1
global current_tr_int1
global c2b_coord_int1
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
global scalar 

tri1 = TriRep( source_tri.Triangulation,mesh0_X(1:size( source_tri.X,1),:)); % define the new mesh for the axial
tri2 = TriRep( source_tri.Triangulation,mesh0_X(size( source_tri.X,1)+1:end,:)); % define the new mesh for the sagittal


b2c_ncoord1 = baryToCart(tri1,  current_tr_int1,  c2b_coord_int1);
b2c_ncoord2 = baryToCart(tri2,  current_tr_int1,  c2b_coord_int1);

lambda = .1;
sqrt_lambda = sqrt(lambda);

n = size(var_array1,1);

incr_row =  n + size(source_tri.X,1);
incr_col = 3 * size(source_tri.X,1);

F  = zeros(size(var_array1,1) + size(source_tri.X,1) + size(source_tri.X,1),1);
Js = zeros(size(var_array1,1) + size(source_tri.X,1) + size(source_tri.X,1), 6 * size(source_tri.X,1));

rows = size(vol_ax,1);
cols = size(vol_ax,2);

for i = 1: n
    
    %% axial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tmp_v1 = axial_m1{sub_3_int1(i)} * [b2c_ncoord1(i,:) 1]'; % 3D point to 2D point in the frame coordinates
    
    fl  = floor(tmp_v1(1) + 1);
    fl2 = floor(tmp_v1(2) + 1);
    cl  = ceil(tmp_v1(1) + 1);
    cl2 = ceil(tmp_v1(2) + 1);
    
    min_max_r1ax  = min(max(fl2,1),rows);
    min_max_r2ax = min(max(cl2,1),rows);
    min_max_c1ax  = min(max(fl,1), cols);
    min_max_c2ax = min(max(cl,1), cols);
    
    neig = [vol_ax(min_max_r1ax, min_max_c1ax, sub_3_int1(i))   vol_ax(min_max_r1ax, min_max_c2ax, sub_3_int1(i));...
            vol_ax(min_max_r2ax,min_max_c1ax, sub_3_int1(i))   vol_ax(min_max_r2ax,min_max_c2ax, sub_3_int1(i))];
    
    new_im_ax = bilinear_interpolation(tmp_v1(2),tmp_v1(1),double(neig));
    
    
    
    %% sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tmp_v1 = sag_m1{sub_2_int1(i)} * [b2c_ncoord2(i,:) 1]'; % 3D point to 2D point in the frame coordinates var_cell{i,j,k}

    fl  = floor(tmp_v1(1) + 1);
    fl2 = floor(tmp_v1(2) + 1);
    cl  = ceil(tmp_v1(1) + 1);
    cl2 = ceil(tmp_v1(2) + 1);

    min_max_r1sg  = min(max(fl2,1),rows);
    min_max_r2sg = min(max(cl2,1),rows);
    min_max_c1sg  = min(max(fl,1),cols);
    min_max_c2sg = min(max(cl,1),cols);
    
    neig = [vol_sag(min_max_r1sg, min_max_c1sg,sub_2_int1(i))   vol_sag(min_max_r1sg, min_max_c2sg,sub_2_int1(i));...
            vol_sag(min_max_r2sg,min_max_c1sg,sub_2_int1(i))   vol_sag(min_max_r2sg,min_max_c2sg,sub_2_int1(i))];
    
    new_im_sag = bilinear_interpolation(tmp_v1(2),tmp_v1(1),double(neig));
    
    
    %% Function
    F(i) = new_im_ax - new_im_sag;
    
    %% Gradient
    w  = zeros(3, 2 * incr_col);
    w2 = zeros(3, 2 * incr_col);
    
    %% Get the index of the points in the mesh (tetrahedron) that contains point i
    %% We multiply by 3 because we have 3 coordinates per each point
    ind_w1 = (tri1.Triangulation(current_tr_int1(i),:)-1)*3 + 1;
    
    w(1,ind_w1 )            = c2b_coord_int1(i,:);
    w2(1,:)   = circshift(w(1,:)',incr_col);

    w(2,:)  = circshift(w(1,:)', 1);
    w(3,:)  = circshift(w(1,:)', 2);
    w2(2,:) = circshift(w2(1,:)',1);
    w2(3,:) = circshift(w2(1,:)',2);
 
    grad_ax  = [gradx_ax( min_max_r2ax, min_max_c2ax, sub_3_int1(i)) grady_ax( min_max_r2ax, min_max_c2ax, sub_3_int1(i))];
    grad_sag = [gradx_sag(min_max_r2sg, min_max_c2sg, sub_2_int1(i)) grady_sag(min_max_r2sg, min_max_c2sg, sub_2_int1(i))];

    Js(i,:) = (grad_ax * axial_m1{sub_3_int1(i)}(1:2,1:3) * w) - (grad_sag * sag_m1{sub_2_int1(i)}(1:2,1:3) * w2);
    
end


for i = 1:size( source_tri.X,1)
    
    mean_tmp1 = mean(tri1.X(list_edges{i},:));
    mean_tmp2 = mean(tri2.X(list_edges{i},:));
    
    lapl_tri1 = abs(tri1.X(i,:) - mean_tmp1);
    lapl_tri2 = abs(tri2.X(i,:) - mean_tmp2);
    
    norm2_lapl_tri1 = sum( lapl_tri1 );
    norm2_lapl_tri2 = sum( lapl_tri2 );
    
    % TODO: bring out the tri1 - mean_tmp1
    F(i + n)        = sqrt_lambda * norm2_lapl_tri1;
    F(i + incr_row) = sqrt_lambda * norm2_lapl_tri2;
    
    jsi1 = (i - 1) * 3 + 1;
    jsi2 = jsi1 + incr_col;
    Js(i +    n    , jsi1:jsi1+2) = sqrt_lambda;%.* lapl_tri1; %(in.sqrt_lambda * ( 1 / norm_lapl_tri1 ) ) .* lapl_tri1; 
    Js(i + incr_row, jsi2:jsi2+2) = sqrt_lambda;%.* lapl_tri2; %(in.sqrt_lambda * ( 1 / norm_lapl_tri2 ) ) .* lapl_tri2; 
    
    for j = 1:length(list_edges{i})
        
        edgein1 = (list_edges{i}(j)-1) * 3 + 1;
        edgein2 = edgein1 + incr_col;
        Js(i + n,     edgein1:edgein1+2) = Js(i + n,     jsi1:jsi1+2).*[scalar(i) scalar(i) scalar(i)];
        Js(i + incr_row, edgein2:edgein2+2) = Js(i + incr_row, jsi2:jsi2+2).*[scalar(i) scalar(i) scalar(i)];
        
    end
    
    
%     F(i +  n)   =  sqrt_lambda*(sum( (tri1.X(i,:) - mean(tri1.X( list_edges{i},:))).^2 ));
%     F(i + incr) =  sqrt_lambda*(sum( (tri2.X(i,:) - mean(tri2.X( list_edges{i},:))).^2 ));
end

F_obj = sum(F.^2);
J = sum(Js);
