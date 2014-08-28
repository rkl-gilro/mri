function [straight_ax_int1, straight_sag_int1, straight_ax_int2, straight_cor_int2, straight_sag_int3, straight_cor_int3] = eval_straight( )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% 
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


global target_tri_ax
global target_tri_sag
global target_tri_cor

global vol_ax_eval
global vol_sag_eval
global vol_cor_eval

global var_cell1_v
global var_cell2_v
global var_cell3_v

global axial_M
global axial_M1
global sag_M
global sag_M1
global cor_M
global cor_M1

global source_tri_v
global crop_rectangle
global crop_volume

global optimizer
global deformation_type

%% First intersection %%
rows_ax = size(vol_ax_eval,1);
cols_ax = size(vol_ax_eval,2);

rows_sag = size(vol_sag_eval,1);
cols_sag = size(vol_sag_eval,2);

tic
disp('------- First Intersection: Axial & Sagittal ----------------')

for k_ax = 1:size(vol_ax_eval,3)
    
    for k_sag = 1:size(vol_sag_eval,3)
        
    tmp  = [var_cell1_v{k_ax, k_sag, :}];
    tmp2 = reshape(tmp,3,optimizer.t+1);
    
    %% Axial 
    current_tr = tsearchn(target_tri_ax.X, target_tri_ax.Triangulation,[tmp2(1,:)' tmp2(2,:)' tmp2(3,:)']); % calculate the tetrahedron where p_3d belongs
    ind_nan1 = find(~isnan(current_tr));
    %% Sagittal
    current_tr = tsearchn(target_tri_sag.X, target_tri_sag.Triangulation,[tmp2(1,:)' tmp2(2,:)' tmp2(3,:)']); % calculate the tetrahedron where p_3d belongs
    ind_nan2 = find(~isnan(current_tr));
    
    %% Make sure the number of points is the same for both references
    %% If not, the smallest is taken
    if length(ind_nan1) < length(ind_nan2)
        ind_nan = ind_nan1;
    else 
        ind_nan = ind_nan2;
    end
    
    num_points = length(ind_nan);
    
    c2b_coord  = cartToBary(target_tri_ax, current_tr(ind_nan),[tmp2(1,ind_nan)' tmp2(2,ind_nan)' tmp2(3,ind_nan)']); %(ind_nan), ind_nan get the barycentric coordinates
    b2c_ncoord_ax1 = baryToCart(source_tri_v, current_tr(ind_nan), c2b_coord); % (ind_nan) get the cartesian coordinates
    
    tmp_v1_ax = axial_M1{k_ax} * [b2c_ncoord_ax1 ones(num_points,1)]'; %   3D point to 2D point in the frame coordinates axial_m1{sub_3(i)}

    %% Make sure the indexes are correct and use bilinear interpolation
    if crop_volume
        fl1ax = floor(tmp_v1_ax(1,:)  - crop_rectangle(1) + 1); % Cols
        fl2ax = floor(tmp_v1_ax(2,:) - crop_rectangle(2) + 1); % Rows
        cl1ax = ceil(tmp_v1_ax(1,:)  - crop_rectangle(1)   + 1); % Cols 
        cl2ax = ceil(tmp_v1_ax(2,:) - crop_rectangle(2) + 1); % Rows  - crop_rectangle(2) 
    else
        fl1ax = floor(tmp_v1_ax(1,:) + 1); % Cols
        fl2ax = floor(tmp_v1_ax(2,:) + 1); % Rows
        cl1ax = ceil(tmp_v1_ax(1,:)  + 1); % Cols
        cl2ax = ceil(tmp_v1_ax(2,:)  + 1); % Rows
    end
    
    min_max_r1ax = min(max(fl2ax,1),rows_ax);%fl2ax;%
    min_max_r2ax = min(max(cl2ax,1),rows_ax);%cl2ax;%
    min_max_c1ax = min(max(fl1ax,1),cols_ax);%fl1ax;%
    min_max_c2ax = min(max(cl1ax,1),cols_ax);%cl1ax;%
    
    new_coordinates_ax_int1{k_ax,k_sag} = [min_max_r2ax' min_max_c2ax'];
    
    res  = deformation_type.T_ax' * [new_coordinates_ax_int1{k_ax,k_sag}(:,1) new_coordinates_ax_int1{k_ax,k_sag}(:,2) k_ax.*ones(num_points,1) ones(num_points,1)]';

    aa  = axial_M{k_ax} * [res(2,:)'+ crop_rectangle(1) - 1  res(1,:)'+ crop_rectangle(2) - 1 ones(num_points,1)]'; 

    new_coordinates_ax2_int1{k_ax,k_sag} = aa;


    %% Sagittal
    
    c2b_coord  = cartToBary(target_tri_sag,current_tr(ind_nan),[tmp2(1,ind_nan)' tmp2(2,ind_nan)' tmp2(3,ind_nan)']); % get the barycentric coordinates
    b2c_ncoord_sg1 = baryToCart(source_tri_v, current_tr(ind_nan), c2b_coord); % get the cartesian coordinates
    
    
    tmp_v1_sag  = sag_M1{k_sag} * [b2c_ncoord_sg1 ones(num_points,1)]'; %    3D point to 2D point in the frame coordinates axial_m1{sub_3(i)}

    %% Make sure the indexes are correct and use bilinear interpolation
    
    if crop_volume
        fl1sag = floor(tmp_v1_sag(1,:) - crop_rectangle(2)  + 1); % 
        fl2sag = floor(tmp_v1_sag(2,:) + 1);
        cl1sag = ceil(tmp_v1_sag(1,:)  - crop_rectangle(2) + 1); % 
        cl2sag = ceil(tmp_v1_sag(2,:)  + 1);
    else
        fl1sag = floor(tmp_v1_sag(1,:) + 1);
        fl2sag = floor(tmp_v1_sag(2,:) + 1);
        cl1sag = ceil(tmp_v1_sag(1,:)  + 1);
        cl2sag = ceil(tmp_v1_sag(2,:)  + 1);
    end
    
    min_max_r1sag = min(max(fl2sag,1),rows_sag);
    min_max_r2sag = min(max(cl2sag,1),rows_sag);
    min_max_c1sag = min(max(fl1sag,1),cols_sag);
    min_max_c2sag = min(max(cl1sag,1),cols_sag);
    
    new_coordinates_sag_int1{k_ax,k_sag} = [min_max_r2sag' min_max_c2sag'];
    
    res  = deformation_type.T_sag' * [new_coordinates_sag_int1{k_ax,k_sag}(:,1) new_coordinates_sag_int1{k_ax,k_sag}(:,2) k_sag.*ones(num_points,1) ones(num_points,1)]';

    ss  = sag_M{k_sag} * [res(2,:)' + crop_rectangle(2) - 1  res(1,:)'- 1 ones(num_points,1)]'; % 1. 

    new_coordinates_sag2_int1{k_ax,k_sag} = ss;


    %% Compute the straightness of the line   
    [~, ~, s] = best_fit_line(aa(1,:)',aa(2,:)',aa(3,:)');
    straight_ax_int1(k_ax, k_sag) = s;


    [~, ~, s] = best_fit_line(ss(1,:)',ss(2,:)',ss(3,:)');
    straight_sag_int1(k_ax, k_sag) = s;
    

    end
end


t = toc



%% Second intersection %%
rows_cor = size(vol_cor_eval,1);
cols_cor = size(vol_cor_eval,2);

tic
disp('------- Second Intersection: Axial & Coronal ----------------')

for k_ax = 1:size(vol_ax_eval,3)
    
    for k_cor = 1:size(vol_cor_eval,3)
        
    tmp  = [var_cell2_v{k_ax,k_cor,:}];
    tmp2 = reshape(tmp,3,optimizer.t+1);
    
    %% Axial 
    current_tr = tsearchn(target_tri_ax.X, target_tri_ax.Triangulation,[tmp2(1,:)' tmp2(2,:)' tmp2(3,:)']); % calculate the tetrahedron where p_3d belongs
    ind_nan1 = find(~isnan(current_tr));
    %% Coronal
    current_tr = tsearchn(target_tri_cor.X, target_tri_cor.Triangulation,[tmp2(1,:)' tmp2(2,:)' tmp2(3,:)']); % calculate the tetrahedron where p_3d belongs
    ind_nan2 = find(~isnan(current_tr));
    
    %% Make sure the number of points is the same for both references
    %% If not, the smallest is taken
    if length(ind_nan1) < length(ind_nan2)
        ind_nan = ind_nan1;
    else 
        ind_nan = ind_nan2;
    end
    
    num_points = length(ind_nan);
    
    c2b_coord  = cartToBary(target_tri_ax, current_tr(ind_nan),[tmp2(1,ind_nan)' tmp2(2,ind_nan)' tmp2(3,ind_nan)']); %(ind_nan), ind_nan get the barycentric coordinates
    b2c_ncoord_ax1 = baryToCart(source_tri_v, current_tr(ind_nan), c2b_coord); % (ind_nan) get the cartesian coordinates
    
    tmp_v1_ax = axial_M1{k_ax} * [b2c_ncoord_ax1 ones(num_points,1)]'; % 3D point to 2D point in the frame coordinates axial_m1{sub_3(i)}
    
    %% Make sure the indexes are correct and use bilinear interpolation
    if crop_volume
        fl1ax = floor(tmp_v1_ax(1,:) - crop_rectangle(1) + 1); % Cols
        fl2ax = floor(tmp_v1_ax(2,:) - crop_rectangle(2) + 1); % Rows
        cl1ax = ceil(tmp_v1_ax(1,:)  - crop_rectangle(1) + 1); % Cols
        cl2ax = ceil(tmp_v1_ax(2,:)  - crop_rectangle(2) + 1); % Rows
    else
        fl1ax = floor(tmp_v1_ax(1,:) + 1); % Cols
        fl2ax = floor(tmp_v1_ax(2,:) + 1); % Rows
        cl1ax = ceil(tmp_v1_ax(1,:)  + 1); % Cols
        cl2ax = ceil(tmp_v1_ax(2,:)  + 1); % Rows
    end
    
    min_max_r1ax = min(max(fl2ax,1),rows_ax);
    min_max_r2ax = min(max(cl2ax,1),rows_ax);
    min_max_c1ax = min(max(fl1ax,1),cols_ax);
    min_max_c2ax = min(max(cl1ax,1),cols_ax);
    
    new_coordinates_ax_int2{k_ax,k_cor} = [min_max_r1ax' min_max_c1ax'];
    res =  deformation_type.T_ax' *  [new_coordinates_ax_int2{k_ax,k_cor}(:,1) new_coordinates_ax_int2{k_ax,k_cor}(:,2) k_ax.*ones(num_points,1) ones(num_points,1)]';
%     res = new_coordinates_ax{k_ax,k_sag}';
    aa = axial_M{k_ax} * [res(2,:)' + crop_rectangle(1) - 1  res(1,:)' + crop_rectangle(2) - 1 ones(num_points,1)]'; % 1. , 2. 
    
    new_coordinates_ax2_int2{k_ax,k_cor} = aa;
    
    
    %% Coronal
 
    c2b_coord  = cartToBary(target_tri_cor,current_tr(ind_nan),[tmp2(1,ind_nan)' tmp2(2,ind_nan)' tmp2(3,ind_nan)']); % get the barycentric coordinates
    b2c_ncoord_cr1 = baryToCart(source_tri_v, current_tr(ind_nan), c2b_coord); % get the cartesian coordinates
    
%     plot3(b2c_ncoord_sg1(1:10:end,1), b2c_ncoord_sg1(1:10:end,2), b2c_ncoord_sg1(1:10:end,3), 'b+');hold on;
    
    tmp_v1_cor = cor_M1{k_cor} * [b2c_ncoord_cr1 ones(num_points,1)]'; % 3D point to 2D point in the frame coordinates axial_m1{sub_3(i)}
    
    %% Make sure the indexes are correct and use bilinear interpolation
    
    if crop_volume
        fl1cor = floor(tmp_v1_cor(1,:) - crop_rectangle(1) + 1);
        fl2cor = floor(tmp_v1_cor(2,:) + 1);
        cl1cor = ceil(tmp_v1_cor(1,:)  - crop_rectangle(1) + 1);
        cl2cor = ceil(tmp_v1_cor(2,:)  + 1);
    else
        fl1cor = floor(tmp_v1_cor(1,:) + 1);
        fl2cor = floor(tmp_v1_cor(2,:) + 1);
        cl1cor = ceil(tmp_v1_cor(1,:)  + 1);
        cl2cor = ceil(tmp_v1_cor(2,:)  + 1);
    end
    
    min_max_r1cor = min(max(fl2cor,1),rows_cor);
    min_max_r2cor = min(max(cl2cor,1),rows_cor);
    min_max_c1cor = min(max(fl1cor,1),cols_cor);
    min_max_c2cor = min(max(cl1cor,1),cols_cor);
    
    new_coordinates_cor_int2{k_ax,k_cor} = [min_max_r1cor' min_max_c1cor'];
    
    res =  deformation_type.T_cor' *  [new_coordinates_cor_int2{k_ax,k_cor}(:,1) new_coordinates_cor_int2{k_ax,k_cor}(:,2) k_cor.*ones(num_points,1) ones(num_points,1)]';
    %     res = new_coordinates_ax{k_ax,k_sag}';
    cc = cor_M{k_cor} * [res(2,:)' + crop_rectangle(1) - 1  res(1,:)'- 1 ones(num_points,1)]'; % 1. 
    
    new_coordinates_cor2_int2{k_ax,k_cor} = cc;

    
    %% Compute the straightness of the line   
    [~, ~, s] = best_fit_line(aa(1,:)',aa(2,:)',aa(3,:)');
    straight_ax_int2(k_ax, k_cor) = s;


    [~, ~, s] = best_fit_line(cc(1,:)',cc(2,:)',cc(3,:)');
    straight_cor_int2(k_ax, k_cor) = s;
    

    end
end

t = toc



%% Third intersection %%

tic
disp('------- Third Intersection: Coronal & Sagittal ----------------')


for k_cor = 1:size(vol_cor_eval,3)
    
    for k_sag = 1:size(vol_sag_eval,3)
        
    tmp  = [var_cell3_v{k_cor,k_sag,:}];
    tmp2 = reshape(tmp,3,optimizer.t+1);
    
    %% Coronal 
    current_tr = tsearchn(target_tri_cor.X, target_tri_cor.Triangulation,[tmp2(1,:)' tmp2(2,:)' tmp2(3,:)']); % calculate the tetrahedron where p_3d belongs
    ind_nan1 = find(~isnan(current_tr));
    %% Sagittal
    current_tr = tsearchn(target_tri_sag.X, target_tri_sag.Triangulation,[tmp2(1,:)' tmp2(2,:)' tmp2(3,:)']); % calculate the tetrahedron where p_3d belongs
    ind_nan2 = find(~isnan(current_tr));
    
    %% Make sure the number of points is the same for both references
    %% If not, the smallest is taken
    if length(ind_nan1) < length(ind_nan2)
        ind_nan = ind_nan1;
    else 
        ind_nan = ind_nan2;
    end
    
    num_points = length(ind_nan);
    
    c2b_coord  = cartToBary(target_tri_cor,current_tr(ind_nan),[tmp2(1,ind_nan)' tmp2(2,ind_nan)' tmp2(3,ind_nan)']); % get the barycentric coordinates
    b2c_ncoord_cr1 = baryToCart(source_tri_v, current_tr(ind_nan), c2b_coord); % get the cartesian coordinates
       
    tmp_v1_cor = cor_M1{k_cor} * [b2c_ncoord_cr1 ones(num_points,1)]'; % 3D point to 2D point in the frame coordinates axial_m1{sub_3(i)}
    
    %% Make sure the indexes are correct and use bilinear interpolation
    
    if crop_volume
        fl1cor = floor(tmp_v1_cor(1,:) - crop_rectangle(1) + 1);
        fl2cor = floor(tmp_v1_cor(2,:) + 1);
        cl1cor = ceil(tmp_v1_cor(1,:)  - crop_rectangle(1) + 1);
        cl2cor = ceil(tmp_v1_cor(2,:)  + 1);
    else
        fl1cor = floor(tmp_v1_cor(1,:) + 1);
        fl2cor = floor(tmp_v1_cor(2,:) + 1);
        cl1cor = ceil(tmp_v1_cor(1,:)  + 1);
        cl2cor = ceil(tmp_v1_cor(2,:)  + 1);
    end
    
    min_max_r1cor = min(max(fl2cor,1),rows_cor);
    min_max_r2cor = min(max(cl2cor,1),rows_cor);
    min_max_c1cor = min(max(fl1cor,1),cols_cor);
    min_max_c2cor = min(max(cl1cor,1),cols_cor);
    
    new_coordinates_cor_int3{k_cor,k_sag} = [min_max_r1cor' min_max_c1cor'];
    
    res =  deformation_type.T_cor' *  [new_coordinates_cor_int3{k_cor,k_sag}(:,1) new_coordinates_cor_int3{k_cor,k_sag}(:,2) k_cor.*ones(num_points,1) ones(num_points,1)]';
    %     res = new_coordinates_ax{k_ax,k_sag}';
    cc = cor_M{k_cor} * [res(2,:)' + crop_rectangle(1) - 1  res(1,:)'- 1 ones(num_points,1)]'; % 1. 
    
    new_coordinates_cor2_int3{k_cor,k_sag} = cc;

    
    %% Sagittal
    
    c2b_coord  = cartToBary(target_tri_sag,current_tr(ind_nan),[tmp2(1,ind_nan)' tmp2(2,ind_nan)' tmp2(3,ind_nan)']); % get the barycentric coordinates
    b2c_ncoord_sg1 = baryToCart(source_tri_v, current_tr(ind_nan), c2b_coord); % get the cartesian coordinates
        
    tmp_v1_sag = sag_M1{k_sag} * [b2c_ncoord_sg1 ones(num_points,1)]'; % 3D point to 2D point in the frame coordinates axial_m1{sub_3(i)}
    
    %% Make sure the indexes are correct and use bilinear interpolation
    
    if crop_volume
        fl1sag = floor(tmp_v1_sag(1,:) - crop_rectangle(2) + 1);
        fl2sag = floor(tmp_v1_sag(2,:) + 1);
        cl1sag = ceil(tmp_v1_sag(1,:)  - crop_rectangle(2) + 1);
        cl2sag = ceil(tmp_v1_sag(2,:)  + 1);
    else
        fl1sag = floor(tmp_v1_sag(1,:) + 1);
        fl2sag = floor(tmp_v1_sag(2,:) + 1);
        cl1sag = ceil(tmp_v1_sag(1,:)  + 1);
        cl2sag = ceil(tmp_v1_sag(2,:)  + 1);
    end
    
    min_max_r1sag = min(max(fl2sag,1),rows_sag);
    min_max_r2sag = min(max(cl2sag,1),rows_sag);
    min_max_c1sag = min(max(fl1sag,1),cols_sag);
    min_max_c2sag = min(max(cl1sag,1),cols_sag);
    
    new_coordinates_sag_int3{k_cor,k_sag} = [min_max_r1sag' min_max_c1sag'];
    
    res = deformation_type.T_sag' *  [new_coordinates_sag_int3{k_cor,k_sag}(:,1) new_coordinates_sag_int3{k_cor,k_sag}(:,2) k_sag.*ones(num_points,1) ones(num_points,1)]';
    %     res = new_coordinates_ax{k_ax,k_sag}';
    ss = sag_M{k_sag} * [res(2,:)' + crop_rectangle(2) - 1  res(1,:)'- 1 ones(num_points,1)]'; % 1. 
    
    new_coordinates_sag2_int3{k_cor,k_sag} = ss;
    
%     plot3(ss(1,1:10:end), ss(2,1:10:end), ss(3,1:10:end), 'm+');hold on
    
        %% Compute the straightness of the line   
    [~, ~, s] = best_fit_line(ss(1,:)',ss(2,:)',ss(3,:)');
    straight_sag_int3(k_cor, k_sag) = s;
    

    [~, ~, s] = best_fit_line(cc(1,:)',cc(2,:)',cc(3,:)');
    straight_cor_int3(k_cor, k_sag) = s;
    

    end
end

t = toc



