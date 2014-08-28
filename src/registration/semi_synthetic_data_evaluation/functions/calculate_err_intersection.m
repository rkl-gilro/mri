function error = calculate_err_intersection(num_slices1, num_slices2, size1, size2, inter_cell, discr_t, source, target_tri1, target_tri2, M1, M2, invM1, invM2, crop1, crop2, show)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Calculate the intersection error in the evaluation given two views
%%
%% Inputs:  1. num_slices1, num_slices2 ->
%%          2. size1, size2 ->
%%          3. inter_cell ->
%%          4. discr_t ->
%%          5. source, target_tri1, target_tri2 ->  
%%          6. M1, M2, invM1, invM2 ->
%%          7. crop1, crop2 -> 
%%          8. show ->
%%
%% Outputs: 1. error ->
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


error = zeros(num_slices1, num_slices2);

for k_slices1 = 1:num_slices1
    
    for k_slices2 = 1:num_slices2
        
    tmp  = [inter_cell{k_slices1, k_slices2, :}];
    tmp2 = reshape(tmp, 3, discr_t+1);
    
    %% First view 
    current_tr = tsearchn(target_tri1.X, target_tri1.Triangulation,[tmp2(1,:)' tmp2(2,:)' tmp2(3,:)']); % calculate the tetrahedron where p_3d belongs
    ind_nan1 = find(~isnan(current_tr));
    
    %% Second view
    current_tr = tsearchn(target_tri2.X, target_tri2.Triangulation,[tmp2(1,:)' tmp2(2,:)' tmp2(3,:)']); % calculate the tetrahedron where p_3d belongs
    ind_nan2 = find(~isnan(current_tr));
    
    %% Make sure the number of points is the same for both references
    %% If not, the smallest is taken
    if length(ind_nan1) < length(ind_nan2)
        ind_nan = ind_nan1;
    else 
        ind_nan = ind_nan2;
    end
    
    num_points = length(ind_nan);
    
    %% Convert the points into the reference coordinate 'source triangulation'
    c2b_coord  = cartToBary(target_tri1, current_tr(ind_nan),[tmp2(1,ind_nan)' tmp2(2,ind_nan)' tmp2(3,ind_nan)']); %(ind_nan), ind_nan get the barycentric coordinates
    b2c_ncoord_ax1 = baryToCart(source, current_tr(ind_nan), c2b_coord); % (ind_nan) get the cartesian coordinates
    
    %% Convert the 3D point into pixel coordinates
    tmp_v1_ax = invM1{k_slices1} * [b2c_ncoord_ax1 ones(num_points,1)]'; %   3D point to 2D point in the frame coordinates invM1{sub_3(i)}

    %% Make sure the indexes are correct and use bilinear interpolation
    if crop_volume
        fl1ax = floor(tmp_v1_ax(1,:) - crop1(1) + 1); % Cols
        fl2ax = floor(tmp_v1_ax(2,:) - crop1(2) + 1); % Rows
        cl1ax = ceil(tmp_v1_ax(1,:)  - crop1(1) + 1); % Cols  - crop_rectangle(1) 
        cl2ax = ceil(tmp_v1_ax(2,:)  - crop1(2) + 1); % Rows  - crop_rectangle(2) 
    else
        fl1ax = floor(tmp_v1_ax(1,:) + 1); % Cols
        fl2ax = floor(tmp_v1_ax(2,:) + 1); % Rows
        cl1ax = ceil(tmp_v1_ax(1,:)  + 1); % Cols
        cl2ax = ceil(tmp_v1_ax(2,:)  + 1); % Rows
    end
    
    min_max_r1ax = min(max(fl2ax,1), size1(1));
    min_max_r2ax = min(max(cl2ax,1), size1(1));
    min_max_c1ax = min(max(fl1ax,1), size1(2));
    min_max_c2ax = min(max(cl1ax,1), size1(2));
    
    new_coordinates_ax_int1{k_slices1,k_slices2} = [min_max_r2ax' min_max_c2ax'];
    
    res  = deformation_type.T_ax' * [new_coordinates_ax_int1{k_slices1,k_slices2}(:,1) new_coordinates_ax_int1{k_slices1,k_slices2}(:,2) k_slices1.*ones(num_points,1) ones(num_points,1)]';

    aa  = M1{k_slices1} * [res(2,:)'+ crop_rectangle(1) - 1  res(1,:)'+ crop_rectangle(2) - 1 ones(num_points,1)]'; 

    new_coordinates_ax2_int1{k_slices1,k_slices2} = aa;

    

    %% Sagittal
    current_tr = tsearchn(target_tri2.X, target_tri2.Triangulation,[tmp2(1,:)' tmp2(2,:)' tmp2(3,:)']); % calculate the tetrahedron where p_3d belongs
    
    ind_nan = find(~isnan(current_tr));
    
    c2b_coord  = cartToBary(target_tri2,current_tr(ind_nan),[tmp2(1,ind_nan)' tmp2(2,ind_nan)' tmp2(3,ind_nan)']); % get the barycentric coordinates
    b2c_ncoord_sg1 = baryToCart(source, current_tr(ind_nan), c2b_coord); % get the cartesian coordinates
    
    
    tmp_v1_sag  = invM2{k_slices2} * [b2c_ncoord_sg1 ones(num_points,1)]'; %    3D point to 2D point in the frame coordinates invM1{sub_3(i)}

    %% Make sure the indexes are correct and use bilinear interpolation
    
    if crop_volume
        fl1sag = floor(tmp_v1_sag(1,:) - crop2(1) + 1); % 
        fl2sag = floor(tmp_v1_sag(2,:) - crop2(2) + 1);
        cl1sag = ceil(tmp_v1_sag(1,:)  - crop2(1) + 1); % 
        cl2sag = ceil(tmp_v1_sag(2,:)  - crop2(2) + 1);
    else
        fl1sag = floor(tmp_v1_sag(1,:) + 1);
        fl2sag = floor(tmp_v1_sag(2,:) + 1);
        cl1sag = ceil(tmp_v1_sag(1,:)  + 1);
        cl2sag = ceil(tmp_v1_sag(2,:)  + 1);
    end
    
    min_max_r1sag = min(max(fl2sag,1), size2(1));
    min_max_r2sag = min(max(cl2sag,1), size2(1));
    min_max_c1sag = min(max(fl1sag,1), size2(2));
    min_max_c2sag = min(max(cl1sag,1), size2(2));
    
    new_coordinates_sag_int1{k_slices1,k_slices2} = [min_max_r2sag' min_max_c2sag'];
    
    res  = deformation_type.T_sag' * [new_coordinates_sag_int1{k_slices1,k_slices2}(:,1) new_coordinates_sag_int1{k_slices1,k_slices2}(:,2) k_slices2.*ones(num_points,1) ones(num_points,1)]';

    ss  = M2{k_slices2} * [res(2,:)'+ crop2(1) - 1  res(1,:)'- 1 ones(num_points,1)]'; % 1. 

    new_coordinates_sag2_int1{k_slices1,k_slices2} = ss;
    
    if show
        plot3(ss(1,1:10:end), ss(2,1:10:end) , ss(3,1:10:end), 'm+');hold on %  + crop_rectangle(2)
        plot3(aa(1,1:10:end), aa(2,1:10:end) , aa(3,1:10:end), 'g+');hold on
    end
    
    error(k_slices1, k_slices2) =  norm(ss-aa) / num_points;

    end
end