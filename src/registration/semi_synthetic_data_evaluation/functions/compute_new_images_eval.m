%% VISUALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global target_tri_ax
global target_tri_sag
global target_tri_cor

global vol_ax_eval
global vol_sag_eval
global vol_cor_eval

global new_axial
global new_sagittal
global new_coronal

global axial_M
global axial_M1
global sag_M
global sag_M1
global cor_M
global cor_M1

global source_tri_v
global crop_rectangle
global crop_volume

size_source = size(source_tri_v,1);

new_axial    = zeros(size(vol_ax_eval));
new_sagittal = zeros(size(vol_sag_eval));
new_coronal  = zeros(size(vol_cor_eval));

%% Axial %%
rows_ax = size(vol_ax_eval,1);
cols_ax = size(vol_ax_eval,2);

tic
disp(['Axial number: ', num2str(size(vol_ax_eval,3))])

if crop_volume
    [j i] = meshgrid(crop_rectangle(1):crop_rectangle(3),crop_rectangle(2):crop_rectangle(4));
else
    [j i] = meshgrid(1:size(vol_ax_eval,2),1:size(vol_ax_eval,1));
end

tmp_i = i(:);
tmp_j = j(:);
    
for k_ax = 1:size(vol_ax_eval,3)
    
    p_3d = axial_M{k_ax} * [j(:)-1 i(:)-1 ones(length(j(:)),1)]'; %zeros(length(j(:)),1)
    
    current_tr = tsearchn(source_tri_v.X, source_tri_v.Triangulation,[p_3d(1,:)' p_3d(2,:)' p_3d(3,:)']); % calculate the tetrahedron where p_3d belongs
    
    ind_nan = find(~isnan(current_tr));
    
    c2b_coord  = cartToBary(source_tri_v, current_tr(ind_nan),[p_3d(1,ind_nan)' p_3d(2,ind_nan)' p_3d(3,ind_nan)']); %(ind_nan), ind_nan get the barycentric coordinates
    b2c_ncoord = baryToCart(target_tri_ax, current_tr(ind_nan), c2b_coord); % (ind_nan) get the cartesian coordinates
    
    tmp_v1_ax = axial_M1{k_ax} * [b2c_ncoord ones(length(ind_nan),1)]'; % 3D point to 2D point in the frame coordinates axial_m1{sub_3(i)}
    
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
    
    
    for in = 1:length(ind_nan) % ind_nan
        
        neigax = [vol_ax_eval(min_max_r1ax(in), min_max_c1ax(in), k_ax) vol_ax_eval(min_max_r1ax(in), min_max_c2ax(in), k_ax);...
                  vol_ax_eval(min_max_r2ax(in), min_max_c1ax(in), k_ax) vol_ax_eval(min_max_r2ax(in), min_max_c2ax(in), k_ax)];
        
%         new_axial(round(tmp_i(in)),round(tmp_j(in)),k_ax) %in -> ind_nan(in)
        if crop_volume      
            new_axial(round(i(ind_nan(in)) - crop_rectangle(2) + 1),round(j(ind_nan(in)) - crop_rectangle(1) + 1),k_ax) = bilinear_interpolation(tmp_v1_ax(2,in), tmp_v1_ax(1,in), double(neigax));
        else
            new_axial(round(i(ind_nan(in)) + 1),round(j(ind_nan(in)) + 1),k_ax) = bilinear_interpolation(tmp_v1_ax(2,in), tmp_v1_ax(1,in), double(neigax));
        end
    end
        

end
t = toc

%% Sagittal %%

rows_sag = size(vol_sag_eval,1);
cols_sag = size(vol_sag_eval,2);

tic
disp(['Sagittal number: ', num2str(size(vol_sag_eval,3))])

    
if crop_volume
    [j i] = meshgrid(crop_rectangle(2):crop_rectangle(4),1:size(vol_sag_eval,1));
    
else
    [j i] = meshgrid(1:size(vol_sag_eval,2),1:size(vol_sag_eval,1));
end
    

for k_sag = 1:size(vol_sag_eval,3)

    
    p_3d = sag_M{k_sag} * [j(:)-1 i(:)-1 ones(length(j(:)),1)]'; %zeros(length(j(:)),1)
    
    current_tr = tsearchn(source_tri_v.X, source_tri_v.Triangulation,[p_3d(1,:)' p_3d(2,:)' p_3d(3,:)']); % calculate the tetrahedron where p_3d belongs
    
    ind_nan = find(~isnan(current_tr));
    
    c2b_coord  = cartToBary(source_tri_v,current_tr(ind_nan),[p_3d(1,ind_nan)' p_3d(2,ind_nan)' p_3d(3,ind_nan)']); % get the barycentric coordinates
    b2c_ncoord = baryToCart(target_tri_sag, current_tr(ind_nan), c2b_coord); % get the cartesian coordinates
    
    tmp_v1_sag = sag_M1{k_sag} * [b2c_ncoord ones(length(ind_nan),1)]'; % 3D point to 2D point in the frame coordinates axial_m1{sub_3(i)}
    
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
    
    for in = 1:length(ind_nan)
        
        neigsag = [vol_sag_eval(min_max_r1sag(in), min_max_c1sag(in), k_sag) vol_sag_eval(min_max_r1sag(in), min_max_c2sag(in), k_sag);...
                   vol_sag_eval(min_max_r2sag(in), min_max_c1sag(in), k_sag) vol_sag_eval(min_max_r2sag(in), min_max_c2sag(in), k_sag)];
        if crop_volume
            new_sagittal(round(i(ind_nan(in))),round(j(ind_nan(in))- crop_rectangle(2) + 1),k_sag) = bilinear_interpolation(tmp_v1_sag(2,in), tmp_v1_sag(1,in), double(neigsag));
        else
            new_sagittal(round(i(ind_nan(in))),round(j(ind_nan(in)) + 1),k_sag) = bilinear_interpolation(tmp_v1_sag(2,in), tmp_v1_sag(1,in), double(neigsag));
        end
    end
    
end

t = toc

%% Coronal %%

rows_cor = size(vol_cor_eval,1);
cols_cor = size(vol_cor_eval,2);

tic
disp(['Coronal number: ', num2str(size(vol_cor_eval,3))])

if crop_volume
    [j i] = meshgrid(crop_rectangle(1):crop_rectangle(3),1:size(vol_cor_eval,1));
else
    [j i] = meshgrid(1:size(vol_cor_eval,2),1:size(vol_cor_eval,1));
end
    
for k_cor = 1:size(vol_cor_eval,3)
    
    p_3d = cor_M{k_cor} * [j(:)-1 i(:)-1 ones(length(j(:)),1)]'; %zeros(length(j(:)),1) 
    
    current_tr = tsearchn(source_tri_v.X,source_tri_v.Triangulation,[p_3d(1,:)' p_3d(2,:)' p_3d(3,:)']); % calculate the tetrahedron where p_3d belongs
    
    ind_nan = find(~isnan(current_tr));
    
    c2b_coord  = cartToBary(source_tri_v,current_tr(ind_nan),[p_3d(1,ind_nan)' p_3d(2,ind_nan)' p_3d(3,ind_nan)']); % get the barycentric coordinates
    b2c_ncoord = baryToCart(target_tri_cor, current_tr(ind_nan), c2b_coord); % get the cartesian coordinates
    
    tmp_v1_cor = cor_M1{k_cor} * [b2c_ncoord ones(length(ind_nan),1)]'; % 3D point to 2D point in the frame coordinates axial_m1{sub_3(i)}
    
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
    
    for in = 1:length(ind_nan)
        
        neigcor = [vol_cor_eval(min_max_r1cor(in), min_max_c1cor(in), k_cor) vol_cor_eval(min_max_r1cor(in), min_max_c2cor(in), k_cor);...
                   vol_cor_eval(min_max_r2cor(in), min_max_c1cor(in), k_cor) vol_cor_eval(min_max_r2cor(in), min_max_c2cor(in), k_cor)];
        if crop_volume
            new_coronal(round(i(ind_nan(in))),round(j(ind_nan(in)) - crop_rectangle(1) + 1),k_cor) = bilinear_interpolation(tmp_v1_cor(2,in), tmp_v1_cor(1,in), double(neigcor));
        else
            new_coronal(round(i(ind_nan(in))),round(j(ind_nan(in)) + 1),k_cor) = bilinear_interpolation(tmp_v1_cor(2,in), tmp_v1_cor(1,in), double(neigcor));
        end
    end
    
end

t = toc
