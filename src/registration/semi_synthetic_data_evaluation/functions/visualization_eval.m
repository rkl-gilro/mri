function visualization_eval( slices )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  Visualisation of the registration results with the original images,
%%  and using the view3d GUI
%%
%%  
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global vol_ax_eval
global vol_sag_eval
global vol_cor_eval

global new_axial
global new_sagittal
global new_coronal

global axial_M
global sag_M
global cor_M

global crop_rectangle
global crop_volume

rows_ax = size(vol_ax_eval,1);
cols_ax = size(vol_ax_eval,2);

rows_sag = size(vol_sag_eval,1);
cols_sag = size(vol_sag_eval,2);

rows_cor = size(vol_cor_eval,1);
cols_cor = size(vol_cor_eval,2);

%% Visualise the new images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
% Deformation Axial and Sagittal
subplot(131);
aa(:,:,1) = convert2u8(new_axial(:,:,slices(1)));
aa(:,:,2) = convert2u8(new_axial(:,:,slices(1)));
aa(:,:,3) = convert2u8(new_axial(:,:,slices(1)));

if crop_volume
    [X, Y, Z, triTexture] = compute_triTexture(axial_M{slices(1)},aa,[crop_rectangle(2)-1 crop_rectangle(4)-1], [crop_rectangle(1)-1 crop_rectangle(3)-1]);
else
    [X, Y, Z, triTexture] = compute_triTexture(axial_M{slices(1)},aa,[0 rows_ax-1],[0 cols_ax-1]);
end

surf(X,Y,Z,triTexture,'FaceColor','texturemap','EdgeColor','none');
hold on
axis equal

ss(:,:,2) = convert2u8(new_sagittal(:,:,slices(2)));
ss(:,:,1) = convert2u8(new_sagittal(:,:,slices(2)));
ss(:,:,3) = convert2u8(new_sagittal(:,:,slices(2)));

if crop_volume
    [X, Y, Z, triTexture] = compute_triTexture(sag_M{slices(2)},ss,[0 rows_sag-1], [crop_rectangle(2)-1 crop_rectangle(4)-1]);
else
    [X, Y, Z, triTexture] = compute_triTexture(sag_M{slices(2)},ss,[0 rows_sag-1],[0 cols_sag-1]);
end

surf(X,Y,Z,triTexture,'FaceColor','texturemap','EdgeColor','none');
hold on
axis equal
axis off; grid off

% Deformation Axial and Coronal
subplot(132);

aa(:,:,1) = convert2u8(new_axial(:,:,slices(1)));
aa(:,:,2) = convert2u8(new_axial(:,:,slices(1)));
aa(:,:,3) = convert2u8(new_axial(:,:,slices(1)));

if crop_volume
    [X, Y, Z, triTexture] = compute_triTexture(axial_M{slices(1)},aa,[crop_rectangle(2)-1 crop_rectangle(4)-1], [crop_rectangle(1)-1 crop_rectangle(3)-1]);
else
    [X, Y, Z, triTexture] = compute_triTexture(axial_M{slices(1)},aa,[0 rows_ax-1],[0 cols_ax-1]);
end

surf(X,Y,Z,triTexture,'FaceColor','texturemap','EdgeColor','none');
hold on
axis equal

cc(:,:,2) = convert2u8(new_coronal(:,:,slices(3)));
cc(:,:,1) = convert2u8(new_coronal(:,:,slices(3)));
cc(:,:,3) = convert2u8(new_coronal(:,:,slices(3)));

if crop_volume
    [X, Y, Z, triTexture] = compute_triTexture(cor_M{slices(3)},cc,[0 rows_cor-1], [crop_rectangle(1)-1 crop_rectangle(3)-1]);
else
    [X, Y, Z, triTexture] = compute_triTexture(cor_M{slices(3)},cc,[0 rows_cor-1], [0 cols_cor-1]);

end

surf(X,Y,Z,triTexture,'FaceColor','texturemap','EdgeColor','none');
hold on
axis equal
axis off; grid off

% Deformation Sagittal and Coronal
subplot(133);

cc(:,:,2) = convert2u8(new_coronal(:,:,slices(3)));
cc(:,:,1) = convert2u8(new_coronal(:,:,slices(3)));
cc(:,:,3) = convert2u8(new_coronal(:,:,slices(3)));

if crop_volume
    [X, Y, Z, triTexture] = compute_triTexture(cor_M{slices(3)},cc,[0 rows_cor-1], [crop_rectangle(1)-1 crop_rectangle(3)-1]);
else
    [X, Y, Z, triTexture] = compute_triTexture(cor_M{slices(3)},cc,[0 rows_cor-1], [0 cols_cor-1]);

end

surf(X,Y,Z,triTexture,'FaceColor','texturemap','EdgeColor','none');
hold on
axis equal


ss(:,:,2) = convert2u8(new_sagittal(:,:,slices(2)));
ss(:,:,1) = convert2u8(new_sagittal(:,:,slices(2)));
ss(:,:,3) = convert2u8(new_sagittal(:,:,slices(2)));

if crop_volume
    [X, Y, Z, triTexture] = compute_triTexture(sag_M{slices(2)},ss,[0 rows_sag-1], [crop_rectangle(2)-1 crop_rectangle(4)-1]);
else
    [X, Y, Z, triTexture] = compute_triTexture(sag_M{slices(2)},ss,[0 rows_sag-1],[0 cols_sag-1]);
end

surf(X,Y,Z,triTexture,'FaceColor','texturemap','EdgeColor','none');
hold on
axis equal
axis off; grid off

%% Visualize Original %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Original Axial and Sagittal
figure;
subplot(131);
aa(:,:,1) = convert2u8(vol_ax_eval(:,:,slices(1)));
aa(:,:,2) = convert2u8(vol_ax_eval(:,:,slices(1)));
aa(:,:,3) = convert2u8(vol_ax_eval(:,:,slices(1)));

if crop_volume
    [X, Y, Z, triTexture] = compute_triTexture(axial_M{slices(1)},aa,[crop_rectangle(2)-1 crop_rectangle(4)-1], [crop_rectangle(1)-1 crop_rectangle(3)-1]);
else
    [X, Y, Z, triTexture] = compute_triTexture(axial_M{slices(1)},aa,[0 rows_ax-1],[0 cols_ax-1]);
end

surf(X,Y,Z,triTexture,'FaceColor','texturemap','EdgeColor','none');
hold on
axis equal

ss(:,:,2) = convert2u8(vol_sag_eval(:,:,slices(2)));
ss(:,:,1) = convert2u8(vol_sag_eval(:,:,slices(2)));
ss(:,:,3) = convert2u8(vol_sag_eval(:,:,slices(2)));

if crop_volume
    [X, Y, Z, triTexture] = compute_triTexture(sag_M{slices(2)},ss,[0 rows_sag-1], [crop_rectangle(2)-1 crop_rectangle(4)-1]);
else
    [X, Y, Z, triTexture] = compute_triTexture(sag_M{slices(2)},ss,[0 rows_sag-1],[0 cols_sag-1]);
end

surf(X,Y,Z,triTexture,'FaceColor','texturemap','EdgeColor','none');
hold on
axis equal
axis off; grid off


% Original Coronal and Sagittal
subplot(132);
cc(:,:,2) = convert2u8(vol_cor_eval(:,:,slices(3)));
cc(:,:,1) = convert2u8(vol_cor_eval(:,:,slices(3)));
cc(:,:,3) = convert2u8(vol_cor_eval(:,:,slices(3)));

if crop_volume
    [X, Y, Z, triTexture] = compute_triTexture(cor_M{slices(3)},cc,[0 size(vol_cor_eval,1)-1], [crop_rectangle(1)-1 crop_rectangle(3)-1]);
else
    [X, Y, Z, triTexture] = compute_triTexture(cor_M{slices(3)},cc,[0 rows_cor-1], [0 cols_cor-1]);

end

surf(X,Y,Z,triTexture,'FaceColor','texturemap','EdgeColor','none');
hold on
axis equal


ss(:,:,2) = convert2u8(vol_sag_eval(:,:,slices(2)));
ss(:,:,1) = convert2u8(vol_sag_eval(:,:,slices(2)));
ss(:,:,3) = convert2u8(vol_sag_eval(:,:,slices(2)));

if crop_volume
    [X, Y, Z, triTexture] = compute_triTexture(sag_M{slices(2)},ss,[0 rows_sag-1], [crop_rectangle(2)-1 crop_rectangle(4)-1]);
else
    [X, Y, Z, triTexture] = compute_triTexture(sag_M{slices(2)},ss,[0 rows_sag-1],[0 cols_sag-1]);
end

surf(X,Y,Z,triTexture,'FaceColor','texturemap','EdgeColor','none');
hold on
axis equal
axis off; grid off

% Original Coronal and Axial
subplot(133);
aa(:,:,1) = convert2u8(vol_ax_eval(:,:,slices(1)));
aa(:,:,2) = convert2u8(vol_ax_eval(:,:,slices(1)));
aa(:,:,3) = convert2u8(vol_ax_eval(:,:,slices(1)));

if crop_volume
    [X, Y, Z, triTexture] = compute_triTexture(axial_M{slices(1)},aa,[crop_rectangle(2)-1 crop_rectangle(4)-1], [crop_rectangle(1)-1 crop_rectangle(3)-1]);
else
    [X, Y, Z, triTexture] = compute_triTexture(axial_M{slices(1)},aa,[0 rows_ax-1],[0 cols_ax-1]);
end

surf(X,Y,Z,triTexture,'FaceColor','texturemap','EdgeColor','none');
hold on
axis equal

cc(:,:,2) = convert2u8(vol_cor_eval(:,:,slices(3)));
cc(:,:,1) = convert2u8(vol_cor_eval(:,:,slices(3)));
cc(:,:,3) = convert2u8(vol_cor_eval(:,:,slices(3)));

if crop_volume
    [X, Y, Z, triTexture] = compute_triTexture(cor_M{slices(3)},cc,[0 size(vol_cor_eval,1)-1], [crop_rectangle(1)-1 crop_rectangle(3)-1]);
else
    [X, Y, Z, triTexture] = compute_triTexture(cor_M{slices(3)},cc,[0 rows_cor-1], [0 cols_cor-1]);

end
surf(X,Y,Z,triTexture,'FaceColor','texturemap','EdgeColor','none');
hold on
axis equal
axis off; grid off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Visualize the new images using the 3D viewer
% images.axial      = new_axial;
% images.axial_info = axial_M;
% 
% images.sagittal      = new_sagittal;
% images.sagittal_info = sag_M;
% 
% images.coronal      = new_coronal;
% images.coronal_info = cor_M;
% 
% view3d(images);


